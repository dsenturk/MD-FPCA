% MD_FPCA_Tutorial.m
%============
% Description:
%============
% Step-by-step implementation of MD-FPCA algorithm using the
% MultilevelFuncLong.m function. 
% 
% MultilevelFuncLong.m assumes data is densely observed on a regular grid  
% in the functional domain and observed either with or without sparsity on  
% a regular grid in the longitudinal domain.
%============
% Functions Implemented: 
%============
% MultilevelFuncLong.m, MeanSmooth2D.m, mean2DKH.m, CovarianceSmooth.m, 
% getRawCovKH.m, gcv_mullwlsn_MFPCA.m, ComputeScores.m,
% and PACE v. 2.17 (http://www.stat.ucdavis.edu/PACE/).
%
% The construction of the MD-FPCA algorithm relies on functions and 
% subroutines from the PACE package (Version 2.17; 2015).
%============
% Tutorial Outline:
%============
% 0. Load .m files and data into working directory.
% 1. Import multidimensional functional data and define file paths.
% 2. Format data for MultilevelFuncLong.m.
% 3. Define input parameters for MultilevelFuncLong.m.
% 4. Implement MultilevelFuncLong.m.
% 5. Visualization of MD-FPCA results.


% Please refer to the MD-FPCA algorithm in Web Appendix A for annotation of 
% MultilievelFuncLong.m steps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Load .m files and data into working directory

% Place MD-FPCA matlab .m files and PACE v. 2.17 (http://www.stat.ucdavis.edu/PACE/)
% into working directory, along with the multilevel_func_data.mat file. 
% These will be used to implement and demonstrate the MD-FPCA algorithm. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Import multidimensional functional data

% Set file path
cd('.../Biometrics Code');
addpath(genpath('.../Biometrics Code'));

% Load example dataset 
load('multilevel_func_data.mat');

% The example dataset is simulated as described in Web Appendix D for N=30
% under a sparse design and high SNR with the following components:

    % Create array of 50 equispaced points on [0,1]
    out1=linspace(0,1,50);
    out2=linspace(0,1,50);
    
    % Create equispaced grid on [0,1]x[0,1]
    [S,T]= meshgrid(linspace(0,1,50));

    % Calculate symmetric mean surface
    mutrue = 10*sqrt(1-(S-.5).^2-(T-.5).^2);

    % Set electrode specific means to zero for the moment
    for j=1:4
    etatrue{j}=zeros(50,50);
    end

    % Stage-one \phi eigensurfaces
    % \phi^(1)_k(t,s)
    varphi1B=sqrt(2).*cos(pi.*(T-S)); % k=1
    varphi2B=sqrt(2).*cos(3*pi.*(T-S)); % k=2

    % \phi^(1)_p(t,s)
    varphi1W=sqrt(2).*sin(pi.*(T-S)); % p=1
    varphi2W=sqrt(2).*sin(3*pi.*(T-S)); % p=2

    % Stage-two \psi eigenfunctions
    % \psi^(1)_kk'(s)
    varpsi11B=sqrt(2).*sin(2*pi.*(out1)); %k=1,k'=1
    varpsi12B=sqrt(2).*cos(2*pi.*(out1)); %k=1,k'=2

    varpsi21B=sqrt(2).*cos(4*pi.*(out1)); %k=2,k'=1
    varpsi22B=sqrt(2).*sin(4*pi.*(out1)); %k=2,k'=2

   % \psi^(1)_pp'(s)
    varpsi11W=sqrt(2).*sin(2*pi.*(out1)); %p=1,p'=1
    varpsi12W=sqrt(2).*cos(2*pi.*(out1)); %p=1,p'=2

    varpsi21W=sqrt(2).*cos(4*pi.*(out1)); %p=2,p'=1
    varpsi22W=sqrt(2).*sin(4*pi.*(out1)); %p=2,p'=2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Format data for MultilevelFuncLong.m

% Store longitudinal time s domain
out1 = unique(multilevel_func_data.longitudinal)';
% Store functional time t domain
out2 = unique(multilevel_func_data.functional);
% Store unique id numbers of subject reference
ids = unique(multilevel_func_data.subject);

% Input variable descriptions:
%
%      y:           1 x N cell containing data on each subject and subunit.
%                   The ith cell entry contains the trajectory across both
%                   functional and longitudinal dimensions for the ith
%                   subject at each subunit j.
%      t:           1 x N cell containing the functional-longitudinal time 
%                   pairs corresponding to the response trajectories in X_ij(t|s) 
%                   for each subject and subunit. Each cell entry is a 
%                   matrix with 2 rows. The number of columns is the total
%                   number of observations for that subject. The first row 
%                   of the matrix corresponds to longitudinal time s and 
%                   second to functional time t.
%      level:       1 x N cell containing the subunit index of functional
%                   repetitions on each subject. Entries correspond to the
%                   entries of y and t pairs.
%      id:          1 x N cell containing the subject ids. Entries
%                   correspond to the entries of y and t pairs.

clearvars y t s level id

% Package input data 

for sub_ind = 1:length(ids);
    bool_id = (multilevel_func_data.subject == ids(sub_ind));
    dat_id = multilevel_func_data(bool_id,:);
    y{sub_ind} = dat_id.y';
    t{sub_ind} = [dat_id.longitudinal, dat_id.functional]';
    level{sub_ind} = dat_id.subunit';
    id{sub_ind} = dat_id.subject';
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Define input parameters for MultilevelFuncLong.m

% (a) Parameter settings for user defined bandwidths

    % Set up first-stage MFPCA parameters for fits across t
    param_step1 = struct('bwmu_Candidates', [], ...
    'bweta_Candidates', [], ...
    'bwxcov_total', [.1 .1], 'bwxcov_Candidates_total', [], ...
    'bwxcov_between', [.1 .1], 'bwxcov_Candidates_between', [], ...
    'FVE_threshold', [0.90 1],'error', 1,...
    'kernel', 'gauss', 'verbose', 'off');

    % Set up second-stage FPCA parameters for fits across s
    param_step2 = setOptions('bwmu', .05, 'bwxcov', [.05 .05], 'selection_k', 'FVE', ...
    'FVE_threshold', 0.90, 'verbose', 'off', 'error', 1);

    % Set bandwidth parameters for global mean smooth
    bwmu = [.2, .2]; % bwmu(1) is s time, bwmu(2) is t time
    % Set up bandwith parameters for subunit-specific deviation smooths
    bweta = [.2, .2]; % bweta(1) is s time, bweta(2) is t time

    % Perform subunit smooths
    eta_switch=0;


% (b) Parameter settings for bandwidth selection by GCV

    % Set up first-stage MFPCA parameters for fits across t
    param_step1 = struct('bwmu_Candidates', [.1 .1; .2 .2; .3 .3], ...
    'bweta_Candidates', [.1 .1; .2 .2; .3 .3], ...
    'bwxcov_total', [0 0], 'bwxcov_Candidates_total', [.1 .1; .2 .2; .3 .3], ...
    'bwxcov_between', [0 0], 'bwxcov_Candidates_between', [.1 .1; .2 .2; .3 .3], ...
    'FVE_threshold', [0.90 1],'error', 1,...
    'kernel', 'gauss', 'verbose', 'off');


    % Set up second-stage FPCA parameters for fits across s
    param_step2 = setOptions('bwmu', .05, 'bwxcov', [.05 .05], 'selection_k', 'FVE', ...
    'FVE_threshold', 0.90, 'verbose', 'off', 'error', 1);

    % Set bandwidth parameters for global mean smooth \mu(t,s)
    bwmu = [0 0]; % bwmu(1) is s time, bwmu(2) is t time
    % Set up bandwith parameters for subunit-specific deviation smooths
    % \eta_j(t,s)
    bweta = [0 0]; % bweta(1) is s time, bweta(2) is t time

    % Perform subunit smooths
    eta_switch=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Implement MultilevelFuncLong.m 

% MultilevelFuncLong.m algorithm
fit30 = MultilevelFuncLong(y, t, level, id, ...
    param_step1, param_step2, bwmu, bweta,eta_switch);

% Note: The MD-FPCA algorithm for the sample dataset takes approximately 15
% minutes to run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Visualization of MD-FPCA results

% (a) Extract data from fit30

    % Store model components
    MD_FPCA_components=getVal(fit30,'names'); % Model components 
    mu=getVal(fit30,'mu'); % \mu(t,s)
    eta=getVal(fit30,'eta'); % \eta_j(t,s)
    MFPCA_step1_30 =getVal(fit30,'MFPCA_step1'); % Stage-one fits
    FPCA_xi_step2_30 = getVal(fit30,'FPCA_xi_stage2'); % Stage-two subject fits
    FPCA_zeta_step2_30 = getVal(fit30,'FPCA_zeta_stage2'); % \Stage-two subunit fits

    % Extract \phi surfaces
    phi_all_30_B = zeros(50, 50, 2); 
    for s_ind =1:50
    phi_all_30_B(:,s_ind,1) = MFPCA_step1_30{s_ind}{11}(:,1);
    if size(MFPCA_step1_30{s_ind}{11}, 2) > 1
        phi_all_30_B(:,s_ind,2) = MFPCA_step1_30{s_ind}{11}(:,2);
    end
    end

    phi_all_30_W = zeros(50, 50, 2);
    for s_ind =1:50
    phi_all_30_W(:,s_ind,1) = MFPCA_step1_30{s_ind}{12}(:,1);
    if size(MFPCA_step1_30{s_ind}{12}, 2) > 1
        phi_all_30_W(:,s_ind,2) = MFPCA_step1_30{s_ind}{12}(:,2);
    end
    end

    % Smooth \phi surfaces
    kernel = 'epan';
    xin = GetMesh(out1, out2); win =  ones(1, 50*50); bw_mu = [.05, 0.05];
    for k = 1:2
        yin = reshape(phi_all_30_B(:,:,k), 50*50, 1);
        [invalid phi_all_30_B(:,:,k)] = mullwlsk_M(bw_mu, kernel, xin, yin, win, out1, out2);
    end

    xin = GetMesh(out1, out2); win =  ones(1, 50*50); bw_mu = [.1, 0.1];
    for k = 1:2
        yin = reshape(phi_all_30_W(:,:,k), 50*50, 1);
        [invalid phi_all_30_W(:,:,k)] = mullwlsk_M(bw_mu, kernel, xin, yin, win, out1, out2);
    end

% (b) Plot model components

    % \mu(t,s)
    % True
    surf(mutrue);
    % Estimated
    surf(mu);

    % \eta
    % True
    surf(etatrue{1});
    surf(etatrue{2});
    surf(etatrue{3});
    surf(etatrue{4});
    % Estimated
    surf(eta{1});
    surf(eta{2});
    surf(eta{3});
    surf(eta{4});

    % \phi^(1)_k(t,s)
    % True
    surf(varphi1B); % k=1
    surf(varphi2B); % k=2
    % Estimated
    surf(-1*phi_all_30_B(:,:,2));% k=1
    surf(-1*phi_all_30_B(:,:,2));% k=2

    % \phi^(2)_p(t,s)
    % True
    surf(varphi1W); % p=1
    surf(varphi2W); % p=2
    % Estimated
    surf(phi_all_30_W(:,:,1)) % p=1
    surf(phi_all_30_W(:,:,2)) % p=2

    % \psi^(1)_kk'(s)
    % True (solid) and estimated (dashed)
    hold on; 
    plot(varpsi11B,'k');
    plot(FPCA_xi_step2_30{1}{4}(:,1),'k--'); %k=1,k'=1
    hold off;

    hold on;
    plot(varpsi12B,'k');
    plot(-1*FPCA_xi_step2_30{1}{4}(:,2),'k--'); %k=1,k'=2 
    hold off;

    hold on;
    plot(varpsi21B,'k');
    plot(-1*FPCA_xi_step2_30{2}{4}(:,1),'k--'); %k=2,k'=1
    hold off;

    hold on;
    plot(varpsi22B,'k');
    plot(FPCA_xi_step2_30{2}{4}(:,2),'k--'); %k=2,k'=2
    hold off;


    % \psi^(2)_pp'(s) 
    % True (solid) and estimated (dashed)
    hold on; 
    plot(varpsi11W,'k');
    plot(FPCA_zeta_step2_30{1}{4}(:,1),'k--'); %p=1,p'=1
    hold off;

    hold on;
    plot(varpsi12W,'k');
    plot(FPCA_zeta_step2_30{1}{4}(:,2),'k--'); %p=1,p'=2
    hold off;

    hold on;
    plot(varpsi21W,'k');
    plot(-1*FPCA_zeta_step2_30{2}{4}(:,1),'k--'); %p=2,p'=1
    hold off;

    hold on;
    plot(varpsi22W,'k');
    plot(FPCA_zeta_step2_30{2}{4}(:,2),'k--'); %p=2,p'=2
    hold off;