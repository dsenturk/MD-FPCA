%MultilevelFuncLong.m
%============
%Description:
%============
% Fits a MD-FPCA model to repeated multilevel functional data.
% MultilevelFuncLong.m assumes data is densely observed on a regular grid  
% in the functional domain and observed either with or without sparsity on  
% a regular grid in the longitudinal domain.
%
% Please refer to the MD-FPCA algorith in Web Appendix A for annotation of 
% MultilievelFuncLong.m steps.
%============
% Functions Implemented: 
%============
% MultilevelFPCA.m, MeanSmooth2D.m,mean2DKH.m,
% CovarianceSmooth.m, getRawCovKH.m, gcv_mullwlsn_MFPCA.m, ComputeScores.m,
% and PACE v. 2.17 (http://www.stat.ucdavis.edu/PACE/).
%
% The construction of the MD-FPCA algorithm relies on functions and 
% subroutines from the PACE package (Version 2.17; 2015).
%============
% Fits a MD-FPCA model to repeated multilevel functional data
%======
%Input:
%======
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
%      param_step1: A structure array containing the multilevel MFPCA parameters 
%                   for the first-stage fits. Parameters are described
%                   below.
%                              bwmu_Candidates:  C x 2 vector of bandwidth candidates, where C is number
%                                                of bandwidth candidates for consideration. The first entry 
%                                                for each row corresponds to longitudinal time s 
%                                                and second entry to functional time t. The final
%                                                bandwidth is selected by GCV from the rows of
%                                                candidates.
%                              bweta_Candidates: C x 2 vector of bandwidth candidates, where C is number
%                                                of bandwidth candidates for consideration. The first entry 
%                                                for each row corresponds to longitudinal time s 
%                                                and second entry to functional time t. The final
%                                                bandwidth is selected by GCV from the rows of
%                                                candidates.
%                              bwxcov_total:     1 x 2 vector, bandwidths for smoothing of the total covariance surface                                                 
%                                                \Sigma_T(t,t'|s)=cov(X_ij(t|s),X_ij(t'|s))
%                                                bwxcov_total(i): ith coordinate of bandwidth vector, 
%                                                i=1,2.
%                                                [0 0]: use GCV for selection of bandwidth candidates.                                                                                          
%                                                bwxcov_total(1)>0 & bwxcov_total(2)>0: user-specified 
%                                                bandwidths.
%                              bwxcov_Candidates_total: C x 2 vector of bandwidth candidates, where C is number
%                                                of bandwidth candidates for consideration. The final
%                                                bandwidth is selected by GCV from the rows of
%                                                candidates.
%                              bwxcov_between:   1*2 vector, bandwidths for covariance surface used 
%                                                for smoothing of the first level covariance surface                                                 
%                                                \Sigma^(1)(t,t'|s)=cov(X_ij(t|s),X_ij'(t'|s)), j not equal to j'                                               
%                                                bwxcov_between(i): ith coordinate of bandwidth 
%                                                vector, i=1,2.
%                                                [0 0]: use GCV for selection of bandwidth candidates.                                                                                           
%                                                bwxcov_total(1)>0 & bwxcov_total(2)>0: user-specified 
%                                                bandwidths.
%                              bwxcov_Candidates_between: C x 2 vector of bandwidth candidates, where C is number
%                                                of bandwidth candidates for consideration. The final
%                                                bandwidth is selected by GCV from the rows of
%                                                candidates.
%                              out1:            Functional time domain, the interval which defines t.
%                              ngrid:           integer, number of support points for the covariance 
%                                               surface in the GCV procedure (selecting bandwidths of 
%                                               covariance). Used for both the total and first level 
%                                               covariance smooths. Note, ngrid is automatically defined by
%                                               t by MultilevelFuncLong.m but can be changed within
%                                               MultilevelFPCA.m                                                                                                                                             
%                              FVE_threshold:   1 x 2 vector of positive numbers between 0 and 1
%                                               used to select the number of first-stage
%                                               eigenfunctions, K and P, that explain
%                                               at least FVE_threshold(1) percentage of 
%                                               variation and accounting for no greater than
%                                               FVE_threshold(2) percentage
%                                               of variation.
%                              no_opt_B:        User defined number of principal components for 
%                                               truncations of the first level covariance, K. Overrides
%                                               FVE_threshold when defined but is [] by default.%                                               
%                              no_opt_W:        User defined number of principal components for 
%                                               truncations of the second level covariance, P. Overrides
%                                               FVE_threshold when defined but is [] by default.
%                              kernel:          a character string to define the kernel to be used 
%                                               in the 1-D or 2-D
%                                               smoothing.
%                                               kernel = 'epan'  ==> Epanechnikov kernel
%                                               'rect'  ==> Rectangular kernel
%                                               'gauss'  ==> Gaussian kernel
%                                                Note: The Gaussian kernel is overall best for sparse 
%                                                  designs but is slower than the other kernels 
%                                                  and if computational speed is of concern then 
%                                                  one may wish to use the Epanechnikov kernel 
%                                                  also for the case of sparse designs.
%                              calc_scores:     0, will not calculate scores within FPCA procedures.                                                   
%                                               1, calculate scores within FPCA procedures..
%                                               Default is 0 and must be changed within the relevant
%                                               FPCA procedure.   
%                              calc_predictions: 0, will not calculate predictions                               
%                                                1, will calculate predictions using \xi_ij(s) and \zeta_ij(s)
%                                                   scores if calc_scores=1 within FPCA procedures.  
%                                               Default is 0 and must be changed within the relevant
%                                               FPCA procedure.                                               
%                              sort_data:       0, data is already sorted by user, do not sort data
%                                               within the MultilevelFPCA function. Default is 0.                                               
%                                               1, sort data within the MultilevelFPCA function.                                              
%                              verbose:         a character string
%                                               'on': display diagnostic messages.                                               
%                                               'off': suppress diagnostic messages.
%                              error:           0, no additional measurement error \sigma(s) assumed.
%                                               1, additional measurement error \sigma(s) assumed.
%       
%      param_step2: A structure array containing the FPCA parameters for the FPCA 
%                   function in the second-stage fits. Parameters are described
%                   below.
%
%                              bwmu:            scalar bwmu>0, bandwidth for mean curve mu.
%                              bwxcov:          1*2 vector, bandwidths for the smooth covariance surface. 
%                              selection_k:     the method of choosing the number of principal components K.
%                                               'AIC': use AIC criterion with pseudo-likelihood of
%                                               measurements (marginal likelihood).
%                                               'BIC': use BIC criterion with pseudo-likelihood of
%                                               measurements (marginal likelihood).                                                                                  
%                                               'FVE' (fraction of variance explained) : use scree plot
%                                               approach to select number of principal
%                                               components), see "FVE_threshold" below.
%                                               positive integer K: user-specified number of principal 
%                                               components.
%                                               Note: BIC and FVE produce the most parsimonious models.
%                              FVE_threshold:   a scalar between 0 and 1.
%                                               It is used to select the number of components k',p' that 
%                                               explain at least FVE_threshold percentage of variation. 
%                              error:           0, no additional measurement error assumed.
%                                               1, additional measurement error is assumed.
%                              verbose:         a character string
%                                               'on': display diagnostic messages.                                               
%                                               'off': suppress diagnostic messages.
%
%      bwmu:        1 x 2 vector used as bandwidths in the mean smooth for 
%                   the response surface. The first entry corresponds to longitudinal
%                   time s and second entry to functional time t. 
%                   [0 0]: implement selection of bandwidth candidates
%                   using GCV.
%                   else: perform two dimensional smoothing using supplied
%                   bandwidths.
%      bweta:       1 x 2 vector used as bandwidths in the subunit specific
%                   mean deviation smooth. The first entry corresponds to 
%                   longitudinal time s and second entry to functional 
%                   time t.
%                   [0 0]: implement selection of bandwidth candidates
%                   using GCV.
%                   else: perform two dimensional smoothing using supplied
%                   bandwidths.
%      eta0:        0, perform subunit mean smooths (default). 
%                   1, suppress subunit mean smooths.
%=======
%Output:
%=======
%      A 1 x 9 cell containing MDFPCA output:
%      mu:              an T x S matrix containing the global mean surface 
%                       \mu(t,s) across both functional time t and longitudinal 
%                       time s.
%      eta:             a 1 x J cell of the subunit-specific mean deviation 
%                       smooths \eta_j(t,s). Each cell entry is a T x S matrix.
%      MFPCA_step1:     a 1 x S cell of first-stage MFPCA fits. Each cell 
%                       contains the MFPCA output described in the 
%                       MultilevelFPCA function at the sth time point.
%                               - First and second level eigenfunctions.
%                               - First and second level eigenvalues.
%                               - First and second level eigenscores.
%                               - Measurement error variance.
%      FPCA_xi_step2:   a 1 x K cell (for each first-stage subject level component)
%                       containing the second-stage FPCA fits. Each cell contains
%                       FPCA output described in the FPCA function of the 
%                       PACE package.
%      FPCA_zeta_step2: a 1 x P cell (for each first-stage subunit level component)
%                       containing the second-stage FPCA fits. Each cell contains
%                       FPCA output described in the FPCA function of the 
%                       PACE package.
%      ypred_step1:     First-stage MFPCA model predictions for each of the 
%                       longitudinal time points s.
%                       Cell entries ypred_step1{s}{i}{j} contains 
%                       predictions across functional time t for subject i 
%                       at longitudinal time s and subunit j.
%      ypred_step2:     Second-stage model predictions for the entire dataset
%                       Cell entries ypred_step2{s}{i}{j} contains 
%                       predictions across functional time t for subject i 
%                       at longitudinal time s and subunit j using first-stage
%                       and second-stage parameter estimates.
%      bw_all:          Bandwidths used or selected by GCV for all
%                       first-stage and second-stage smooths.
%      OUTPUTnames:  entry names of output for reference
%
function [OUTPUT] = MultilevelFuncLong(y, t, level, id, param_step1, ...
                                       param_step2, ...
                                       bwmu, bweta, eta0)
                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Format input data

% Pool t, y, level and id cells into single vectors or matrices
tt = cell2mat(t);
yy = cell2mat(y);
levelmat = cell2mat(level);
idmat = cell2mat(id);

% For bandwidth output
bw_all = {};

% Store longitudinal and functional domains
out1 = unique(tt(1,:)); % Time grids for s
out2 = unique(tt(2,:)); % Time grids for t

% Store variation explained threshold for step1
FVE_threshold = param_step1.FVE_threshold(1);
param_step1.FVE_threshold = [1, 1]; % Accept all components in first-stage fit
param_step1.calc_scores = 0; % Do not calculate scores in first-stage fits
param_step1.calc_predictions = 0; % Do not calculate predictions from first-stage fits

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Estimate global mean \mu(t,s) across s and t with a two-dimensional smoother

fprintf(1,'Mean Smoothing: Global Mean Smooth \n');

% 2-D mean smooth
[mu, bw_mu] = MeanSmooth2D(y, t, out1, out2, bwmu, param_step1.bwmu_Candidates);

% Center trajectories by \mu(t,s) (used for calculation of sub-unit specific shifts)
ycenter1 = {};
for sub_ind = 1:length(y)
    [~, tobs1] = ismember(t{sub_ind}(1,:), out1);
    [~, tobs2] = ismember(t{sub_ind}(2,:), out2);
    tobs_ind = (tobs1-1)*size(mu, 1) + tobs2;
    ycenter1{sub_ind} = y{sub_ind} - mu(tobs_ind);  
end
bw_all{1} = bw_mu; % store bandwidth used for global mean smooths


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Estimate sub-unit specific shift \eta_j(t,s) across s and t with a two-dimensional smoother

fprintf(1,'Mean Smoothing: Subunit Surface Smooth \n');

% Subset centered data at each level for smoothing
uniq_level = unique(cell2mat(level));
for level_ind = 1:length(uniq_level)
    
    % suppress estimation of \eta_j(t,s)
    if eta0 == 1    
        bw_eta = [];
        eta{level_ind} = zeros(length(out2), length(out1));
       
    % default, estimate \eta_j(t,s)    
    else 
        for sub_ind = 1:length(y)
            bool_ind = (level{sub_ind} == uniq_level(level_ind));
            ylevel{level_ind}{sub_ind} = ycenter1{sub_ind}(bool_ind);
            tlevel{level_ind}{sub_ind} = t{sub_ind}(:,bool_ind);
            idlevel{level_ind}{sub_ind} = id{sub_ind}(:,bool_ind);
        end  % subject loop
        
        % Smooth \eta_j(t,s) surface
        [eta{level_ind}, bw_eta{level_ind}] = MeanSmooth2D(ylevel{level_ind}, tlevel{level_ind}, ...
            out1, out2, bweta, param_step1.bweta_Candidates);
        
    end
    
end  % level loop
bw_all{2} = bw_eta; % store bandwidth used for sub-unit specific smooths


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. First-stage Karhunen-Loeve decomposition 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3a-c. Fit multilevel FPCA (MFPCA) models at each longitudinal time s across t and estimate model eigencomponents

fprintf(1,'First-stage MFPCA: Model Fits Across t for each s \n');

% Set up MFPCA parameters
param_step1.out1 = out2; % define functional time domain for MFPCA
param_step1.ngrid = length(out2);
param_step1.sort_data = 0; % suppress sorting of data
% begin loop that performs MFPCA at each longitudinal time s across t
for s_ind = 1:length(out1) 
    
    % Set up MFPCA means
    param_step1.xmu = mu(:,s_ind)'; % extract \mu(t|s)
    for level_ind = 1:length(eta)
        param_step1.xeta(:,level_ind) = eta{level_ind}(:,s_ind); % extract \eta_j(t|s)
    end    
    
    % Set up MFPCA input data observed at a fixed longitudinal time
    y_step1 = {}; t_step1 = {}; level_step1 = {}; id_step1 = {};
    for sub_ind = 1:length(y)
        bool_ind = ismember(t{sub_ind}(1,:), out1(s_ind));
        if sum(bool_ind) > 0
            y_step1{length(y_step1) + 1} = y{sub_ind}(bool_ind);
            t_step1{length(t_step1) + 1} = t{sub_ind}(2, bool_ind);
            level_step1{length(level_step1) + 1} = level{sub_ind}(bool_ind);
            id_step1{length(id_step1) + 1} = id{sub_ind}(bool_ind);
        end
    end  % end subject loop
    id_all_step1{s_ind} = unique(cell2mat(id_step1));

    % Fit MFPCA for each s time across t time
    MFPCA_step1{s_ind} = MultilevelFPCA(y_step1, t_step1, level_step1, id_step1, param_step1);
    % Store fits for all components
    MFPCA_step1_allcomp{s_ind} = MFPCA_step1{s_ind};
    
    bw_step1_total{s_ind} = MFPCA_step1{s_ind}{18};
    bw_step1_between{s_ind} = MFPCA_step1{s_ind}{19};
    
end  % end first-stage loop over s time

% store bandwidths used for first-stage covariance smooths (records
% bandwidths selected by GCV if automatic bandwidth selection implemented)
bw_all{3} = bw_step1_total; 
bw_all{4} = bw_step1_between;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3d-e. Select number of first-stage subject (K) and subunit (P) eigenfunctions \phi^(1)_k(t|s) and \phi^(2)_p(t|s)

fprintf(1,'First-stage MFPCA: Selecting Number of Components \n');

% Select number of first level components K
for s_ind = 1:length(out1)
    lambda_tmp = MFPCA_step1{s_ind}{9}; % Extract First-stage subject eigenvalues \lambda^(1)_k(s)
    FVE = cumsum(lambda_tmp) / sum(lambda_tmp);
    K_tmp(s_ind) = find(FVE >= FVE_threshold, 1);
end
K = max(K_tmp);  % Number of first level components

% Select number of second level components P
for s_ind = 1:length(out1)
    lambda_tmp = MFPCA_step1{s_ind}{10}; % Extract First-stage subunit eigenvalues \lambda^(2)_p(s)
    FVE = cumsum(lambda_tmp) / sum(lambda_tmp);
    L_tmp(s_ind) = find(FVE >= FVE_threshold, 1);
end
L = max(L_tmp);  % Number of second level components P

for s_ind = 1:length(out1)
    num_compb(s_ind) = size(MFPCA_step1{s_ind}{11}, 2);
    num_compw(s_ind) = size(MFPCA_step1{s_ind}{12}, 2);
end

if sum(num_compb < K) > 0 % Check that K does not exceed number of components estimated at each MFPCA fit
    K = min(num_compb);
end
if sum(num_compw < L) > 0 % Check that P does not exceed number of components estimated at each MFPCA fit
    L = min(num_compw);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3f. Determine the sign of first-stage eigenfunctions \phi^(1)_k(t|s) and \phi^(2)_p(t|s) 
fprintf(1,'First-stage MFPCA: Sign of Eigenfunctions \n');

% Smoothness is insured by minimizing the distance between \phi(t|s and
% the preceding five eigenfunctions, {\phi^(1)_k(t|s')|s'=1,2,3,4,5}.

% Determine sign of first-stage subject eigenfunctions \phi^(1)_k(t|s) to insure smoothness in s
tn = length(out2);
m = length(out1);
phi_all_B = zeros(tn, m, K);
lambda_all_B = zeros(m, K);

% Subject level
for k=1:K
    for i = 1:m
        phi(:,i) = MFPCA_step1{i}{11}(:,k);
        lambda(i) = MFPCA_step1{i}{9}(k); 
    end
    
    for i=1:m
        
        if i==1
            phi(:,i)=phi(:,i);

        elseif i==2;

        if (norm(phi(:,i-1) - phi(:,i)) > norm(phi(:,i-1) + phi(:,i)))  
            phi(:,i)=-1*phi(:,i);
            
        end        


        elseif i==3;

        int1=[phi(:,i-2)-phi(:,i),phi(:,i-1)-phi(:,i)];
        int2=[phi(:,i-2)+phi(:,i),phi(:,i-1)+phi(:,i)];

        if (trapz(trapz(power(int1',2))) >  trapz(trapz(power(int2',2))))
            phi(:,i)=-1*phi(:,i);
            
        end 

        elseif i==4;

        int1=[phi(:,i-3)-phi(:,i),phi(:,i-2)-phi(:,i),phi(:,i-1)-phi(:,i)];
        int2=[phi(:,i-3)+phi(:,i),phi(:,i-2)+phi(:,i),phi(:,i-1)+phi(:,i)];

        if (trapz(trapz(power(int1',2))) >  trapz(trapz(power(int2',2))))
            phi(:,i)=-1*phi(:,i);
            
        end 

        elseif i==5;

        int1=[phi(:,i-4)-phi(:,i),phi(:,i-3)-phi(:,i),phi(:,i-2)-phi(:,i),phi(:,i-1)-phi(:,i)];
        int2=[phi(:,i-4)+phi(:,i),phi(:,i-3)+phi(:,i),phi(:,i-2)+phi(:,i),phi(:,i-1)+phi(:,i)];

        if (trapz(trapz(power(int1',2))) >  trapz(trapz(power(int2',2))))
            phi(:,i)=-1*phi(:,i);
            
        end 

        else

            int1=[phi(:,i-5)-phi(:,i),phi(:,i-4)-phi(:,i),phi(:,i-3)-phi(:,i),phi(:,i-2)-phi(:,i),phi(:,i-1)-phi(:,i)];
            int2=[phi(:,i-5)+phi(:,i),phi(:,i-4)+phi(:,i),phi(:,i-3)+phi(:,i),phi(:,i-2)+phi(:,i),phi(:,i-1)+phi(:,i)];

        if (trapz(trapz(power(int1',2))) >  trapz(trapz(power(int2',2))))
            phi(:,i)=-1*phi(:,i);
            
        end

        end
        
       
    end
    phi_all_B(:, :, k)=phi;
    lambda_all_B(:,k)=lambda;
end



% Determine sign of first-stage subunit eigenfunctions \phi^(2)_p(t|s) to insure smoothness in s
tn = length(out2);
m = length(out1);
phi_all_W = zeros(tn, m, L);
lambda_all_W = zeros(m, L);
% Subunit level
for k=1:L
    for i = 1:m
        phi(:,i) = MFPCA_step1{i}{12}(:,k);
        lambda(i) = MFPCA_step1{i}{10}(k); 
    end
    
    for i=1:m
        
        if i==1
            phi(:,i)=phi(:,i);

        elseif i==2;

        if (norm(phi(:,i-1) - phi(:,i)) > norm(phi(:,i-1) + phi(:,i)))  
            phi(:,i)=-1*phi(:,i);
            
        end        


        elseif i==3;

        int1=[phi(:,i-2)-phi(:,i),phi(:,i-1)-phi(:,i)];
        int2=[phi(:,i-2)+phi(:,i),phi(:,i-1)+phi(:,i)];

        if (trapz(trapz(power(int1',2))) >  trapz(trapz(power(int2',2))))
            phi(:,i)=-1*phi(:,i);
            
        end 

        elseif i==4;

        int1=[phi(:,i-3)-phi(:,i),phi(:,i-2)-phi(:,i),phi(:,i-1)-phi(:,i)];
        int2=[phi(:,i-3)+phi(:,i),phi(:,i-2)+phi(:,i),phi(:,i-1)+phi(:,i)];

        if (trapz(trapz(power(int1',2))) >  trapz(trapz(power(int2',2))))
            phi(:,i)=-1*phi(:,i);
            
        end 

        elseif i==5;

        int1=[phi(:,i-4)-phi(:,i),phi(:,i-3)-phi(:,i),phi(:,i-2)-phi(:,i),phi(:,i-1)-phi(:,i)];
        int2=[phi(:,i-4)+phi(:,i),phi(:,i-3)+phi(:,i),phi(:,i-2)+phi(:,i),phi(:,i-1)+phi(:,i)];

        if (trapz(trapz(power(int1',2))) >  trapz(trapz(power(int2',2))))
            phi(:,i)=-1*phi(:,i);
            
        end 

        else

            int1=[phi(:,i-5)-phi(:,i),phi(:,i-4)-phi(:,i),phi(:,i-3)-phi(:,i),phi(:,i-2)-phi(:,i),phi(:,i-1)-phi(:,i)];
            int2=[phi(:,i-5)+phi(:,i),phi(:,i-4)+phi(:,i),phi(:,i-3)+phi(:,i),phi(:,i-2)+phi(:,i),phi(:,i-1)+phi(:,i)];

        if (trapz(trapz(power(int1',2))) >  trapz(trapz(power(int2',2))))
            phi(:,i)=-1*phi(:,i);
            
        end

        end
        
       
    end
    phi_all_W(:, :, k)=phi;
    lambda_all_W(:,k)=lambda;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3g. Calculate first-stage subject and subunit eigenscores, \xi_{ik}(s) and \zeta}_{ijp}(s), using multilevel BLUP
fprintf(1,'First-stage MFPCA: Calculating Scores \n');

for s_ind = 1:length(out1) % begin loop to estimate \xi_{ik}(s) and \zeta}_{ijp}(s) at each longitudinal time s across t
    
    % Set up MFPCA means
    mu_s = mu(:,s_ind)';
    eta_s = [];
    for level_ind = 1:length(eta)
        eta_s(:,level_ind) = eta{level_ind}(:,s_ind);
    end
    
    % Set up input data
    y_step1 = {}; t_step1 = {}; level_step1 = {}; id_step1 = {};
    for sub_ind = 1:length(y)
        bool_ind = ismember(t{sub_ind}(1,:), out1(s_ind));
        if sum(bool_ind) > 0
            y_step1{length(y_step1) + 1} = y{sub_ind}(bool_ind);
            t_step1{length(t_step1) + 1} = t{sub_ind}(2, bool_ind);
            level_step1{length(level_step1) + 1} = level{sub_ind}(bool_ind);
            id_step1{length(id_step1) + 1} = id{sub_ind}(bool_ind);
        end
    end  % end subject loop
    

    % Convert to vector/matrix form
    tt_step1 = cell2mat(t_step1);
    yy_step1 = cell2mat(y_step1);
    levelmat_step1 = cell2mat(level_step1);
    idmat_step1 = cell2mat(id_step1);
    
    % Center trajectories by \mu(t|s)
    [~, tobs] = ismember(tt_step1, out2);
    yycenter1 = yy_step1 - mu_s(tobs);
    
    % Center trajectories by \mu(t|s) and \eta_j(t|s)
    uniq_level = unique(levelmat_step1);
    yycenter2 = [];
    for i = 1:length(uniq_level)
        bool_ind = (levelmat_step1 == uniq_level(i));
        [~, tobs2] = ismember(tt_step1(bool_ind), out2);
        yycenter2(bool_ind) = yycenter1(bool_ind) - eta_s(tobs2,i)';
    end
    
    % Organize input for ComputeScores function
    lambda_B_s = lambda_all_B(s_ind, :); % First-stage subject eigenvalues \lambda^(1)_k(s)
    lambda_W_s = lambda_all_W(s_ind, :); % First-stage subunit eigenvalues \lambda^(2)_p(s)
    phi_B_s = reshape(phi_all_B(:,s_ind,:), length(out2), K); % First-stage subject eigenfunctions \phi^(1)_k(t|s)
    phi_W_s = reshape(phi_all_W(:,s_ind,:), length(out2), L); % First-stage subunit eigenfunctions \phi^(2)_p(t|s)
    sigma_s = MFPCA_step1{s_ind}{8}; % Measurement error \sigma(s)
    
    % Compute BLUP estimates of first and second level eigenscores,\xi_{ik}(s) and \zeta}_{ijp}(s), using multilevel BLUP
    [xi_step1 zeta_step1] = ComputeScores(yycenter2, tt_step1, idmat_step1, levelmat_step1, ...
        lambda_B_s, lambda_W_s, phi_B_s, phi_W_s, sigma_s);
    
    % Calculate first-stage model predictions using multilevel BLUP score estimates
    uniq_id = unique(idmat_step1);
    uniq_level = unique(levelmat_step1);
    for i = 1:length(uniq_id)
        for j = 1:length(uniq_level)
            pred{j} = mu_s' + eta_s(:,j) + phi_B_s * xi_step1(i,:)' + phi_W_s * zeta_step1{j}(i,:)';
        end
        ypred{i} = pred;
    end
    
    % Store first-stage model predictions for each longitudinal time s
    ypred_step1{s_ind} = ypred;

    % Store all first-stage MFPCA output
    MFPCA_step1{s_ind}{9} = lambda_B_s; % First-stage subject eigenvalues \lambda^(1)_k(s)
    MFPCA_step1{s_ind}{10} = lambda_W_s; % First-stage subunit eigenvalues \lambda^(2)_p(s)
    MFPCA_step1{s_ind}{11} = phi_B_s; % First-stage subject eigenfunctions \phi^(1)_k(t|s)
    MFPCA_step1{s_ind}{12} = phi_W_s; % First-stage subunit eigenfunctions \phi^(2)_p(t|s)
    MFPCA_step1{s_ind}{13} = xi_step1; % First-stage subject eigenscores \xi_{ik}(s)
    MFPCA_step1{s_ind}{14} = zeta_step1; % First-stage subunit eigenscores \zeta}_{ijp}(s)

end  % end first-stage loop over s time


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Second-stage Karhunen-Loeve decomposition 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4a,b,d. Fit FPCA to first level eigenscores \xi_{ik}(s) across s
fprintf(1,'Second-stage FPCA: Model First-stage Subject Eigenscores Across s \n');

% Set up data for First-stage subject eigenscore fits
param_step2.ngrid = length(out1);
param_step2.xmu = zeros(1, length(out1));
param_step2.maxk = param_step2.ngrid - 2;
for K_ind = 1:K
    
    y_xi_step2 = cell(1, length(y)); 
    t_xi_step2 = cell(1, length(y));
    for s_ind = 1:length(out1)
        
        xi_tmp = getVal(MFPCA_step1{s_ind}, 'xi');        
        for sub_ind = 1:size(xi_tmp, 1)
            bool_ind = find(ismember(unique(idmat), id_all_step1{s_ind}(sub_ind)));
            y_xi_step2{bool_ind}(end + 1) = xi_tmp(sub_ind, K_ind);
            t_xi_step2{bool_ind}(end + 1) = out1(s_ind);
        end  % end sub_ind loop
        
    end  % end s time loop
    
    % Fit First-stage subject FPCA for each component k
    [T, FPCA_xi_step2{K_ind}] = evalc('FPCA(y_xi_step2, t_xi_step2, param_step2);');
    bw_step2_xi{K_ind} = FPCA_xi_step2{K_ind}{12};
    
end  % end K_ind loop

bw_all{5} = bw_step2_xi; % store second-stage smoothing bandwidth used for covariance surface of \xi_{ik}(s) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4a,c,d. Fit FPCA to second level eigenscores scores \zeta_{ijp}(s) across s
fprintf(1,'Second-stage FPCA: Model First-stage Subunit Eigenscores Across s \n');

% Set up data for first-stage subunit eigenscore fits
uniq_level = unique(levelmat);
J = length(uniq_level);
for L_ind = 1:L
    
    y_zeta_step2 = cell(1, length(y)*J);
    t_zeta_step2 = cell(1, length(y)*J);
    for s_ind = 1:length(out1)
        
        zeta_tmp = getVal(MFPCA_step1{s_ind}, 'zeta');
        for sub_ind = 1:size(zeta_tmp{1}, 1)
            for level_ind = 1:length(uniq_level)
                
                bool_ind = find(ismember(unique(idmat), id_all_step1{s_ind}(sub_ind)));
                y_zeta_step2{J*(bool_ind - 1) + level_ind}(end + 1) = zeta_tmp{level_ind}(sub_ind, L_ind);
                t_zeta_step2{J*(bool_ind - 1) + level_ind}(end + 1) = out1(s_ind);
            end  % end level_ind loop
        end  % end sub_ind loop
        
    end  % end s time loop
    
    % Omit data that contains a zero. Zeros means the subunit is missing
    % in the data and this will skew the fits
    for i = 1:length(y_zeta_step2)
        y_zeta_step2{i} = y_zeta_step2{i}(~(y_zeta_step2{i} == 0));
        t_zeta_step2{i} = t_zeta_step2{i}(~(y_zeta_step2{i} == 0));
    end
        
    % Fit First-stage subunit FPCA for each component p
    [T,FPCA_zeta_step2{L_ind}] = evalc('FPCA(y_zeta_step2, t_zeta_step2, param_step2);');
    bw_step2_zeta{L_ind} = FPCA_zeta_step2{L_ind}{12};
    
end  % end L_ind loop
bw_all{6} = bw_step2_zeta; % store second-stage smoothing bandwidth used for covariance surface of \zeta}_{ijp}(s) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4e. Estimate second-stage eigenscores \xi'_{ikk'} and \zeta'_{ijpp'}

fprintf(1,'Second-stage FPCA: Predictions \n');

for s_ind = 1:length(out1)
    
    
    mu_s = mu(:,s_ind)';  % extract \mu(t|s)
    eta_s = []; 
    for level_ind = 1:length(eta) % extract \eta_j(t|s) 
        eta_s(:,level_ind) = eta{level_ind}(:, s_ind);
    end
    phi_B_s = reshape(phi_all_B(:,s_ind,:), length(out2), K);
    phi_W_s = reshape(phi_all_W(:,s_ind,:), length(out2), L);
    
    % Calculate second-stage predictions
    for i = 1:length(y)
        for j = 1:length(uniq_level)
            
            % Extract second-stage score predictions
            xi_pred = [];
            for K_ind = 1:K
                xi_pred(end + 1) = FPCA_xi_step2{K_ind}{17}{i}(s_ind);
            end
            zeta_pred = [];
            for L_ind = 1:L
                zeta_pred(end + 1) = FPCA_zeta_step2{L_ind}{17}{J*(i - 1) + j}(s_ind);
            end
            
            % Calculate prediction
            pred{j} = mu_s' + eta_s(:,j) + phi_B_s * xi_pred' + phi_W_s * zeta_pred';
        end
        ypred_s{i} = pred;
    end
    
    % Store second-stage prediction
    ypred_step2{s_ind} = ypred_s;
    
end % end s_ind scores loop


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Organize output

OUTPUTnames = {'mu', 'eta', 'MFPCA_step1', 'FPCA_xi_stage2', ...
               'FPCA_zeta_stage2', 'ypred_stage1', ...
               'ypred_stage2','bw_all','names'};
bwnames = {'bw_mu', 'bw_eta', 'bw_stage1_total', 'bw_stage1_level_one', ...
           'bw_stage2_xi','bw_stage2_zeta'};
bw_all{7} = bwnames;
OUTPUT = {mu, eta, MFPCA_step1, FPCA_xi_step2, FPCA_zeta_step2, ...
          ypred_step1, ypred_step2, bw_all, OUTPUTnames};



end % end MultilevelFuncLong function