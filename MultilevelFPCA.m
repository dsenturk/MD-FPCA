%MultilevelFPCA.m
%============
%Description:
%============
%Fits a multilevel FPCA model to multilevel functional data 
%============
% Functions Implemented: 
%============
% MeanSmooth2D.m,mean2DKH.m, CovarianceSmooth.m, getRawCovKH.m,
% gcv_mullwlsn_MFPCA.m, ComputeScores.m,
% and PACE v. 2.17 (http://www.stat.ucdavis.edu/PACE/).
%
% The construction of the MD-FPCA algorithm relies on functions and 
% subroutines from the PACE package (Version 2.17; 2015).
%======
%Input:
%======
%      y:           1 x N cell containing data on each subject and level.
%                   The ith cell entry contains the trajectory for the ith
%                   subject at each repetition/level.
%      t:           1 x N cell containing the time points corresponding to
%                   the response trajectories in y for each subject and
%                   level.
%      level:       1 x N cell containing the level index of functional
%                   repetitions on each subject. Entries correspond to the
%                   entries of y.
%      id:          1 x N cell containing the subject ids. Entries
%                   correspond to the entries of y.
%      param:       A structure array containing the FPCA parameters for the 
%                   MultilevelFPCA function.  These are described below:
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
%                                                \Sigma_T(t,t'|s)=cov(X_ij(t|s),X_ij(t'|s)).
%                                                bwxcov_total(i): ith coordinate of bandwidth vector, 
%                                                i=1,2.
%                                                [0 0]: use GCV for selection of bandwidth candidates                                                                                             
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
%                                               surface in the CV procedure (selecting bandwidths of 
%                                               covariance). Used for both the total and first level 
%                                               covariance smooths.
%                              FVE_threshold:   1 x 2 vector of positive numbers between 0 and 1
%                                               used to select the number of first-stage
%                                               eigenfunctions, K and P, that explain
%                                               at least FVE_threshold(1) percentage of 
%                                               variation and accounting for no greater than
%                                               FVE_threshold(2) percentage of variation.
%                              no_opt_B:        User defined number of principal components for 
%                                               truncations of the first level covariance, K.            
%                              no_opt_W:        User defined number of principal components for 
%                                               truncations of the second level covariance, P.
%                              kernel:          a character string to define the kernel to be used 
%                                               in the 1-D or 2-D smoothing.
%                                               kernel = 'epan'  ==> Epanechnikov kernel
%                                               'rect'  ==> Rectangular kernel
%                                               'gauss'  ==> Gaussian kernel
%                                                Note: The Gaussian kernel is overall best for sparse 
%                                                  designs but is slower than the other kernels 
%                                                  and if computational speed is of concern then 
%                                                  one may wish to use the Epanechnikov kernel 
%                                                  also for the case of sparse designs.
%                              calc_scores:     0, will not calculate scores.
%                                               1, calculate scores.
%                              calc_predictions: 0, will not calculate predictions.
%                                                1, will calculate predictions using \xi_ij(s) and \zeta_ij(s)
%                                                   scores if calc_scores=1.
%                              xmu:             1 x T vector, user-defined mean function \mu(t|s) valued at 
%                                               distinct input time points, in ascending order from 
%                                               all subjects, corresponding to unique time points. 
%                              xeta:            T x J vector, user-defined \eta_j(t|s) function valued at 
%                                               distinct input time points, in ascending order from 
%                                               all subjects, corresponding to unique time points. 
%                              sort_data:       0, data is already sorted by user, do not sort data
%                                               withing the MultilevelFPCA function
%                                               1, sort data within the MultilevelFPCA function
%                              verbose:         a character string
%                                               'on': display diagnostic messages
%                                               'off': suppress diagnostic messages
%                              error:           0, no additional measurement error \sigma(s) assumed.
%                                               1, additional measurement error \sigma(s) assumed.%=======
%Output:
%=======
%      A 1 x 20 cell containing multilevel FPCA output:
%      mu:           T x 1 vector, mean function \mu(t|s) used to
%                    center observations for estimation of covariance
%                    surfaces.
%      eta:          T x J matrix, sub-unit shift \eta_j(t|s) used to 
%                    center observations for estimation of covariance
%                    surfaces.
%      xcov_total:   T x T matrix, smoothed total covariance surface
%                    \Sigma_T(t,t'| s).
%      xcov_between: T x T matrix, smoothed first level covariance surface
%                    \Sigma^(1)(t,t'| s).
%      xcov_within:  T x T matrix, smoothed second level covariance surface
%                    \Sigma^(2)(t,t'| s).
%      no_opt_B:     User defined number of principal components for 
%                    truncations of the first level covariance, K.            
%      no_opt_W:     User defined number of principal components for 
%                    truncations of the second level covariance, P.
%      sigma:        scalar, estimate of measurement error variance
%                    \sigma(s) if error=1, while it is [] if error=0. 
%      lambda_B:     1 x K vector, first level eigenvalues \lambda^(1)_k(s).
%      lambda_W:     1 x P vector, second level eigenvalues \lambda^(2)_p(s).
%      phi_B:        T x K matrix, first level eigenfunctions \phi^(1)_k(t|s).
%      phi_W:        T x P matrix, second level eigenfunctions \phi^(2)_p(t|s).
%      xi:           N x K matrix of first level eigenscores \xi_{ik}(s) if 
%                    calc_scores=1, else [].
%      zeta:         1 x J cell each containing a N x P matrix of second 
%                    level eigenscores \zeta}_{ijp}(s) if calc_scores=1, else [].                    
%      ypred:        1 x N cell array, ypred{i} is the 1 x 4 cell array of predictions 
%                    for the ith subject for all subunits evaluated at all 
%                    unique functional time points defined by out if
%                    calc_predictions=1, else [].
%      bw_xcov_total: 1 x 2 vector of bandwidths used for total covariance
%                     smooth.
%      bw_xcov_between: 1 x 2 vector of bandwidths used for first level
%                       covariance smooth.
%      OUTPUTnames:  entry names of output for reference.
%
function [OUTPUT] = MultilevelFPCA(y, t, level, id, param) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store all functional parameters from param
bwxcov_total = param.bwxcov_total;
bwxcov_between = param.bwxcov_between;
out1=param.out1; ngrid = param.ngrid;
FVE_threshold = param.FVE_threshold;
error = param.error; kernel = param.kernel;
verbose = param.verbose;
calc_scores=param.calc_scores;
calc_predictions=param.calc_predictions;
xmu=param.xmu; xeta=param.xeta;
bwmu_Candidates=param.bwmu_Candidates;
bweta_Candidates=param.bweta_Candidates;
bwxcov_Candidates_total=param.bwxcov_Candidates_total;
bwxcov_Candidates_between=param.bwxcov_Candidates_between;
no_opt_B=[];
no_opt_W=[];
% Sort data
if param.sort_data == 1
    yy = cell2mat(y); tt = cell2mat(t); 
    levellevel = cell2mat(level); idid = cell2mat(id);
    dat_tmp = sortrows([idid' levellevel' tt' yy'], [1 2 3]);
    yy = dat_tmp(:,4); tt = dat_tmp(:,3); 
    levellevel = dat_tmp(:,2); idid = dat_tmp(:,1);
    
    % Create new data cells 
    ids = unique(idid);
    for sub_ind = 1:length(ids);
        bool_id = (idid == ids(sub_ind));
        y{sub_ind} = transpose(yy(bool_id));
        t{sub_ind} = transpose(tt(bool_id));
        level{sub_ind} = transpose(levellevel(bool_id));
        id{sub_ind} = transpose(idid(bool_id));
    end;
    clearvars dat_tmp;  
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform mean smoothing \mu(t|s)
tt = cell2mat(t);
yy = cell2mat(y);
levelmat = cell2mat(level);
idmat = cell2mat(id);

% Import estimate of \mu(t|s) calculated in MultilevelFuncLong.m
    mu = xmu;
    bw_mu = [];
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find subunit specific deviations from \mu(t|s), \eta_j(t|s)
% Center trajectories by global mean
[~, tobs] = ismember(tt, out1);
yycenter1 = yy - mu(tobs);

% 
uniq_level = unique(levelmat);
yycenter2 = [];
for i = 1:length(uniq_level)
    bool_ind = (levelmat == uniq_level(i));
    
% Import estimate of \eta_j(t|s) calculated in MultilevelFuncLong.m
        eta(:,i) = xeta(:,i);
        bw_eta = [];
  
    
    % Center all responses by \mu(t|s) and \eta_j(t|s) 
    [~, tobs2] = ismember(tt(bool_ind), out1);
    yycenter2(bool_ind) = yycenter1(bool_ind) - eta(tobs2,i)';
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate and smooth covariance surfaces
if strcmp(verbose, 'on') == 1
    fprintf(1,'Part II: Perform covariance smoothing \n');
end

if strcmp(verbose, 'on') == 1
    fprintf(1,'     IIa: Obtain Raw Covariance\n');
end
% Find raw total covariance and first level covariance 
[rcov_same, rcov_diff] = getRawCovKH(yycenter2, tt, idmat, levelmat, 0);

% Smooth raw total covariance, \Sigma_T(t,t'|s) 
if strcmp(verbose, 'on') == 1
    fprintf(1,'     IIb: Total Covariance Smooth \n');
end
[xcov_total bw_xcov_total] = CovarianceSmooth(t, rcov_same, bwxcov_total, ...
    out1, ngrid, kernel, error, verbose, bwxcov_Candidates_total);

% Smooth raw first level covariance, \Sigma^(1)(t,t'|s)
if strcmp(verbose, 'on') == 1
    fprintf(1,'     IIc: first level Covariance Smooth \n');
end
[xcov_between bw_xcov_between] = CovarianceSmooth(t, rcov_diff, bwxcov_between, ...
    out1, ngrid, kernel, 0, verbose, bwxcov_Candidates_between);

% Calculate smooth second level covariance, \Sigma^(2)(t,t'|s)
xcov_within = xcov_total - xcov_between;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform separate FPCA on first level and second level covariance surfaces, \Sigma^(1)(t,t'|s) and \Sigma^(2)(t,t'|s)
if strcmp(verbose, 'on') == 1
    fprintf(1,'Part III: Perform principal components analysis\n');
end


% Find number of principal components to be used for first level FPCA
if isempty(no_opt_B)
    [~, ~, lambda_B] = no_FVE(xcov_between, FVE_threshold(1));
    FVEsum = cumsum(lambda_B)./sum(lambda_B);
    FVE = lambda_B ./ sum(lambda_B);
    no_opt_B = find(FVEsum >= FVE_threshold(1) & ...
        FVE <= FVE_threshold(2) , 1, 'first');
end
     
% Find number of principal components to be used for second level FPCA
if isempty(no_opt_W)
    [~, ~, lambda_W] = no_FVE(xcov_within, FVE_threshold(1));
    FVEsum = cumsum(lambda_W)./sum(lambda_W);
    FVE = lambda_W ./ sum(lambda_W);
    no_opt_W = find(FVEsum >= FVE_threshold(1) & ...
        FVE <= FVE_threshold(2) , 1, 'first');
end

% Estimate measurement error \sigma(s)
if error==1
    [invalid, sigma]=pc_covE(t,out1,bw_xcov_total,ngrid,1,kernel,rcov_same);
else
    sigma=[];
end

% Perform FPCA on \Sigma^(1)(t,t'|s) and \Sigma^(2)(t,t'|s)
[T,lambda_W, phi_W, eigen_W, noeig_W] = evalc('getEigens(xcov_within,out1,out1,no_opt_W,1);');
[T,lambda_B, phi_B, eigen_B, noeig_B] =  evalc('getEigens(xcov_between,out1,out1,no_opt_B,1);');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Order eigenvalues and eigenfunctions 

% Order eigenvalues
[lambda_B,idx_b]=sort(lambda_B,'descend');
[lambda_W,idx_w]=sort(lambda_W,'descend');

% Order eigenfunctions
phi_B=phi_B(:,idx_b);
phi_W=phi_W(:,idx_w);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate first and second level eigenscores, \xi_{ik}(s) and \zeta}_{ijp}(s) 
if calc_scores == 1
    
    if strcmp(verbose, 'on') == 1
        fprintf(1,'Part IV: Compute FPC scores\n');
    end

    [xi zeta] = ComputeScores(yycenter2, tt, idmat, levelmat, ... % Compute eigenscores using multilevel BLUP of Di et al. (2014)
        lambda_B, lambda_W, phi_B, phi_W, sigma);
else
    xi=[]; zeta=[];
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate subject and subunit level predictions
if (calc_scores == 1) & (calc_predictions == 1)
    
    if strcmp(verbose, 'on') == 1
        fprintf(1,'Part V: Compute Predictions \n');
    end

    uniq_id = unique(idmat);
    uniq_level = unique(levelmat);
    for i = 1:length(uniq_id)
        for j = 1:length(uniq_level)
            pred{j} = mu' + eta(:,j) + phi_B * xi(i,:)' + phi_W * zeta{j}(i,:)'; % Calculate prediction using estimated components
        end
        ypred{i} = pred;
    end
else
    ypred = [];
end

% Define output names and assign output

OUTPUTnames = {'mu', 'eta', 'xcov_total', 'xcov_between', 'xcov_within', ...
    'no_opt_B', 'no_opt_W', 'sigma', 'lambda_B', 'lambda_W', 'phi_B', ...
    'phi_W', 'xi', 'zeta', 'ypred', 'bw_mu', 'bw_eta', 'bw_xcov_total', ...
    'bw_xcov_between'};
OUTPUT = {mu, eta, xcov_total, xcov_between, xcov_within, no_opt_B, ...
    no_opt_W, sigma, lambda_B, lambda_W, phi_B, phi_W, xi, zeta, ypred, ...
    bw_mu, bw_eta, bw_xcov_total, bw_xcov_between, OUTPUTnames};



end