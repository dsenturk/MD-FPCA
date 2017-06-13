%ComputeScores.m
%============
%Description:
%============
%Computes scores for multilevel data using BLUP of Di et al. (2014).
%======
%Input:
%======
%      ycentered:  vector of responses centered by \mu(t|s) and \eta_j(t|s).
%      tt:         vector of time points corresponding to vector ycentered.
%      idmat:      vector of id numbers corresponding to vector ycentered
%      levelmat:   vector of j-level repetitions indices corresponding to
%                  vector ycentered.
%      lambda_B:   1 x K vector of eigenvalues corresponding to the
%                  first level covariance truncated eigen decomposition,
%                  \lambda^(1)_k(s).
%      lambda_W:   1 x P vector of eigenvalues corresponding to the
%                  second level covariance truncated eigen decomposition,
%                  \lambda^(2)_p(s).
%      phi_B:      T x K matrix of eigenfunctions for the first level
%                  covariance, \phi^(1)_k(t|s). Each column is a separate 
%                  eigenfunction .
%      phi_W:      T x P matrix of eigenfunctions for the second level
%                  covariance, \phi^(2)_p(t|s). Each column is a separate 
%                  eigenfunction.
%      sigma:      Estimate measurement error variance, \sigma(s). If data 
%                  is not assumed to contain error, set sigma = [].
%=======
%Output:
%=======
%      xi:         N x K matrix containing first level eigenscores 
%                  \xi_{ik}(s), where each row corresponds to a subject 
%                  and each column corresponds to a first level 
%                  eigenvalue-eigenfunction pair.
%      zeta:       1 x J cell of second level eigenscores \zeta_{ijp}(s).
%                  Each cell contains a N x P matrix of the jth second 
%                  level eigenscores, where each row corresponds to a 
%                  subject and each column corresponds to a second level 
%                  eigenvalue-eigenfunction pair.
%
function [xi, zeta] = ComputeScores(ycentered, tt, idmat, levelmat, ...
    lambda_B, lambda_W, phi_B, phi_W, sigma)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store mean and sub-unit centered data
yycenter2 = ycentered;
out1 = unique(tt);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate covariance matrices for BLUP estimation

% Reconsruct first and and second level covariance surfaces eigencomponents
est_B = phi_B * diag(lambda_B) * phi_B'; % First level covariance surface
est_W = phi_W * diag(lambda_W) * phi_W'; % Second level covariance surface


% Measurement error variance
if ~isempty(sigma)
    est_sigma = diag(sigma*ones(1, length( est_B(:,1) )));
else
    est_sigma = diag(0*ones(1, length( est_B(:,1) )));
end

% Store covariance surfaces
within_level_tmp = est_W + est_sigma; 
between_level = est_B;

% Calculate SIGMA covariance matrix for estimation of eigenscores using BLUP
J = length(unique(levelmat)); % Extract number of subunits
SIGMA = repmat(between_level, J); 
SIGMA = kron(eye(J), within_level_tmp) + SIGMA;

% Calculate A and B matrices
LAMBDA_W = diag(lambda_W);  
LAMBDA_B = diag(lambda_B);
A = repmat(LAMBDA_B * phi_B', 1, J); % cov(\xi_{ik}(s), Y_i)
B = kron(eye(J), LAMBDA_W * phi_W'); % cov(\zeta_{ijp}(s), Y_i)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate subject eigenscores

% Matrix of complete times and levels
t_comp = repmat(out1, 1, J);
tmp = repmat(unique(levelmat), length(out1), 1);
t_comp = [t_comp' tmp(:)]';
if length(tt) == 1
    tmp_comp = mat2dataset([1, t_comp'], ...
        'VarNames',{'index', 't', 'j'});
else
    tmp_comp = mat2dataset([(1:length(t_comp))', t_comp'], ...
        'VarNames',{'index', 't', 'j'});
end

uniq_id = unique(idmat);

% Calculate subject subject-subunit specific eigenscores

for i = 1:length(uniq_id) 
    % Subset data for ith subject
    bool_ind = (idmat == uniq_id(i));
    ysub_center = yycenter2(bool_ind);
    ttsub = tt(bool_ind);
    levelsub = levelmat(bool_ind);
    
    % Observed times for ith subject
    t_obs = [ttsub' levelsub']';
    tmp_obs = mat2dataset(t_obs', 'VarNames',{'t', 'j'});
    tmp = join(tmp_obs, tmp_comp);
    obs_index = tmp.index;
    
    % Subject terms for observed data
    pinvSIGMA = inv(SIGMA(obs_index, obs_index));
     
    % Estimated first-stage eigenscores using BLUP
    xi(i,:) = A(:,obs_index) * pinvSIGMA * ysub_center';
    zeta_tmp(i,:) = B(:,obs_index) * pinvSIGMA * ysub_center';
end


% Reorganize zeta scores into 1 x J cell each containing N x P matrix 
no_opt_W = length(lambda_W);
extract_index = 1:no_opt_W;
for j = 1:J
    zeta{j} = zeta_tmp(:, extract_index);
    extract_index = extract_index + no_opt_W;
end


end


