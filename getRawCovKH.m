%getRawCovKH.m
%============
%Description:
%============
%Finds the total raw covariance and the level one raw covariance for 
%multilevel FPCA implementation.
%============
% Functions Implemented: 
%============
% PACE v. 2.17 (http://www.stat.ucdavis.edu/PACE/).
%
% The construction of the MD-FPCA algorithm relies on functions and 
% subroutines from the PACE package (Version 2.17; 2015).
%======
%Input:
%======
%      yycentered:  vector of responses centered by \mu(t|s) and \eta_j(t|s)      
%      tt:          vector of time points corresponding to vector ycentered. 
%      idmat:       vector of id numbers corresponding to vector ycentered.
%      levelmat:    vector of j-level repetitions indices corresponding to
%                   vector ycentered.
%      error:      0, no additional measurement error assumed.
%                  1, additional measurement error is assumed.
%=======
%Output:
%=======
%      same_struct: a structure array that contains tpairn, cxxn, indx, win 
%                   and cyy for the total raw covariance \Sigma_T(t,t'|s).
%      diff_struct: a structure array object that contains tpairn, cxxn, indx, win 
%                   and cyy for the level one raw covariance
%                   \Sigma^(1)(t,t'|s).
%
%         Elements of the above structs are described below:
%         tpairn:      2 * N  vector denotes the  pairs of time points for 
%                      subject concatenating as two vectors 
%                      if error = 1, all (t_ij, t_ij) will be removed.  
%         cxxn:        1 * N vector of raw covariance corresponding to tpairn      
%         indx:        1 * N vector of indices for each subject.
%         win:         1 * N weight matrix for the 2-D smoother for covariance
%                       function.
%         cyy:         1 * M vector of raw covariance corresponding to all 
%                       pairs of time points, i.e., it is the same as cxxn if 
%                       error = 0.


function [same_struct diff_struct] = getRawCovKH(yycentered, tt, idmat, ...
    levelmat, error)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define matrices for raw covariance estimation

t1_all_same = []; t2_all_same = [];
t1_all_diff = []; t2_all_diff = [];
y1_all_same = []; y2_all_same = [];
y1_all_diff = []; y2_all_diff = [];
indx_same = []; indx_diff = [];
count = [];

uniq_id = unique(idmat);
ncohort = length(uniq_id);
out1  = unique(tt);

for i = 1:ncohort
    
    % Subset subject i
    bool_ind = (idmat == uniq_id(i));
    yysub = yycentered(bool_ind);
    ttsub = tt(bool_ind);
    levelsub = levelmat(bool_ind);
    
    % Create all combinations of within subject observations
    tvec = 1:length(ttsub);
    [x1, x2] = meshgrid(tvec, tvec);
    x1 = x1(:); x2 = x2(:);
    
    yy1 = yysub(x1); yy2 = yysub(x2);
    tt1 = ttsub(x1); tt2 = ttsub(x2);
    levelsub1 = levelsub(x1); levelsub2 = levelsub(x2);
    
    % Index for same subunit or different subunit observations
    same_level = (levelsub1 == levelsub2);
    diff_level = (levelsub1 > levelsub2);
    tt1same = tt1(same_level);    tt2same = tt2(same_level);
    tt1diff = tt1(diff_level);    tt2diff = tt2(diff_level);
    yy1same = yy1(same_level);    yy2same = yy2(same_level);
    yy1diff = yy1(diff_level);    yy2diff = yy2(diff_level);
    
    % Store for each subject
    t1_all_same((length(t1_all_same) + 1):(length(t1_all_same) + ...
        length(tt1same))) = tt1same;
    t2_all_same((length(t2_all_same) + 1):(length(t2_all_same) + ...
        length(tt2same))) = tt2same;
    t1_all_diff((length(t1_all_diff) + 1):(length(t1_all_diff) + ...
        length(tt1diff))) = tt1diff;
    t2_all_diff((length(t2_all_diff) + 1):(length(t2_all_diff) + ...
        length(tt2diff))) = tt2diff;
    
    y1_all_same((length(y1_all_same) + 1):(length(y1_all_same) + ...
        length(yy1same))) = yy1same;
    y2_all_same((length(y2_all_same) + 1):(length(y2_all_same) + ...
        length(yy2same))) = yy2same;
    y1_all_diff((length(y1_all_diff) + 1):(length(y1_all_diff) + ...
        length(yy1diff))) = yy1diff;
    y2_all_diff((length(y2_all_diff) + 1):(length(y2_all_diff) + ...
        length(yy2diff))) = yy2diff;
        

end

% Store output
tpairn_same = [t1_all_same' t2_all_same']';
tpairn_diff = [t1_all_diff' t2_all_diff']';
cyy_same = y1_all_same .* y2_all_same;
cyy_diff = y1_all_diff .* y2_all_diff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average across subunits and subjects to speed up algorithm
% For regular dense or regular with missing data only

% Create unique time pairs and average for each pair
uniq_same = unique(tpairn_same', 'rows')';
tmp_uniq = mat2dataset([(1:length(uniq_same))', uniq_same'], ...
    'VarNames',{'index', 't1', 't2'});
tmp_all = mat2dataset(tpairn_same', 'VarNames',{'t1', 't2'});
tmp = join(tmp_all, tmp_uniq);

% Store averaged data
cyy_same = accumarray(tmp.index, cyy_same',[], @(x)mean(x))';
tpairn_same = uniq_same;

% For diff
uniq_diff = unique(tpairn_diff', 'rows')';
tmp_uniq = mat2dataset([(1:length(uniq_diff))', uniq_diff'], ...
    'VarNames',{'index', 't1', 't2'});
tmp_all = mat2dataset(tpairn_diff', 'VarNames',{'t1', 't2'});
tmp = join(tmp_all, tmp_uniq);
 
% Store averaged data
cyy_diff = accumarray(tmp.index, cyy_diff',[], @(x)mean(x))';
tpairn_diff = uniq_diff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if error == 1
    tneq_same = find(tpairn_same(1,:) ~= tpairn_same(2,:));
%     indx_same = indx_same(tneq_same);
    tpairn_same = tpairn_same(:,tneq_same);
    cxxn_same = cyy_same(tneq_same);     
elseif error == 0
    cxxn_same = cyy_same;     
end

cxxn_diff = cyy_diff;  

win_same = ones(1, length(cxxn_same));
win_diff = ones(1, length(cxxn_diff));

same_struct = struct('tpairn',tpairn_same, 'cxxn',cxxn_same, ...
    'indx',indx_same, 'win',win_same,'cyy',cyy_same,'count',count);
diff_struct = struct('tpairn',tpairn_diff, 'cxxn',cxxn_diff, ...
    'indx',indx_diff, 'win',win_diff,'cyy',cyy_diff,'count',count);
     

end
