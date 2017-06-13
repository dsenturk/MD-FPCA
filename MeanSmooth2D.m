%MeanSmooth2D.m
%============
%Description:
%============
%Performs local weighted least squares kernel smoothing on a
%two-dimensional with user specified bandwidths or a bandwidth selected
%from a list of candidates by GCV.
%============
% Functions Implemented: 
%============
% mean2DKH.m and PACE v. 2.17 (http://www.stat.ucdavis.edu/PACE/).
%
% The construction of the MD-FPCA algorithm relies on functions and 
% subroutines from the PACE package (Version 2.17; 2015).
%======
%Input:
%======
%      y:           1 x N cell containing data on each subject or subunit.
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
%      out1:        Time grid for longitudinal time, s.
%      out2:        Time grid for functional time, t.
%      bwmu:        1 x 2 vector used as bandwidths in the mean smooth for 
%                   the response surface. The first entry corresponds to longitudinal
%                   time s and second entry to functional time t. 
%                   [0 0]: implement selection of bandwidth candidates
%                   using GCV.
%                   else: perform two dimensional smoothing using supplied
%                   bandwidths.
%      bwmu_Candidates1: C x 2 vector of bandwidth candidates, where C is number
%                       of bandwidth candidates for consideration. The first entry 
%                       for each row corresponds to longitudinal time s 
%                       and second entry to functional time t. The final
%                       bandwidth is selected by GCV from the C rows of
%                       candidates.
%=======
%Output:
%=======
%      mu:          Two dimensional smoothed response surface.
%      bwmu:        Bandwidth used for final smoothing. If bandwidth was specificed 
%                   by user, this will equal bwmu, otherwise bwmu will equal 
%                   the bandwidth selected from the C rows of candidates by
%                   GCV.
%                   


function [mu, bw_mu] = MeanSmooth2D(y, t, out1, out2, bwmu, bwmu_Candidates1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pool data into vectors/matrices
tt = cell2mat(t);
yy = cell2mat(y);

% Obtain smoothed mean surface
rawmean = getRaw2dMean(tt, yy);
outmu = GetMesh(out1, out2);
win1 = ones(1,size(rawmean.tt,2));

% Perform two dimensional smoothing with user specificed bandwidth
if bwmu(1) > 0 & bwmu(2) > 0
    
    % Smooth data
    bw_mu = bwmu;
    [invalid, mu] = mean2dKH(bw_mu,'epan',rawmean.tt,rawmean.yy',win1,outmu,rawmean.count);
  
% Select bandwidth from list of candidates using GCV   
elseif bwmu(1) == 0 || bwmu(2) == 0
   
    bw = bwmu_Candidates1; 
    N = length(rawmean.yy);   %No. of pairs
    r1 = range(out1);
    r2 = range(out2);
    k0 = mykernel(0, 'gauss');
    
    for bw_ind = 1:size(bw, 1)
        
        bwmu = bw(bw_ind, :);
        [invalid, mutmp] = mean2dKH(bwmu,'epan',rawmean.tt,rawmean.yy',win1,outmu,rawmean.count);
        
        
        musmooth = mat2dataset([mutmp', rawmean.tt'], 'VarNames',{'musmooth', 't1', 't2'});
        datobs = mat2dataset([yy', tt'], 'VarNames',{'datobs', 't1', 't2'});
        dat_merge=join(musmooth, datobs,'type','rightouter','mergekeys',true);
                
        cvsum = (dat_merge.datobs - dat_merge.musmooth)' * (dat_merge.datobs - dat_merge.musmooth);
        gcv(bw_ind) = cvsum / ( (1-(r1*k0)/(N*bwmu(1))) * (1-(r2*k0)/(N*bwmu(2))) ); % Calculate GCV statistic
        
        if all(gcv == Inf)
            fprintf(1,'Error: All GCV values are Inf, setting to minimum bandwidth candidate!\n');
            bw_mu = min(bw);
        else
            bw_mu = bw(find(gcv == min(gcv),1,'first'),:); % Select bw_mu to be minimized of GCV statistic
        end

    end  % end bw_ind loop
    
    [invalid, mu] = mean2dKH(bw_mu,'epan',rawmean.tt,rawmean.yy',win1,outmu,rawmean.count); % Perform final smooth with bandwidth selected by GCV
        
end


% Convert to matrix form 
mu = reshape(mu, length(out2), length(out1));
clear rawmean win1 tt yy;

end