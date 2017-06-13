%Mean2dKH.m
%============
%Description:
%============
%Performs local weighted least squares kernel smoothing on a
%two-dimensional with user specified bandwidths. 
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
%      bwmu:        1 x 2 vector used as bandwidths in the mean smooth for 
%                   the response surface. The first entry corresponds to longitudinal
%                   time s and second entry to functional time t. 
%      kernel:      Kernel used for smoothing. Default is 'epan'
%                   (Epanechnikov).                            
%      xin:         Matrix of predictors.
%      yin:         Vector of responses.
%      win:         Vector of case weights.
%      out:         Rectangular 2-D grid, domain for the response surface
%      count:       Number of observations at a point in the rectangular 2-D
%                   grid defind by out.

%=======
%Output:
%=======
%      mu:          Two dimensional smoothed response surface.
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [invalid,mu]=mean2dKH(bw,kernel,xin,yin,win, out, count)


active=find(win~=0);
xin=xin(:,active);
yin=yin(active);
win=win(active);
invalid=0;
%gap=[];

    mat = num2cell(out, 1);
    [invalid, mu] = cellfun(@(x) GetEstimate(xin, yin, win, bw, kernel, x, count), mat);
    mu(invalid == 1) = [];
    return;
end






