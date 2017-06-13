%gcv_mullwlsn_MFPCA.m
%============
%Description:
%============
%Performs covariance smoothing and bandwidth candidate selection.
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
%      t:          Vector of subjects' observed times.
%      rcov:       Raw covariance as outputted by the getRawCovKH matlab
%                  function.
%      regular:    Default 0.
%      bwxcov_Candidates: C x 2 vector of bandwidth candidates, where C is number
%                         of bandwidth candidates for consideration. The final
%                         bandwidth is selected by GCV from the C rows of
%                         candidates.          
%      out1:       1 x T vector of observable time points.
%      ngrid1:     integer, number of support points for the covariance 
%                  surface.
%      kernel:     a character string to define the kernel to be used in 
%                  the 1-D or 2-D smoothing.
%                  kernel = 'epan'  ==> Epanechnikov kernel
%                          'rect'  ==> Rectangular kernel
%                          'gauss'  ==> Gaussian kernel
%                  Note: The Gaussian kernel is overall best for sparse 
%                        designs but is slower than the other kernels and 
%                        if computational speed is of concern then one may 
%                        wish to use the Epanechnikov kernel also for the 
%                        case of sparse designs.
%      error:      0, no additional measurement error assumed.
%                  1, additional measurement error is assumed.
%      verbose:    a character string.
%                  'on': display diagnostic messages.
%                  'off': suppress diagnostic messages.
%=======
%Output:
%=======
%      bw_xcov:     Bandwidth used for final smoothing, bw_xcov will equal 
%                   the bandwidth selected from the C rows of candidates by
%                   GCV.
%      gcv:         1 x C vector of generalized cross validation statistics  
%                   calculated for each C pairs of bandwidth candidates.


function [bw_xcov,gcv]=gcv_mullwlsn_MFPCA(t,ngrid,regular,error,kernel,rcov, verbose, bwxcov_Candidates)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform covariance smooth and bandwidth selection

tt = cell2mat(t);
out1=unique(tt);
clear tt;
a0=min(out1);
b0=max(out1);

clear t;

rcovcount = rcov.count;
if error == 1
    tpairn = rcov.tpairn;
    tneq=find(tpairn(1,:)~=tpairn(2,:));
    cyy = rcov.cyy;
    tpairn = tpairn(:,tneq);
    cxxn=cyy(tneq);
    win=ones(1,length(cxxn));
    if regular == 1
        rcovcount = rcovcount(tneq);
    end
else
    tpairn = rcov.tpairn;
    cxxn = rcov.cxxn;
    win =  rcov.win;
end


clear rcov cyy tneq;

N = length(cxxn);   %No. of pairs for the covariance
r = range(out1);
clear out1;

% Bandwidth Candidates
if isempty(bwxcov_Candidates)
    fprintf(1,'Error: Must input bandwidth candidates!\n');
else
    bw = bwxcov_Candidates;
end

k0 = mykernel(0, kernel);

out21 = linspace(a0,b0,ngrid); 

gcv = Inf*ones(size(bw,1),1);
% Perform covariance smoothing with bandwidth candidate and calculate cross
% validation statistic
for k = 1:size(bw,1)
    
    if regular == 1
        [invalid, xcov]= mullwlsk(bw(k,:), kernel, tpairn, cxxn', win, out21, out21,rcovcount); 
    else
        [invalid, xcov]= mullwlsk(bw(k,:), kernel, tpairn, cxxn', win, out21, out21);
    end
    %interpolate the smooth covariance from (out21,out21) to (tpairn(1,:), tpairn(2,:))
    if invalid ~= 1
        newxcov = interp2(out21,out21,xcov, tpairn(1,:),tpairn(2,:),'spline');
        clear xcov;
        if regular == 1
            cvsum = (cxxn./(rcovcount')-newxcov)*(cxxn./(rcovcount')-newxcov)';
        else
            cvsum = (cxxn-newxcov)*(cxxn-newxcov)';
        end
        clear newxcov;
        bottom = 1-(1/N)*((r*k0)/bw(k))^2;
        gcv(k) = cvsum/(bottom)^2;
        tmp = gcv(~isinf(gcv));
        if length(tmp) > 1 && gcv(k) > gcv(k-1)
            leave = 1;
            break;
        end
    end
end


if all(gcv == Inf)
    fprintf(1,'Error: All GCV values are Inf, setting to minimum bandwidth candidate!\n');
    bw_xcov = min(bw);
else
    bw_xcov = bw(find(gcv == min(gcv),1,'first'),:); % Select bandwidth that is associated with minimum GCV statistic
end


if strcmp(kernel,'gauss') == 0 && strcmp(verbose, 'on') == 1
    fprintf(1,['GCV bandwidth choice for COV function: (' num2str(bw_xcov(1)) ',' num2str(bw_xcov(2)) ')\n']);
end
end
