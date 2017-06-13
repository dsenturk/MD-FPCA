%CovarianceSmooth.m
%============
%Description:
%============
%Performs covariance smoothing with bandwidth candidate selection or user
%specified bandwidths.
%============
% Functions Implemented: 
%============
% getRawCovKH.m and PACE v. 2.17 (http://www.stat.ucdavis.edu/PACE/).
%
% The construction of the MD-FPCA algorithm relies on functions and 
% subroutines from the PACE package (Version 2.17; 2015).
%======
%Input:
%======
%      t:          Vector of subjects' observed times.
%      rcov:       Raw covariance as outputted by the getRawCovKH matlab
%                  function.
%      bwxcov:     1*2 vector, bandwidths for covariance surface used for
%                  smoothing of raw covariance, rov.
%                  bwxcov(i): ith coordinate of bandwidth vector, i=1,2.
%                  [0 0]: implement selection of bandwidth candidates
%                  using GCV.
%                  else: perform two dimensional smoothing using supplied
%                  bandwidths.
%      bwxcov_Candidates: C x 2 vector of bandwidth candidates, where C is number
%                         of bandwidth candidates for consideration. The final
%                         bandwidth is selected by GCV from the C rows of
%                         candidates.          
%      out1:       1 x T vector of observable time points.
%      ngrid1:     integer, number of support points for the covariance 
%                  surface.
%      kernel:     a character string to define the kernel to be used in 
%                  the 1-D or 2-D smoothing
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
%      xcov:        T x T smoothed covariance matrix where T is the length
%                   of vector out1.
%      bw_xcov:     Bandwidth used for final covariance smoothing. If bandwidth 
%                   was specificed  by user, this will equal bwxcov, otherwise  
%                   bw_xcov will equal the bandwidth selected from the C 
%                   rows of candidates by GCV.


function [xcov bw_xcov] = CovarianceSmooth(t, rcov, bwxcov, out1, ...
    ngrid1, kernel, error, verbose, bwxcov_Candidates)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform covariance smooth and bandwidth selection

% Select bandwidth from list of candidates using GCV
if bwxcov(1)==0 || bwxcov(2)==0
        
        [bw_xcov] = gcv_mullwlsn_MFPCA(t,ngrid1,0,error, kernel,rcov,verbose, ...
            bwxcov_Candidates);
        
        if isempty(bw_xcov)
            fprintf(1,'Error: FPCA is aborted because the observed data is too sparse to estimate the covariance function!');
            return;
        end
        
        
        if strcmp(verbose, 'on') == 1
            fprintf(1,['GCV bandwidth choice for COV function : (' num2str(bw_xcov(1)) ',' num2str(bw_xcov(2)) ')\n']);
        end
        
% Select user specified bandwidth        
elseif bwxcov > 0
    bw_xcov = bwxcov;
elseif bwxcov(1) < 0 || bwxcov(2) < 0
    fprintf(1,'Error: Bandwidth choice for the covariance function must be positive!\n');
    return;
end
    
% If measurement error is assumed, smooth over diagonal entries   
rcov1 = rcov;
if error == 1
    tpairn = rcov1.tpairn;
    tneq=find(tpairn(1,:)~=tpairn(2,:));
    cyy = rcov1.cyy;
    rcov1.tpairn = tpairn(:,tneq);
    rcov1.cxxn=cyy(tneq);
    rcov1.win=ones(1,length(rcov1.cxxn));
end

% Perform covariance smooth with user specifed bandwidth or GCV selected
% bandwidth
[invalid,xcov]=mullwlsk(bw_xcov,kernel,rcov1.tpairn,rcov1.cxxn',rcov1.win,out1,out1);  

   
xcov = (xcov+xcov')/2;   %transform the smoothed covariance matrix to guarantee it is a symmetric matrix.


end
