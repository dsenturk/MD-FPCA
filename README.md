CONTENTS OF THIS FOLDER
——————————————	

* MD_FPCA_Tutorial.m : 		Step-by-step implementation of MD-FPCA algorithm using the
					MultilevelFuncLong.m function. MultilevelFuncLong.m assumes data is densely 
					observed on a regular grid in the functional domain and observed either with or 
					without sparsity on a regular grid in the longitudinal domain.

* multilevel_func_data.mat : 	Sample dataset for performing MD_FPCA_Tutorial.m

* Multilevelfunctong.m : 		Fits a MD-FPCA model to repeated multilevel functional data.
					MultilevelFuncLong.m assumes data is densely observed on a regular grid  
					in the functional domain and observed either with or without sparsity on  
					a regular grid in the longitudinal domain. Please refer to the MD-FPCA algorithm 
					in Web Appendix A for annotation of MultilievelFuncLong.m steps.

* MultilevelFPCA.m :		Fits a multilevel FPCA model to multilevel functional data. 

* MeanSmooth2D.m :		Performs local weighted least squares kernel smoothing on a
					two-dimensional with user specified bandwidths or a bandwidth selected
					from a list of candidates by GCV.

* mean2dKH.m :			Performs local weighted least squares kernel smoothing on a
					two-dimensional with user specified bandwidths. 

* getRawCovKH.m :			Finds the total raw covariance and the level one raw covariance for 
					multilevel FPCA implementation.

* gcv_mullwlsn_MFPCA.m :	Performs covariance smoothing and bandwidth candidate selection.

* CovarianceSmooth.m :		Performs covariance smoothing with bandwidth candidate selection or user
					specified bandwidths.

* ComputeScores.m :		Computes scores for multilevel data using BLUP of Di et al. (2014).			
		 
INTRODUCTION
——————————————	

The contents of this folder allow for implementation of the MD-FPCA  decomposition described in  "A multi-dimensional functional principal components analysis of EEG data
" by Hasenstab and Scheffler et al. (2017) (http://onlinelibrary.wiley.com/doi/10.1111/biom.12635/epdf). Users may apply the proposed MD-FPCA decomposition (Multilevelfunctong.m) 
to the sample dataset multilevel_func_data.mat. Detailed instructions on how to perform the aforementioned procedure and visualize results are included in MD_FPCA_Tutorial.m. 
The remaining functions are dependencies that are called by Multilevelfunctong.m.

REQUIREMENTS
——————————————	

The included Matlab programs were developed using Matlab 2015b and require the Matlab package PACE v. 2.17 (http://www.stat.ucdavis.edu/PACE/).

INSTALLATION
——————————————	

Place Matlab .m files and PACE v. 2.17 (http://www.stat.ucdavis.edu/PACE/)
into working directory, along with the multilevel_func_data.mat file. 
These will be used to implement and demonstrate the MD-FPCA algorithm through steps detailed in MD_FPCA_Tutorial.m.
