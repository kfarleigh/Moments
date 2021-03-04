# Moments
This repository contains scripts for 3D and 4D demographic models using the program Moments, models were published in [Farleigh et al., 2021](). See the Citing this repository subheading for the references in this repository. 
## Purpose
Perform demographic model optimization and comparisons with the python package [moments](https://bitbucket.org/simongravel/moments/src/master/). 
## Overview
In the main repoitory there are two scripts (moments_Run_Optimizations.py & Optimize_Function.py) that must be in the working directory to run properly. These scripts were developed by [Dr. Daniel Portik](https://github.com/dportik/moments_pipeline) and have been modified for our purposes. The optimization script will perform the optimization routine proposed by [Portik et al., 2017](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14266). The routine was originally written for [dadi](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695) but will also work for moments. 
### What you'll need
1. A site frequency spectrum (SFS), this can be created in moments or with a program like [easySFS](https://github.com/isaacovercast/easySFS)
2. The moments_Run_Optimizations.py and Optimize_Function.py scripts, they must be in the same directory as your SFS and Model.py script. 
## Workflow
1. Create your SFS and place in your working directory
2. Download the moments_Run_Optimizations.py and Optimize_Function.py scripts into your working directory with your SFS and Models.py script
3. Run moments
4. Inspect your results 
## Citing this repository
If you use any of the resources from this repository please consider citing the following publications. The models were developed as a part of [Farleigh et al., 2021](), the optimization routines were developed by [Dr. Daniel Portik](https://github.com/dportik/moments_pipeline) and were originally published in [2017](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14266). 

* Farleigh et al., citation
* Gutenkunst, R.N., Hernandez, R.D., Williamson, S.H., and C.D. Bustamante. 2009. Inferring the joint demographic history of multiple populations from multidimensional SNP frequency data. PLoS Genetics 5: e1000695.
* Jouganous, J., Long, W., Ragsdale, A. P., and S. Gravel. 2017. Inferring the joint demographic history of multiple populations: Beyond the diffusion approximation. Genetics 117: 1549-1567.
* Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O., Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017. Evaluating mechanisms of diversification in a Guineo-Congolian forest frog using demographic model selection. Molecular Ecology 26: 5245-5263. https://doi.org/10.1111/mec.14266


