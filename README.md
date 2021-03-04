# Moments
This repository contains scripts for 3D and 4D demographic models using the program Moments, models were published in [Farleigh et al., 2021](). See the Citing this repository subheading for the references in this repository. 

## Purpose
Perform demographic model optimization and comparisons with the python package [moments](https://bitbucket.org/simongravel/moments/src/master/). 
## Overview
In the main repository there are two scripts (moments_Run_Optimizations.py & Optimize_Function.py) that must be in the working directory to run properly. These scripts were developed by [Dr. Daniel Portik](https://github.com/dportik/moments_pipeline) and have been modified for our purposes. The optimization script will perform the optimization routine proposed by [Portik et al., 2017](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14266). The routine was originally written for [dadi](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695) but will also work for moments. 

The 3D and 4D subdirectories contain the python scripts for the respective sets of models, there are also conceptual figures that demonstrate what each model does in these subdirectories.
### What you'll need
1. A site frequency spectrum (SFS), this can be created in moments or with a program like [easySFS](https://github.com/isaacovercast/easySFS)
2. The moments_Run_Optimizations.py and Optimize_Function.py scripts, they must be in the same directory as your SFS and Model.py script. 
3. Summarize_Outputs.py, must also be in your working directory
## Workflow
1. Create your SFS and place in your working directory
2. Download the moments_Run_Optimizations.py and Optimize_Function.py scripts into your working directory with your SFS and Models.py script
3. Run moments
4. Inspect your results 
## Example Usage
Load python
```
module load anaconda-python3
source /software/python/anaconda3/etc/profile.d/conda.sh
```
Import the modules 
```
python3
import moments
import matplotlib
import pylab
```
Run the Models (assuming that your models are named Models_3D.py), make sure you understand the optimization settings before copying mine. See Dr. Daniel Portik's [repository](https://github.com/dportik/moments_pipeline) for an explaination. 
```
import Models_3D
import Optimize_Functions

#Set the number of rounds here
rounds = 4
#define the lists for optional arguments
#you can change these to alter the settings of the optimization routine
reps = [10,20,30,40]
maxiters = [3,5,10,15]
folds = [3,2,2,1]
prefix = "Final"
fs_folded = True

# Run the Models 
Optimize_Functions.Optimize_Routine(fs, prefix, "sim_split_no_mig", Models_3D.sim_split_no_mig, rounds, 4, fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, param_labels = "nu1, nu2, nu3, T1")

# Summarize the outputn (after leaving python)
python ./Summarize_Outputs.py ./
```
## Contact
If you have any questions or issues with this repository please post it on the issues page of this repository or email me at keakafarleigh@gmail.com. Also, if you would like to run moments but do not see the models that you want to run here feel free to reach out. 

## Citing this repository
If you use any of the resources from this repository please consider citing the following publications. The models were developed as a part of [Farleigh et al., 2021](), the optimization routine was developed by [Dr. Daniel Portik](https://github.com/dportik/moments_pipeline) and was originally published in [2017](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14266). 

* Farleigh et al., citation
* Gutenkunst, R.N., Hernandez, R.D., Williamson, S.H., and C.D. Bustamante. 2009. Inferring the joint demographic history of multiple populations from multidimensional SNP frequency data. PLoS Genetics 5: e1000695.
* Jouganous, J., Long, W., Ragsdale, A. P., and S. Gravel. 2017. Inferring the joint demographic history of multiple populations: Beyond the diffusion approximation. Genetics 117: 1549-1567.
* Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O., Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017. Evaluating mechanisms of diversification in a Guineo-Congolian forest frog using demographic model selection. Molecular Ecology 26: 5245-5263. https://doi.org/10.1111/mec.14266


