# Lewis and Vazquez-Grande (2018) r* Estimation

This code uses Bayesian methods to estimate the natural rate of interest (r*) using data through 2018-Q3.  This specification was introduced in [Lewis and Vazquez-Grande (2018)](https://onlinelibrary.wiley.com/doi/abs/10.1002/jae.2671), which is included in the `paper` folder for reference, along with an appendix which outlines the state space structure of the model and some other estimation details.

## Executing the code
Cloning the repository and running the MATLAB script `estimateModel.m` will produce estimates of the unobserved series for the natural rate of interest ( *r*\* ), the longer-run natural rate of interest ( *r*\*<sub>*LR*</sub> ), potential output ( *y*<sup>*P*</sup> ), the growth rate of potential output ( *g* ) and the non-growth component of the natural rate ( *z* ).  Relevant percentiles from the posterior estimates of these objects will be in the structure `outVals` when the estimation completes, and two files will be generated in the `data\output` folder containing the filtered and smoothed estimates from `outVals`.  The folder contains two files which should be replicable by running the code as-is, though you will need to specify a directory for the larger .mat files (see below).  As is stated in the comments of the code itself, users should remove the first 6 lines of `estimateModel.m` if they do not want to use the random number generator seed used for replication purposes.

## Updating the data
The data are provided through 2018-Q3 and will not be updated further.  However, the code draws its underlying data from the file `lvg.data.mat`.  Users can update the data in the spreadsheet `data\input\lvg.data.xls` as desired and run `dataXLS2MAT.m` to create a new version of `lvg.data.mat` which they can use to generate new model output.

## .mat file output
The `estimateModel.m` script will generate two .mat files which will both be sizeable.  The location of these files needs to be specified using the code in lines 14 to 22 of `estimateModel.m` or else these output files will not be saved.  The first will be called `AltModel_YYYYMMDD_HHMMSS.mat` (if the `baseAlt` variable is left set to 1) where the date and time in the filename are specified by the time the script is executed.  This file will have a structure in it called `res`, which will contain the data and the posterior draws of the model's parameters.  The second file will have the same name, with `_outVals` added to the end.  That file will contain the data and parameter draws, but also all of the draws of the unobserved variables implied by the filter, as well as time series paths of those unobserved objects based on percentiles specified in the code.

## .csv file output
After the code is run, two files generated in the `data\output` folder will contain percentiles of the unobserved objects of interest: the natural rate ( *r*\* ), the longer-run version ( *r*\*<sub>*LR*</sub> ), potential output ( *y*<sup>*P*</sup> ), the growth rate of potential output ( *g* ) and the non-growth component of the natural rate ( *z* ).  One file contains the one-sided "filtered" estimates, and the other contains the two-sided "smoothed" estimates.  The files `lvgEstimatesFiltered.replicate.csv` and `lvgEstimatesSmoothed.replicate.csv` contain the values that should result from running the code as originally released.

## Description of `outVals` output
The structure `outVals` has the following components, several of which are duplicates from the input data.  Note that below there are *T* periods, *M* draws from the posterior, *K* parameters, *N*<sub>*y*</sub> observable series, *N*<sub>*u*</sub> exogenous series, *N*<sub>*S*</sub> elements in the state vector, and *PP* is the number of percentiles shown in the posterior paths which aggregate up in the information from the posterior draws (these are originally the 2.5, 5, 16, 50, 84, 95, and 97.5 percentiles). 

|Field Name | Description | 
|:--- |:---- |
|`.y`| Observed Data in the model (*N*<sub>*y*</sub> by *T*) |
|`.u`| Exogenous Data in the model (*N*<sub>*u*</sub> by *T*) |
|`.dates`| The corresponding dates for the observed and exogenous data (*T* by 1) |
|`.vintage`| The date when the data was pulled |
|`.varNames`| The names of the parameters in text |
|`.varNamesLatex`| The names of the parameters in LaTeX |
|`.theta`| The draws of the parameters from the joint posterior (*M* by *K*) |
|`.LL`| The log-likelihood value based on each draw from the joint posterior (1 by *M*) |
|`.LLt`| The period log-likelihood based on each draw from the joint posterior (T by *M*) |
|`.BSDraw`| The **backward-sampled** (*smoothed*) draws of the state variable vectors (*N*<sub>*S*</sub> by *T* by *M*) | 
|`.FFDraw`| The **forward-filtered** (*filtered*) draws of the state variable vectors (*N*<sub>*S*</sub> by *T* by *M*) | 
|`.yPBSDraws`| The **backward-sampled** (*smoothed*) draws of potential output (*T* by *M*) |
|`.yPFFDraws`| The **forward-filtered** (*filtered*) draws of potential output (*T* by *M*) |
|`.zBSDraws`| The **backward-sampled** (*smoothed*) draws of the non-growth component of *r*\* (*T* by *M*) |
|`.zFFDraws`| The **forward-filtered** (*filtered*) draws of the non-growth component of *r*\* (*T* by *M*) |
|`.gBSDraws`| The **backward-sampled** (*smoothed*) draws of potential output growth (*T* by *M*) |
|`.gFFDraws`| The **forward-filtered** (*filtered*) draws of potential output growth (*T* by *M*) |
|`.rSBSDraws`| The **backward-sampled** (*smoothed*) draws of *r*\* (*T* by *M*) |
|`.rSFFDraws`| The **forward-filtered** (*filtered*) draws of *r*\* (*T* by *M*) |
|`.rSLRBSDraws`| The **backward-sampled** (*smoothed*) draws of longer-run *r*\* (*T* by *M*) |
|`.rSLRFFDraws`| The **forward-filtered** (*filtered*) draws of longer-run *r*\* (*T* by *M*) |
|`.rSBSpath`| The **backward-sampled** (*smoothed*) path of *r*\* (*T* by *PP*) |
|`.rSFFpath`| The **forward-filtered** (*filtered*) path of *r*\* (*T* by *PP*) |
|`.rSLRBSpath`| The **backward-sampled** (*smoothed*) path of longer-run *r*\* (*T* by *PP*) |
|`.rSLRFFpath`| The **forward-filtered** (*filtered*) path of longer-run *r*\* (*T* by *PP*) |
|`.yPBSpath`| The **backward-sampled** (*smoothed*) path of potential output (*T* by *PP*) |
|`.yPFFpath`| The **forward-filtered** (*filtered*) path of potential output (*T* by *PP*) |
|`.zBSpath`| The **backward-sampled** (*smoothed*) path of the non-growth component of *r*\* (*T* by *PP*) |
|`.zFFpath`| The **forward-filtered** (*filtered*) path of the non-growth component of *r*\* (*T* by *PP*) |
|`.gBSpath`| The **backward-sampled** (*smoothed*) path of the growth rate of potential output (*T* by *PP*) |
|`.gFFpath`| The **forward-filtered** (*filtered*) path of the growth rate of potential output (*T* by *PP*) |
|`.BF`| The Bayes Factor in favor of the model with transitory shocks, as discussed in the LVG (2018) paper |


