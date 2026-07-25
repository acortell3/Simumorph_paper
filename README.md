# REPRODUCIBILITY FOR

Supplementary material and code for "Simulating shape variation in material culture: The Simumorph R-package", by Alfredo Cortell-Nicolau and Anne Kandler, currently under review.

## Contents and structure
The present project contains the following six folders:

* Figures: Figures used in the manuscript.
* Parametric\_space\_outlines: Contains the outlines of the shapes used for GMM.
* R: Scripts necessary for the full reproducibility of the paper (see 'usage' section).
* SI3\_Results: Contains the .rds objects used in the SI.
* SI\_figures: Contains two folders (SI2 and SI3) with all the figures of the supplementary material.
* Utilities: Contains the necessary .rds objects for the main script (e.g. shapes, covariance matrices, etc.) 

## Usage
To fully reproduce the article, the different scripts must be run in the order shown in the folder R, and below. Paths are local, so they should work out of the box, but please have a look if you need to modify them in your setup. Here's a brief summary of what each script does. Please, bear in mind that parallel computing has been used for this paper. 35 cores for the main script (reproducible script) and 20 for the sensitivity analysis, and has taken several days in both cases. Do not attempt to replicate with few cores. If you want to first try it, reduce the number of simulations where convenient.

1. 01\_Reproducible\_script.R: It reproduces in full the contents show in the paper.
2. 02\_SI3\_script.R: It reproduces the analysis for SI3.
3. 03\_Sensitivity\_analysis.R: Only the sensitivity analysis of the parameters alpha, delta, epsilon and s, as required.
4. 04\_SI3\_plots.R: If the scripts above have been run, all objects will have been generated and this produces the plots shown in SI3. The rds objects, as used in the paper, can also be found at the folder SI3\_results, so this script should work as a stand-alone.
5. 05\_Fig\_1\_mosquito.R: This produces the figure 1 of the paper. Since it doesn't take part of the analysis per se, it has been saved separately.

## The Simumorph R package
Additional information can be found at the dedicated Simumorph github [here](https://github.com/acortell3/Simumorph). 

## Additional safety
This folder has also been saved as a zenodo file here, to preserve the exact content of the paper.

## Wrap up
And I guess that's all you need to know, but please do reach out if you have any doubt on how to implement this!




