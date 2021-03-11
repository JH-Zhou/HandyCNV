# HandyCNV
![R-CMD-check](https://github.com/JH-Zhou/HandyCNV/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JH-Zhou/HandyCNV/actions/workflows/R-CMD-check.yaml)
An R package for Standardized Summary, Annotation, Comparison, and Visualization of CNV, CNVR and ROH

![Fig.1 Working Flow](https://github.com/JH-Zhou/HandyCNV/blob/master/vignettes/High_resolution_pipeline_HandyCNV.png)

# Introduction
This package was originally designed for the Post-analysis of CNV results inferred from PennCNV and CNVPartition (GenomeStudio). But now it has now been expanded to accept input files in standard formats for a wider range of applications. Our motivation is to provide a standard, reproducible and time-saving pipeline for the post-analysis of CNVs and ROHs detected from SNP genotyping data for the majority of diploid Species. The functions provided in this package can be categorised into five sectors: Conversion, Summary, Annotation, Comparison and Visualization. The most useful features provided are: integrating summarized results, generating lists of CNVR, annotating the results with known gene positions, plotting CNVR distribution maps, and producing customised visualisations of CNVs and ROHs with gene and other related information on one plot. This package also supports a range of customisations, including the colour, size of high resolution figures, and output folder, avoiding conflict between the results of different runs. Running through all functions detailed in the vignette could help us to identify and explore the most interesting genomic regions more easily. In the following sections, we will present how to use the functions provided in this package to solve these problems.

# Prerequisite
First, to run this package, we need to make sure that R (Version >= 3.5.2) is installed in your computer (R download link: https://www.r-project.org/).<br> Once R is installed, the 'HandyCNV' package can be installed from Github repository by running the following script:
```{r}
#library(remotes)
#install_github(repo = "JH-Zhou/HandyCNV", auth_token = "3d2ac98e4c297bab332f1e68b3b2d49f3a17d6aa")
```
Then, we need to load the 'HandyCNV' package in order to run the following examples. This can be done using the `library` function as shown below.
```{r setup}
library(HandyCNV)
```
