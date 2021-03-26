# HandyCNV
[![R-CMD-check](https://github.com/JH-Zhou/HandyCNV/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JH-Zhou/HandyCNV/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/JH-Zhou/HandyCNV/branch/master/graph/badge.svg?token=5ERXZH0GFW)](https://codecov.io/gh/JH-Zhou/HandyCNV)

An R package for Standardized Summary, Annotation, Comparison, and Visualization of CNV, CNVR and ROH

![Fig.1 Working Flow](https://github.com/JH-Zhou/HandyCNV/blob/master/vignettes/High_resolution_pipeline_HandyCNV.png)

# Introduction
This package was originally designed for the Post-analysis of CNV results inferred from PennCNV and CNVPartition (GenomeStudio). However, it has now been expanded to accept input files in standard formats for a wider range of applications. Our motivation is to provide a standard, reproducible and time-saving pipeline for the post-analysis of CNVs and ROHs detected from SNP genotyping data for the majority of diploid Species. The functions provided in this package can be categorised into five sections: Conversion, Summary, Annotation, Comparison and Visualization. The most useful features provided are: integrating summarized results, generating lists of CNVR, annotating the results with known gene positions, plotting CNVR distribution maps, and producing customised visualisations of CNVs and ROHs with gene and other related information on one plot. This package also supports a range of customisations, including the colour, size of high resolution figures, and output folder, avoiding conflict between the results of different runs. Running through all functions detailed in the vignette could help us to identify and explore the most interesting genomic regions more easily. In the following sections, we will present how to use the functions provided in this package to solve these problems.

## Installation and Prerequisites
First, to run this package, we need to make sure that R (Version >= 3.5.2) is installed in your computer (R download link: https://www.r-project.org/). Once R is installed, the 'HandyCNV' package can be installed from Github repository by running the following script:
```{r}
#library(remotes)
#install_github(repo = "JH-Zhou/HandyCNV", auth_token = "3d2ac98e4c297bab332f1e68b3b2d49f3a17d6aa")
```
Then, we need to load the 'HandyCNV' package in order to run the following examples. This can be done using the `library` function as shown below.
```{r setup}
library(HandyCNV)
```

To start playing with this package, we first need to prepare at least one CNV result list. With only CNV results as an input file, we can explore the functions from section 1 to section 10 as below. But to get more interesting and potentially valuable results, we will have to prepare additional input files, including the map file and signal intensity file for the SNP chip used to generate the CNV list, the pedigree of the samples, and Plink format (bed, bim and fam) genotype files. Many of these files are already required during CNV discovery, and by generating them eaarly in the pipeline, you will find that the rest of workflow will be much easier and faster. We will be introducing how to prepare each of the input files in the relevant chunks below.

## Demo data
We have provided some internal demo data, which should installed with package, in order to demonstrate how to use this package. We can access the demo data using the built-in R function `system.file` (see examples below). In the demo data, the first CNV results file is the default output from PennCNV with the ARS-UCD 1.2 (ARS) cattle reference genome. The second CNV results file contains the default output from CNVPartition with UMD 3.1 (UMD) cattle reference genome, and the third CNV list is an example of a standard input file that could be prepared from CNV results inferred by other tools.<br>

Note: If you prefer not to test with this demo data, just skip this step to section 1 and use your own data instead.
```{r}
cnv_penn_ars <- system.file("extdata", "Demo_data/cnv_results/ARS_PennCNV_WGC.goodcnv", package = "HandyCNV")
cnv_part_umd <- system.file("extdata", "Demo_data/cnv_results/UMD_CNVPartition.txt", package = "HandyCNV")
cnv_other_umd <- system.file("extdata", "Demo_data/cnv_results/UMD_Standard_CNV.txt", package = "HandyCNV")
```

## Setting up a working directory
To start the analysis, let's first set up the working directory. This will help ensure that all results will be saved in the same directory.
```{r eval=FALSE}
setwd(dir = "C:/Users/handy_test") #remember to replace the path with your own 
```

# What issues can this package solve?
As above we mentioned, the functions in the package have been categorised into five sections: Conversion, Summary, Annotation, Comparison and Visualization.

Functions in the Conversion section can, for example, be used to convert coordinates between the default and target mapfiles provided by users; this would be useful if we are detecting CNVs with two different reference genomes. The converted map files can be output in formats suitable for use in Plink, PennCNV and our package 'HandyCNV'. It is also possible to convert coordinates for intervals (CNV, CNVR, QTL et al) between the defaults and target maps, which might useful when comparing the results to those of other researchers. Functions included in the Summarization section allow the formatting and plotting of outputs, such as producing detailed summaries of CNVs and ROHs, making CNV summary plots, generating CNV regions (CNVRs), and plotting CNVR distribution maps.

The Annotation section allows the downloading of reference gene lists from the UCSC website, enabling the annotation of CNV, CNVR, ROH or other intervals with gene locations. The Comparison section facilitates the comparison of CNVs, CNVRs, Gene Frequent Lists and any other intervals. The highlight of this section is that comparison of CNVs and CNVRs will produce reports detailing comparison results on both individual and population levels, and all comparison results are relative to both input files. Finally, the Visualization section contains functions to support customised CNV and ROH plotting by chromosome, specific sample, regions of interest, or target genes (plotting by gene is not yet available for ROH). It is also possible to annotate the plots with other information, such as gene locations, log R ratio, B Allele Frequency, SNP genotype, LD or the source of CNVs. These functions can also plot genes separately from reference gene lists, which might useful for comparing outputs with plots from other studies.
