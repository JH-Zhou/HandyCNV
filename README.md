# Welcome to HandyCNV
[![R-CMD-check](https://github.com/JH-Zhou/HandyCNV/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JH-Zhou/HandyCNV/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/JH-Zhou/HandyCNV/branch/master/graph/badge.svg?token=5ERXZH0GFW)](https://codecov.io/gh/JH-Zhou/HandyCNV)

An R package for Standardized Summary, Annotation, Comparison, and Visualization of CNV, CNVR and ROH

<p align = "center">
<img src = "https://github.com/JH-Zhou/HandyCNV/blob/master/man/figures/High_resolution_pipeline_HandyCNV.jpg">
</p>
<p align = "center">
Main functions and outputs form HandyCNV
</p>

# Introduction
This package was originally designed for the Post-analysis of CNV results inferred from PennCNV and CNVPartition (GenomeStudio). However, it has now been expanded to accept input files in standard formats for a wider range of applications. Our motivation is to provide a standard, reproducible and time-saving pipeline for the post-analysis of CNVs and ROHs detected from SNP genotyping data for the majority of diploid Species. The functions provided in this package can be categorised into five sections: Conversion, Summary, Annotation, Comparison and Visualization. The most useful features provided are: integrating summarized results, generating lists of CNVR, annotating the results with known gene positions, plotting CNVR distribution maps, and producing customised visualisations of CNVs and ROHs with gene and other related information on one plot. This package also supports a range of customisations, including the colour, size of high resolution figures, and output folder, avoiding conflict between the results of different runs. Running through all functions detailed in the vignette could help us to identify and explore the most interesting genomic regions more easily.

# Vignettes and Manual
The details examples please visit our Github pages: https://jh-zhou.github.io/HandyCNV/

# Installation and Prerequisites
First, to run this package, we need to make sure that R (Version >= 3.5.2) is installed in your computer (R download link: https://www.r-project.org/). Once R is installed, the 'HandyCNV' package can be installed from Github repository by running the following script. If you rarely used R, it may take more time to install the 'HandyCNV' for the first time.
## 1. Method one, install from Github Repo directly
```{r}
install.packages("remotes") # Run this code if you haven't install 'remotes' package before 
remotes::install_github(repo = "JH-Zhou/HandyCNV@v1.1.5")
```

## 2. Method two, install manually 

If the first method cannnot work well for some reasons, we can manually download the 'Source code (Zip)' from the newly released tag at here: [Download Source Code](https://github.com/JH-Zhou/HandyCNV/releases/tag/v1.1.5)

Then install the Source Code from the local path by following code:
```{r}
install.packages("remotes") # Run this code if you haven't install 'remotes' package before 
remotes::install_local(path = "C:/Users/HandyCNV-1.1.5.zip") # Repalce 'C:/Users/' to your local path where you downloaded the Source Code
```
Then, we need to load the 'HandyCNV' package in order to run the following examples. This can be done using the `library` function as shown below.
```{r setup}
library(HandyCNV)
```

# What issues can this package solve?
Click the following link to browse the output in examples.

[1. How do we prepare the standard cnv input file for HandyCNV and get a quick summary?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#1-how-do-we-prepare-the-standard-cnv-input-file-for-handycnv-and-get-a-quick-summary-)

[2. How do we visualise CNVs?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#1-how-do-we-prepare-the-standard-cnv-input-file-for-handycnv-and-get-a-quick-summary-)

[3. What types of summary plots are available?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#3-what-types-of-summary-plots-are-available-)

[4. How do we generate CNVRs (CNV Regions) from CNV results?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#4-how-do-we-generate-cnvrs-cnv-regions-from-cnv-results-)

[5. How do we annotate genes for CNVRs and CNVs?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#5-how-do-we-annotate-genes-for-cnvrs-and-cnvs-)

[6. How to plot CNVR distribution map?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#6-how-to-plot-cnvr-distribution-map-)

[7. How can we plot all high frequency CNVRs at once?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#7-how-can-we-plot-all-high-frequency-cnvrs-at-once-)

[8. Can we compare CNVs between different result sets?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#8-can-we-compare-cnvs-between-different-result-sets-)

[9. What about CNVRs?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#9-what-about-cnvrs-)

[10. How do we find the consensus set of genes common to multiple CNV result sets?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#10-how-do-we-find-the-consensus-set-of-genes-common-to-multiple-cnv-result-sets-)

[11. How can we create the map file used in HandyCNV, to allow comparing CNVs between different reference genomes?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#11-how-can-we-create-the-map-file-used-in-handycnv-to-allow-comparing-cnvs-between-different-reference-genomes-)

[12. How can we integrate the CNVs and annotated genes with additional information, such as Log R ratio, B Allele Frequency, call rate, heterozygosity, missing value rate and Linkage Disequilibrium, and plot it as one figure?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#12-how-can-we-integrate-the-cnvs-and-annotated-genes-with-additional-information-such-as-log-r-ratio-b-allele-frequency-call-rate-heterozygosity-missing-value-rate-and-linkage-disequilibrium-and-plot-it-as-one-figure-)

[13. How can we make a plot to show the source of the CNVs?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#13-how-can-we-make-a-plot-to-show-the-source-of-the-cnvs-)

[14. How can we plot just the genes in a specific region, to save as a seperate figure?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#14-how-can-we-plot-just-the-genes-in-a-specific-region-to-save-as-a-seperate-figure-)

[15. How can we find regions with high frequencies of runs of homozygosity?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#15-how-can-we-find-regions-with-high-frequencies-of-runs-of-homozygosity-)

[16. How do we visualize ROHs?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#16-how-do-we-visualise-rohs-)

[17. How do we get haplotype for ROH region?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#17-how-do-we-get-haplotype-for-roh-region-)

[18. How do we convert coordinates for CNV, CNVR, ROH, or any other intervals?](https://jh-zhou.github.io/HandyCNV/articles/HandyCNV.html#18-how-do-we-convert-coordinates-for-cnv-cnvr-roh-or-any-other-intervals-)

# Feature request
If you have any special requirements for this package, please feel free to sumbit your demands via this link: [Submit Requirments](https://github.com/JH-Zhou/HandyCNV/issues/new?assignees=&labels=&template=feature-request.md), we are happy to add the new features to meet your needs. 

# Bug report
If you find any errors while using this package, please tell us via this link: [Bug Report](https://github.com/JH-Zhou/HandyCNV/issues/new?assignees=&labels=&template=bug_report.md), we will fix it as soon as possible.

# Citation
If this tool is useful for your academic research, please cite our publication: [Browse publication](https://www.biorxiv.org/content/10.1101/2021.04.05.438403v1)

J. Zhou, L. Liu, T. J. Lopdell, D. J. Garrick, and Y. Shi, “HandyCNV: Standardized Summary, Annotation, Comparison, and Visualization of CNV, CNVR and ROH,” doi: 10.1101/2021.04.05.438403.

# Current release: HandyCNV v1.1.5 Release Date: 2021/08/29

# What's new
Minor modifications, such as unifying input file formats and correcting spelling errors.

# Previous release: HandyCNV v1.1.4 Release Date: 2021/07/23

# What's new

## Major improvements:
1. Most functions now support reading variable object as input files;
2. Most functions now support returning the main output as object to R environment for the further operation;
3. New function 'get_samples' to extract samples ID by searching interested gene from CNV annotation list.
## Minor changes:
1. The 'call_cnvr' funtion now support generating CNVRs from CNV list that contains Chromosomes without CNVs information;
2. Add links of Horse_quCab2.0 genome reference and sheep 'oviAri3' reference genome into 'get_refgene' function;
3. Setup a standard table to support present comparison plot with empty group in 'compare_cnvr' function;
4. Add '-' as separator between the two recoded haplotypes in 'get_haplotype' function.   

# Previous release: HandyCNV v1.1.3 Release Date: 2021/05/26

# What's new

1. New function to plot SNP density from SNP genotyping map.
```{r, warning=FALSE}
plot_snp_density(map = "convert_map/target_plink.map", 
                 max_chr = 24, #optional
                 top_density = 60, #optional 
                 low_density = 20, #optional
                 color_top = "red", #optional
                 color_low = "blue", #optional
                 color_mid = "black", #optional
                 legend_position = c(0.9, 0.1), #optional 
                 x_label = "Physical position\n物理位置", #optional
                 y_label = "SNPs/Mb\n每1Mb区间的SNP数",#optional
                 ncol_1 = 5) 
#save the plot by 'ggsave'
#ggsave(filename = "snp_density.png", width = 26, height = 18, units = "cm", dpi = 350)
```
![Fig.Demo SNP density distribution](https://github.com/JH-Zhou/HandyCNV/blob/master/vignettes/snp_density.png)

2. Revised CNVs status distribution plot in 'cnv_summary_plot' function, force to appear the boxplot and line on chromosome that has no CNVs.


# Previous release: HandyCNV v1.1.2 Release Date: 2021/04/18

# What's new

Corrected the version number.

# Released Version: HandyCNV v1.1.1 Release Date: 2021/04/14

# What's new

## 1. Update call_cnvr.R …
1. The CNV list could only be loaded from the local directory through a 'Path' before, now supports to read data from working environment by checking the type of input file;
2. Support to return the CNVR list to working environment.

## 2. Update cnv_clean.R …
1. Support to return Clean CNV List to working environment.

## 3. Update cnv_visualising.R …
1. Support to load CNV List from working environment.

## 4. Update compare_cnvr.R …
1. Add new function to generate the Unique and Mutual CNVRs by uniting the overlapped CNVRs between two results. There are two purposes of this work, one is to better understand the overlapping CNVRs, the other is to mark the common regions on CNVR distribution map in 'cnvr_plot' function;

## 5. Update 'cnvr_plot' …
1. Add 'overlap_cnvr' argument to support to mark overlapped region on CNVR distribution map;
2. Add 'label_prop' argument to show the proportion of CNVRs length to total length of relative chromosome on CNVR map;
3. Add 'chr_col' argument to customize the color of Chromosome;
4. Add 'overlap_col' argument to customize color of overlapped CNVRs
5. Reduce the margin of final CNVR distribution map.

### New feature demo
```{r}
# Demo code:
cnvr_plot(cnvr = "./cnvr_combine_part_penn_umd/cnvr.txt", assembly = "UMD", 
          sample_size =  388, common_cnv_threshold = 0.05, 
          overlap_cnvr = "./compare_cnvr_Penn_UMD_Vs_Part_UMD/common_cnvr.txt", 
          gain_col = "deeppink1", loss_col = "deepskyblue3", mixed_col = "springgreen3", 
          folder = "./cnvr_combine_part_penn_umd/cnvr_map_common")
```       
![Fig.2 CNVR Map](https://github.com/JH-Zhou/HandyCNV/blob/master/vignettes/cnvr_plot.png)

## 6. Update compare_gene.R …
1. Add 'color_label' argument to display the color of genes that passed common threshold in the Heatmap;
