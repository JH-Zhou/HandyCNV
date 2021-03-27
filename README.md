# Welcome to HandyCNV
[![R-CMD-check](https://github.com/JH-Zhou/HandyCNV/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JH-Zhou/HandyCNV/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/JH-Zhou/HandyCNV/branch/master/graph/badge.svg?token=5ERXZH0GFW)](https://codecov.io/gh/JH-Zhou/HandyCNV)

An R package for Standardized Summary, Annotation, Comparison, and Visualization of CNV, CNVR and ROH

![Fig.1 Working Flow](https://github.com/JH-Zhou/HandyCNV/raw/master/man/figures/High_resolution_pipeline_HandyCNV.png)

# Introduction
This package was originally designed for the Post-analysis of CNV results inferred from PennCNV and CNVPartition (GenomeStudio). However, it has now been expanded to accept input files in standard formats for a wider range of applications. Our motivation is to provide a standard, reproducible and time-saving pipeline for the post-analysis of CNVs and ROHs detected from SNP genotyping data for the majority of diploid Species. The functions provided in this package can be categorised into five sections: Conversion, Summary, Annotation, Comparison and Visualization. The most useful features provided are: integrating summarized results, generating lists of CNVR, annotating the results with known gene positions, plotting CNVR distribution maps, and producing customised visualisations of CNVs and ROHs with gene and other related information on one plot. This package also supports a range of customisations, including the colour, size of high resolution figures, and output folder, avoiding conflict between the results of different runs. Running through all functions detailed in the vignette could help us to identify and explore the most interesting genomic regions more easily.

# Vignettes and Manual
The details examples please visit our Github pages: https://jh-zhou.github.io/HandyCNV/

# Installation and Prerequisites
First, to run this package, we need to make sure that R (Version >= 3.5.2) is installed in your computer (R download link: https://www.r-project.org/). Once R is installed, the 'HandyCNV' package can be installed from Github repository by running the following script. If you rarely used R, it may take more time to install the 'HandyCNV' for the first time.
```{r}
#install.packages("remotes")
#library(remotes)
#install_github(repo = "JH-Zhou/HandyCNV", auth_token = "3d2ac98e4c297bab332f1e68b3b2d49f3a17d6aa")
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

