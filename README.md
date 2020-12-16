# HandyCNV
R package for fast and normalized summarizing CNV results from SNP intensity data 

# All files we need is CNV results, SNP maps, Reference gene file, Pedigree,Plink bim, bed and fam files, Snp signal intensity
#1. CNV results, now only support PennCNV and CNVPartition output files (will adding a standard format to fit all kinds of demands)
#2. SNP maps (If user need to convert the coordinates, should provide at least two assembly map file) (need to add the four fixed colums requirments)
#3. rederence gene files, this files can download from UCSC website
#4. Pedigree (Not necessary, only use in function 'cnv_visual' to plot the source of CNV, its useful when check some high heritability CNVs)
#5. Plink bim, bed and fam files (Not necessary, only use in the tenth function 'plot_all')
#6. Snp signal intensity (Not necessary, only use in the tenth function 'plot_all')


#here is an example of how to use HandyCNV
#Copy this code in Rstudio to install HandyCNV: remotes::install_github(repo = "JH-Zhou/HandyCNV", auth_token = "3d2ac98e4c297bab332f1e68b3b2d49f3a17d6aa")


setwd("C:/Users/Jinghang/Desktop/HandyCNV_test_file_1") # Users should set the working directory by themself
library(HandyCNV)
library(data.table)
library(gaston)
library(graphics)
library(base2grob)
library(cowplot)
library(tidyr)
library(dplyr)
library(ggrepel)
library(scales)
library(rgl)
#1. First function: convert_map
HandyCNV::convert_map(umd_map = "map/final_403_XJB_CNV.map", 
                      standard_ARS_map = "map/9913_ARS1.2_139977_GGPHDV3_marker_name_180910.map")

#2. Second function:cnv_clean
#2.1 Clean CNV results for PennCNV, need add a seperate penn_id-sep
HandyCNV::cnv_clean(penncnv = "cnv_results/ARS_PennCNV_WGC.goodcnv", penn_id_sep = "cnv/") #for ARS PennCNV results
HandyCNV::cnv_clean(penncnv = "cnv_results/UMD_PennCNV_WGC.goodcnv", penn_id_sep = "cnv/") #for UMD PennCNV results
#2.2 Clean CNV for CNVParttition
HandyCNV::cnv_clean(cnvpartition = "cnv_results/UMD_CNVPartition.txt")

#3. Third function: cnv_summary_plot
HandyCNV::cnv_summary_plot(clean_cnv = "2_1_clean_cnv_ARS_Penn/penncnv_clean.cnv", 
                           plot_sum_1 = "yes",
                           plot_sum_2 = "yes") #for ARS PennCNV results

HandyCNV::cnv_summary_plot(clean_cnv = "2_2_clean_cnv_UMD_Penn/penncnv_clean.cnv", 
                           plot_sum_1 = "yes",
                           plot_sum_2 = "yes") #for UMD PennCNV results

HandyCNV::cnv_summary_plot(clean_cnv = "2_3_clean_cnv_UMD_Part/cnvpart_clean.cnv", 
                           plot_sum_1 = "yes",
                           plot_sum_2 = "yes") #for UMD CNVPartition results

#4 Forth function:call_cnvr
HandyCNV::call_cnvr(clean_cnv = "2_1_clean_cnv_ARS_Penn/penncnv_clean.cnv") #for ARS PennCNV results
HandyCNV::call_cnvr(clean_cnv = "2_2_clean_cnv_UMD_Penn/penncnv_clean.cnv") #for UMD PennCNV results
HandyCNV::call_cnvr(clean_cnv = "2_3_clean_cnv_UMD_Part/cnvpart_clean.cnv") #for UMD CNVPartition results

#5 Fifth function:call_gene, notice the version of reference gene should correspond to the CNV verison
#which means the ARS CNV and CNVR results should use ARS reference genome
#The UMD CNV and CNVR results should use UMD reference genome
HandyCNV::call_gene(refgene = "gene_annotation/refGene_ars1.2.txt", 
          cnvr = "4_1_call_cnvr_ARS_Penn/cnvr.txt", 
          clean_cnv = "2_1_clean_cnv_ARS_Penn/penncnv_clean.cnv") #for ARS PennCNV results

HandyCNV::call_gene(refgene = "gene_annotation/refGene_umd3.1.txt", 
                    cnvr = "4_2_call_cnvr_UMD_Penn/cnvr.txt", 
                    clean_cnv = "2_2_clean_cnv_UMD_Penn/penncnv_clean.cnv") #for UMD PennCNV results

HandyCNV::call_gene(refgene = "gene_annotation/refGene_umd3.1.txt",
                    cnvr = "4_3_call_cnvr_UMD_Part/cnvr.txt",
                    clean_cnv = "2_3_clean_cnv_UMD_Part/cnvpart_clean.cnv") #for UMD CNVPartition results

#6 Sixth function: here we only take ARS PennCNV results as example 
#6.1 plot all CNV distribution on all Chromsomes
HandyCNV::cnv_visual(clean_cnv = "2_1_clean_cnv_ARS_Penn/penncnv_clean.cnv", max_chr_length = 160)

#6.2 plot one of your interest chromosome 
HandyCNV::cnv_visual(clean_cnv = "2_1_clean_cnv_ARS_Penn/penncnv_clean.cnv", chr_id = 17)

#6.3 plot one of your interest CNV region 
HandyCNV::cnv_visual(clean_cnv = "2_1_clean_cnv_ARS_Penn/penncnv_clean.cnv", chr_id = 17, start_position = 20, end_position = 25)

#6.4 plot one of your interest individual
HandyCNV::cnv_visual(clean_cnv = "2_1_clean_cnv_ARS_Penn/penncnv_clean.cnv", individual_id = "202031420083_R05C01")

#6.5 plot the cnvr sourse information
#user should provide a pedigree at least include the Sample_ID, Sire_ID two columns
#Now only support three source information, Sire_ID, Herd and Source(The sire's original country)
HandyCNV::cnv_visual(clean_cnv = "2_1_clean_cnv_ARS_Penn/penncnv_clean.cnv", chr_id = 17, start_position = 20, end_position = 25, 
                     report_id = "yes", 
                     pedigree = "Pedigree.csv") 

#7 Seventh function:cnvr_plot, we only take the ARS PennCNV results as example
#7.1 plot the cnvr distribution on each chromosome, output figture name called:cnvr_plot.png
#Uesr should define the assembly argument specific the ARS or UMD results will been plot
HandyCNV::cnvr_plot(cnvr = "4_1_call_cnvr_ARS_Penn/cnvr.txt", assembly = "ARS")

#7.2 plot all the high frequent CNVR, input file have to include the annotated CNV results 'cnv_annotation'.
#user could custom the frequency threshold by define the sample_size and common_cnv_threshold arguments
#in this example we set sample_size = 388 and common_cnvThreshold = 0.05
#which means only plot the cnvr with frequecy higher than 388 * 0.05 = 19.4
HandyCNV::cnvr_plot(cnvr = "4_1_call_cnvr_ARS_Penn/cnvr.txt", 
                    cnv_annotation = "5_1_call_gene_ARS_Penn/cnv_annotation.txt", 
                    sample_size = 388, 
                    common_cnv_threshold = 0.05)

#8 Eighth function: compare_CNV
#8.1 compare the CNV results with same reference assembly only need to input two clean CNV results
HandyCNV::compare_cnv(cnv_umd = "2_2_clean_cnv_UMD_Penn/penncnv_clean.cnv", 
                      cnv_ars = "2_3_clean_cnv_UMD_Part/cnvpart_clean.cnv")

#8.2 compare the CNV results with different reference assembly, need to provide the two verion's map extra
#this two version assembly map file was generated in the first step, called 'umd_ars_map.map'
#this function will convert the coordinate for bith cnv results at first
#then make comparison, the first principal is that the comparison will execute by the lastest ARS coordinate
HandyCNV::compare_cnv(cnv_umd = "2_2_clean_cnv_UMD_Penn/penncnv_clean.cnv", 
                      cnv_ars = "2_1_clean_cnv_ARS_Penn/penncnv_clean.cnv",
                      umd_ars_map = "1_covert_map/umd_ars_map.map")


#9 Ninth function: compare_cnvr
#9.1 compare the CNVR results with same reference assembly only need to input two CNVR results
HandyCNV::compare_cnvr(cnvr_umd = "4_2_call_cnvr_UMD_Penn/cnvr.txt", 
                       cnvr_ars = "4_3_call_cnvr_UMD_Part/cnvr.txt")

#9.2 compare the CNVR results with different reference assembly, need to provide the two verion's map extra
#this two version assembly map file was generated in the first step, called 'umd_ars_map.map'
#this function will convert the coordinate for bith cnvr results at first then make comparison #
#the first principal is that the comparison will execute by the lastest ARS coordinate
#the difference between compare CNVR to CNV is that CNVR more concern about the overlapping length
#but compare CNV we more concern about the overlapping number of CNVs
#and CNV have the individual level and population level
#but CNVR only on the population level
HandyCNV::compare_cnvr(cnvr_umd = "4_2_call_cnvr_UMD_Penn/cnvr.txt", 
                       cnvr_ars = "4_1_call_cnvr_ARS_Penn/cnvr.txt", 
                       umd_ars_map = "1_covert_map/umd_ars_map.map")
