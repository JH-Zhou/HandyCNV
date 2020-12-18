require(data.table, quietly = TRUE)
require(dplyr, quietly = TRUE)
call_gene <- function(refgene, cnvr= NULL, clean_cnv = NULL){
  if(!file.exists("call_gene")){
    dir.create("call_gene")
    print("A new folder 'call_gene' was created in working directory.")
  }

  gene <- fread(file = refgene, header = FALSE)
  names(gene) <- c("bin", "name", "Chr", "strand", "Start", "End",
                   "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
                   "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  gene$Chr <- as.integer(sub("chr", "", gene$Chr)) #convert Chr to integer in order to use foverlap in next step
  cnvr <- fread(file = cnvr)
  setkey(gene, Chr, Start, End)
  cnvr_gene <- foverlaps(cnvr, gene)

  print("Starting to check gene annotation status in CNVR file....")
  #Summarise how many CNVRs got gene annotated
  cnvr_gene$Check_gene <- ""  #add a empty column, to assign if this CNVR has a gene
  for (i in 1:nrow(cnvr_gene)){
    if (is.na(cnvr_gene$Start[i]) & is.na(cnvr_gene$End[i])) {
      cnvr_gene$Check_gene[i] = "Without-gene"
    } else{
      cnvr_gene$Check_gene[i] = "Has-gene"
    }
  }

  cnvr_has_gene <- subset(cnvr_gene, Check_gene == "Has-gene") #extract a subset which contain the CNVR has gene only
  gene_number <- length(unique(cnvr_has_gene$name2)) #count the number of unique genes
  num_CNVR_has_gene <- length(unique(cnvr_has_gene$CNVR_ID)) #count the number of unique CNVR which got gene annotation
  num_CNVR_no_gene <- length(unique(cnvr_gene$CNVR_ID)) - num_CNVR_has_gene

  gene_summary <- data.frame(matrix(nrow = 1, ncol = 3)) #creat a data table to save summary results
  names(gene_summary) <- c("Num of CNVR Has Gene", "Num of CNVR Without Gene", "Total Number of Genes") #assign columns name
  gene_summary[1, 1:3] <- c(num_CNVR_has_gene, num_CNVR_no_gene, gene_number) #assign relative value into table
  print("The summary of annotation results for CNVR as shown below:")
  print(gene_summary)

  fwrite(cnvr_gene, file = "call_gene/cnvr_annotation.txt", sep = "\t", quote = FALSE)
  fwrite(gene_summary, file = "call_gene/call_gene_summary.txt", sep = "\t", quote = FALSE)

  cnv <- fread(file = clean_cnv)
  setkey(gene, Chr, Start, End)
  cnv_annotation <- foverlaps(cnv, gene)
  names(cnv_annotation)[c(5, 6, 18, 19)] <- c("g_Start", "g_End", "CNV_Start", "CNV_End")
  print("Starting to check gene annotation status in CNV file....")

  cnv_annotation$Check_gene <- ""
  for (i in 1:nrow(cnv_annotation)){
    if (is.na(cnv_annotation$g_Start[i]) & is.na(cnv_annotation$g_End[i])){
      cnv_annotation$Check_gene[i] <- "Witout-gene"
    } else {
      cnv_annotation$Check_gene[i] <- "Has-gene"
    }
  }

  cnv_gene <- subset(cnv_annotation, Check_gene == "Has-gene") #extract CNV which has gene
  gene_freq <- cnv_gene %>% group_by(name2) %>% count(name2, name = "Frequency", sort = TRUE)
  print(paste0(nrow(gene_freq), " genes were matched in the CNV and CNVR, top 10 frequent genes as below: "))
  print(gene_freq[1:10, ])

  fwrite(cnv_annotation, file = "call_gene/cnv_annotation.txt", sep = "\t", quote = FALSE)
  fwrite(gene_freq, file = "call_gene/gene_freq_cnv.txt", sep = "\t", quote = FALSE)
  if(file.exists("call_gene/cnvr_annotation.txt") & file.exists("call_gene/call_gene_summary.txt") & file.exists("call_gene/cnv_annotation.txt") & file.exists("call_gene/gene_freq_cnv.txt")) {
    print("Task done. CNV, CNVR Annotation, summary files and the Frequency of Gene in CNV Region have been saved in Working directory.")
  } else {
    print("Task failed, please check your input file format carefully!!")
  }
}


