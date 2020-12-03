require(data.table, quietly = TRUE)
require(dplyr, quietly = TRUE)
call_gene <- function(refgene = NULL, cnvr= NULL, clean_cnv = NULL){
  if(!file.exists("call_gene")){
    dir.create("call_gene")
  }

  gene <- fread(file = refgene)
  names(gene) <- c("bin", "name", "Chr", "strand", "Start", "End",
                   "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
                   "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  gene$Chr <- as.integer(sub("chr", "", gene$Chr))
  cnvr <- fread(file = cnvr)
  setkey(gene, Chr, Start, End)
  cnvr_gene <- foverlaps(cnvr, gene)
  fwrite(cnvr_gene, file = "call_gene/cnvr_annotation.txt", sep = "\t", quote = FALSE)
  if(file.exists("call_gene/cnvr_annotation.txt")) {
    print("Task done, CNVR Annotation file saved in your Working directory.")
  } else {
    print("Task failed, please check your input file format carefully!!")
  }

  cnv <- fread(file = clean_cnv)
  setkey(gene, Chr, Start, End)
  cnv_annotation <- foverlaps(cnv, gene)
  names(cnv_annotation)[c(5, 6, 18, 19)] <- c("g_Start", "g_End", "CNV_Start", "CNV_End")
  gene_freq <- cnv_annotation %>% group_by(name2) %>% count(name2, name = "Frequent", sort = TRUE)
  print(paste0(nrow(gene_freq)-1, " genes were matched in the CNV and CNVR, top 10 frequent genes as below: "))
  print(gene_freq[1:11, ])

  fwrite(cnv_annotation, file = "call_gene/cnv_annotation.txt", sep = "\t", quote = FALSE)
  fwrite(gene_freq, file = "call_gene/gene_freq_cnv.txt", sep = "\t", quote = FALSE)
  if(file.exists("call_gene/cnv_annotation.txt") & file.exists("call_gene/gene_freq_cnv.txt")) {
    print("Task done. CNV, CNVR Annotation files and the Frequency of Gene in CNV Region have been saved in your Working directory.")
  } else {
    print("Task failed, please check your input file format carefully!!")
  }
}


