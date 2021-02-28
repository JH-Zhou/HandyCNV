#' Title call_gene
#'
#' This function used for annotating genes for any interval by giving the interval list and reference gene list.
#' The interval file requires at least consists of four columns, they are ID, Chr, Start and End.
#' The first column must be the ID column, the column of Chr should only contain the number of chromosome, the units of Start and End columns are the basepair.
#'
#' @param refgene internal reference gene-list are "ARS_ens", "ARS_UCSC" and "UMD_UCSC". Or provide the reference genes list corresponding to your data
#' @param interval Could be CNV, ROH, QTL or any interval list. At least comprised by four columns Interval_ID, Chr, Start and End
#' @param clean_cnv The output data from cnv_clean function
#' @param folder set the name of folder to save results
#' @import dplyr
#'
#' @importFrom  data.table fread fwrite setkey foverlaps
#'
#' @return Annotated file and brief summary
#' @export call_gene
#'
call_gene <- function(refgene = "ARS_ens", interval = NULL, clean_cnv = NULL, folder = "ARS"){
  if(!file.exists(paste0("call_gene_", folder))){
    dir.create(paste0("call_gene_", folder))
    print(paste0("A new folder ",  "call_gene_", folder, " was created in working directory."))
  }

  #if(missing(refgene)){
    gene <- fread(file = refgene, header = TRUE)
  #} else{
  #  gene <- fread(file = refgene, header = FALSE)
  #  names(gene) <- c("bin", "name", "Chr", "strand", "Start", "End",
  #                   "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
  #                   "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  #}

  #if(refgene == "ARS_ens"){
  #  refgene = system.file("extdata", "Demo_data/gene_annotation/ensGene_ars_210202.txt", package = "HandyCNV")
  #  gene <- fread(file = refgene, header = TRUE)
  #} else if(refgene == "ARS_UCSC"){
  #  refgene = system.file("extdata", "Demo_data/gene_annotation/refGene_ars1.2.txt", package = "HandyCNV")
  #  gene <- fread(file = refgene, header = FALSE)
  #  names(gene) <- c("bin", "name", "Chr", "strand", "Start", "End",
  #                   "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
  #                   "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  #} else if(refgene == "UMD_UCSC"){
  #  refgene = system.file("extdata", "Demo_data/gene_annotation/refGene_umd3.1.txt", package = "HandyCNV")
  #  gene <- fread(file = refgene, header = FALSE)
  #  names(gene) <- c("bin", "name", "Chr", "strand", "Start", "End",
  #                   "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
  #                   "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  #} else{
  #  gene <- fread(file = refgene, header = FALSE)
  #    names(gene) <- c("bin", "name", "Chr", "strand", "Start", "End",
  #                     "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
  #                     "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  #}

  #else{
  # print("Wrong argument input. Defaults refgene can only be 'ARS_ens' or 'ARS_UCSC' or 'UMD_UCSC'")
  #}

  gene$Chr <- suppressWarnings(as.integer(sub("chr", "", gene$Chr))) #convert Chr to integer in order to use foverlap in next step
  interval <- fread(file = interval, header = TRUE)
  if(ncol(interval) < 4){
    cat("Not enought columns provide in 'interval' list, \nat least four columns required:\n'Interval_ID' 'Chr' 'Start' 'End'
        ")
  }
  names(interval)[1] = "ID"    #Replace names of table in order to the summarize in the final step,
                               #It required provide the Interval ID in the first column
                               #And at least four columns, Interval_ID (CNVR or ROH or QTL), Chr, Start and End
  setkey(gene, Chr, Start, End)
  cnvr_gene <- foverlaps(interval, gene, by.x = c("Chr", "Start", "End"), type = "any")

  print("Starting to check gene annotation status in interval file....")
  #Summarize how many CNVRs got gene annotated
  cnvr_gene <- cnvr_gene %>%
               mutate(Check_gene = if_else(is.na(Start) & is.na(End), true = "Non_gene", false = "Has_gene"))

  cnvr_has_gene <- subset(cnvr_gene, Check_gene == "Has_gene") #extract a subset which contain the CNVR has gene only
  gene_number <- length(unique(cnvr_has_gene$name2)) #count the number of unique genes
  num_CNVR_has_gene <- nrow(unique(cnvr_has_gene[,c("Chr", "i.Start", "i.End")])) #count the number of unique CNVR which got gene annotation
  num_CNVR_no_gene <- nrow(unique(cnvr_gene[,c("Chr", "i.Start", "i.End")])) - num_CNVR_has_gene #"i.Start" and "i.End" are the location of Chr

  gene_summary <- data.frame(matrix(nrow = 1, ncol = 3)) #create a data table to save summary results
  names(gene_summary) <- c("Interval_Has_Gene", "Invertal_Without_Gene", "Total_Number_of_Genes") #assign columns name
  gene_summary[1, 1:3] <- c(num_CNVR_has_gene, num_CNVR_no_gene, gene_number) #assign relative value into table
  print("The summary of annotation results as shown below:")
  print(gene_summary)

  #report gene as list format for annotation use in David Annotation
  gene_list <- cnvr_gene %>%
               select(c("ID", "Chr", "name2"))
  fwrite(gene_list, file = paste0("call_gene_", folder, "/gene_list.annotation"), sep = "\t", quote = F)

  #summarize results group by ID then paste all gene into one cell
  window_gene <- cnvr_gene %>%
    group_by(ID) %>%
    summarise(gene_name = paste(unique(name2), collapse = ",")) %>% #collapse the gene name into windows
    left_join(unique(cnvr_gene[,c("ID", "Chr", "i.Start", "i.End")]), by = "ID") %>% #matching the rest information to window_id
    arrange(Chr, i.Start) %>% #descending order by number_roh
    mutate(num_gene = lengths(strsplit(gene_name, split = ","))) %>%#count the number of genes in window
    select(ID, Chr, i.Start, i.End, num_gene, gene_name)

  fwrite(window_gene, file = paste0("call_gene_", folder, "/interval_gene_summarise_table.txt"), sep = "\t", quote = F)

  fwrite(cnvr_gene, file = paste0("call_gene_", folder, "/interval_annotation.txt"), sep = "\t", quote = FALSE)
  fwrite(gene_summary, file = paste0("call_gene_", folder, "/call_gene_summary.txt"), sep = "\t", quote = FALSE)

  ###########################################################3
  #when the input of clean_cnv provided, run the rest of codes
  if (!is.null(clean_cnv)){
    cnv <- fread(file = clean_cnv)
    setkey(gene, Chr, Start, End)
    cnv_annotation <- foverlaps(cnv, gene)
    cnv_annotation <- cnv_annotation %>%
                      rename(g_Start = Start,
                             g_End = End,
                             CNV_Start = i.Start,
                             CNV_End = i.End)
    print("Starting to check gene annotation status in CNV file....")

    cnv_annotation <- cnv_annotation %>%
                      mutate(Check_gene = if_else(is.na(g_Start) & is.na(g_End), true = "Non_gene", false = "Has_gene"))


    cnv_gene <- subset(cnv_annotation, Check_gene == "Has_gene") #extract CNV which has gene
    gene_freq <- cnv_gene %>% group_by(name2) %>% count(name2, name = "Frequency", sort = TRUE)
    print(paste0(nrow(gene_freq), " genes were matched in the CNV and CNVR, top 10 frequent genes as below: "))
    print(gene_freq[1:10, ])

    # assign CNVR ID and positions into gene frequency list
    cnvr_has_gene_unique <- unique(subset(cnvr_has_gene, select = c("name2", "ID", "Chr", "i.Start", "i.End")))
    # in this file still has some duplicated genes which located in multiple CNVRs, here we need to report its out
    # the function check_dup_gene is used to report all duplicated genes
    check_dup_gene <- function(cnvr_has_gene_unique) {
      from_top <- which(duplicated(cnvr_has_gene_unique$name2))
      from_last <- which(duplicated(cnvr_has_gene_unique$name2, fromLast = TRUE))
      all_dup <- sort(unique(c(from_top, from_last)))
    }

    dup_gene <- cnvr_has_gene_unique[check_dup_gene(cnvr_has_gene_unique), ] # extract all duplicated genes
    print("These genes are duplicated in multiple CNVR, will only remain the first CNVR_ID in final gene frequency list report!")
    print(dup_gene)

    cnvr_has_gene_unique_pure <- cnvr_has_gene_unique[-which(duplicated(cnvr_has_gene_unique$name2)), ] # exclude the rows with duplicated gene name, only remained the first value
    gene_freq_location <- merge(gene_freq, cnvr_has_gene_unique_pure, by = "name2", all.x = TRUE)
    names(gene_freq_location)[names(gene_freq_location) == "i.Start"] = "CNVR_Start"
    names(gene_freq_location)[names(gene_freq_location) == "i.End"] = "CNVR_End"
    gene_freq_location <- gene_freq_location %>%
                          arrange(-Frequency)

    fwrite(cnv_annotation, file = paste0("call_gene_", folder, "/cnv_annotation.txt"), sep = "\t", quote = FALSE)
    fwrite(gene_freq_location, file = paste0("call_gene_", folder, "/gene_freq_cnv.txt"), sep = "\t", quote = FALSE)
    if(file.exists(paste0("call_gene_", folder, "/interval_annotation.txt")) & file.exists(paste0("call_gene_", folder, "/call_gene_summary.txt")) & file.exists(paste0("call_gene_", folder, "/cnv_annotation.txt")) & file.exists(paste0("call_gene_", folder, "/gene_freq_cnv.txt"))) {
      print("Task done. CNV, CNVR Annotation, summary files and the Frequency of Gene in CNV Region have been saved in Working directory.")
    } else {
      print("Task failed, please check your input file format carefully!!")
    }
  } else {
    if(file.exists(paste0("call_gene_", folder, "/interval_annotation.txt")) & file.exists(paste0("call_gene_", folder, "/call_gene_summary.txt"))) {
      print("Task done. Annotation and summary files have been saved in Working directory.")
    } else {
      print("Task failed, please check your input file format carefully!!")
    }
  }

}


