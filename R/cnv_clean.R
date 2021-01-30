#' Title  cnv_clean
#' Used for cleaning the default CNV results into a standard format for the further use.
#' Now support to read the results from PennCNV and CNVPartition
#'
#' @param cnvpartition
#' @param penncnv
#' @param penn_id_sep
#'
#'
#' @import data.table tidyr
#'
#' @return
#' @export
#'
#' @examples
cnv_clean <- function(cnvpartition = NULL, penncnv = NULL, penn_id_sep = "cnv/") {
  #creat a directory to store output files
  if (!file.exists("clean_cnv")){
    dir.create("clean_cnv")
  }

  if(!is.null(cnvpartition)) {
    cnvpart <- fread(file = cnvpartition, skip = 7)
    names(cnvpart) <- c("Sample_ID", "Chr", "Start", "End", "CNV_Value", "CNV_Conf", "Comment", "Empty")
    #cnvpart_pure <- cnvpart_pure[-c(grep("X", cnvpart_pure$Chr)), ] #delete cnv on the chr X
    cnvpart$Chr <- sub("X", "30", cnvpart$Chr)
    cnvpart$Start <- as.numeric(cnvpart$Start)
    cnvpart$End <- as.numeric(cnvpart$End)
    cnvpart$CNV_Conf <- as.numeric(cnvpart$CNV_Conf)
    cnvpart$Length <- cnvpart$End - cnvpart$Start + 1 #add a new column as Length
    cnvpart <- cnvpart[cnvpart$CNV_Conf >= 35, ] #delete the confedence score lesser than 35
    cnvpart_roh <- cnvpart[c(grep("2", cnvpart$CNV_Value)),]
    cnvpart_pure <- cnvpart[-c(grep("2", cnvpart$CNV_Value)), ] #delete 2 copy cnv
    cnvpart_pure <- cnvpart_pure[cnvpart_pure$Length <= 5000000, ] #delete CNV larger than 5 Mb

    average_indiv_cnv <- round(nrow(cnvpart_pure)/length(unique(cnvpart_pure$Sample_ID)), 2)
    print(paste0("There are ", length(unique(cnvpart_pure$Sample_ID))," individuals with ",  nrow(cnvpart_pure), " CNVs in total."))
    print(paste0("The average number of CNV on each Individual is ", average_indiv_cnv))

    summary_cnvpart <- cnvpart_pure %>% group_by(CNV_Value) %>% summarise("N" = n(), "Average Length" = mean(Length), "Min Length" = min(Length), "Max Length" = max(Length))
    print("The basic summary of each CNV type as following:")
    print(summary_cnvpart)

    fwrite(summary_cnvpart, file = "cnvpart_summary", sep = "\t", quote = FALSE, col.names = TRUE)

    fwrite(cnvpart_pure, file = "clean_cnv/cnvpart_clean.cnv", sep = "\t", quote = FALSE)
    fwrite(cnvpart_roh, file = "clean_cnv/cnvpart_roh.cnv", sep = "\t", quote = FALSE)

    if (file.exists("clean_cnv/cnvpart_clean.cnv")){
       print("Task finished, Clean CNV, ROH, CNV Summary results were saved in your working directory.")
     }

    else {
       print("Task failed, please check input cnv data carefully!!!")
     }
  }

  else{
    penn <- fread(file = penncnv, sep = " ", header = FALSE)
    penn <- separate(penn, V1, into = c("V1", "Start"), sep =":")
    penn <- separate(penn, Start, into = c("Start", "End"), sep = "-")
    penn <- separate(penn, V2, into = c("V2", "Num_SNP"), sep = "=")
    penn <- separate(penn, V3, into = c("V3", "Length"), sep = "=")
    penn <- separate(penn, V4, into = c("V4", "CNV_Number"), sep = ",")
    penn <- separate(penn, CNV_Number, into = c("CNV_Number", "CNV_Value"), sep = "=")
    penn <- separate(penn, V5, into = c("V5", "Sample_ID"), sep = penn_id_sep) # **this area need to define by userself
    penn <- separate(penn, Sample_ID, into = c("Sample_ID", "Format"), sep = ".txt")
    penn <- separate(penn, V6, into = c("V6", "Start_SNP"), sep = "=")
    penn <- separate(penn, V7, into = c("V7", "End_SNP"), sep = "=")
    penn <- separate(penn, V1, into = c("V1", "Chr"), sep = "chr")
    penn <- separate(penn, V8, into = c("V8", "Confedence_Score"), sep = "=") #If you did't use -confedence function while running penn detect.pl, don't need this step
    penn$Start <- as.numeric(penn$Start)
    penn$End <- as.numeric(penn$End)
    penn$Confedence_Score <- as.numeric(penn$Confedence_Score)
    penn$Length <- penn$End - penn$Start + 1
    penn <- penn[, c("Sample_ID", "Chr", "Start", "End", "Num_SNP", "Length", "CNV_Value", "Start_SNP", "End_SNP")] #"Confedence_Score")]
    #penn <- penn[penn$Confedence_Score >= 10] #delete the cnv which confedence lesser than 10
    penn <- penn[penn$Len <= 5000000, ]

    average_indiv_cnv <- round(nrow(penn)/length(unique(penn$Sample_ID)), 2)
    print(paste0("There are ", length(unique(penn$Sample_ID)), " individuals with ", nrow(penn), " CNVs in total."))
    print(paste0("The average number of CNVs on each Individual are around ", average_indiv_cnv))

    summary_cnv <- penn %>% group_by(CNV_Value) %>% summarise("N" = n(), "Average Length" = mean(Length), "Min Length" = min(Length), "Max Length" = max(Length))
    print("The basic summary of each CNV type as following:")
    print(summary_cnv)
    fwrite(penn, file = "clean_cnv/penncnv_clean.cnv", sep ="\t", quote = FALSE)
    fwrite(summary_cnv, file = "clean_cnv/penncnv_summary", sep = "\t", quote = FALSE, col.names = TRUE)
    if (file.exists("clean_cnv/penncnv_clean.cnv")){
      print("Task finished, please check Clean CNV results in your working directory.")
    }

    else {
      print("Task failed, please check input cnv data carefully!!!")
    }
  }
}
