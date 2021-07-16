#' Clean CNV
#'
#' Import CNV call results produced by the software packages PennCNV and CNVPartition, and converts them into a standard format for use in other functions in the `HandyCNV` package.
#' It now support to accept a CNV list in the standard format, the standard format should have at least five columns with header: Sample_ID, Chr, Start, End, CNV_Value
#'
#' @param cnvpartition load CNV results from CNVPartition
#' @param penncnv load CNV results from PennCNV
#' @param penn_id_sep the separator in the `Sample ID` column of PennCNV results. Useful if the ID is bound to the path
#' @param standard_cnv Load a user-generated CNV input file. The following columns must be present: Sample_ID, Chr, Start, End, CNV_Value
#' @param drop_length exclude CNVs longer than this threshold, unit is "Mb"
#' @param folder set the name of the output folder
#'
#' @import dplyr
#' @importFrom data.table fread fwrite
#' @importFrom tidyr separate
#'
#' @return Formatted CNV results and brief summary files.
#' @export cnv_clean
#'
cnv_clean <- function(cnvpartition = NULL, penncnv = NULL, standard_cnv = NULL, drop_length = 5, penn_id_sep = "cnv/", folder = "cnv_clean") {
  #create a directory to store output files
  if (!file.exists(folder)){
    dir.create(folder)
    cat(paste0("Output folder '", folder, "' has been created in the working directory.\n"))
  }

  if(!is.null(cnvpartition)) {
    cnvpart <- fread(file = cnvpartition, skip = 7, header = FALSE)
    names(cnvpart) <- c("Sample_ID", "Chr", "Start", "End", "CNV_Value", "CNV_Conf", "Comment", "Empty")
    #cnvpart_pure <- cnvpart_pure[-c(grep("X", cnvpart_pure$Chr)), ] #delete cnv on the chr X
    cnvpart$Chr <- sub("X", "99", cnvpart$Chr)
    cnvpart$Start <- as.numeric(cnvpart$Start)
    cnvpart$End <- as.numeric(cnvpart$End)
    cnvpart$CNV_Conf <- as.numeric(cnvpart$CNV_Conf)
    cnvpart$Length <- cnvpart$End - cnvpart$Start + 1 #add a new column as Length
    cnvpart <- cnvpart[cnvpart$CNV_Conf >= 35, ] #delete the confidence score lesser than 35
    cnvpart_roh <- cnvpart[c(grep("2", cnvpart$CNV_Value)),]
    cnvpart_pure <- cnvpart[-c(grep("2", cnvpart$CNV_Value)), ] #delete 2 copy cnv
    cnvpart_pure <- cnvpart_pure[cnvpart_pure$Length <= drop_length*1000000, ] #delete CNV larger than 5 Mb

    average_indiv_cnv <- round(nrow(cnvpart_pure)/length(unique(cnvpart_pure$Sample_ID)), 2)
    cat(paste0("There are ", length(unique(cnvpart_pure$Sample_ID))," individuals with ",  nrow(cnvpart_pure), " CNVs in total.\n"))
    cat(paste0("The average number of CNV on each Individual is ", average_indiv_cnv, "\n"))

    summary_cnvpart <- cnvpart_pure %>% group_by(CNV_Value) %>% summarise("N" = n(), "Average Length" = round(mean(Length),digits = 0), "Min Length" = min(Length), "Max Length" = max(Length))
    cat("Basic summary stats by CNV type:\n")
    print(summary_cnvpart)

    fwrite(summary_cnvpart, file = paste0(folder, "/cnvpart_summary.txt"), sep = "\t", quote = FALSE, col.names = TRUE)

    fwrite(cnvpart_pure, file = paste0(folder, "/cnvpart_clean.cnv"), sep = "\t", quote = FALSE)
    fwrite(cnvpart_roh, file = paste0(folder, "/cnvpart_roh.cnv"), sep = "\t", quote = FALSE)

    return(cnvpart_pure)

    if (file.exists(paste0(folder,"/cnvpart_clean.cnv"))){
       cat(paste0("Task finished. Clean CNV, ROH, and CNV summary results have been saved in the '", folder, "' directory.\n"))
     }

    else {
       stop("Task failed: please check the format of your input CNV data carefully!!!")
     }
  }

    else if (!(is.null(penncnv))){
    penn <- fread(file = penncnv, sep = " ", header = FALSE)
    penn <- separate(penn, V1, into = c("V1", "Start"), sep =":")
    penn <- separate(penn, Start, into = c("Start", "End"), sep = "-")
    penn <- separate(penn, V2, into = c("V2", "Num_SNP"), sep = "=")
    penn <- separate(penn, V3, into = c("V3", "Length"), sep = "=")
    penn <- separate(penn, V4, into = c("V4", "CNV_Number"), sep = ",")
    penn <- separate(penn, CNV_Number, into = c("CNV_Number", "CNV_Value"), sep = "=")
    penn <- separate(penn, V5, into = c("V5", "Sample_ID"), sep = penn_id_sep) # **this area need to define by user
    penn <- separate(penn, Sample_ID, into = c("Sample_ID", "Format"), sep = ".txt")
    penn <- separate(penn, V6, into = c("V6", "Start_SNP"), sep = "=")
    penn <- separate(penn, V7, into = c("V7", "End_SNP"), sep = "=")
    penn <- separate(penn, V1, into = c("V1", "Chr"), sep = "chr")
    penn <- separate(penn, V8, into = c("V8", "Confidence_Score"), sep = "=") #If you didn't use -confidence function while running penn detect.pl, don't need this step
    penn$Start <- as.numeric(penn$Start)
    penn$End <- as.numeric(penn$End)
    penn$Confedence_Score <- as.numeric(penn$Confedence_Score)
    penn$Length <- penn$End - penn$Start + 1
    penn <- penn[, c("Sample_ID", "Chr", "Start", "End", "Num_SNP", "Length", "CNV_Value", "Start_SNP", "End_SNP")] #"Confidence_Score")]
    #penn <- penn[penn$Confidence_Score >= 10] #delete the cnv which confidence lesser than 10
    penn <- penn[penn$Length <= drop_length*1000000, ]

    average_indiv_cnv <- round(nrow(penn)/length(unique(penn$Sample_ID)), 2)
    cat(paste0("There are ", length(unique(penn$Sample_ID)), " individuals with ", nrow(penn), " CNVs in total.\n"))
    cat(paste0("The average number of CNVs present per individual is ", average_indiv_cnv, "\n"))

    summary_cnv <- penn %>%
                   group_by(CNV_Value) %>%
                   summarise("N" = n(),
                            "Average Length" = round(mean(Length), digits = 0),
                            "Min Length" = min(Length),
                            "Max Length" = max(Length))
    cat("Basic summary stats by CNV type:\n")
    print(summary_cnv)

    fwrite(penn, file = paste0(folder, "/penncnv_clean.cnv"), sep ="\t", quote = FALSE)
    fwrite(summary_cnv, file = paste0(folder, "/penncnv_summary.txt"), sep = "\t", quote = FALSE, col.names = TRUE)

    return(penn)
    if (file.exists(paste0(folder, "/penncnv_clean.cnv"))){
      cat(paste0("Task finished. Please check clean CNV results in the '", folder, "' directory.\n"))
    }
    else {
      stop("Task failed: please check input CNV data carefully!!!")
    }
    }
    else{

      if(typeof(standard_cnv) == "character"){
        user_cnv = fread(standard_cnv, header = TRUE)
      } else {
        user_cnv = standard_cnv
      }

      default_title <- c("Sample_ID", "Chr", "Start", "End", "CNV_Value")
      #check input format
      if(!(all(colnames(user_cnv) %in% default_title))){
        stop(paste0("Incorrect columns provided: the following columns are required: ", paste0(default_title, collapse=', ') ))
      }

      #calculate length and QC
      user_cnv <- user_cnv %>%
                  mutate(Length = End - Start + 1) %>%
                  filter(Length < drop_length * 1000000)

      #brief summary
      average_indiv_cnv <- round(nrow(user_cnv)/length(unique(user_cnv$Sample_ID)), 2)
      cat(paste0("There are ", length(unique(user_cnv$Sample_ID))," individuals with ",  nrow(user_cnv), " CNVs in total.\n"))
      cat(paste0("The average number of CNV per individual is ", average_indiv_cnv, "\n"))

      #summary CNV
      summary_cnv <- user_cnv %>%
                     group_by(CNV_Value) %>%
                     summarise("N" = n(),
                              "Average Length" = round(mean(Length),digits = 0),
                              "Min Length" = min(Length),
                              "Max Length" = max(Length))
      cat("Basic summary stats by CNV type:\n")
      print(summary_cnv)

      #write output
      fwrite(user_cnv, file = paste0(folder, "/cleancnv.cnv"), sep = "\t", quote = FALSE, col.names = TRUE)
      fwrite(summary_cnv, file = paste0(folder, "/cnv_summary.txt"), sep = "\t", quote = FALSE, col.names = TRUE)

      return(user_cnv)
      if (file.exists(paste0(folder, "/cleancnv.cnv"))){
        cat(paste0("Task finished. Please check clean CNV results in the '", folder, "' directory.\n"))
      }
      else {
        stop("Task failed: please check input CNV data carefully!!!")
      }
    }
}
