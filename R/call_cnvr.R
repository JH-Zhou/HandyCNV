#' Title clean CNV
#' Formatting CNV results from PennCNV and CNVPartition to the standard format
#' Call CNVR, summary individual CNVR type, report CNVR frequency and type
#'
#' @param clean_cnv the clean cnv file was generate by clean_cnv function
#' @param roh roh file from clean_cnv function, only CNVPartition results will generate the roh results
#' @param chr_set set the maximum number of chromosome in CNV list
#' @param folder text, set the name of folder
#' @import dplyr
#' @importFrom data.table fread fwrite setkey foverlaps
#' @importFrom reshape2 dcast
#' @importFrom grDevices dev.off pdf png
#'
#' @return Four or five tables. Including CNVR list, brief summary and individual cnv and CNVR ID.
#' @export call_cnvr
#'
call_cnvr <- function(clean_cnv, roh = NULL, chr_set = 29, folder = "UMD") {
  if(!file.exists(paste0("call_cnvr_", folder))){
    dir.create(paste0("call_cnvr_", folder))
  }

  clean_cnv <- data.table::fread(file = clean_cnv, sep = "\t", header = TRUE)

  merge_cnvr <- function(cnv) {
    if (nrow(cnv) == 1) {
      return(cnv)
    }

    cnv <- cnv[order(cnv$Chr, cnv$Start),]
    cnvr_union = cnv[1, ]

    for (i in 2:nrow(cnv)) {
      rest_cnv <- cnv[i, ]

      if (cnvr_union$End[nrow(cnvr_union)] < rest_cnv$Start) {
        cnvr_union <- dplyr::bind_rows(cnvr_union, rest_cnv)
      } else if (cnvr_union$End[nrow(cnvr_union)] == rest_cnv$Start) {
        cnvr_union$End[nrow(cnvr_union)] <- rest_cnv$End
      }
      if (rest_cnv$End > cnvr_union$End[nrow(cnvr_union)]) {
        cnvr_union$End[nrow(cnvr_union)] <- rest_cnv$End
      }
    }
    return(cnvr_union)
  }

  cnvr <- data.frame()

  for ( i in 1:chr_set){
    cnv_chr <- clean_cnv[which(clean_cnv$Chr == i), ]
    cnvr_chr <- merge_cnvr(cnv = cnv_chr)
    cnvr <- rbind(cnvr, cnvr_chr)
    print(paste0("Choromsome ", i, " was processed."))
  }

  #extract CNVR information, recode for CNVR
  cnvr_union <- cnvr[, c("Chr", "Start", "End")]
  cnvr_union$CNVR_ID <- paste0("CNVR_", seq(1, nrow(cnvr_union), 1))
  data.table::setkey(cnvr_union, Chr, Start, End)
  cnv_cnvr <- data.table::foverlaps(clean_cnv, cnvr_union)
  names(cnv_cnvr)[names(cnv_cnvr) == "i.Start"] <- "CNV_Start"
  names(cnv_cnvr)[names(cnv_cnvr) == "i.End"] <- "CNV_End"

  #add frequent of CNVR
  cnvr_frequent <- cnv_cnvr %>% group_by(CNVR_ID) %>% count(CNVR_ID, name = "Frequent")
  cnvr_union_f <- merge(cnvr_union, cnvr_frequent, by = "CNVR_ID", sort = F)
  fwrite(cnvr_union_f, file = paste0("call_cnvr_", folder, "/no_freq_cnvr.txt"), sep = "\t", quote = FALSE)

  if (is.null(roh)) {
    #add type of CNVR
    cnvr_type <- cnv_cnvr %>% group_by(CNVR_ID) %>% count(CNV_Value, name = "Count")
    cnvr_type2 <- reshape2::dcast(cnvr_type, CNVR_ID ~ CNV_Value, value.var = "Count")
    cnvr_type2$Type <- NA
    for (i in 1:nrow(cnvr_type2)) {
      if (is.na(cnvr_type2[i, c("0")]) & is.na(cnvr_type2[i, c("1")])){
        cnvr_type2$Type[i] <- "Gain"
      }

      else if (is.na(cnvr_type2[i, c("3")]) & is.na(cnvr_type2[i, c("4")])) {
        cnvr_type2$Type[i] <- "Loss"
      }

      else{cnvr_type2$Type[i] <- "Mixed"}
    }

    cnvr_f_type <- merge(cnvr_union_f, cnvr_type2, sort = F)
    cnvr_f_type$Length <- cnvr_f_type$End - cnvr_f_type$Start + 1

    print(paste0(nrow(cnvr_f_type), " CNVR generated in total."))

    #The overall summary
    cnvr_summary <- cnvr_f_type %>%
      group_by(Type) %>%
      summarise(N = n(), "Average Length" = mean(Length), "Min Legnth" = min(Length), "Max Length" = max(Length), "Total Length" = sum(Length))

    #Partial summary
    cnvr_chr_summary <- cnvr_f_type %>%
      group_by(Chr) %>%
      summarise("Total Length" = sum(Length), "Number of CNVR" = n())

    print("Overall summary of CNVR as following:")
    print(cnvr_summary)
    print("Partial summary of CNVR on each Chromosome as following:")
    print(cnvr_chr_summary)

    fwrite(cnv_cnvr, file = paste0("call_cnvr_", folder, "/individual_cnv_cnvr.txt"), sep = "\t", quote = FALSE)
    fwrite(cnvr_f_type, file = paste0("call_cnvr_", folder, "/cnvr.txt"), sep = "\t", quote = FALSE)
    fwrite(cnvr_summary, file = paste0("call_cnvr_", folder, "/cnvr_summary.txt"), sep = "\t", quote = FALSE, col.names = TRUE)
    fwrite(cnvr_chr_summary, file = paste0("call_cnvr_", folder, "/cnvr_chr_summary.txt"), sep = "\t", quote = FALSE, col.names = TRUE)

    if(file.exists(paste0("call_cnvr_", folder,"/cnvr.txt")) & file.exists(paste0("call_cnvr_", folder, "/individual_cnv_cnvr.txt"))) {
      print("Task done, CNVR results saved in the working directory.")
    } else {print("WARNING, lack of output file, please check format of your input file!!")}
  }

  else {
    cnvr_union_f$length <- cnvr_union_f$End - cnvr_union_f$Start + 1
    fwrite(cnvr_union_f, file = paste0("call_cnvr_", folder, "/roh.txt"), sep = "\t", quote = FALSE)

    if(file.exists(paste0("call_cnvr_", folder, "/roh.txt"))) {
      print("Task done, ROH results saved in the working directory.")
    } else {print("WARNING, lack of output file, please check format of your input file!!")}
  }
}


