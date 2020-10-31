cnv_clean <- function(cnvpartition = NULL, penncnv = NULL, penn_id_sep = "cnv/") {
  if(!is.null(cnvpartition)) {
    cnvpart <- fread(file = cnvpartition, skip = 7)
    names(cnvpart) <- c("Sample_ID", "Chr", "Start", "End", "CNV_Value", "CNV_Conf", "Comment", "Empty")
    cnvpart_pure <- cnvpart[-c(grep("2", cnvpart$CNV_Value)), ] #delete 2 copy cnv
    #cnvpart_pure <- cnvpart_pure[-c(grep("X", cnvpart_pure$Chr)), ] #delete cnv on the chr X
    cnvpart$Start <- as.numeric(cnvpart$Start)
    cnvpart$End <- as.numeric(cnvpart$End)
    cnvpart$CNV_Conf <- as.numeric(cnvpart$CNV_Conf)
    cnvpart_pure$Length <- cnvpart_pure$End - cnvpart_pure$Start + 1 #add a new column as Length
    cnvpart_pure <- cnvpart_pure[cnvpart_pure$Length <= 5000000, ] #delete CNV larger than 5 Mb
    cnvpart_pure <- cnvpart_pure[cnvpart_pure$CNV_Conf >= 35, ] #delete the confedence score lesser than 35
    fwrite(cnvpart_pure, file = "cnvpart_clean.cnv", sep = "\t", quote = FALSE)

    if (file.exists("cnvpart_clean.cnv")){
       print("Task finished, please check Clean CNV results in your working directory.")
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
    fwrite(penn, file = "penncnv_clean.cnv", sep ="\t", quote = FALSE)

    if (file.exists("penncnv_clean.cnv")){
      print("Task finished, please check Clean CNV results in your working directory.")
    }

    else {
      print("Task failed, please check input cnv data carefully!!!")
    }
  }
}
