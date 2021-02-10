#' Title compare_cnv
#'The idea to compare CNV between ars and umd are find out how many differents are there
#' we defined 6 comparison standards of CNV
#' 1) overlaped
#' 1.same start and end, same SNP inside, fully overlap
#' 2.same start and end, diffrent snp number, fully overlap
#' 3.different start or end, overlaped, partial overlap
#' 2) non-overlap
#' 4.missing start or end position
#' 5.End <= start
#' 6.different start or end, non-overlap
#' according to the codintion, the first thing is to match coordinates for both version
#' then find overlap cnv and non overlapcnv
#' then summarize how many CNVs are in above standards
#'
#' @param cnv_umd
#' @param cnv_ars
#' @param umd_ars_map
#' @param width_1
#' @param height_1
#' @import dplyr data.table scales ggplot2
#'
#' @return
#' @export compare_cnv
#'
#' @examples
compare_cnv <- function(cnv_umd, cnv_ars, umd_ars_map = NULL, width_1 = 14, height_1 = 11) {

  #default plot function
  plot_comparison <- function(cnv_checkover, title_fig, width_1 =14, height_1 =11) {
    cnv <- cnv_checkover
    title_f = title_fig
    cnv_cal <- cnv %>% count(CNV_Value, Check_overlap) %>%
      mutate(percent_total = n / sum(n)) %>%
      group_by(CNV_Value) %>%
      mutate(percent_group = n /sum(n))

    cnv_freq <- cnv_cal %>% group_by(CNV_Value) %>% summarise(percent_total  = sum(percent_total), num = sum(n))

    png(res = 300, filename = paste0(title_f, ".png"), width = width_1, height = height_1, bg = "transparent", units = "cm")
    compare_plot <- ggplot(cnv_cal, aes(x = CNV_Value, y = percent_total, fill = Check_overlap)) +
      geom_col() +
      geom_text(data = subset(cnv_cal, Check_overlap == "Overlap"), aes(label = scales::percent(percent_group)), color = "blue", position = position_stack(0.5)) +
      geom_text(inherit.aes = FALSE, data = cnv_freq, aes(x = CNV_Value, y = percent_total, label = num), vjust = -0.5, hjust = 1.3) +
      geom_text(inherit.aes = FALSE, data = cnv_freq, aes(x = CNV_Value, y = percent_total, label = scales::percent(percent_total)), vjust = -0.5, hjust  =-0.5) +
      scale_y_continuous(labels=scales::percent) +
      theme_classic() +
      theme(legend.position = c(0.95, 0.9), legend.title = element_blank()) +
      labs(x = "CNV Value", y = "Percentage of CNV Number", caption = "Note: The integer and percentage on the top of bar plot indicates the total number and percentage of CNV.\nPercentage with blue color indicates the percent of Overlapped CNVs within a CNV Value group.")
    print(compare_plot)
    dev.off()
    fwrite(cnv_cal, file = paste0(title_f, ".summary"), sep = "\t", quote = FALSE)
  }

  #dealing with data
  if (is.null(umd_ars_map)) {
    cnv_ars <- fread(cnv_ars)
    cnv_ars$version <- "Verision_ARS" # add verision in dataframe
    colnames(cnv_ars) <- paste(colnames(cnv_ars), "ARS", sep = "_") #add suffix to all colnames
    ars_colnames <- colnames(cnv_ars) #set original column names use for extarting columns after matching

    cnv_umd <- fread(cnv_umd)
    cnv_umd$version <- "Version_UMD"
    colnames(cnv_umd) <- paste(colnames(cnv_umd), "UMD", sep = "_")
    umd_colnames <- colnames(cnv_umd)

    #######compare results in UMD at first-------------------------------------------------------------------------
    #1. find overlaped CNV between UCD and ARS
    setkey(cnv_ars, Sample_ID_ARS, Chr_ARS, Start_ARS, End_ARS)
    overlap_umd <- foverlaps(cnv_umd, cnv_ars, by.x = c("Sample_ID_UMD", "Chr_UMD", "Start_UMD", "End_UMD"), type = "any", nomatch = NULL)

    #might have some duplicated rows after find overlap, because of some CNV in ARS larger than UMD
    # or some CNV in UMD larger than ARS caused by the SNP position and density
    uniq_overlap_umd <- subset(overlap_umd, select = umd_colnames)

    #find out non-overlap CNV, original CNV - overlaped CNV
    final_overlap_umd <- unique(uniq_overlap_umd)

    #2. find non-overlap CNV between UND and ARS
    non_overlap_umd <- dplyr::setdiff(cnv_umd, final_overlap_umd)

    if (nrow(non_overlap_umd) == 0) {
      print("These two input data are completely same, comparison will stoped.")
    }

    #3. cheking if overlap and non-overlap results are correct by counting the total number
    if (nrow(non_overlap_umd) + nrow(final_overlap_umd) == nrow(cnv_umd)) {
      print(paste0("Comparison Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total ", nrow(cnv_umd)))
    } else {print("Comparison failed in validation, please use the original output files from HandyCNV as the input files")}

    #merge all cnv results then make compasion plot
    final_overlap_umd$Check_overlap <- "Overlap"
    non_overlap_umd$Check_overlap <- "Non-Overlap"
    checkover_indiv_1 <- rbind(final_overlap_umd, non_overlap_umd)
    colnames(checkover_indiv_1) <- sub("_UMD", "", colnames(checkover_indiv_1))

    plot_comparison(cnv_checkover = checkover_indiv_1, title_fig = "checkover_indiv_1", width_1 = width_1, height_1 = height_1)

    #4. write out the overlap and non-overlap CNVs
    fwrite(final_overlap_umd, file = "overlap_cnv_1.indiv", sep = "\t", quote = FALSE)
    fwrite(non_overlap_umd, file = "non_overlap_cnv_1.indiv", sep = "\t", quote = FALSE)
    fwrite(checkover_indiv_1, file = "cnv_all_inidv_1.checkoverlap", sep = "\t", quote = FALSE)

    #5.summarize difference
    overlap_percent_indiv_1 <- round(nrow(final_overlap_umd) / nrow(cnv_umd), 3) * 100
    non_overlap_percent_indiv_1 <- round(nrow(non_overlap_umd) / nrow(cnv_umd), 3) * 100
    print(paste0("The number of overlaped CNVs on individual level is ", nrow(final_overlap_umd), ", which is around ", overlap_percent_indiv_1, " percent in the first file"))
    print(paste0("The number of Non-overlaped CNVs on individual level is ", nrow(non_overlap_umd), ", which is around ", non_overlap_percent_indiv_1, " percent in the first file"))


    #setkey(cnv_umd_ars, Chr_UMD, Start_UMD, End_UMD)
    #find overlap on population level
    setkey(cnv_ars, Chr_ARS, Start_ARS, End_ARS)
    pop_overlap <- foverlaps(cnv_umd, cnv_ars, by.x = c("Chr_UMD", "Start_UMD", "End_UMD"), type = "any", nomatch = NULL)
    pop_overlap_umd <- subset(pop_overlap, select = umd_colnames)
    final_pop_overlap_umd <- unique(pop_overlap_umd)
    non_overlap_pop_umd <- dplyr::setdiff(cnv_umd, final_pop_overlap_umd)
    if (nrow(non_overlap_pop_umd) + nrow(final_pop_overlap_umd) == nrow(cnv_umd)) {
      print(paste0("Comparison to the first file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total of ", nrow(cnv_umd)))
    } else {print("Comparison to the first file was failed in validation, please use the original output files from HandyCNV as the input files")}

    fwrite(final_pop_overlap_umd, file = "overlap_cnv_1.popu", sep = "\t", quote = FALSE)
    fwrite(non_overlap_pop_umd, file = "non_overlap_cnv_1.popu", sep = "\t", quote = FALSE)

    #make comparison plot
    final_pop_overlap_umd$Check_overlap <- "Overlap"
    non_overlap_pop_umd$Check_overlap <- "Non-Overlap"
    checkover_pop_1 <- rbind(final_pop_overlap_umd, non_overlap_pop_umd)
    colnames(checkover_pop_1) <- sub("_UMD", "", colnames(checkover_pop_1))


    plot_comparison(cnv_checkover = checkover_pop_1, title_fig = "checkover_pop_1", width_1 = width_1, height_1 = height_1)

    fwrite(checkover_pop_1, file = "cnv_all_population_1.checkoverlap", sep = "\t", quote = FALSE)

    #5.summarize difference
    #UMD
    overlap_percent_pop <- round(nrow(final_pop_overlap_umd) / nrow(cnv_umd), 3) * 100
    non_overlap_percent_pop <- round(nrow(non_overlap_pop_umd) / nrow(cnv_umd), 3) * 100
    print(paste0("The number of overlaped CNVs on population level is ", nrow(final_pop_overlap_umd), ", which is around ", overlap_percent_pop, " percent in first file."))
    print(paste0("The number of Non-overlaped CNVs on population level is ", nrow(non_overlap_pop_umd), ", which is around ", non_overlap_percent_pop, " percent in first file"))

    overlap_summary <- data.frame("Item" = c("Overlapped CNVs", "Non-overlapped CNVs"),
                                  "Individual Level" = c(overlap_percent_indiv_1, non_overlap_percent_indiv_1),
                                  "Population Level" = c(overlap_percent_pop, non_overlap_percent_pop))
    print("The final comparison results of the first file as follows: ")
    print(overlap_summary)
    fwrite(overlap_summary, file = "overlap_cnv_1.summary", sep = "\t", quote = FALSE)
    print("Comparison to the first file was finished.")


    ##########compare results in ARS at second-------------------------------------------------------------------
    #1. find overlaped CNV between UCD and ARS
    setkey(cnv_umd, Sample_ID_UMD, Chr_UMD, Start_UMD, End_UMD)
    #because some snp cannot find in ARS map, so we removed these CNV used right_umd instead
    overlap_ars <- foverlaps(cnv_ars, cnv_umd, by.x = c("Sample_ID_ARS", "Chr_ARS", "Start_ARS", "End_ARS"), type = "any", nomatch = NULL)

    #might have some duplicated rows after find overlap, because of some CNV in ARS larger than UMD
    # or some CNV in UMD larger than ARS caused by the SNP position and density
    uniq_overlap_ars <- subset(overlap_ars, select = ars_colnames)

    #find out non-overlap CNV, original CNV - overlaped CNV
    final_overlap_ars <- unique(uniq_overlap_ars)

    #2. find non-overlap CNV between UND and ARS
    non_overlap_ars <- dplyr::setdiff(cnv_ars, final_overlap_ars)

    #3. cheking if overlap and non-overlap results are correct by counting the total number
    if (nrow(non_overlap_ars) + nrow(final_overlap_ars) == nrow(cnv_ars)) {
      print(paste0("Comparison to second file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total ", nrow(cnv_ars)))
    } else {print("Comparison to second file was failed in validation, please use the original output files from HandyCNV as the input files")}

    #4. write out the overlap and non-overlap CNVs
    fwrite(final_overlap_ars, file = "overlap_cnv_ars_2.indiv", sep = "\t", quote = FALSE)
    fwrite(non_overlap_ars, file = "non_overlap_cnv_ars_2.indiv", sep = "\t", quote = FALSE)

    final_overlap_ars$Check_overlap <- "Overlap"
    non_overlap_ars$Check_overlap <- "Non-Overlap"
    checkover_indiv_2 <- rbind(final_overlap_ars, non_overlap_ars)
    colnames(checkover_indiv_2) <- sub("_ARS", "", colnames(checkover_indiv_2))

    fwrite(checkover_indiv_2, file = "cnv_all_indiv_2.checkoverlap", sep = "\t", quote = FALSE)

    plot_comparison(cnv_checkover = checkover_indiv_2, title_fig = "checkover_indiv_2", width_1 = width_1, height_1 = height_1)

    #5.summarize difference
    overlap_percent_indiv_2 <- round(nrow(final_overlap_ars) / nrow(cnv_ars), 3) * 100
    non_overlap_percent_indiv_2 <- round(nrow(non_overlap_ars) / nrow(cnv_ars), 3) * 100
    print(paste0("The number of overlaped CNVs on individual level is ", nrow(final_overlap_ars), ", which is around ", overlap_percent_indiv_2, " percent in second file."))
    print(paste0("The number of Non-overlaped CNVs on individual level is ", nrow(non_overlap_ars), ", which is around ", non_overlap_percent_indiv_2, " percent in second file."))


    #setkey(cnv_umd_ars, Chr_UMD, Start_UMD, End_UMD)
    #find overlap on population level
    setkey(cnv_umd, Chr_UMD, Start_UMD, End_UMD)
    pop_overlap_2 <- foverlaps(cnv_ars, cnv_umd, by.x = c("Chr_ARS", "Start_ARS", "End_ARS"), type = "any", nomatch = NULL)
    pop_overlap_ars <- subset(pop_overlap_2, select = ars_colnames)
    final_pop_overlap_ars <- unique(pop_overlap_ars)
    non_overlap_pop_ars <- dplyr::setdiff(cnv_ars, final_pop_overlap_ars)
    if (nrow(non_overlap_pop_ars) + nrow(final_pop_overlap_ars) == nrow(cnv_ars)) {
      print(paste0("Comparison to the second file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total of ", nrow(cnv_ars)))
    } else {print("Comparison to the second file was failed in validation, please use the original output files from HandyCNV as the input files")}

    fwrite(final_pop_overlap_ars, file = "overlap_cnv_ars.popu", sep = "\t", quote = FALSE)
    fwrite(non_overlap_pop_ars, file = "non_overlap_cnv_ars.popu", sep = "\t", quote = FALSE)

    final_pop_overlap_ars$Check_overlap <- "Overlap"
    non_overlap_pop_ars$Check_overlap <- "Non-Overlap"
    checkover_pop_2 <- rbind(final_pop_overlap_ars, non_overlap_pop_ars)
    colnames(checkover_pop_2) <- sub("_ARS", "", colnames(checkover_pop_2))


    plot_comparison(cnv_checkover = checkover_pop_2, title_fig = "checkover_pop_2", width_1 = width_1, height_1 = height_1)
    fwrite(checkover_pop_2, file = "checkover_pop_2.txt", sep = "\t", quote = FALSE)

    #5.summarize difference
    #ARS
    overlap_percent_pop_2 <- round(nrow(final_pop_overlap_ars) / nrow(cnv_ars), 3) * 100
    non_overlap_percent_pop_2 <- round(nrow(non_overlap_pop_ars) / nrow(cnv_ars), 3) * 100
    print(paste0("The number of overlaped CNVs in seceond file on population level are ", nrow(final_pop_overlap_ars), ", which is around ", overlap_percent_pop_2, " percent"))
    print(paste0("The number of Non-overlaped CNVs in second file on population level are ", nrow(non_overlap_pop_ars), ", which is around ", non_overlap_percent_pop_2, " percent"))

    overlap_summary_2 <- data.frame("Item" = c("Overlapped CNVs", "Non-overlapped CNVs"),
                                    "Individual Level" = c(overlap_percent_indiv_2, non_overlap_percent_indiv_2),
                                    "Population Level" = c(overlap_percent_pop_2, non_overlap_percent_pop_2))
    print("The final comparison results of the second file as follows: ")
    print(overlap_summary_2)
    fwrite(overlap_summary_2, file = "overlap_cnv_2.summary", sep = "\t", quote = FALSE)
    print("Task done. Comparison results were saved in your working directory")
  }

  else {
    #convert coordinate for CNV and CNVR
    cnv_ars <- fread(cnv_ars)
    cnv_ars$version <- "Verision_ARS" # add verision in dataframe
    colnames(cnv_ars) <- paste(colnames(cnv_ars), "ARS", sep = "_") #add suffix to all colnames
    ars_colnames <- colnames(cnv_ars) #set original column names use for extarting columns after matching

    cnv_umd <- fread(cnv_umd)
    cnv_umd$version <- "Version_UMD"
    colnames(cnv_umd) <- paste(colnames(cnv_umd), "UMD", sep = "_")
    umd_colnames <- colnames(cnv_umd)

    two_map <- fread(umd_ars_map)
    colnames(two_map) <- paste(colnames(two_map), "Map", sep = "_")

    #cnv_umd_ars <- merge(cnv_umd, two_map, by.x = c("Chr", "Start"), by.y = c("Chr_UMD", "Position_UMD"), all.x = TRUE)

    #1. convert the umd result to ars
    #matching the start position by snp name, validating by matching the original position and new mathced position
    #cnv_umd_ars <- merge(cnv_umd, two_map, by.x = "Start_SNP_UMD", by.y = "Name_Map", all.x = TRUE)

    #There are some duplicated rows after merge progress, because of there are some different SNP with same location in the map file
    #To sovle this problem, we should check duplicated row after each merge step and remove the duplicates
    cnv_umd_ars <- merge(cnv_umd, two_map, by.x = c("Chr_UMD", "Start_UMD"), by.y = c("Chr_def_Map", "Position_def_Map"), all.x = TRUE)
    dup_index <- grep("TRUE", duplicated(cnv_umd_ars[, c("Chr_UMD", "Start_UMD", "Sample_ID_UMD")]))
    if (length(dup_index) > 0){
      cnv_umd_ars <- cnv_umd_ars[-dup_index, ]
      print(paste0("There are ", length(dup_index), " duplicated SNPs after converting the Start Position for the first file, they were all successfully removed."))
    } else {
      print("Non duplicated SNPs were detected after converting the Start Position for the first file.")
    }

    if (all(cnv_umd_ars$Start_UMD == cnv_umd_ars$Position_def_Map)){
      print("Matching results by the Start position passed the validation.")
    } else {print("matching results didn't pass the validation, please check your input data, make sure the input file all generate by HandyCNV")}

    cnv_umd_ars <- subset(cnv_umd_ars, select = c(umd_colnames, "Position_tar_Map"))
    names(cnv_umd_ars)[names(cnv_umd_ars) == "Position_tar_Map"] <- "Start_ARS_Map"
    umd_temp_colnames <- colnames(cnv_umd_ars)

    #matching the end position of CNV
    #cnv_umd_ars <- merge(cnv_umd_ars, two_map, by.x = "End_SNP_UMD", by.y = "Name_Map", all.x = TRUE)
    cnv_umd_ars <- merge(cnv_umd_ars, two_map, by.x = c("Chr_UMD", "End_UMD"), by.y = c("Chr_def_Map", "Position_def_Map"), all.x = TRUE)
    dup_index <- grep("TRUE", duplicated(cnv_umd_ars[, c("Chr_UMD", "End_UMD", "Sample_ID_UMD")]))
    if (length(dup_index) > 0){
      cnv_umd_ars <- cnv_umd_ars[-dup_index, ]
      print(paste0("There are ", length(dup_index), " duplicated SNPs after converting the End Position for the first file, they were all successfully removed."))
    } else {
      print("Non duplicated SNPs were detected after converting the End Position for the first file.")
    }

    if (all(cnv_umd_ars$End_UMD == cnv_umd_ars$Position_def_Map)){
      print("Matching results by the End position passed the validation.")
    } else {print("matching results didn't pass the validation, please check your input data, make sure the input file all generate by HandyCNV")}

    cnv_umd_ars <- subset(cnv_umd_ars, select = c(umd_temp_colnames, "Position_tar_Map"))
    names(cnv_umd_ars)[names(cnv_umd_ars) == "Position_tar_Map"] <- "End_ARS_Map"
    umd_convert_colnames <- colnames(cnv_umd_ars)


    #2.convert ars result to umd
    #cnv_ars_umd <- merge(cnv_ars, two_map, by.x = "Start_SNP_ARS", by.y = "Name_Map", all.x = TRUE)
    cnv_ars_umd <- merge(cnv_ars, two_map, by.x = c("Chr_ARS", "Start_ARS"), by.y = c("Chr_tar_Map", "Position_tar_Map"), all.x = TRUE)
    dup_index <- grep("TRUE", duplicated(cnv_ars_umd[, c("Chr_ARS", "Start_ARS", "Sample_ID_ARS")]))
    if (length(dup_index) > 0){
      cnv_ars_umd <- cnv_ars_umd[-dup_index, ]
      print(paste0("There are ", length(dup_index), " duplicated SNPs after converting the Start Position for the Second file, they were all successfully removed."))
    } else {
      print("Non duplicated SNPs were detected after converting the Start Position for the Second file.")
    }

    #checking if the matching results correct?
    if (all(cnv_ars_umd$Start_ARS == cnv_ars_umd$Position_tar_Map)){
      print("Matching results by the Start position passed the validation.")
    } else {print("matching results didn't pass the validation, please check your input data, make sure the input file all generate by HandyCNV")}
    cnv_ars_umd <- subset(cnv_ars_umd, select = c(ars_colnames, "Position_def_Map"))
    names(cnv_ars_umd)[names(cnv_ars_umd) == "Position_def_Map"] <- "Start_UMD_Map"
    ars_temp_colnames <- colnames(cnv_ars_umd)

    #matching end position
    #cnv_ars_umd <- merge(cnv_ars_umd, two_map, by.x = "End_SNP_ARS", by.y = "Name_Map", all.x = TRUE)
    cnv_ars_umd <- merge(cnv_ars_umd, two_map, by.x = c("Chr_ARS", "End_ARS"), by.y = c("Chr_tar_Map", "Position_tar_Map"), all.x = TRUE)
    dup_index <- grep("TRUE", duplicated(cnv_ars_umd[, c("Chr_ARS", "End_ARS", "Sample_ID_ARS")]))
    if (length(dup_index) > 0){
      cnv_ars_umd <- cnv_ars_umd[-dup_index, ]
      print(paste0("There are ", length(dup_index), " duplicated SNPs after converting the End Postion for the Second File, they were all successfully removed."))
    } else {
      print("Non duplicated SNPs were detected after converting the End Postion for the Second File.")
    }

    if (all(cnv_ars_umd$End_ARS == cnv_ars_umd$Position_tar_Map)){
      print("Matching results by the End position passed the validation.")
    } else {print("matching results didn't pass the validation, please check your input data, make sure the input file all generate by HandyCNV")}
    cnv_ars_umd <- subset(cnv_ars_umd, select = c(ars_temp_colnames, "Position_def_Map"))
    names(cnv_ars_umd)[names(cnv_ars_umd) == "Position_def_Map"] <- "End_UMD_Map"
    ars_convert_colnames <- colnames(cnv_ars_umd)

    #write out CNV with converted coordinates
    fwrite(cnv_umd_ars, file = "cleancnv_UtoA.coord", sep = "\t", quote = FALSE)
    fwrite(cnv_ars_umd, file = "cleancnv_AtoU.coord", sep = "\t", quote = FALSE)

    #two requrments for foverlap: 1) start <= end, 2) no NA in both start and end
    cnv_umd_ars$End_ARS_Map[is.na(cnv_umd_ars$End_ARS_Map)] <- 0
    cnv_umd_ars$Start_ARS_Map[is.na(cnv_umd_ars$Start_ARS_Map)] <- 0
    wrong_umd <- cnv_umd_ars[which((cnv_umd_ars$End_ARS_Map - cnv_umd_ars$Start_ARS_Map) <= 0), ]
    right_umd <- cnv_umd_ars[-c(which((cnv_umd_ars$End_ARS_Map - cnv_umd_ars$Start_ARS_Map) <= 0)), ]
    wrong_ars <- cnv_ars_umd[which((cnv_ars_umd$End_UMD_Map - cnv_ars_umd$Start_UMD_Map) <= 0), ]
    wrong_ars <- wrong_ars[-c(which(wrong_ars$End_UMD_Map == 0)), ]
    right_ars <- setdiff(cnv_ars_umd, wrong_ars)

    #######compare results in UMD at first-------------------------------------------------------------------------
    #1. find overlaped CNV between UCD and ARS
    setkey(cnv_ars_umd, Sample_ID_ARS, Chr_ARS, Start_ARS, End_ARS)
    #because some snp cannot find in ARS map, so we removed these CNV used right_umd instead
    overlap_umd <- foverlaps(right_umd, cnv_ars_umd, by.x = c("Sample_ID_UMD", "Chr_UMD", "Start_ARS_Map", "End_ARS_Map"), type = "any", nomatch = NULL)

    #might have some duplicated rows after find overlap, because of some CNV in ARS larger than UMD
    # or some CNV in UMD larger than ARS caused by the SNP position and density
    uniq_overlap_umd <- subset(overlap_umd, select = umd_convert_colnames)

    #find out non-overlap CNV, original CNV - overlaped CNV
    teamp_overlap_umd <- subset(overlap_umd, select = umd_convert_colnames)
    final_overlap_umd <- unique(teamp_overlap_umd)

    #2. find non-overlap CNV between UND and ARS
    non_overlap_umd <- dplyr::setdiff(cnv_umd_ars, final_overlap_umd)

    #3. cheking if overlap and non-overlap results are correct by counting the total number
    if (nrow(non_overlap_umd) + nrow(final_overlap_umd) == nrow(cnv_umd_ars)) {
      print(paste0("Comparison Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total ", nrow(cnv_umd_ars)))
    } else {print("Comparison failed in validation, please use the original output files from HandyCNV as the input files")}

    #4. write out the overlap and non-overlap CNVs
    fwrite(final_overlap_umd, file = "overlap_cnv.indiv", sep = "\t", quote = FALSE)
    fwrite(non_overlap_umd, file = "non_overlap_cnv.indiv", sep = "\t", quote = FALSE)


    final_overlap_umd$Check_overlap <- "Overlap"
    non_overlap_umd$Check_overlap <- "Non-Overlap"
    checkover_indiv_1 <- rbind(final_overlap_umd, non_overlap_umd)
    colnames(checkover_indiv_1) <- sub("_UMD", "", colnames(checkover_indiv_1))

    plot_comparison(cnv_checkover = checkover_indiv_1, title_fig = "checkover_indiv_1", width_1 = width_1, height_1 = height_1)

    #5.summarize difference
    overlap_percent_indiv <- round(nrow(final_overlap_umd) / nrow(cnv_umd_ars), 3) * 100
    non_overlap_percent_indiv <- round(nrow(non_overlap_umd) / nrow(cnv_umd_ars), 3) * 100
    print(paste0("The number of overlaped CNVs on individual level is ", nrow(final_overlap_umd), ", which is around ", overlap_percent_indiv, " percent"))
    print(paste0("The number of Non-overlaped CNVs on individual level is ", nrow(non_overlap_umd), ", which is around ", non_overlap_percent_indiv, " percent"))


    #setkey(cnv_umd_ars, Chr_UMD, Start_UMD, End_UMD)
    #find overlap on population level
    setkey(cnv_ars_umd, Chr_ARS, Start_ARS, End_ARS)
    pop_overlap <- foverlaps(right_umd, cnv_ars_umd, by.x = c("Chr_UMD", "Start_ARS_Map", "End_ARS_Map"), type = "any", nomatch = NULL)
    pop_overlap_umd <- subset(pop_overlap, select = umd_convert_colnames)
    final_pop_overlap_umd <- unique(pop_overlap_umd)
    non_overlap_pop_umd <- dplyr::setdiff(cnv_umd_ars, final_pop_overlap_umd)
    if (nrow(non_overlap_pop_umd) + nrow(final_pop_overlap_umd) == nrow(cnv_umd_ars)) {
      print(paste0("Comparison to the first file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total of ", nrow(cnv_umd_ars)))
    } else {print("Comparison to the first file was failed in validation, please use the original output files from HandyCNV as the input files")}

    fwrite(final_pop_overlap_umd, file = "overlap_cnv.popu", sep = "\t", quote = FALSE)
    fwrite(non_overlap_pop_umd, file = "non_overlap_cnv.popu", sep = "\t", quote = FALSE)

    #plot comparison
    final_pop_overlap_umd$Check_overlap <- "Overlap"
    non_overlap_pop_umd$Check_overlap <- "Non-Overlap"
    checkover_pop_1 <- rbind(final_pop_overlap_umd, non_overlap_pop_umd)
    colnames(checkover_pop_1) <- sub("_UMD", "", colnames(checkover_pop_1))


    plot_comparison(cnv_checkover = checkover_pop_1, title_fig = "checkover_pop_1", width_1 = width_1, height_1 = height_1)

    fwrite(checkover_pop_1, file = "cnv_all_population_1.checkoverlap", sep = "\t", quote = FALSE)

    #5.summarize difference
    #UMD
    overlap_percent_pop <- round(nrow(final_pop_overlap_umd) / nrow(cnv_umd_ars), 3) * 100
    non_overlap_percent_pop <- round(nrow(non_overlap_pop_umd) / nrow(cnv_umd_ars), 3) * 100
    print(paste0("The number of overlaped CNVs on population level is ", nrow(final_pop_overlap_umd), ", which is around ", overlap_percent_pop, " percent in first file."))
    print(paste0("The number of Non-overlaped CNVs on population level is ", nrow(non_overlap_pop_umd), ", which is around ", non_overlap_percent_pop, " percent in first file"))

    overlap_summary <- data.frame("Item" = c("Overlapped CNVs", "Non-overlapped CNVs"),
                                  "Individual Level" = c(overlap_percent_indiv, non_overlap_percent_indiv),
                                  "Population Level" = c(overlap_percent_pop, non_overlap_percent_pop))
    print("The final comparison results of the first file as follows: ")
    print(overlap_summary)
    fwrite(overlap_summary, file = "overlap_cnv.summary", sep = "\t", quote = FALSE)
    print("Comparison to the first file was finished.")


    ##########compare results in ARS at first-------------------------------------------------------------------
    #1. find overlaped CNV between UCD and ARS
    setkey(right_umd, Sample_ID_UMD, Chr_UMD, Start_ARS_Map, End_ARS_Map)
    #because some snp cannot find in ARS map, so we removed these CNV used right_umd instead
    overlap_ars <- foverlaps(cnv_ars_umd, right_umd, by.x = c("Sample_ID_ARS", "Chr_ARS", "Start_ARS", "End_ARS"), type = "any", nomatch = NULL)

    #might have some duplicated rows after find overlap, because of some CNV in ARS larger than UMD
    # or some CNV in UMD larger than ARS caused by the SNP position and density
    uniq_overlap_ars <- subset(overlap_ars, select = ars_convert_colnames)

    #find out non-overlap CNV, original CNV - overlaped CNV
    teamp_overlap_ars <- subset(overlap_ars, select = ars_convert_colnames)
    final_overlap_ars <- unique(teamp_overlap_ars)

    #2. find non-overlap CNV between UND and ARS
    non_overlap_ars <- dplyr::setdiff(cnv_ars_umd, final_overlap_ars)

    #3. cheking if overlap and non-overlap results are correct by counting the total number
    if (nrow(non_overlap_ars) + nrow(final_overlap_ars) == nrow(cnv_ars_umd)) {
      print(paste0("Comparison to second file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total ", nrow(cnv_ars_umd)))
    } else {print("Comparison to second file was failed in validation, please use the original output files from HandyCNV as the input files")}

    #4. write out the overlap and non-overlap CNVs
    fwrite(final_overlap_ars, file = "overlap_cnv_ars.indiv", sep = "\t", quote = FALSE)
    fwrite(non_overlap_ars, file = "non_overlap_cnv_ars.indiv", sep = "\t", quote = FALSE)

    final_overlap_ars$Check_overlap <- "Overlap"
    non_overlap_ars$Check_overlap <- "Non-Overlap"
    checkover_indiv_2 <- rbind(final_overlap_ars, non_overlap_ars)
    colnames(checkover_indiv_2) <- sub("_ARS", "", colnames(checkover_indiv_2))

    fwrite(checkover_indiv_2, file = "cnv_all_indiv_2.checkoverlap", sep = "\t", quote = FALSE)

    plot_comparison(cnv_checkover = checkover_indiv_2, title_fig = "checkover_indiv_2", width_1 = width_1, height_1 = height_1)

    #5.summarize difference
    overlap_percent_indiv_2 <- round(nrow(final_overlap_ars) / nrow(cnv_ars_umd), 3) * 100
    non_overlap_percent_indiv_2 <- round(nrow(non_overlap_ars) / nrow(cnv_ars_umd), 3) * 100
    print(paste0("The number of overlaped CNVs on individual level is ", nrow(final_overlap_ars), ", which is around ", overlap_percent_indiv_2, " percent in second file."))
    print(paste0("The number of Non-overlaped CNVs on individual level is ", nrow(non_overlap_ars), ", which is around ", non_overlap_percent_indiv_2, " percent in second file."))


    #setkey(cnv_umd_ars, Chr_UMD, Start_UMD, End_UMD)
    #find overlap on population level
    setkey(right_umd, Chr_UMD, Start_ARS_Map, End_ARS_Map)
    pop_overlap_2 <- foverlaps(cnv_ars_umd, right_umd, by.x = c("Chr_ARS", "Start_ARS", "End_ARS"), type = "any", nomatch = NULL)
    pop_overlap_ars <- subset(pop_overlap_2, select = ars_convert_colnames)
    final_pop_overlap_ars <- unique(pop_overlap_ars)
    non_overlap_pop_ars <- dplyr::setdiff(cnv_ars_umd, final_pop_overlap_ars)
    if (nrow(non_overlap_pop_ars) + nrow(final_pop_overlap_ars) == nrow(cnv_ars_umd)) {
      print(paste0("Comparison to the second file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total of ", nrow(cnv_ars_umd)))
    } else {print("Comparison to the second file was failed in validation, please use the original output files from HandyCNV as the input files")}

    fwrite(final_pop_overlap_ars, file = "overlap_cnv_ars.popu", sep = "\t", quote = FALSE)
    fwrite(non_overlap_pop_ars, file = "non_overlap_cnv_ars.popu", sep = "\t", quote = FALSE)


    final_pop_overlap_ars$Check_overlap <- "Overlap"
    non_overlap_pop_ars$Check_overlap <- "Non-Overlap"
    checkover_pop_2 <- rbind(final_pop_overlap_ars, non_overlap_pop_ars)
    colnames(checkover_pop_2) <- sub("_ARS", "", colnames(checkover_pop_2))


    plot_comparison(cnv_checkover = checkover_pop_2, title_fig = "checkover_pop_2", width_1 = width_1, height_1 = height_1)
    fwrite(checkover_pop_2, file = "checkover_pop_2.txt", sep = "\t", quote = FALSE)

    #5.summarize difference
    #UMD
    overlap_percent_pop_2 <- round(nrow(final_pop_overlap_ars) / nrow(cnv_ars_umd), 3) * 100
    non_overlap_percent_pop_2 <- round(nrow(non_overlap_pop_ars) / nrow(cnv_ars_umd), 3) * 100
    print(paste0("The number of overlaped CNVs in seceond file on population level are ", nrow(final_pop_overlap_ars), ", which is around ", overlap_percent_pop_2, " percent"))
    print(paste0("The number of Non-overlaped CNVs in second file on population level are ", nrow(non_overlap_pop_ars), ", which is around ", non_overlap_percent_pop_2, " percent"))

    overlap_summary_2 <- data.frame("Item" = c("Overlapped CNVs", "Non-overlapped CNVs"),
                                    "Individual Level" = c(overlap_percent_indiv_2, non_overlap_percent_indiv_2),
                                    "Population Level" = c(overlap_percent_pop_2, non_overlap_percent_pop_2))
    print("The final comparison results of the second file as follows: ")
    print(overlap_summary_2)
    fwrite(overlap_summary_2, file = "overlap_cnv_2.summary", sep = "\t", quote = FALSE)
    print("Task done. Comparison results were saved in your working directory")
  }
}

