
#' Title compare_cnvr
#'#The idea to compare CNV between ars and umd are find out how many differents are there
#we defined 6 comparison standards of CNV
#1) overlaped
#1.same start and end, same SNP inside, fully overlap
#2.same start and end, diffrent snp number, fully overlap
#3.different start or end, overlaped, partial overlap
#2) non-overlap
#4.missing start or end position
#5.End <= start
#6.different start or end, non-overlap
#according to the codintion, the first thing is to match coordinates for both version
#then find overlap cnv and non overlapcnv
#then summarize how many CNVRs are in above standards
#'
#' @param cnvr_umd first cnvr list, not limited on umd or ars, default file is the result from call_cnvr function
#' @param cnvr_ars second cnvr list, not limited on umd or ars, default file is the result from call_cnvr function
#' @param umd_ars_map map file contains coordinats in both version of map. only need in comparison between the results from different versions. default file is generated from convert_map function
#' @param width_1 number to set the width of final plot size, unit is 'cm'
#' @param height_1 number to set the height of final plot size, unit is 'cm'
#' @param hjust_prop default value is 0.0. used to adjust horizontal position of the number of overlapped CNVR in the plot
#' @param hjust_num default value is 1.5. used to adjust horizontal position of the number of overlapped CNVR in the plot
#'
#' @import dplyr ggplot2
#' @importFrom data.table fread fwrite setkey foverlaps
#'
#' @return
#' @export
#'
#' @examples
compare_cnvr <- function(cnvr_umd, cnvr_ars, umd_ars_map = NULL, width_1 = 15, height_1 = 15, hjust_prop = 0.0, hjust_num = 1.5) {

  #default plot function
  plot_comparison <- function(cnv_checkover, title_fig, width_1 = 15, height_1 = 15, hjust_prop = 0.0, hjust_num = 1.5) {
    cnv <- cnv_checkover
    title_f = title_fig
    drop_name = "Overlap_length"
     cnvr_cal = unique(subset(cnv, select = !(colnames(cnv) %in% drop_name))) %>% group_by(Type, Check_overlap) %>%
      summarise(origi_length = sum(Length), num_CNVR = n_distinct(CNVR_ID))
    cnvr_over <- cnv %>% group_by(Type, Check_overlap) %>%
      summarise(overlap_len = sum(Overlap_length), num_overlap_oppsite = n_distinct(Overlap_length))

    cnvr_over[1,3] = cnvr_cal[1,3] + cnvr_cal[2,3] - cnvr_over[2,3]
    cnvr_over[3,3] = cnvr_cal[3,3] + cnvr_cal[4,3] - cnvr_over[4,3]
    cnvr_over[5,3] = cnvr_cal[5,3] + cnvr_cal[6,3] - cnvr_over[6,3]

    cnvr_cal <- merge(cnvr_cal, cnvr_over, by = c("Type", "Check_overlap"), all.x = TRUE)

    cnvr_cal = cnvr_cal %>% group_by(Type) %>%
      mutate(prop_overlap_len = overlap_len / sum(origi_length), prop_num = num_CNVR / sum(num_CNVR))
    #cnvr_cal %>% add_row(Total = "")

    print(paste0("CNVR comparison summary results in ", title_f," as following:"))
    print(cnvr_cal)

    #prepare title for plot
    overlap_sum <- cnvr_cal %>%
                   group_by(Check_overlap) %>%
                   summarise(num_CNVR = sum(num_CNVR),
                             overlap_len = sum(overlap_len))

    summary_title <- paste0(overlap_sum$num_CNVR[2]," CNVRs ",
                            round(overlap_sum$overlap_len[2]/sum(cnvr_cal$origi_length), 4) * 100, "% overlap Length")

    png(res = 300, filename = paste0(title_f, ".png"), width = width_1, height = height_1, bg = "transparent", units = "cm")
    compare_plot <- ggplot(cnvr_cal) +
      geom_col(aes(x = Type, y = overlap_len, fill = Check_overlap)) +
      scale_fill_brewer(palette = 'Set1') +
      geom_text(data = subset(cnvr_cal, Check_overlap == "Overlap"), aes(x= Type, y = overlap_len,label = scales::percent(prop_overlap_len)), color = "black", position = position_stack(0.5), hjust = hjust_prop) +
      geom_text(data = subset(cnvr_cal, Check_overlap == "Overlap"), aes(x= Type, y = overlap_len,label = num_CNVR), color = "orange", position = position_stack(0.5), hjust = hjust_num) +
      #geom_text(aes(x = Type, y = percent_total, label = num), vjust = -0.5, hjust = 1.3) +
      #geom_text(inherit.aes = FALSE, data = cnv_freq, aes(x = Type, y = percent_total, label = scales::percent(percent_total)), vjust = -0.5, hjust  =-0.5) +
      scale_y_continuous(labels=scales::unit_format(unit = "", scale = 1e-6)) +
      theme_classic() +
      theme(legend.position = "top", legend.title = element_blank()) +
      labs(x = summary_title, y = "Total Length (Mb)")
    print(compare_plot)
    dev.off()
    fwrite(cnvr_cal, file = paste0(title_f, ".summary"), sep = "\t", quote = FALSE)
  }

  #dealing with data
  if (is.null(umd_ars_map)) {
    cnv_ars <- fread(cnvr_ars)
    cnv_ars$version <- "Verision_ARS" # add verision in dataframe
    colnames(cnv_ars) <- paste(colnames(cnv_ars), "ARS", sep = "_") #add suffix to all colnames
    ars_colnames <- colnames(cnv_ars) #set original column names use for extarting columns after matching

    cnv_umd <- fread(cnvr_umd)
    cnv_umd$version <- "Version_UMD"
    colnames(cnv_umd) <- paste(colnames(cnv_umd), "UMD", sep = "_")
    umd_colnames <- colnames(cnv_umd)

    #find overlap on population level
    setkey(cnv_ars, Chr_ARS, Start_ARS, End_ARS)
    pop_overlap <- data.table::foverlaps(cnv_umd, cnv_ars, by.x = c("Chr_UMD", "Start_UMD", "End_UMD"), type = "any", nomatch = NULL)
    pop_overlap$Overlap_length <- pmin(pop_overlap$End_ARS, pop_overlap$End_UMD) - pmax(pop_overlap$Start_ARS, pop_overlap$Start_UMD) + 1
    pop_overlap_umd <- subset(pop_overlap, select = umd_colnames)
    final_pop_overlap_umd <- unique(pop_overlap_umd)
    non_overlap_pop_umd <- dplyr::setdiff(cnv_umd, final_pop_overlap_umd)
    if (nrow(non_overlap_pop_umd) + nrow(final_pop_overlap_umd) == nrow(cnv_umd)) {
      print(paste0("Comparison to the first file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVRs which are in total of ", nrow(cnv_umd)))
    } else {print("Comparison to the first file was failed in validation, please use the original output files from HandyCNV as the input files")}

    fwrite(final_pop_overlap_umd, file = "overlap_cnv_1.popu", sep = "\t", quote = FALSE)
    fwrite(non_overlap_pop_umd, file = "non_overlap_cnv_1.popu", sep = "\t", quote = FALSE)

    #make comparison plot
    final_pop_overlap_umd$Check_overlap <- "Overlap"
    non_overlap_pop_umd$Check_overlap <- "Non-Overlap"
    checkover_pop_1 <- rbind(final_pop_overlap_umd, non_overlap_pop_umd)
    colnames(checkover_pop_1) <- sub("_UMD", "", colnames(checkover_pop_1))
    checkover_pop_length <- merge(checkover_pop_1, pop_overlap[, c("CNVR_ID_UMD", "Overlap_length")], by.x = "CNVR_ID", by.y = "CNVR_ID_UMD", all.x = TRUE)
    #checkover_pop_length <- unique(checkover_pop_length)
    drop_name = "Overlap_length"
    checkover_pop_uniqe_1 = unique(subset(checkover_pop_length, select = !(colnames(checkover_pop_length) %in% drop_name))) #used for calculate overlapping


    try(plot_comparison(cnv_checkover = checkover_pop_length, title_fig = "checkover_pop_1", width_1 = width_1, height_1 = height_1, hjust_prop = hjust_prop, hjust_num = hjust_num), silent = TRUE)

    fwrite(checkover_pop_length, file = "checkover_pop_1.txt", sep = "\t", quote = FALSE)

    #5.summarize difference
    #UMD overlap number
    overlap_percent_pop <- round(nrow(final_pop_overlap_umd) / nrow(cnv_umd), 3) * 100
    non_overlap_percent_pop <- round(nrow(non_overlap_pop_umd) / nrow(cnv_umd), 3) * 100
    print(paste0("The number of overlaped CNVRs on population level is ", nrow(final_pop_overlap_umd), ", which is around ", overlap_percent_pop, " percent in first file."))
    print(paste0("The number of Non-overlaped CNVRs on population level is ", nrow(non_overlap_pop_umd), ", which is around ", non_overlap_percent_pop, " percent in first file"))

    #overlap length
    overlap_length_1 <- sum(checkover_pop_length$Overlap_length, na.rm = TRUE)
    cnvr_length_1 <- sum(checkover_pop_uniqe_1$Length, na.rm = TRUE)
    overlap_length_prop <- round((overlap_length_1 / cnvr_length_1),3) * 100
    print(paste0("The overlapping length is ", overlap_length_1, " bp, which is around ", overlap_length_prop, " percent in first file"))

    overlap_summary <- data.frame("Item" = c("Overlapped CNVRs", "Non-overlapped CNVRs", "In Total"),
                                  "Number of CNVR" = c(nrow(final_pop_overlap_umd), nrow(non_overlap_pop_umd), nrow(final_pop_overlap_umd) + nrow(non_overlap_pop_umd)),
                                  "Proportion of Number (%)" = c(overlap_percent_pop, non_overlap_percent_pop, 100),
                                  "Length(bp)" = c(overlap_length_1, cnvr_length_1- overlap_length_1, cnvr_length_1),
                                  "Proportion of Length (%)" = c(overlap_length_prop, 100 - overlap_length_prop, 100))
    print("The final comparison results of the first file as follows: ")
    print(overlap_summary)
    fwrite(overlap_summary, file = "overlap_cnv_1.summary", sep = "\t", quote = FALSE)
    print("Comparison to the first file was finished.")


  ##########compare results in ARS at second-------------------------------------------------------------------
   #setkey(cnv_umd_ars, Chr_UMD, Start_UMD, End_UMD)
    #find overlap on population level
    setkey(cnv_umd, Chr_UMD, Start_UMD, End_UMD)
    pop_overlap_2 <- data.table::foverlaps(cnv_ars, cnv_umd, by.x = c("Chr_ARS", "Start_ARS", "End_ARS"), type = "any", nomatch = NULL)
    pop_overlap_2$Overlap_length <- pmin(pop_overlap_2$End_ARS, pop_overlap_2$End_UMD) - pmax(pop_overlap_2$Start_ARS, pop_overlap_2$Start_UMD) + 1
    pop_overlap_ars <- subset(pop_overlap_2, select = ars_colnames)
    final_pop_overlap_ars <- unique(pop_overlap_ars)
    non_overlap_pop_ars <- dplyr::setdiff(cnv_ars, final_pop_overlap_ars)
    if (nrow(non_overlap_pop_ars) + nrow(final_pop_overlap_ars) == nrow(cnv_ars)) {
      print(paste0("Comparison to the second file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVRs which are in total of ", nrow(cnv_ars)))
    } else {print("Comparison to the second file was failed in validation, please use the original output files from HandyCNV as the input files")}

    fwrite(final_pop_overlap_ars, file = "overlap_cnv_ars.popu", sep = "\t", quote = FALSE)
    fwrite(non_overlap_pop_ars, file = "non_overlap_cnv_ars.popu", sep = "\t", quote = FALSE)

    final_pop_overlap_ars$Check_overlap <- "Overlap"
    non_overlap_pop_ars$Check_overlap <- "Non-Overlap"
    checkover_pop_2 <- rbind(final_pop_overlap_ars, non_overlap_pop_ars)
    colnames(checkover_pop_2) <- sub("_ARS", "", colnames(checkover_pop_2))
    checkover_pop_length_2 <- merge(checkover_pop_2, pop_overlap_2[, c("CNVR_ID_ARS", "Overlap_length")], by.x = "CNVR_ID", by.y = "CNVR_ID_ARS", all.x = TRUE)
    #checkover_pop_length_2 <- unique(checkover_pop_length_2)
    drop_name = "Overlap_length"
    checkover_pop_uniqe_2 = unique(subset(checkover_pop_length_2, select = !(colnames(checkover_pop_length_2) %in% drop_name))) #used for calculate overlapping


    try(plot_comparison(cnv_checkover = checkover_pop_length_2, title_fig = "checkover_pop_2", width_1 = width_1, height_1 = height_1, hjust_prop = hjust_prop, hjust_num = hjust_num), silent = TRUE)
    fwrite(checkover_pop_length_2, file = "checkover_pop_2.txt", sep = "\t", quote = FALSE)

    #5.summarize difference
    #ARS
    overlap_percent_pop_2 <- round(nrow(final_pop_overlap_ars) / nrow(cnv_ars), 3) * 100
    non_overlap_percent_pop_2 <- round(nrow(non_overlap_pop_ars) / nrow(cnv_ars), 3) * 100
    print(paste0("The number of overlaped CNVRs in seceond file on population level are ", nrow(final_pop_overlap_ars), ", which is around ", overlap_percent_pop_2, " percent"))
    print(paste0("The number of Non-overlaped CNVRs in second file on population level are ", nrow(non_overlap_pop_ars), ", which is around ", non_overlap_percent_pop_2, " percent"))

    #overlap length
    overlap_length_2 <- sum(checkover_pop_length_2$Overlap_length, na.rm = TRUE)
    cnvr_length_2 <- sum(checkover_pop_uniqe_2$Length, na.rm = TRUE)
    overlap_length_prop_2 <- round((overlap_length_2 / cnvr_length_2),3) * 100
    print(paste0("The overlapping length is ", overlap_length_2, " bp, which is around ", overlap_length_prop_2, " percent in first file"))

    overlap_summary_2 <- data.frame("Item" = c("Overlapped CNVRs", "Non-overlapped CNVRs", "In Total"),
                                    "Number of CNVR" = c(nrow(final_pop_overlap_ars), nrow(non_overlap_pop_ars), nrow(final_pop_overlap_ars) + nrow(non_overlap_pop_ars)),
                                    "Proportion of Number (%)" = c(overlap_percent_pop_2, non_overlap_percent_pop_2, 100),
                                    "Length(bp)" = c(overlap_length_2, cnvr_length_2- overlap_length_2, cnvr_length_2),
                                    "Proportion of Length (%)" = c(overlap_length_prop_2, 100 - overlap_length_prop_2, 100))
    print("The final comparison results of the second file as follows: ")
    print(overlap_summary_2)
    fwrite(overlap_summary_2, file = "overlap_cnv_2.summary", sep = "\t", quote = FALSE)
    print("Task done. Comparison results were saved in the working directory")
  }

  else {
    #convert coordinate for CNV and CNVR
    cnv_ars <- fread(cnvr_ars)
    cnv_ars$version <- "Verision_ARS" # add verision in dataframe
    colnames(cnv_ars) <- paste(colnames(cnv_ars), "ARS", sep = "_") #add suffix to all colnames
    ars_colnames <- colnames(cnv_ars) #set original column names use for extarting columns after matching

    cnv_umd <- fread(cnvr_umd)
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
    print("Starting to convert the coordinates for the first input file...")
    cnv_umd_ars <- merge(cnv_umd, two_map, by.x = c("Chr_UMD", "Start_UMD"), by.y = c("Chr_def_Map", "Position_def_Map"), all.x = TRUE)
    dup_index <- grep("TRUE", duplicated(cnv_umd_ars[, c("Chr_UMD", "Start_UMD", "CNVR_ID_UMD")]))
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
    dup_index <- grep("TRUE", duplicated(cnv_umd_ars[, c("Chr_UMD", "End_UMD", "CNVR_ID_UMD")]))
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
    print("Starting to convert the coordinates for the second input file...")
    cnv_ars_umd <- merge(cnv_ars, two_map, by.x = c("Chr_ARS", "Start_ARS"), by.y = c("Chr_tar_Map", "Position_tar_Map"), all.x = TRUE)
    dup_index <- grep("TRUE", duplicated(cnv_ars_umd[, c("Chr_ARS", "Start_ARS", "CNVR_ID_ARS")]))
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
    dup_index <- grep("TRUE", duplicated(cnv_ars_umd[, c("Chr_ARS", "End_ARS", "CNVR_ID_ARS")]))
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
    fwrite(cnv_umd_ars, file = "cnvr_UtoA.coord", sep = "\t", quote = FALSE)
    fwrite(cnv_ars_umd, file = "cnvr_AtoU.coord", sep = "\t", quote = FALSE)

    #two requrments for foverlap: 1) start <= end, 2) no NA in both start and end
    #cnv_umd_ars$End_ARS_Map[is.na(cnv_umd_ars$End_ARS_Map)] <- 0
    #cnv_umd_ars$Start_ARS_Map[is.na(cnv_umd_ars$Start_ARS_Map)] <- 0
    #After conversion of coordinate, we should check how many CNVRs got wrong position
    #The types of wrong position caused by the difference between the two version of map
    #1, NA, unkown position in start or end
    #2, 0 position in start or end
    #3, Start position lager than End Position
    #So we need to find out these CNVR with wrong position, then extarct the right CNVR for comparison
    wrong_umd <- subset(cnv_umd_ars, subset = (Start_ARS_Map == 0 | End_ARS_Map == 0 | is.na(Start_ARS_Map) | is.na(End_ARS_Map) | End_ARS_Map - Start_ARS_Map <= 0))
    right_umd <- subset(cnv_umd_ars, subset = !(Start_ARS_Map == 0 | End_ARS_Map == 0 | is.na(Start_ARS_Map) | is.na(End_ARS_Map) | End_ARS_Map - Start_ARS_Map <= 0))
    #right_umd <- setdiff(x = cnv_umd_ars, y = wrong_umd)
    #which(cnv_umd_ars$Start_ARS_Map == 0)
    #which(cnv_umd_ars$End_ARS_Map == 0)
    #grep("TRUE", is.na(cnv_umd_ars$Start_ARS_Map))
    #grep("TRUE", is.na(cnv_umd_ars$End_ARS_Map))
    #which((cnv_umd_ars$End_ARS_Map - cnv_umd_ars$Start_ARS_Map) <= 0)
    #right_umd <- cnv_umd_ars[-c(which((cnv_umd_ars$End_ARS_Map - cnv_umd_ars$Start_ARS_Map) <= 0)), ]
    #wrong_ars <- cnv_ars_umd[which((cnv_ars_umd$End_UMD_Map - cnv_ars_umd$Start_UMD_Map) <= 0), ]
    #wrong_ars <- wrong_ars[-c(which(wrong_ars$End_UMD_Map == 0)), ]
    #right_ars <- setdiff(cnv_ars_umd, wrong_ars)

    #######compare results in UMD at first-------------------------------------------------------------------------
    #setkey(cnv_umd_ars, Chr_UMD, Start_UMD, End_UMD)
    #find overlap on population level
    print("Starting to find overlapping CNVR...")
    setkey(cnv_ars_umd, Chr_ARS, Start_ARS, End_ARS)
    pop_overlap <- data.table::foverlaps(right_umd, cnv_ars_umd, by.x = c("Chr_UMD", "Start_ARS_Map", "End_ARS_Map"), type = "any", nomatch = NULL)
    pop_overlap$Overlap_length <- pmin(pop_overlap$End_ARS, pop_overlap$End_ARS_Map) - pmax(pop_overlap$Start_ARS, pop_overlap$Start_ARS_Map) + 1
    pop_overlap <- unique(pop_overlap)
    pop_overlap_umd <- subset(pop_overlap, select = umd_convert_colnames)
    final_pop_overlap_umd <- unique(pop_overlap_umd)
    non_overlap_pop_umd <- dplyr::setdiff(cnv_umd_ars, final_pop_overlap_umd)
    if (nrow(non_overlap_pop_umd) + nrow(final_pop_overlap_umd) == nrow(cnv_umd_ars)) {
      print(paste0("Comparison to the first file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVRs which are in total of ", nrow(cnv_umd_ars)))
    } else {print("Comparison to the first file was failed in validation, please use the original output files from HandyCNV as the input files")}

    fwrite(final_pop_overlap_umd, file = "overlap_cnvr_1.popu", sep = "\t", quote = FALSE)
    fwrite(non_overlap_pop_umd, file = "non_overlap_cnvr_1.popu", sep = "\t", quote = FALSE)

    #plot comparison
    final_pop_overlap_umd$Check_overlap <- "Overlap"
    non_overlap_pop_umd$Check_overlap <- "Non-Overlap"
    checkover_pop_1 <- rbind(final_pop_overlap_umd, non_overlap_pop_umd)
    colnames(checkover_pop_1) <- sub("_UMD", "", colnames(checkover_pop_1))
    checkover_pop_length <- merge(checkover_pop_1, pop_overlap[, c("CNVR_ID_UMD", "Overlap_length")], by.x = "CNVR_ID", by.y = "CNVR_ID_UMD", all.x = TRUE)
    #checkover_pop_length <- unique(checkover_pop_length)
    drop_name = "Overlap_length"
    checkover_pop_uniqe = unique(subset(checkover_pop_length, select = !(colnames(checkover_pop_length) %in% drop_name))) #used for calculate overlapping


    plot_comparison(cnv_checkover = checkover_pop_length, title_fig = "checkover_pop_1", width_1 = width_1, height_1 = height_1, hjust_prop = hjust_prop, hjust_num = hjust_num)

    fwrite(checkover_pop_length, file = "checkover_pop_1.txt", sep = "\t", quote = FALSE)

    #5.summarize difference
    #UMD
    overlap_percent_pop <- round(nrow(final_pop_overlap_umd) / nrow(cnv_umd_ars), 3) * 100
    non_overlap_percent_pop <- round(nrow(non_overlap_pop_umd) / nrow(cnv_umd_ars), 3) * 100
    print(paste0("The number of overlaped CNVRs on population level is ", nrow(final_pop_overlap_umd), ", which is around ", overlap_percent_pop, " percent in first file."))
    print(paste0("The number of Non-overlaped CNVRs on population level is ", nrow(non_overlap_pop_umd), ", which is around ", non_overlap_percent_pop, " percent in first file"))

    #overlap length
    overlap_length_1 <- sum(checkover_pop_length$Overlap_length, na.rm = TRUE)
    cnvr_length_1 <- sum(checkover_pop_uniqe$Length, na.rm = TRUE)
    overlap_length_prop <- round((overlap_length_1 / cnvr_length_1),3) * 100
    print(paste0("The overlapping length is ", overlap_length_1, " bp, which is around ", overlap_length_prop, " percent in first file"))

    overlap_summary <- data.frame("Item" = c("Overlapped CNVRs", "Non-overlapped CNVRs", "In Total"),
                                  "Number of CNVR" = c(nrow(final_pop_overlap_umd), nrow(non_overlap_pop_umd), nrow(final_pop_overlap_umd) + nrow(non_overlap_pop_umd)),
                                  "Proportion of Number (%)" = c(overlap_percent_pop, non_overlap_percent_pop, 100),
                                  "Length(bp)" = c(overlap_length_1, cnvr_length_1- overlap_length_1, cnvr_length_1),
                                  "Proportion of Length (%)" = c(overlap_length_prop, 100 - overlap_length_prop, 100))


    #overlap_summary <- data.frame("Item" = c("Overlapped CNVRs", "Non-overlapped CNVRs"),
    #                              "Population Level" = c(overlap_percent_pop, non_overlap_percent_pop))
    print("The final comparison results of the first file as follows: ")
    print(overlap_summary)
    fwrite(overlap_summary, file = "overlap_cnv.summary", sep = "\t", quote = FALSE)
    print("Comparison to the first file was finished.")


    ##########compare results in ARS-------------------------------------------------------------------
    #find overlap on population level
    setkey(right_umd, Chr_UMD, Start_ARS_Map, End_ARS_Map)
    pop_overlap_2 <- data.table::foverlaps(cnv_ars_umd, right_umd, by.x = c("Chr_ARS", "Start_ARS", "End_ARS"), type = "any", nomatch = NULL)
    pop_overlap_2$Overlap_length <- pmin(pop_overlap_2$End_ARS, pop_overlap_2$End_ARS_Map) - pmax(pop_overlap_2$Start_ARS, pop_overlap_2$Start_ARS_Map) + 1
    pop_overlap_ars <- subset(pop_overlap_2, select = ars_convert_colnames)
    final_pop_overlap_ars <- unique(pop_overlap_ars)
    non_overlap_pop_ars <- dplyr::setdiff(cnv_ars_umd, final_pop_overlap_ars)
    if (nrow(non_overlap_pop_ars) + nrow(final_pop_overlap_ars) == nrow(cnv_ars_umd)) {
      print(paste0("Comparison to the second file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVRs which are in total of ", nrow(cnv_ars_umd)))
    } else {print("Comparison to the second file was failed in validation, please use the original output files from HandyCNV as the input files")}

    fwrite(final_pop_overlap_ars, file = "overlap_cnvr_2.popu", sep = "\t", quote = FALSE)
    fwrite(non_overlap_pop_ars, file = "non_overlap_cnvr_2.popu", sep = "\t", quote = FALSE)


    final_pop_overlap_ars$Check_overlap <- "Overlap"
    non_overlap_pop_ars$Check_overlap <- "Non-Overlap"
    checkover_pop_2 <- rbind(final_pop_overlap_ars, non_overlap_pop_ars)
    colnames(checkover_pop_2) <- sub("_ARS", "", colnames(checkover_pop_2))
    checkover_pop_length_2 <- merge(checkover_pop_2, pop_overlap_2[, c("CNVR_ID_ARS", "Overlap_length")], by.x = "CNVR_ID", by.y = "CNVR_ID_ARS", all.x = TRUE)
    #checkover_pop_length_2 <- unique(checkover_pop_length_2)

    drop_name = "Overlap_length"
    checkover_pop_uniqe_2 = unique(subset(checkover_pop_length_2, select = !(colnames(checkover_pop_length_2) %in% drop_name))) #used for calculate overlapping


    plot_comparison(cnv_checkover = checkover_pop_length_2, title_fig = "checkover_pop_2", width_1 = width_1, height_1 = height_1, hjust_prop = hjust_prop, hjust_num = hjust_num)
    fwrite(checkover_pop_length_2, file = "checkover_pop_2.txt", sep = "\t", quote = FALSE)

    #5.summarize difference
    #UMD
    overlap_percent_pop_2 <- round(nrow(final_pop_overlap_ars) / nrow(cnv_ars_umd), 3) * 100
    non_overlap_percent_pop_2 <- round(nrow(non_overlap_pop_ars) / nrow(cnv_ars_umd), 3) * 100
    print(paste0("The number of overlaped CNVRs in seceond file on population level are ", nrow(final_pop_overlap_ars), ", which is around ", overlap_percent_pop_2, " percent"))
    print(paste0("The number of Non-overlaped CNVRs in second file on population level are ", nrow(non_overlap_pop_ars), ", which is around ", non_overlap_percent_pop_2, " percent"))

    #overlap length
    overlap_length_2 <- sum(checkover_pop_length_2$Overlap_length, na.rm = TRUE)
    cnvr_length_2 <- sum(checkover_pop_uniqe_2$Length, na.rm = TRUE)
    overlap_length_prop_2 <- round((overlap_length_2 / cnvr_length_2),3) * 100
    print(paste0("The overlapping length is ", overlap_length_2, " bp, which is around ", overlap_length_prop_2, " percent in first file"))

    overlap_summary_2 <- data.frame("Item" = c("Overlapped CNVRs", "Non-overlapped CNVRs", "In Total"),
                                    "Number of CNVR" = c(nrow(final_pop_overlap_ars), nrow(non_overlap_pop_ars), nrow(final_pop_overlap_ars) + nrow(non_overlap_pop_ars)),
                                    "Proportion of Number (%)" = c(overlap_percent_pop_2, non_overlap_percent_pop_2, 100),
                                    "Length(bp)" = c(overlap_length_2, cnvr_length_2- overlap_length_2, cnvr_length_2),
                                    "Proportion of Length (%)" = c(overlap_length_prop_2, 100 - overlap_length_prop_2, 100))
    #overlap_summary_2 <- data.frame("Item" = c("Overlapped CNVRs", "Non-overlapped CNVRs"),
    #                                "Population Level" = c(overlap_percent_pop_2, non_overlap_percent_pop_2))
    print("The final comparison results of the second file as follows: ")
    print(overlap_summary_2)
    fwrite(overlap_summary_2, file = "overlap_cnv_2.summary", sep = "\t", quote = FALSE)
    print("Task done. Comparison results were saved in the working directory")
  }
}

