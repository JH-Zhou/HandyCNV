#' Title compare_cnvr
#'#The idea to compare CNVR between different results are find out how many differences are there
#we defined 6 comparison standards of CNV
#1) overlapped
#1.same start and end, same SNP inside, fully overlap
#2.same start and end, different snp number, fully overlap
#3.different start or end, overlapped, partial overlap
#2) non-overlap
#4.missing start or end position
#5.End <= start
#6.different start or end, non-overlap
#according to the condition, the first thing is to match coordinates for both version
#then find overlap cnv and non overlap cnv
#then summarize how many CNVRs are in above standards
#'
#' @param cnvr_def first cnvr list, not limited on def or tar, default file is the result from call_cnvr function
#' @param cnvr_tar second cnvr list, not limited on def or tar, default file is the result from call_cnvr function
#' @param def_tar_map map file contains coordinates in both version of map. only need in comparison between the results from different versions. default file is generated from convert_map function
#' @param width_1 number to set the width of final plot size, unit is 'cm'
#' @param height_1 number to set the height of final plot size, unit is 'cm'
#' @param hjust_prop default value is 0.0. used to adjust horizontal position of the number of overlapped CNVR in the plot
#' @param hjust_num default value is 1.5. used to adjust horizontal position of the number of overlapped CNVR in the plot
#' @param col_1 set color for overlapped bar
#' @param col_2 set color for non-overlapped bar
#' @param folder set name of folder to save results
#' @import dplyr ggplot2
#' @importFrom data.table fread fwrite setkey foverlaps
#'
#' @return Details comparison results of CNVRs between input lists
#' @export compare_cnvr

compare_cnvr <- function(cnvr_def, cnvr_tar, def_tar_map = NULL, width_1 = 15, height_1 = 15, hjust_prop = 0.0, hjust_num = 1.5, folder = "compare_cnvr", col_1 = "gray20", col_2 = "springgreen4") {
  if(!file.exists(folder)){
    dir.create(folder)
    print(paste0("A new folder ", folder, " was created in working directory."))
  }

  #default plot function
  plot_comparison <- function(cnv_checkover, title_fig, width_1 = 15, height_1 = 15, hjust_prop = 0.0, hjust_num = 1.5) {
    cnv <- cnv_checkover
    title_f = title_fig
    drop_name = "Overlap_length"
    cnvr_cal = unique(subset(cnv, select = !(colnames(cnv) %in% drop_name))) %>%
               group_by(Type, Check_overlap) %>%
               summarise(origi_length = sum(Length),
                         num_CNVR = n_distinct(CNVR_ID))
    cnvr_over <- cnv %>%
                 group_by(Type, Check_overlap) %>%
                 summarise(overlap_len = sum(Overlap_length),
                           num_overlap_oppsite = n_distinct(Overlap_length, na.rm = TRUE))

    cnvr_over[1,3] = cnvr_cal[1,3] + cnvr_cal[2,3] - cnvr_over[2,3]
    cnvr_over[3,3] = cnvr_cal[3,3] + cnvr_cal[4,3] - cnvr_over[4,3]
    cnvr_over[5,3] = cnvr_cal[5,3] + cnvr_cal[6,3] - cnvr_over[6,3]

    cnvr_cal <- merge(cnvr_cal, cnvr_over, by = c("Type", "Check_overlap"), all.x = TRUE)

    cnvr_cal = cnvr_cal %>%
               group_by(Type) %>%
               mutate(prop_overlap_len = round(overlap_len / sum(origi_length),3),
                      prop_num = round(num_CNVR / sum(num_CNVR),3))
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

    #add manual color to btar
    color_bar <- c("Non-Overlap" = col_1,
                   "Overlap" = col_2)
    png(res = 300, filename = paste0(folder, "/", title_f, ".png"), width = width_1, height = height_1, bg = "transparent", units = "cm")
    compare_plot <- ggplot(cnvr_cal) +
      geom_col(aes(x = Type, y = overlap_len, fill = Check_overlap)) +
      scale_fill_manual(values = color_bar) +
      geom_text(data = subset(cnvr_cal, Check_overlap == "Overlap"), aes(x= Type, y = overlap_len,label = scales::percent(prop_overlap_len, accuracy = 0.1)), color = "orange", position = position_stack(0.5), hjust = hjust_prop) +
      geom_text(data = subset(cnvr_cal, Check_overlap == "Overlap"), aes(x= Type, y = overlap_len,label = num_CNVR), color = "black", position = position_stack(0.5), hjust = hjust_num) +
      #geom_text(aes(x = Type, y = percent_total, label = num), vjust = -0.5, hjust = 1.3) +
      #geom_text(inherit.aes = FALSE, data = cnv_freq, aes(x = Type, y = percent_total, label = scales::percent(percent_total)), vjust = -0.5, hjust  =-0.5) +
      scale_y_continuous(labels=scales::unit_format(unit = "", scale = 1e-6)) +
      theme_classic() +
      theme(legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.5,"line"),
            legend.text  = element_text(size = 6)) +
      labs(x = summary_title, y = "Total Length (Mb)")
    print(compare_plot)
    dev.off()
    fwrite(cnvr_cal, file = paste0(folder, "/", title_f, ".summary"), sep = "\t", quote = FALSE)
  }

  #dealing with data
  if (is.null(def_tar_map)) {
    cnv_tar <- fread(cnvr_tar)
    cnv_tar$version <- "Verision_TAR" # add verision in dataframe
    colnames(cnv_tar) <- paste(colnames(cnv_tar), "TAR", sep = "_") #add suffix to all colnames
    tar_colnames <- colnames(cnv_tar) #set original column names use for extarting columns after matching

    cnv_def <- fread(cnvr_def)
    cnv_def$version <- "Version_DEF"
    colnames(cnv_def) <- paste(colnames(cnv_def), "DEF", sep = "_")
    def_colnames <- colnames(cnv_def)

    #find overlap on population level
    setkey(cnv_tar, Chr_TAR, Start_TAR, End_TAR)
    pop_overlap <- data.table::foverlaps(cnv_def, cnv_tar, by.x = c("Chr_DEF", "Start_DEF", "End_DEF"), type = "any", nomatch = NULL)
    pop_overlap$Overlap_length <- pmin(pop_overlap$End_TAR, pop_overlap$End_DEF) - pmax(pop_overlap$Start_TAR, pop_overlap$Start_DEF) + 1
    pop_overlap_def <- subset(pop_overlap, select = def_colnames)
    final_pop_overlap_def <- unique(pop_overlap_def)
    non_overlap_pop_def <- dplyr::setdiff(cnv_def, final_pop_overlap_def)
    if (nrow(non_overlap_pop_def) + nrow(final_pop_overlap_def) == nrow(cnv_def)) {
      print(paste0("Comparison to the first file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVRs which are in total of ", nrow(cnv_def)))
    } else {print("Comparison to the first file was failed in validation, please use the original output files from HandyCNV as the input files")}

    fwrite(final_pop_overlap_def, file = paste0(folder, "/overlap_cnv_1.popu"), sep = "\t", quote = FALSE)
    fwrite(non_overlap_pop_def, file = paste0(folder, "/non_overlap_cnv_1.popu"), sep = "\t", quote = FALSE)

    #make comparison plot
    final_pop_overlap_def$Check_overlap <- "Overlap"
    non_overlap_pop_def$Check_overlap <- "Non-Overlap"
    checkover_pop_1 <- rbind(final_pop_overlap_def, non_overlap_pop_def)
    colnames(checkover_pop_1) <- sub("_DEF", "", colnames(checkover_pop_1))
    checkover_pop_length <- merge(checkover_pop_1, pop_overlap[, c("CNVR_ID_DEF", "Overlap_length")], by.x = "CNVR_ID", by.y = "CNVR_ID_DEF", all.x = TRUE)
    #checkover_pop_length <- unique(checkover_pop_length)
    drop_name = "Overlap_length"
    checkover_pop_uniqe_1 = unique(subset(checkover_pop_length, select = !(colnames(checkover_pop_length) %in% drop_name))) #used for calculate overlapping


    try(plot_comparison(cnv_checkover = checkover_pop_length, title_fig = "checkover_pop_1", width_1 = width_1, height_1 = height_1, hjust_prop = hjust_prop, hjust_num = hjust_num), silent = TRUE)

    fwrite(checkover_pop_length, file = paste0(folder, "/checkover_pop_1.txt"), sep = "\t", quote = FALSE)

    #5.summarize difference
    #DEF overlap number
    overlap_percent_pop <- round(nrow(final_pop_overlap_def) / nrow(cnv_def), 3) * 100
    non_overlap_percent_pop <- round(nrow(non_overlap_pop_def) / nrow(cnv_def), 3) * 100
    print(paste0("The number of overlaped CNVRs on population level is ", nrow(final_pop_overlap_def), ", which is around ", overlap_percent_pop, " percent in first file."))
    print(paste0("The number of Non-overlaped CNVRs on population level is ", nrow(non_overlap_pop_def), ", which is around ", non_overlap_percent_pop, " percent in first file"))

    #overlap length
    overlap_length_1 <- sum(pop_overlap$Overlap_length, na.rm = TRUE)
    #overlap_length_1 <- sum(checkover_pop_length$Overlap_length, na.rm = TRUE)
    cnvr_length_1 <- sum(checkover_pop_uniqe_1$Length, na.rm = TRUE)
    overlap_length_prop <- round((overlap_length_1 / cnvr_length_1),3) * 100
    print(paste0("The overlapping length is ", overlap_length_1, " bp, which is around ", overlap_length_prop, " percent in first file"))

    overlap_summary <- data.frame("Item" = c("Overlapped CNVRs", "Non-overlapped CNVRs", "In Total"),
                                  "Number of CNVR" = c(nrow(final_pop_overlap_def), nrow(non_overlap_pop_def), nrow(final_pop_overlap_def) + nrow(non_overlap_pop_def)),
                                  "Proportion of Number (%)" = c(overlap_percent_pop, non_overlap_percent_pop, 100),
                                  "Length(bp)" = c(overlap_length_1, cnvr_length_1- overlap_length_1, cnvr_length_1),
                                  "Proportion of Length (%)" = c(overlap_length_prop, 100 - overlap_length_prop, 100))
    print("The final comparison results of the first file as follows: ")
    print(overlap_summary)
    fwrite(overlap_summary, file = paste0(folder, "/overlap_cnv_1.summary"), sep = "\t", quote = FALSE)
    print("Comparison to the first file was finished.")


  ##########compare results in TAR at second-------------------------------------------------------------------
   #setkey(cnv_def_tar, Chr_DEF, Start_DEF, End_DEF)
    #find overlap on population level
    setkey(cnv_def, Chr_DEF, Start_DEF, End_DEF)
    pop_overlap_2 <- data.table::foverlaps(cnv_tar, cnv_def, by.x = c("Chr_TAR", "Start_TAR", "End_TAR"), type = "any", nomatch = NULL)
    pop_overlap_2$Overlap_length <- pmin(pop_overlap_2$End_TAR, pop_overlap_2$End_DEF) - pmax(pop_overlap_2$Start_TAR, pop_overlap_2$Start_DEF) + 1
    pop_overlap_tar <- subset(pop_overlap_2, select = tar_colnames)
    final_pop_overlap_tar <- unique(pop_overlap_tar)
    non_overlap_pop_tar <- dplyr::setdiff(cnv_tar, final_pop_overlap_tar)
    if (nrow(non_overlap_pop_tar) + nrow(final_pop_overlap_tar) == nrow(cnv_tar)) {
      print(paste0("Comparison to the second file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVRs which are in total of ", nrow(cnv_tar)))
    } else {print("Comparison to the second file was failed in validation, please use the original output files from HandyCNV as the input files")}

    fwrite(final_pop_overlap_tar, file = paste0(folder, "/overlap_cnv_tar.popu"), sep = "\t", quote = FALSE)
    fwrite(non_overlap_pop_tar, file = paste0(folder, "/non_overlap_cnv_tar.popu"), sep = "\t", quote = FALSE)

    final_pop_overlap_tar$Check_overlap <- "Overlap"
    non_overlap_pop_tar$Check_overlap <- "Non-Overlap"
    checkover_pop_2 <- rbind(final_pop_overlap_tar, non_overlap_pop_tar)
    colnames(checkover_pop_2) <- sub("_TAR", "", colnames(checkover_pop_2))
    checkover_pop_length_2 <- merge(checkover_pop_2, pop_overlap_2[, c("CNVR_ID_TAR", "Overlap_length")], by.x = "CNVR_ID", by.y = "CNVR_ID_TAR", all.x = TRUE)
    #checkover_pop_length_2 <- unique(checkover_pop_length_2)
    drop_name = "Overlap_length"
    checkover_pop_uniqe_2 = unique(subset(checkover_pop_length_2, select = !(colnames(checkover_pop_length_2) %in% drop_name))) #used for calculate overlapping


    try(plot_comparison(cnv_checkover = checkover_pop_length_2, title_fig = "checkover_pop_2", width_1 = width_1, height_1 = height_1, hjust_prop = hjust_prop, hjust_num = hjust_num), silent = TRUE)
    #plot_comparison(cnv_checkover = checkover_pop_length_2, title_fig = "checkover_pop_2", width_1 = width_1, height_1 = height_1, hjust_prop = hjust_prop, hjust_num = hjust_num)
    fwrite(checkover_pop_length_2, file = paste0(folder, "/checkover_pop_2.txt"), sep = "\t", quote = FALSE)

    #5.summarize difference
    #TAR
    overlap_percent_pop_2 <- round(nrow(final_pop_overlap_tar) / nrow(cnv_tar), 3) * 100
    non_overlap_percent_pop_2 <- round(nrow(non_overlap_pop_tar) / nrow(cnv_tar), 3) * 100
    print(paste0("The number of overlaped CNVRs in seceond file on population level are ", nrow(final_pop_overlap_tar), ", which is around ", overlap_percent_pop_2, " percent"))
    print(paste0("The number of Non-overlaped CNVRs in second file on population level are ", nrow(non_overlap_pop_tar), ", which is around ", non_overlap_percent_pop_2, " percent"))

    #overlap length
    overlap_length_2 <- sum(pop_overlap_2$Overlap_length, na.rm = TRUE)
    #overlap_length_2 <- sum(checkover_pop_length_2$Overlap_length, na.rm = TRUE)
    cnvr_length_2 <- sum(checkover_pop_uniqe_2$Length, na.rm = TRUE)
    overlap_length_prop_2 <- round((overlap_length_2 / cnvr_length_2),3) * 100
    print(paste0("The overlapping length is ", overlap_length_2, " bp, which is around ", overlap_length_prop_2, " percent in first file"))

    overlap_summary_2 <- data.frame("Item" = c("Overlapped CNVRs", "Non-overlapped CNVRs", "In Total"),
                                    "Number of CNVR" = c(nrow(final_pop_overlap_tar), nrow(non_overlap_pop_tar), nrow(final_pop_overlap_tar) + nrow(non_overlap_pop_tar)),
                                    "Proportion of Number (%)" = c(overlap_percent_pop_2, non_overlap_percent_pop_2, 100),
                                    "Length(bp)" = c(overlap_length_2, cnvr_length_2- overlap_length_2, cnvr_length_2),
                                    "Proportion of Length (%)" = c(overlap_length_prop_2, 100 - overlap_length_prop_2, 100))
    print("The final comparison results of the second file as follows: ")
    print(overlap_summary_2)
    fwrite(overlap_summary_2, file = paste0(folder, "/overlap_cnv_2.summary"), sep = "\t", quote = FALSE)
    print("Task done. Comparison results were saved in the working directory")
  }

  # else {
  #   #convert coordinate for CNV and CNVR
  #   cnv_tar <- fread(cnvr_tar)
  #   cnv_tar$version <- "Verision_TAR" # add version in dataframe
  #   colnames(cnv_tar) <- paste(colnames(cnv_tar), "TAR", sep = "_") #add suffix to all colnames
  #   tar_colnames <- colnames(cnv_tar) #set original column names use for extracting columns after matching
  #
  #   cnv_def <- fread(cnvr_def)
  #   cnv_def$version <- "Version_DEF"
  #   colnames(cnv_def) <- paste(colnames(cnv_def), "DEF", sep = "_")
  #   def_colnames <- colnames(cnv_def)
  #
  #   two_map <- fread(def_tar_map)
  #   colnames(two_map) <- paste(colnames(two_map), "Map", sep = "_")
  #
  #   #cnv_def_tar <- merge(cnv_def, two_map, by.x = c("Chr", "Start"), by.y = c("Chr_DEF", "Position_DEF"), all.x = TRUE)
  #
  #   #1. convert the def result to tar
  #   #matching the start position by snp name, validating by matching the original position and new mathced position
  #   #cnv_def_tar <- merge(cnv_def, two_map, by.x = "Start_SNP_DEF", by.y = "Name_Map", all.x = TRUE)
  #
  #   #There are some duplicated rows after merge progress, because of there are some different SNP with same location in the map file
  #   #To solve this problem, we should check duplicated row after each merge step and remove the duplicates
  #   print("Starting to convert the coordinates for the first input file...")
  #   cnv_def_tar <- merge(cnv_def, two_map, by.x = c("Chr_DEF", "Start_DEF"), by.y = c("Chr_def_Map", "Position_def_Map"), all.x = TRUE)
  #   dup_index <- grep("TRUE", duplicated(cnv_def_tar[, c("Chr_DEF", "Start_DEF", "CNVR_ID_DEF")]))
  #   if (length(dup_index) > 0){
  #     cnv_def_tar <- cnv_def_tar[-dup_index, ]
  #     print(paste0("There are ", length(dup_index), " duplicated SNPs after converting the Start Position for the first file, they were all successfully removed."))
  #   } else {
  #     print("Non duplicated SNPs were detected after converting the Start Position for the first file.")
  #   }
  #
  #   if (all(cnv_def_tar$Start_DEF == cnv_def_tar$Position_def_Map)){
  #     print("Matching results by the Start position passed the validation.")
  #   } else {print("matching results didn't pass the validation, please check your input data, make sure the input file all generate by HandyCNV")}
  #
  #   cnv_def_tar <- subset(cnv_def_tar, select = c(def_colnames, "Position_tar_Map"))
  #   names(cnv_def_tar)[names(cnv_def_tar) == "Position_tar_Map"] <- "Start_TAR_Map"
  #   def_temp_colnames <- colnames(cnv_def_tar)
  #
  #   #matching the end position of CNV
  #   #cnv_def_tar <- merge(cnv_def_tar, two_map, by.x = "End_SNP_DEF", by.y = "Name_Map", all.x = TRUE)
  #   cnv_def_tar <- merge(cnv_def_tar, two_map, by.x = c("Chr_DEF", "End_DEF"), by.y = c("Chr_def_Map", "Position_def_Map"), all.x = TRUE)
  #   dup_index <- grep("TRUE", duplicated(cnv_def_tar[, c("Chr_DEF", "End_DEF", "CNVR_ID_DEF")]))
  #   if (length(dup_index) > 0){
  #     cnv_def_tar <- cnv_def_tar[-dup_index, ]
  #     print(paste0("There are ", length(dup_index), " duplicated SNPs after converting the End Position for the first file, they were all successfully removed."))
  #   } else {
  #     print("Non duplicated SNPs were detected after converting the End Position for the first file.")
  #   }
  #
  #   if (all(cnv_def_tar$End_DEF == cnv_def_tar$Position_def_Map)){
  #     print("Matching results by the End position passed the validation.")
  #   } else {print("matching results didn't pass the validation, please check your input data, make sure the input file all generate by HandyCNV")}
  #
  #   cnv_def_tar <- subset(cnv_def_tar, select = c(def_temp_colnames, "Position_tar_Map"))
  #   names(cnv_def_tar)[names(cnv_def_tar) == "Position_tar_Map"] <- "End_TAR_Map"
  #   def_convert_colnames <- colnames(cnv_def_tar)
  #
  #
  #   #2.convert tar result to def
  #   #cnv_tar_def <- merge(cnv_tar, two_map, by.x = "Start_SNP_TAR", by.y = "Name_Map", all.x = TRUE)
  #   print("Starting to convert the coordinates for the second input file...")
  #   cnv_tar_def <- merge(cnv_tar, two_map, by.x = c("Chr_TAR", "Start_TAR"), by.y = c("Chr_tar_Map", "Position_tar_Map"), all.x = TRUE)
  #   dup_index <- grep("TRUE", duplicated(cnv_tar_def[, c("Chr_TAR", "Start_TAR", "CNVR_ID_TAR")]))
  #   if (length(dup_index) > 0){
  #     cnv_tar_def <- cnv_tar_def[-dup_index, ]
  #     print(paste0("There are ", length(dup_index), " duplicated SNPs after converting the Start Position for the Second file, they were all successfully removed."))
  #   } else {
  #     print("Non duplicated SNPs were detected after converting the Start Position for the Second file.")
  #   }
  #
  #   #checking if the matching results correct?
  #   if (all(cnv_tar_def$Start_TAR == cnv_tar_def$Position_tar_Map)){
  #     print("Matching results by the Start position passed the validation.")
  #   } else {print("matching results didn't pass the validation, please check your input data, make sure the input file all generate by HandyCNV")}
  #   cnv_tar_def <- subset(cnv_tar_def, select = c(tar_colnames, "Position_def_Map"))
  #   names(cnv_tar_def)[names(cnv_tar_def) == "Position_def_Map"] <- "Start_DEF_Map"
  #   tar_temp_colnames <- colnames(cnv_tar_def)
  #
  #   #matching end position
  #   #cnv_tar_def <- merge(cnv_tar_def, two_map, by.x = "End_SNP_TAR", by.y = "Name_Map", all.x = TRUE)
  #   cnv_tar_def <- merge(cnv_tar_def, two_map, by.x = c("Chr_TAR", "End_TAR"), by.y = c("Chr_tar_Map", "Position_tar_Map"), all.x = TRUE)
  #   dup_index <- grep("TRUE", duplicated(cnv_tar_def[, c("Chr_TAR", "End_TAR", "CNVR_ID_TAR")]))
  #   if (length(dup_index) > 0){
  #     cnv_tar_def <- cnv_tar_def[-dup_index, ]
  #     print(paste0("There are ", length(dup_index), " duplicated SNPs after converting the End Postion for the Second File, they were all successfully removed."))
  #   } else {
  #     print("Non duplicated SNPs were detected after converting the End Postion for the Second File.")
  #   }
  #
  #   if (all(cnv_tar_def$End_TAR == cnv_tar_def$Position_tar_Map)){
  #     print("Matching results by the End position passed the validation.")
  #   } else {print("matching results didn't pass the validation, please check your input data, make sure the input file all generate by HandyCNV")}
  #   cnv_tar_def <- subset(cnv_tar_def, select = c(tar_temp_colnames, "Position_def_Map"))
  #   names(cnv_tar_def)[names(cnv_tar_def) == "Position_def_Map"] <- "End_DEF_Map"
  #   tar_convert_colnames <- colnames(cnv_tar_def)
  #
  #   #write out CNV with converted coordinates
  #   fwrite(cnv_def_tar, file = paste0(folder, "/cnvr_UtoA.coord"), sep = "\t", quote = FALSE)
  #   fwrite(cnv_tar_def, file = paste0(folder, "/cnvr_AtoU.coord"), sep = "\t", quote = FALSE)
  #
  #   #two requirements for foverlap: 1) start <= end, 2) no NA in both start and end
  #   #cnv_def_tar$End_TAR_Map[is.na(cnv_def_tar$End_TAR_Map)] <- 0
  #   #cnv_def_tar$Start_TAR_Map[is.na(cnv_def_tar$Start_TAR_Map)] <- 0
  #   #After conversion of coordinate, we should check how many CNVRs got wrong position
  #   #The types of wrong position caused by the difference between the two version of map
  #   #1, NA, unknown position in start or end
  #   #2, 0 position in start or end
  #   #3, Start position lager than End Position
  #   #So we need to find out these CNVR with wrong position, then extarct the right CNVR for comparison
  #   wrong_def <- subset(cnv_def_tar, subset = (Start_TAR_Map == 0 | End_TAR_Map == 0 | is.na(Start_TAR_Map) | is.na(End_TAR_Map) | End_TAR_Map - Start_TAR_Map <= 0))
  #   right_def <- subset(cnv_def_tar, subset = !(Start_TAR_Map == 0 | End_TAR_Map == 0 | is.na(Start_TAR_Map) | is.na(End_TAR_Map) | End_TAR_Map - Start_TAR_Map <= 0))
  #   print(paste0("There are ", nrow(right_def), " CNVRs successful coordinates coversion in Input list 1"))
  #   print(paste0("There are ", nrow(wrong_def), " CNVRs failed to coordinates coversion in Input list 1"))
  #   #right_def <- setdiff(x = cnv_def_tar, y = wrong_def)
  #   #which(cnv_def_tar$Start_TAR_Map == 0)
  #   #which(cnv_def_tar$End_TAR_Map == 0)
  #   #grep("TRUE", is.na(cnv_def_tar$Start_TAR_Map))
  #   #grep("TRUE", is.na(cnv_def_tar$End_TAR_Map))
  #   #which((cnv_def_tar$End_TAR_Map - cnv_def_tar$Start_TAR_Map) <= 0)
  #   #right_def <- cnv_def_tar[-c(which((cnv_def_tar$End_TAR_Map - cnv_def_tar$Start_TAR_Map) <= 0)), ]
  #   #wrong_tar <- cnv_tar_def[which((cnv_tar_def$End_DEF_Map - cnv_tar_def$Start_DEF_Map) <= 0), ]
  #   #wrong_tar <- wrong_tar[-c(which(wrong_tar$End_DEF_Map == 0)), ]
  #   #right_tar <- setdiff(cnv_tar_def, wrong_tar)
  #
  #   #######compare results in DEF at first-------------------------------------------------------------------------
  #   #setkey(cnv_def_tar, Chr_DEF, Start_DEF, End_DEF)
  #   #find overlap on population level
  #   print("Starting to find overlapping CNVR...")
  #   setkey(cnv_tar_def, Chr_TAR, Start_TAR, End_TAR)
  #   pop_overlap <- data.table::foverlaps(right_def, cnv_tar_def, by.x = c("Chr_DEF", "Start_TAR_Map", "End_TAR_Map"), type = "any", nomatch = NULL)
  #   pop_overlap$Overlap_length <- pmin(pop_overlap$End_TAR, pop_overlap$End_TAR_Map) - pmax(pop_overlap$Start_TAR, pop_overlap$Start_TAR_Map) + 1
  #   pop_overlap <- unique(pop_overlap)
  #   pop_overlap_def <- subset(pop_overlap, select = def_convert_colnames)
  #   final_pop_overlap_def <- unique(pop_overlap_def)
  #   non_overlap_pop_def <- dplyr::setdiff(cnv_def_tar, final_pop_overlap_def)
  #   if (nrow(non_overlap_pop_def) + nrow(final_pop_overlap_def) == nrow(cnv_def_tar)) {
  #     print(paste0("Comparison to the first file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVRs which are in total of ", nrow(cnv_def_tar)))
  #   } else {print("Comparison to the first file was failed in validation, please use the original output files from HandyCNV as the input files")}
  #
  #   fwrite(final_pop_overlap_def, file = paste0(folder, "/overlap_cnvr_1.popu"), sep = "\t", quote = FALSE)
  #   fwrite(non_overlap_pop_def, file = paste0(folder, "/non_overlap_cnvr_1.popu"), sep = "\t", quote = FALSE)
  #
  #   #plot comparison
  #   final_pop_overlap_def$Check_overlap <- "Overlap"
  #   non_overlap_pop_def$Check_overlap <- "Non-Overlap"
  #   checkover_pop_1 <- rbind(final_pop_overlap_def, non_overlap_pop_def)
  #   colnames(checkover_pop_1) <- sub("_DEF", "", colnames(checkover_pop_1))
  #   checkover_pop_length <- merge(checkover_pop_1, pop_overlap[, c("CNVR_ID_DEF", "Overlap_length")], by.x = "CNVR_ID", by.y = "CNVR_ID_DEF", all.x = TRUE)
  #   #checkover_pop_length <- unique(checkover_pop_length)
  #   drop_name = "Overlap_length"
  #   checkover_pop_uniqe = unique(subset(checkover_pop_length, select = !(colnames(checkover_pop_length) %in% drop_name))) #used for calculate overlapping
  #
  #
  #   plot_comparison(cnv_checkover = checkover_pop_length, title_fig = "checkover_pop_1", width_1 = width_1, height_1 = height_1, hjust_prop = hjust_prop, hjust_num = hjust_num)
  #
  #   fwrite(checkover_pop_length, file = paste0(folder, "/checkover_pop_1.txt"), sep = "\t", quote = FALSE)
  #
  #   #5.summarize difference
  #   #DEF
  #   overlap_percent_pop <- round(nrow(final_pop_overlap_def) / nrow(cnv_def_tar), 3) * 100
  #   non_overlap_percent_pop <- round(nrow(non_overlap_pop_def) / nrow(cnv_def_tar), 3) * 100
  #   print(paste0("The number of overlaped CNVRs on population level is ", nrow(final_pop_overlap_def), ", which is around ", overlap_percent_pop, " percent in first file."))
  #   print(paste0("The number of Non-overlaped CNVRs on population level is ", nrow(non_overlap_pop_def), ", which is around ", non_overlap_percent_pop, " percent in first file"))
  #
  #   #overlap length
  #   overlap_length_1 <- sum(pop_overlap$Overlap_length, na.rm = TRUE)
  #   #overlap_length_1 <- sum(checkover_pop_length$Overlap_length, na.rm = TRUE)
  #   cnvr_length_1 <- sum(checkover_pop_uniqe$Length, na.rm = TRUE)
  #   overlap_length_prop <- round((overlap_length_1 / cnvr_length_1),3) * 100
  #   print(paste0("The overlapping length is ", overlap_length_1, " bp, which is around ", overlap_length_prop, " percent in first file"))
  #
  #   overlap_summary <- data.frame("Item" = c("Overlapped CNVRs", "Non-overlapped CNVRs", "In Total"),
  #                                 "Number of CNVR" = c(nrow(final_pop_overlap_def), nrow(non_overlap_pop_def), nrow(final_pop_overlap_def) + nrow(non_overlap_pop_def)),
  #                                 "Proportion of Number (%)" = c(overlap_percent_pop, non_overlap_percent_pop, 100),
  #                                 "Length(bp)" = c(overlap_length_1, cnvr_length_1- overlap_length_1, cnvr_length_1),
  #                                 "Proportion of Length (%)" = c(overlap_length_prop, 100 - overlap_length_prop, 100))
  #
  #
  #   #overlap_summary <- data.frame("Item" = c("Overlapped CNVRs", "Non-overlapped CNVRs"),
  #   #                              "Population Level" = c(overlap_percent_pop, non_overlap_percent_pop))
  #   print("The final comparison results of the first file as follows: ")
  #   print(overlap_summary)
  #   fwrite(overlap_summary, file = paste0(folder, "/overlap_cnv.summary"), sep = "\t", quote = FALSE)
  #   print("Comparison to the first file was finished.")
  #
  #
  #   ##########compare results in TAR-------------------------------------------------------------------
  #   #find overlap on population level
  #   setkey(right_def, Chr_DEF, Start_TAR_Map, End_TAR_Map)
  #   pop_overlap_2 <- data.table::foverlaps(cnv_tar_def, right_def, by.x = c("Chr_TAR", "Start_TAR", "End_TAR"), type = "any", nomatch = NULL)
  #   pop_overlap_2$Overlap_length <- pmin(pop_overlap_2$End_TAR, pop_overlap_2$End_TAR_Map) - pmax(pop_overlap_2$Start_TAR, pop_overlap_2$Start_TAR_Map) + 1
  #   pop_overlap_tar <- subset(pop_overlap_2, select = tar_convert_colnames)
  #   final_pop_overlap_tar <- unique(pop_overlap_tar)
  #   non_overlap_pop_tar <- dplyr::setdiff(cnv_tar_def, final_pop_overlap_tar)
  #   if (nrow(non_overlap_pop_tar) + nrow(final_pop_overlap_tar) == nrow(cnv_tar_def)) {
  #     print(paste0("Comparison to the second file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVRs which are in total of ", nrow(cnv_tar_def)))
  #   } else {print("Comparison to the second file was failed in validation, please use the original output files from HandyCNV as the input files")}
  #
  #   fwrite(final_pop_overlap_tar, file = paste0(folder, "/overlap_cnvr_2.popu"), sep = "\t", quote = FALSE)
  #   fwrite(non_overlap_pop_tar, file = paste0(folder, "/non_overlap_cnvr_2.popu"), sep = "\t", quote = FALSE)
  #
  #
  #   final_pop_overlap_tar$Check_overlap <- "Overlap"
  #   non_overlap_pop_tar$Check_overlap <- "Non-Overlap"
  #   checkover_pop_2 <- rbind(final_pop_overlap_tar, non_overlap_pop_tar)
  #   colnames(checkover_pop_2) <- sub("_TAR", "", colnames(checkover_pop_2))
  #   checkover_pop_length_2 <- merge(checkover_pop_2, pop_overlap_2[, c("CNVR_ID_TAR", "Overlap_length")], by.x = "CNVR_ID", by.y = "CNVR_ID_TAR", all.x = TRUE)
  #   #checkover_pop_length_2 <- unique(checkover_pop_length_2)
  #
  #   drop_name = "Overlap_length"
  #   checkover_pop_uniqe_2 = unique(subset(checkover_pop_length_2, select = !(colnames(checkover_pop_length_2) %in% drop_name))) #used for calculate overlapping
  #
  #
  #   plot_comparison(cnv_checkover = checkover_pop_length_2, title_fig = "checkover_pop_2", width_1 = width_1, height_1 = height_1, hjust_prop = hjust_prop, hjust_num = hjust_num)
  #   fwrite(checkover_pop_length_2, file = paste0(folder, "/checkover_pop_2.txt"), sep = "\t", quote = FALSE)
  #
  #   #5.summarize difference
  #   #DEF
  #   overlap_percent_pop_2 <- round(nrow(final_pop_overlap_tar) / nrow(cnv_tar_def), 3) * 100
  #   non_overlap_percent_pop_2 <- round(nrow(non_overlap_pop_tar) / nrow(cnv_tar_def), 3) * 100
  #   print(paste0("The number of overlaped CNVRs in seceond file on population level are ", nrow(final_pop_overlap_tar), ", which is around ", overlap_percent_pop_2, " percent"))
  #   print(paste0("The number of Non-overlaped CNVRs in second file on population level are ", nrow(non_overlap_pop_tar), ", which is around ", non_overlap_percent_pop_2, " percent"))
  #
  #   #overlap length
  #   overlap_length_2 <- sum(pop_overlap_2$Overlap_length, na.rm = TRUE)
  #   #overlap_length_2 <- sum(checkover_pop_length_2$Overlap_length, na.rm = TRUE)
  #   cnvr_length_2 <- sum(checkover_pop_uniqe_2$Length, na.rm = TRUE)
  #   overlap_length_prop_2 <- round((overlap_length_2 / cnvr_length_2),3) * 100
  #   print(paste0("The overlapping length is ", overlap_length_2, " bp, which is around ", overlap_length_prop_2, " percent in first file"))
  #
  #   overlap_summary_2 <- data.frame("Item" = c("Overlapped CNVRs", "Non-overlapped CNVRs", "In Total"),
  #                                   "Number of CNVR" = c(nrow(final_pop_overlap_tar), nrow(non_overlap_pop_tar), nrow(final_pop_overlap_tar) + nrow(non_overlap_pop_tar)),
  #                                   "Proportion of Number (%)" = c(overlap_percent_pop_2, non_overlap_percent_pop_2, 100),
  #                                   "Length(bp)" = c(overlap_length_2, cnvr_length_2- overlap_length_2, cnvr_length_2),
  #                                   "Proportion of Length (%)" = c(overlap_length_prop_2, 100 - overlap_length_prop_2, 100))
  #   #overlap_summary_2 <- data.frame("Item" = c("Overlapped CNVRs", "Non-overlapped CNVRs"),
  #   #                                "Population Level" = c(overlap_percent_pop_2, non_overlap_percent_pop_2))
  #   print("The final comparison results of the second file as follows: ")
  #   print(overlap_summary_2)
  #   fwrite(overlap_summary_2, file = paste0(folder, "/overlap_cnv_2.summary"), sep = "\t", quote = FALSE)
  #   print("Task done. Comparison results were saved in the working directory")
  # }
}

