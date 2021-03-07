#' Title compare_cnv
#'The idea to compare CNVs between different results are find out how many different are there
#' we defined 6 comparison standards of CNV
#' 1) overlapped
#' 1.same start and end, same SNP inside, fully overlap
#' 2.same start and end, different snp number, fully overlap
#' 3.different start or end, overlapped, partial overlap
#' 2) non-overlap
#' 4.missing start or end position
#' 5.End <= start
#' 6.different start or end, non-overlap
#' according to the condition, the first thing is to match coordinates for both version
#' then find overlap cnv and non overlap cnv
#' then summarize how many CNVs are in above standards
#'
#' @param cnv_def first cnv list, defaults file is the result from clean_cnv function
#' @param cnv_tar second cnv list, defaults file is the result from clean_cnv function
#' @param def_tar_map map file contains coordinates in both version of map. only need in comparison between the results from different versions. default file is generated from convert_map function
#' @param width_1 integer, default value is 14, number to set the width of final plot size, unit is 'cm'
#' @param height_1 integer, default value is 14,number to set the height of final plot size, unit is 'cm'
#' @param legend_x decimal digit, default value is 0.9, consistent with ggplot manual coordinates of legend
#' @param legend_y decimal digit, default value is 0.9, consistent with ggplot manual coordinates of legend
#' @param col_1 set color for overlapped bar
#' @param col_2 set color for non-overlapped bar
#' @param folder set name of folder to save results
#' @param plot_caption If TRUE, present Note Caption in Comparison plot
#'
#' @import dplyr scales ggplot2
#'
#' @importFrom data.table fread fwrite setkey foverlaps
#'
#' @return Details comparison results of CNVs between input lists.
#' @export compare_cnv
#'
compare_cnv <- function(cnv_def, cnv_tar, def_tar_map = NULL, width_1 = 14, height_1 = 11, legend_x = 0.9, legend_y = 0.9, folder = "compare_cnv", col_1 = "pink", col_2 = "lightblue", plot_caption = TRUE) {
  if(!file.exists(folder)){
    dir.create(folder)
    print(paste0("A new folder ", folder, " was created in working directory."))
  }

  #default plot function
  plot_comparison <- function(cnv_checkover, title_fig, width_1 =14, height_1 =11, legend_x = 0.9, legend_y = 0.9, plot_caption = TRUE) {
    cnv <- cnv_checkover
    title_f = title_fig
    cnv_cal <- cnv %>% count(CNV_Value, Check_overlap) %>%
      mutate(percent_total = round(n / sum(n), 4)) %>%
      group_by(CNV_Value) %>%
      mutate(percent_group = round(n /sum(n), 4))

    cnv_freq <- cnv_cal %>%
                group_by(CNV_Value) %>%
                summarise(percent_total = sum(percent_total),
                          num = sum(n))

    # manual add color
    color_bar <- c("Non-Overlap" = col_1,
                   "Overlap" = col_2)
    png(res = 300, filename = paste0(folder, "/", title_f, ".png"), width = width_1, height = height_1, bg = "transparent", units = "cm")
    compare_plot <- ggplot(cnv_cal, aes(x = CNV_Value, y = percent_total, fill = Check_overlap)) +
      geom_col() +
      scale_fill_manual(values = color_bar) +
      geom_text(data = subset(cnv_cal, Check_overlap == "Overlap"), aes(label = scales::percent(percent_group)), color = "blue", position = position_stack(0.5)) +
      geom_text(inherit.aes = FALSE, data = cnv_freq, aes(x = CNV_Value - 0.4, y = percent_total, label = num), vjust = -0.5) +
      geom_text(inherit.aes = FALSE, data = cnv_freq, aes(x = CNV_Value + 0.2, y = percent_total, label = scales::percent(percent_total,accuracy = 0.1)), vjust = -0.5, color = "red") +
      scale_y_continuous(labels=scales::percent) +
      theme_classic() +
      theme(legend.position = c(legend_x, legend_y),
            legend.title = element_blank(),
            legend.key.size = unit(0.5,"line"),
            legend.text  = element_text(size = 6)) +
      {if(plot_caption == "TRUE") labs(x = "CNV Value", y = "Percentage of CNV Number", caption = "Note: Black integer - Number of CNVs. Red percentage - Proportion of CNVs in total.\nBlue percentage - Proportion of Overlapped CNVs in each type group.")} +
      {if(plot_caption == "FALSE") labs(x = "CNV Value", y = "Percentage of CNV Number")}


    print(compare_plot)
    dev.off()
    fwrite(cnv_cal, file = paste0(folder, "/", title_f, ".summary"), sep = "\t", quote = FALSE)
  }

  #dealing with data
  if (is.null(def_tar_map)) {
    cnv_tar <- fread(cnv_tar)
    cnv_tar$version <- "Verision_TAR" # add version in dataframe
    colnames(cnv_tar) <- paste(colnames(cnv_tar), "TAR", sep = "_") #add suffix to all colnames
    tar_colnames <- colnames(cnv_tar) #set original column names use for extracting columns after matching

    cnv_def <- fread(cnv_def)
    cnv_def$version <- "Version_DEF"
    colnames(cnv_def) <- paste(colnames(cnv_def), "DEF", sep = "_")
    def_colnames <- colnames(cnv_def)

    #######compare results in DEF at first-------------------------------------------------------------------------
    #1. find overlapped CNV between UCD and TAR
    setkey(cnv_tar, Sample_ID_TAR, Chr_TAR, Start_TAR, End_TAR)
    overlap_def <- foverlaps(cnv_def, cnv_tar, by.x = c("Sample_ID_DEF", "Chr_DEF", "Start_DEF", "End_DEF"), type = "any", nomatch = NULL)

    #might have some duplicated rows after find overlap, because of some CNV in TAR larger than DEF
    # or some CNV in DEF larger than TAR caused by the SNP position and density
    uniq_overlap_def <- subset(overlap_def, select = def_colnames)

    #find out non-overlap CNV, original CNV - overlapped CNV
    final_overlap_def <- unique(uniq_overlap_def)

    #2. find non-overlap CNV between UND and TAR
    non_overlap_def <- dplyr::setdiff(cnv_def, final_overlap_def)

    if (nrow(non_overlap_def) == 0) {
      print("These two input data are completely same, comparison will stoped.")
    }

    #3. checking if overlap and non-overlap results are correct by counting the total number
    if (nrow(non_overlap_def) + nrow(final_overlap_def) == nrow(cnv_def)) {
      print(paste0("Comparison Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total ", nrow(cnv_def)))
    } else {print("Comparison failed in validation, please use the original output files from HandyCNV as the input files")}

    #merge all cnv results then make comparison plot
    final_overlap_def$Check_overlap <- "Overlap"
    non_overlap_def$Check_overlap <- "Non-Overlap"
    checkover_indiv_1 <- rbind(final_overlap_def, non_overlap_def)
    colnames(checkover_indiv_1) <- sub("_DEF", "", colnames(checkover_indiv_1))

    plot_comparison(cnv_checkover = checkover_indiv_1, title_fig = "checkover_indiv_1", width_1 = width_1, height_1 = height_1, legend_x = legend_x, legend_y = legend_y, plot_caption = plot_caption)

    #4. write out the overlap and non-overlap CNVs
    fwrite(final_overlap_def, file = paste0(folder, "/overlap_cnv_1.indiv"), sep = "\t", quote = FALSE)
    fwrite(non_overlap_def, file = paste0(folder, "/non_overlap_cnv_1.indiv"), sep = "\t", quote = FALSE)
    fwrite(checkover_indiv_1, file = paste0(folder, "/cnv_all_inidv_1.checkoverlap"), sep = "\t", quote = FALSE)

    #5.summarize difference
    overlap_percent_indiv_1 <- round(nrow(final_overlap_def) / nrow(cnv_def), 3) * 100
    non_overlap_percent_indiv_1 <- round(nrow(non_overlap_def) / nrow(cnv_def), 3) * 100
    print(paste0("The number of overlaped CNVs on individual level is ", nrow(final_overlap_def), ", which is around ", overlap_percent_indiv_1, " percent in the first file"))
    print(paste0("The number of Non-overlaped CNVs on individual level is ", nrow(non_overlap_def), ", which is around ", non_overlap_percent_indiv_1, " percent in the first file"))


    #setkey(cnv_def_tar, Chr_DEF, Start_DEF, End_DEF)
    #find overlap on population level
    setkey(cnv_tar, Chr_TAR, Start_TAR, End_TAR)
    pop_overlap <- foverlaps(cnv_def, cnv_tar, by.x = c("Chr_DEF", "Start_DEF", "End_DEF"), type = "any", nomatch = NULL)

    pop_overlap_def <- subset(pop_overlap, select = def_colnames)
    final_pop_overlap_def <- unique(pop_overlap_def)
    non_overlap_pop_def <- dplyr::setdiff(cnv_def, final_pop_overlap_def)
    if (nrow(non_overlap_pop_def) + nrow(final_pop_overlap_def) == nrow(cnv_def)) {
      print(paste0("Comparison to the first file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total of ", nrow(cnv_def)))
    } else {print("Comparison to the first file was failed in validation, please use the original output files from HandyCNV as the input files")}

    fwrite(final_pop_overlap_def, file = paste0(folder, "/overlap_cnv_1.popu"), sep = "\t", quote = FALSE)
    fwrite(non_overlap_pop_def, file = paste0(folder, "/non_overlap_cnv_1.popu"), sep = "\t", quote = FALSE)

    #make comparison plot
    final_pop_overlap_def$Check_overlap <- "Overlap"
    non_overlap_pop_def$Check_overlap <- "Non-Overlap"
    checkover_pop_1 <- rbind(final_pop_overlap_def, non_overlap_pop_def)
    colnames(checkover_pop_1) <- sub("_DEF", "", colnames(checkover_pop_1))


    plot_comparison(cnv_checkover = checkover_pop_1, title_fig = "checkover_pop_1", width_1 = width_1, height_1 = height_1, legend_x = legend_x, legend_y = legend_y, plot_caption = plot_caption)

    fwrite(checkover_pop_1, file = paste0(folder, "/cnv_all_population_1.checkoverlap"), sep = "\t", quote = FALSE)

    #5.summarize difference
    #DEF
    overlap_percent_pop <- round(nrow(final_pop_overlap_def) / nrow(cnv_def), 3) * 100
    non_overlap_percent_pop <- round(nrow(non_overlap_pop_def) / nrow(cnv_def), 3) * 100
    print(paste0("The number of overlaped CNVs on population level is ", nrow(final_pop_overlap_def), ", which is around ", overlap_percent_pop, " percent in first file."))
    print(paste0("The number of Non-overlaped CNVs on population level is ", nrow(non_overlap_pop_def), ", which is around ", non_overlap_percent_pop, " percent in first file"))

    overlap_summary <- data.frame("Item" = c("Overlapped CNVs", "Non-overlapped CNVs"),
                                  "Individual Level" = c(overlap_percent_indiv_1, non_overlap_percent_indiv_1),
                                  "Population Level" = c(overlap_percent_pop, non_overlap_percent_pop))
    print("The final comparison results of the first file as follows: ")
    print(overlap_summary)
    fwrite(overlap_summary, file = paste0(folder, "/overlap_cnv_1.summary"), sep = "\t", quote = FALSE)
    print("Comparison to the first file was finished.")


    ##########compare results in TAR at second-------------------------------------------------------------------
    #1. find overlapped CNV between UCD and TAR
    setkey(cnv_def, Sample_ID_DEF, Chr_DEF, Start_DEF, End_DEF)
    #because some snp cannot find in TAR map, so we removed these CNV used right_def instead
    overlap_tar <- foverlaps(cnv_tar, cnv_def, by.x = c("Sample_ID_TAR", "Chr_TAR", "Start_TAR", "End_TAR"), type = "any", nomatch = NULL)

    #might have some duplicated rows after find overlap, because of some CNV in TAR larger than DEF
    # or some CNV in DEF larger than TAR caused by the SNP position and density
    uniq_overlap_tar <- subset(overlap_tar, select = tar_colnames)

    #find out non-overlap CNV, original CNV - overlapped CNV
    final_overlap_tar <- unique(uniq_overlap_tar)

    #2. find non-overlap CNV between UND and TAR
    non_overlap_tar <- dplyr::setdiff(cnv_tar, final_overlap_tar)

    #3. checking if overlap and non-overlap results are correct by counting the total number
    if (nrow(non_overlap_tar) + nrow(final_overlap_tar) == nrow(cnv_tar)) {
      print(paste0("Comparison to second file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total ", nrow(cnv_tar)))
    } else {print("Comparison to second file was failed in validation, please use the original output files from HandyCNV as the input files")}

    #4. write out the overlap and non-overlap CNVs
    fwrite(final_overlap_tar, file = paste0(folder, "/overlap_cnv_tar_2.indiv"), sep = "\t", quote = FALSE)
    fwrite(non_overlap_tar, file = paste0(folder, "/non_overlap_cnv_tar_2.indiv"), sep = "\t", quote = FALSE)

    final_overlap_tar$Check_overlap <- "Overlap"
    non_overlap_tar$Check_overlap <- "Non-Overlap"
    checkover_indiv_2 <- rbind(final_overlap_tar, non_overlap_tar)
    colnames(checkover_indiv_2) <- sub("_TAR", "", colnames(checkover_indiv_2))

    fwrite(checkover_indiv_2, file = paste0(folder, "/cnv_all_indiv_2.checkoverlap"), sep = "\t", quote = FALSE)

    plot_comparison(cnv_checkover = checkover_indiv_2, title_fig = "checkover_indiv_2", width_1 = width_1, height_1 = height_1, legend_x = legend_x, legend_y = legend_y, plot_caption = plot_caption)

    #5.summarize difference
    overlap_percent_indiv_2 <- round(nrow(final_overlap_tar) / nrow(cnv_tar), 3) * 100
    non_overlap_percent_indiv_2 <- round(nrow(non_overlap_tar) / nrow(cnv_tar), 3) * 100
    print(paste0("The number of overlaped CNVs on individual level is ", nrow(final_overlap_tar), ", which is around ", overlap_percent_indiv_2, " percent in second file."))
    print(paste0("The number of Non-overlaped CNVs on individual level is ", nrow(non_overlap_tar), ", which is around ", non_overlap_percent_indiv_2, " percent in second file."))


    #setkey(cnv_def_tar, Chr_DEF, Start_DEF, End_DEF)
    #find overlap on population level
    setkey(cnv_def, Chr_DEF, Start_DEF, End_DEF)
    pop_overlap_2 <- foverlaps(cnv_tar, cnv_def, by.x = c("Chr_TAR", "Start_TAR", "End_TAR"), type = "any", nomatch = NULL)
    pop_overlap_tar <- subset(pop_overlap_2, select = tar_colnames)
    final_pop_overlap_tar <- unique(pop_overlap_tar)
    non_overlap_pop_tar <- dplyr::setdiff(cnv_tar, final_pop_overlap_tar)
    if (nrow(non_overlap_pop_tar) + nrow(final_pop_overlap_tar) == nrow(cnv_tar)) {
      print(paste0("Comparison to the second file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total of ", nrow(cnv_tar)))
    } else {print("Comparison to the second file was failed in validation, please use the original output files from HandyCNV as the input files")}

    fwrite(final_pop_overlap_tar, file = paste0(folder, "/overlap_cnv_tar.popu"), sep = "\t", quote = FALSE)
    fwrite(non_overlap_pop_tar, file = paste0(folder, "/non_overlap_cnv_tar.popu"), sep = "\t", quote = FALSE)

    final_pop_overlap_tar$Check_overlap <- "Overlap"
    non_overlap_pop_tar$Check_overlap <- "Non-Overlap"
    checkover_pop_2 <- rbind(final_pop_overlap_tar, non_overlap_pop_tar)
    colnames(checkover_pop_2) <- sub("_TAR", "", colnames(checkover_pop_2))


    plot_comparison(cnv_checkover = checkover_pop_2, title_fig = "checkover_pop_2", width_1 = width_1, height_1 = height_1, legend_x = legend_x, legend_y = legend_y, plot_caption = plot_caption)
    fwrite(checkover_pop_2, file = paste0(folder, "/checkover_pop_2.txt"), sep = "\t", quote = FALSE)

    #5.summarize difference
    #TAR
    overlap_percent_pop_2 <- round(nrow(final_pop_overlap_tar) / nrow(cnv_tar), 3) * 100
    non_overlap_percent_pop_2 <- round(nrow(non_overlap_pop_tar) / nrow(cnv_tar), 3) * 100
    print(paste0("The number of overlaped CNVs in seceond file on population level are ", nrow(final_pop_overlap_tar), ", which is around ", overlap_percent_pop_2, " percent"))
    print(paste0("The number of Non-overlaped CNVs in second file on population level are ", nrow(non_overlap_pop_tar), ", which is around ", non_overlap_percent_pop_2, " percent"))

    overlap_summary_2 <- data.frame("Item" = c("Overlapped CNVs", "Non-overlapped CNVs"),
                                    "Individual Level" = c(overlap_percent_indiv_2, non_overlap_percent_indiv_2),
                                    "Population Level" = c(overlap_percent_pop_2, non_overlap_percent_pop_2))
    print("The final comparison results of the second file as follows: ")
    print(overlap_summary_2)
    fwrite(overlap_summary_2, file = paste0(folder, "/overlap_cnv_2.summary"), sep = "\t", quote = FALSE)
    print("Task done. Comparison results were saved in your working directory")
  }

  # else {
  #   #convert coordinate for CNV and CNVR
  #   cnv_tar <- fread(cnv_tar)
  #   cnv_tar$version <- "Verision_TAR" # add version in dataframe
  #   colnames(cnv_tar) <- paste(colnames(cnv_tar), "TAR", sep = "_") #add suffix to all colnames
  #   tar_colnames <- colnames(cnv_tar) #set original column names use for extracting columns after matching
  #
  #   cnv_def <- fread(cnv_def)
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
  #   #matching the start position by snp name, validating by matching the original position and new matched position
  #   #cnv_def_tar <- merge(cnv_def, two_map, by.x = "Start_SNP_DEF", by.y = "Name_Map", all.x = TRUE)
  #
  #   #There are some duplicated rows after merge progress, because of there are some different SNP with same location in the map file
  #   #To solve this problem, we should check duplicated row after each merge step and remove the duplicates
  #   cnv_def_tar <- merge(cnv_def, two_map, by.x = c("Chr_DEF", "Start_DEF"), by.y = c("Chr_def_Map", "Position_def_Map"), all.x = TRUE)
  #   dup_index <- grep("TRUE", duplicated(cnv_def_tar[, c("Chr_DEF", "Start_DEF", "Sample_ID_DEF")]))
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
  #   dup_index <- grep("TRUE", duplicated(cnv_def_tar[, c("Chr_DEF", "End_DEF", "Sample_ID_DEF")]))
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
  #   cnv_tar_def <- merge(cnv_tar, two_map, by.x = c("Chr_TAR", "Start_TAR"), by.y = c("Chr_tar_Map", "Position_tar_Map"), all.x = TRUE)
  #   dup_index <- grep("TRUE", duplicated(cnv_tar_def[, c("Chr_TAR", "Start_TAR", "Sample_ID_TAR")]))
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
  #   dup_index <- grep("TRUE", duplicated(cnv_tar_def[, c("Chr_TAR", "End_TAR", "Sample_ID_TAR")]))
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
  #   fwrite(cnv_def_tar, file = paste0(folder, "/cleancnv_UtoA.coord"), sep = "\t", quote = FALSE)
  #   fwrite(cnv_tar_def, file = paste0(folder, "/cleancnv_AtoU.coord"), sep = "\t", quote = FALSE)
  #
  #   #two requirements for foverlap: 1) start <= end, 2) no NA in both start and end
  #   cnv_def_tar$End_TAR_Map[is.na(cnv_def_tar$End_TAR_Map)] <- 0
  #   cnv_def_tar$Start_TAR_Map[is.na(cnv_def_tar$Start_TAR_Map)] <- 0
  #   wrong_def <- cnv_def_tar[which((cnv_def_tar$End_TAR_Map - cnv_def_tar$Start_TAR_Map) <= 0), ]
  #   right_def <- cnv_def_tar[-c(which((cnv_def_tar$End_TAR_Map - cnv_def_tar$Start_TAR_Map) <= 0)), ]
  #   wrong_tar <- cnv_tar_def[which((cnv_tar_def$End_DEF_Map - cnv_tar_def$Start_DEF_Map) <= 0), ]
  #   wrong_tar <- wrong_tar[-c(which(wrong_tar$End_DEF_Map == 0)), ]
  #   right_tar <- setdiff(cnv_tar_def, wrong_tar)
  #
  #   #######compare results in DEF at first-------------------------------------------------------------------------
  #   #1. find overlapped CNV between UCD and TAR
  #   setkey(cnv_tar_def, Sample_ID_TAR, Chr_TAR, Start_TAR, End_TAR)
  #   #because some snp cannot find in TAR map, so we removed these CNV used right_def instead
  #   overlap_def <- foverlaps(right_def, cnv_tar_def, by.x = c("Sample_ID_DEF", "Chr_DEF", "Start_TAR_Map", "End_TAR_Map"), type = "any", nomatch = NULL)
  #
  #   #might have some duplicated rows after find overlap, because of some CNV in TAR larger than DEF
  #   # or some CNV in DEF larger than TAR caused by the SNP position and density
  #   uniq_overlap_def <- subset(overlap_def, select = def_convert_colnames)
  #
  #   #find out non-overlap CNV, original CNV - overlapped CNV
  #   teamp_overlap_def <- subset(overlap_def, select = def_convert_colnames)
  #   final_overlap_def <- unique(teamp_overlap_def)
  #
  #   #2. find non-overlap CNV between UND and TAR
  #   non_overlap_def <- dplyr::setdiff(cnv_def_tar, final_overlap_def)
  #
  #   #3. cheking if overlap and non-overlap results are correct by counting the total number
  #   if (nrow(non_overlap_def) + nrow(final_overlap_def) == nrow(cnv_def_tar)) {
  #     print(paste0("Comparison Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total ", nrow(cnv_def_tar)))
  #   } else {print("Comparison failed in validation, please use the original output files from HandyCNV as the input files")}
  #
  #   #4. write out the overlap and non-overlap CNVs
  #   fwrite(final_overlap_def, file = paste0(folder, "/overlap_cnv.indiv"), sep = "\t", quote = FALSE)
  #   fwrite(non_overlap_def, file = paste0(folder, "/non_overlap_cnv.indiv"), sep = "\t", quote = FALSE)
  #
  #
  #   final_overlap_def$Check_overlap <- "Overlap"
  #   non_overlap_def$Check_overlap <- "Non-Overlap"
  #   checkover_indiv_1 <- rbind(final_overlap_def, non_overlap_def)
  #   colnames(checkover_indiv_1) <- sub("_DEF", "", colnames(checkover_indiv_1))
  #
  #   plot_comparison(cnv_checkover = checkover_indiv_1, title_fig = "checkover_indiv_1", width_1 = width_1, height_1 = height_1, legend_x = legend_x, legend_y = legend_y)
  #
  #   #5.summarize difference
  #   overlap_percent_indiv <- round(nrow(final_overlap_def) / nrow(cnv_def_tar), 3) * 100
  #   non_overlap_percent_indiv <- round(nrow(non_overlap_def) / nrow(cnv_def_tar), 3) * 100
  #   print(paste0("The number of overlaped CNVs on individual level is ", nrow(final_overlap_def), ", which is around ", overlap_percent_indiv, " percent"))
  #   print(paste0("The number of Non-overlaped CNVs on individual level is ", nrow(non_overlap_def), ", which is around ", non_overlap_percent_indiv, " percent"))
  #
  #
  #   #setkey(cnv_def_tar, Chr_DEF, Start_DEF, End_DEF)
  #   #find overlap on population level
  #   setkey(cnv_tar_def, Chr_TAR, Start_TAR, End_TAR)
  #   pop_overlap <- foverlaps(right_def, cnv_tar_def, by.x = c("Chr_DEF", "Start_TAR_Map", "End_TAR_Map"), type = "any", nomatch = NULL)
  #   pop_overlap_def <- subset(pop_overlap, select = def_convert_colnames)
  #   final_pop_overlap_def <- unique(pop_overlap_def)
  #   non_overlap_pop_def <- dplyr::setdiff(cnv_def_tar, final_pop_overlap_def)
  #   if (nrow(non_overlap_pop_def) + nrow(final_pop_overlap_def) == nrow(cnv_def_tar)) {
  #     print(paste0("Comparison to the first file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total of ", nrow(cnv_def_tar)))
  #   } else {print("Comparison to the first file was failed in validation, please use the original output files from HandyCNV as the input files")}
  #
  #   fwrite(final_pop_overlap_def, file = paste0(folder, "/overlap_cnv.popu"), sep = "\t", quote = FALSE)
  #   fwrite(non_overlap_pop_def, file = paste0(folder, "/non_overlap_cnv.popu"), sep = "\t", quote = FALSE)
  #
  #   #plot comparison
  #   final_pop_overlap_def$Check_overlap <- "Overlap"
  #   non_overlap_pop_def$Check_overlap <- "Non-Overlap"
  #   checkover_pop_1 <- rbind(final_pop_overlap_def, non_overlap_pop_def)
  #   colnames(checkover_pop_1) <- sub("_DEF", "", colnames(checkover_pop_1))
  #
  #
  #   plot_comparison(cnv_checkover = checkover_pop_1, title_fig = "checkover_pop_1", width_1 = width_1, height_1 = height_1, legend_x = legend_x, legend_y = legend_y)
  #
  #   fwrite(checkover_pop_1, file = paste0(folder, "/cnv_all_population_1.checkoverlap"), sep = "\t", quote = FALSE)
  #
  #   #5.summarize difference
  #   #DEF
  #   overlap_percent_pop <- round(nrow(final_pop_overlap_def) / nrow(cnv_def_tar), 3) * 100
  #   non_overlap_percent_pop <- round(nrow(non_overlap_pop_def) / nrow(cnv_def_tar), 3) * 100
  #   print(paste0("The number of overlaped CNVs on population level is ", nrow(final_pop_overlap_def), ", which is around ", overlap_percent_pop, " percent in first file."))
  #   print(paste0("The number of Non-overlaped CNVs on population level is ", nrow(non_overlap_pop_def), ", which is around ", non_overlap_percent_pop, " percent in first file"))
  #
  #   overlap_summary <- data.frame("Item" = c("Overlapped CNVs", "Non-overlapped CNVs"),
  #                                 "Individual Level" = c(overlap_percent_indiv, non_overlap_percent_indiv),
  #                                 "Population Level" = c(overlap_percent_pop, non_overlap_percent_pop))
  #   print("The final comparison results of the first file as follows: ")
  #   print(overlap_summary)
  #   fwrite(overlap_summary, file = paste0(folder, "/overlap_cnv.summary"), sep = "\t", quote = FALSE)
  #   print("Comparison to the first file was finished.")
  #
  #
  #   ##########compare results in TAR at first-------------------------------------------------------------------
  #   #1. find overlaped CNV between UCD and TAR
  #   setkey(right_def, Sample_ID_DEF, Chr_DEF, Start_TAR_Map, End_TAR_Map)
  #   #because some snp cannot find in TAR map, so we removed these CNV used right_def instead
  #   overlap_tar <- foverlaps(cnv_tar_def, right_def, by.x = c("Sample_ID_TAR", "Chr_TAR", "Start_TAR", "End_TAR"), type = "any", nomatch = NULL)
  #
  #   #might have some duplicated rows after find overlap, because of some CNV in TAR larger than DEF
  #   # or some CNV in DEF larger than TAR caused by the SNP position and density
  #   uniq_overlap_tar <- subset(overlap_tar, select = tar_convert_colnames)
  #
  #   #find out non-overlap CNV, original CNV - overlapped CNV
  #   teamp_overlap_tar <- subset(overlap_tar, select = tar_convert_colnames)
  #   final_overlap_tar <- unique(teamp_overlap_tar)
  #
  #   #2. find non-overlap CNV between UND and TAR
  #   non_overlap_tar <- dplyr::setdiff(cnv_tar_def, final_overlap_tar)
  #
  #   #3. cheking if overlap and non-overlap results are correct by counting the total number
  #   if (nrow(non_overlap_tar) + nrow(final_overlap_tar) == nrow(cnv_tar_def)) {
  #     print(paste0("Comparison to second file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total ", nrow(cnv_tar_def)))
  #   } else {print("Comparison to second file was failed in validation, please use the original output files from HandyCNV as the input files")}
  #
  #   #4. write out the overlap and non-overlap CNVs
  #   fwrite(final_overlap_tar, file = paste0(folder, "/overlap_cnv_tar.indiv"), sep = "\t", quote = FALSE)
  #   fwrite(non_overlap_tar, file = paste0(folder, "/non_overlap_cnv_tar.indiv"), sep = "\t", quote = FALSE)
  #
  #   final_overlap_tar$Check_overlap <- "Overlap"
  #   non_overlap_tar$Check_overlap <- "Non-Overlap"
  #   checkover_indiv_2 <- rbind(final_overlap_tar, non_overlap_tar)
  #   colnames(checkover_indiv_2) <- sub("_TAR", "", colnames(checkover_indiv_2))
  #
  #   fwrite(checkover_indiv_2, file = paste0(folder, "/cnv_all_indiv_2.checkoverlap"), sep = "\t", quote = FALSE)
  #
  #   plot_comparison(cnv_checkover = checkover_indiv_2, title_fig = "checkover_indiv_2", width_1 = width_1, height_1 = height_1, legend_x = legend_x, legend_y = legend_y)
  #
  #   #5.summarize difference
  #   overlap_percent_indiv_2 <- round(nrow(final_overlap_tar) / nrow(cnv_tar_def), 3) * 100
  #   non_overlap_percent_indiv_2 <- round(nrow(non_overlap_tar) / nrow(cnv_tar_def), 3) * 100
  #   print(paste0("The number of overlaped CNVs on individual level is ", nrow(final_overlap_tar), ", which is around ", overlap_percent_indiv_2, " percent in second file."))
  #   print(paste0("The number of Non-overlaped CNVs on individual level is ", nrow(non_overlap_tar), ", which is around ", non_overlap_percent_indiv_2, " percent in second file."))
  #
  #
  #   #setkey(cnv_def_tar, Chr_DEF, Start_DEF, End_DEF)
  #   #find overlap on population level
  #   setkey(right_def, Chr_DEF, Start_TAR_Map, End_TAR_Map)
  #   pop_overlap_2 <- foverlaps(cnv_tar_def, right_def, by.x = c("Chr_TAR", "Start_TAR", "End_TAR"), type = "any", nomatch = NULL)
  #   pop_overlap_tar <- subset(pop_overlap_2, select = tar_convert_colnames)
  #   final_pop_overlap_tar <- unique(pop_overlap_tar)
  #   non_overlap_pop_tar <- dplyr::setdiff(cnv_tar_def, final_pop_overlap_tar)
  #   if (nrow(non_overlap_pop_tar) + nrow(final_pop_overlap_tar) == nrow(cnv_tar_def)) {
  #     print(paste0("Comparison to the second file was Passed Validation, the number of Overlap and Non-overlap CNV equal to the original number of CNVs which are in total of ", nrow(cnv_tar_def)))
  #   } else {print("Comparison to the second file was failed in validation, please use the original output files from HandyCNV as the input files")}
  #
  #   fwrite(final_pop_overlap_tar, file = paste0(folder, "/overlap_cnv_tar.popu"), sep = "\t", quote = FALSE)
  #   fwrite(non_overlap_pop_tar, file = paste0(folder, "/non_overlap_cnv_tar.popu"), sep = "\t", quote = FALSE)
  #
  #
  #   final_pop_overlap_tar$Check_overlap <- "Overlap"
  #   non_overlap_pop_tar$Check_overlap <- "Non-Overlap"
  #   checkover_pop_2 <- rbind(final_pop_overlap_tar, non_overlap_pop_tar)
  #   colnames(checkover_pop_2) <- sub("_TAR", "", colnames(checkover_pop_2))
  #
  #
  #   plot_comparison(cnv_checkover = checkover_pop_2, title_fig = "checkover_pop_2", width_1 = width_1, height_1 = height_1, legend_x = legend_x, legend_y = legend_y)
  #   fwrite(checkover_pop_2, file = paste0(folder, "/checkover_pop_2.txt"), sep = "\t", quote = FALSE)
  #
  #   #5.summarize difference
  #   #DEF
  #   overlap_percent_pop_2 <- round(nrow(final_pop_overlap_tar) / nrow(cnv_tar_def), 3) * 100
  #   non_overlap_percent_pop_2 <- round(nrow(non_overlap_pop_tar) / nrow(cnv_tar_def), 3) * 100
  #   print(paste0("The number of overlaped CNVs in seceond file on population level are ", nrow(final_pop_overlap_tar), ", which is around ", overlap_percent_pop_2, " percent"))
  #   print(paste0("The number of Non-overlaped CNVs in second file on population level are ", nrow(non_overlap_pop_tar), ", which is around ", non_overlap_percent_pop_2, " percent"))
  #
  #   overlap_summary_2 <- data.frame("Item" = c("Overlapped CNVs", "Non-overlapped CNVs"),
  #                                   "Individual Level" = c(overlap_percent_indiv_2, non_overlap_percent_indiv_2),
  #                                   "Population Level" = c(overlap_percent_pop_2, non_overlap_percent_pop_2))
  #   print("The final comparison results of the second file as follows: ")
  #   print(overlap_summary_2)
  #   fwrite(overlap_summary_2, file = paste0(folder, "/overlap_cnv_2.summary"), sep = "\t", quote = FALSE)
  #   print("Task done. Comparison results were saved in your working directory")
  # }
}

