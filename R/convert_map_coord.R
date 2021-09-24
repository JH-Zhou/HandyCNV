#' Convert map file
#'
#' Prepare map files for PennCNV and Plink with Default and Target assembly.
#' Four columns in fixed order are required for both default and target map file, they are Chr, Name, Morgan, Position
#'
#' @param default_map the map file to be converted
#' @param target_map the map file used to convert
#' @param defMap_title customize the title of default map in plot
#' @param tarMap_title customize the title of target map in plot
#' @param chr_set set the number of autosomes
#' @param species the name of species
#' @param col_1 set color for the type of Match in SNP comparison plot
#' @param col_2 set color for the type of Mismatch in SNP comparison plot
#' @param col_3 set color for the bar of Target_Map in SNP density plot
#' @param col_4 set color for line of Target_map in density plot
#' @param col_5 set color for point of Default_Map in SNP density plot
#' @param col_6 set color for line of Default_map in SNP density plot
#' @param legend_position default value = c(0.9, 0.9), used to adjust the legend in map comparison plot
#' @param height_c set height for Comparison plot
#' @param width_c set width for comparison plot
#' @param height_d set height for SNP density plot
#' @param width_d set width for SNP density plot
#' @import ggplot2 dplyr
#' @importFrom data.table fread fwrite setkey foverlaps setDT
#'
#' @return Details of comparison results between given map, and standard input map files used in PennCNV and Plink
#' @export convert_map
#'
convert_map <- function(default_map, target_map, chr_set = 29, defMap_title = "Default", tarMap_title = "Target", species = "Bovine", col_1 = "green4", col_2 = "red1", col_3 = "deeppink2", col_4 = "deeppink2", col_5 ="turquoise3", col_6 = "turquoise3", legend_position = c(0.9, 0.9), height_c = 16, width_c  = 22, height_d = 15, width_d = 25){
  if(!dir.exists(paths = "convert_map")){
    dir.create(path = "convert_map")
    print("New folder convert_map was created.")
  }

  if(typeof(default_map) == "character"){
    def_map <- fread(default_map)
  } else {
    def_map = default_map
  }


  if(typeof(target_map) == "character"){
    tar_map <- fread(target_map)
  } else {
    tar_map = target_map
  }


  if(species == "Bovine"){
    names(def_map) <- c("Chr_def", "Name","Morgan_def", "Position_def")
    #replace X, Y, MT to 30, 31, and 32
    def_map$Chr_def <- sub("X", "30", def_map$Chr_def)
    def_map$Chr_def <- sub("Y", "31", def_map$Chr_def)
    def_map$Chr_def <- sub("MT", "33", def_map$Chr_def)
    def_map$Chr_def <- as.integer(def_map$Chr_def)

    names(tar_map) <- c("Chr_tar","Name", "Mb_tar", "Position_tar")
    #replace X, Y, MT to 30, 31, and 32
    tar_map$Chr_tar <- sub("X", "30", tar_map$Chr_tar)
    tar_map$Chr_tar <- sub("Y", "31", tar_map$Chr_tar)
    tar_map$Chr_tar <- sub("MT", "33", tar_map$Chr_tar)
    tar_map$Chr_tar <- as.integer(tar_map$Chr_tar)

  } else if(!(species == "Bovine")) {
    names(def_map) <- c("Chr_def", "Name","Morgan_def", "Position_def")
    #replace X, Y, MT to 30, 31, and 32
    def_map$Chr_def <- sub("X", "87", def_map$Chr_def)
    def_map$Chr_def <- sub("Y", "88", def_map$Chr_def)
    def_map$Chr_def <- sub("MT", "89", def_map$Chr_def)
    def_map$Chr_def <- as.integer(def_map$Chr_def)


    names(tar_map) <- c("Chr_tar","Name", "Mb_tar", "Position_tar")
    #replace X, Y, MT to 30, 31, and 32
    tar_map$Chr_tar <- sub("X", "87", tar_map$Chr_tar)
    tar_map$Chr_tar <- sub("Y", "88", tar_map$Chr_tar)
    tar_map$Chr_tar <- sub("MT", "89", tar_map$Chr_tar)
    tar_map$Chr_tar <- as.integer(tar_map$Chr_tar)
  }


  #####1.find out the difference between two map
  print("Starting comparing the differences between the two versions")
  def_tar_map <- merge(def_map, tar_map, by = "Name", all = TRUE, sort = FALSE)
  def_tar_map$Chr_def <- as.integer(def_tar_map$Chr_def) # convert chr to numeric, for later plot in order
  def_tar_map$Chr_tar <- as.integer(def_tar_map$Chr_tar) # convert chr to numeric, for later plot in order

  #check how many negative position in new assembly and replace it to positive
  print("Checking if negtive position exsits in new converted map....")
  negtive_position <- table(def_tar_map$Position_tar < 0)
  print(negtive_position)
  print(paste0("There are ", negtive_position[2], " SNPs with negtive postion in new assembly." ))

  #replace the negative position to positive
  def_tar_map$Position_tar <- sub("-", "", def_tar_map$Position_tar)
  def_tar_map$Position_tar <- as.numeric(def_tar_map$Position_tar)

  #check matching results between two maps
  def_tar_map <- def_tar_map %>%
                 replace_na(list(Position_def = 0, Position_tar = 0, Chr_def = 99, Chr_tar = 99)) %>%
                 mutate(Match = if_else(Chr_def == Chr_tar, true = "Match", false = "Mismatch")) %>% #check if SNP on the same Chr between tow version of MAP
                 mutate(Type = case_when(Chr_tar == 99 & !(Chr_def == 99) ~ "Unfound in Target map",
                          !(Chr_tar == 99) & Chr_def == 99 ~ "Unfound in Default map",
                          Chr_tar == 99 & Chr_def == 99 ~ "Both missing",
                          Chr_def == 0 & Chr_tar == 0 ~ "Both unknown Chromosome",
                          Chr_def == 0 & !(Chr_tar == 0) ~ "Unknown Chromosome in Default map",
                          !(Chr_def == 0) & Chr_tar == 0 ~ "Unknown Chromosome in Target map",
                          Chr_def == Chr_tar ~ "Same Chromosome",
                          TRUE ~ "Different Chromosome"))

  check_match <- table(def_tar_map$Match)
  print("The matching results of SNP chromosome ID in the two versions are as follows:")
  print(check_match)
  comparasion_snp <- data.frame(Type = table(def_tar_map$Type))
  print("The details of comparison as follows: ")
  print(comparasion_snp)

  print("Ploting SNP comparison results.....")
  if(species == "Bovine"){
    note = "Note: 30 = X, 31 = Y, 32 = PAR, 33 = MT, 99 = Only found in one map"
  } else if (!(species == "Bovine")){
    note = "Note: 87 = X, 88 = Y, 89 = MT, 99 = Only found in one map"
  }
  #add manual color
  color_point <- c("Match" = col_1,
                   "Mismatch" = col_2)
  png(filename = "convert_map/comparison_snp.png", height = height_c, width = width_c, units = "cm", res = 350)
  comparison_snp <- ggplot(def_tar_map, aes(x = Position_def, y = Position_tar, color = Match)) +
    geom_point() +
    scale_color_manual(values = color_point) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = c(0.90, 0.07)) +
    scale_y_continuous(labels = unit_format(unit = "", scale = 1e-6)) +
    scale_x_continuous(labels = unit_format(unit = "", scale = 1e-6)) +
    facet_wrap(~Chr_tar) +
    labs(x = paste0("Position on ", defMap_title),
         y = paste0("Position on ", tarMap_title),
         color = NULL,
         caption = note)
  print(comparison_snp)
  dev.off()

  if (file.exists("convert_map/comparison_snp.png")) {
    print("Comparison of SNP by Chromosome plot was saved in working directory.")
  }


  #SNP density on each map
  print("Preparing the map files for PennCNV and Plink....")
  def_tar_map_fix <- merge(def_map, tar_map, by = "Name", all.x = TRUE, sort = FALSE) #cannot sort SNP
  # prepare map file for pennCNV and plink
  penncnv_map_def <- def_tar_map_fix %>%
                     select("Name" = Name, "Chr" = Chr_def, "Position" = Position_def)

  penncnv_map_tar <- def_tar_map_fix %>%
                     select("Name" = Name, "Chr" = Chr_tar, "Position" = Position_tar)

  plink_map_tar <- def_tar_map_fix %>%
                   select(Chr_tar, Name, Mb_tar, Position_tar) %>%
                   replace_na(list(Chr_tar = 0, Name = 0, Mb_tar = 0, Position_tar = 0))

  print("Writing converted map for PennCNV and Plink formats...")
  fwrite(def_tar_map, file = "convert_map/def_tar_map.map", sep ="\t", quote = FALSE)
  fwrite(penncnv_map_tar, file = "convert_map/target_penncnv.map", sep = "\t", quote = FALSE, col.names = TRUE)
  fwrite(penncnv_map_def, file = "convert_map/default_penncnv.map", sep = "\t", quote = FALSE, col.names = TRUE)
  fwrite(plink_map_tar, file = "convert_map/target_plink.map", sep = "\t", quote = FALSE, col.names = FALSE)
  fwrite(comparasion_snp, file = "convert_map/diffirence_in_two_map.txt", sep ="\t", quote = FALSE, col.names = TRUE)

  if (file.exists("convert_map/target_penncnv.map") & file.exists("convert_map/default_penncnv.map") & file.exists("convert_map/target_plink.map")) {
    print("Target map was converted to PennCNV map file...")
    print("Converted map for PennCNV format were saved in working directory.")
    print("Converted map for Plink format were saved in working directory.")
  }


  #compare snp difference between two versions of map
  #summarize how many SNP on each Chromosome
  snp_def <- data.frame(table(def_tar_map$Chr_def))
  names(snp_def) <- c("Chr", "SNP_Freq_def")

  snp_tar <- data.frame(table(def_tar_map$Chr_tar))
  names(snp_tar) <- c("Chr", "SNP_Freq_tar")

  snp_tar_def <- merge(snp_tar, snp_def, all = TRUE)
  snp_tar_def[is.na(snp_tar_def)] <- 0
  snp_tar_def$Difference <- abs(snp_tar_def$SNP_Freq_tar - snp_tar_def$SNP_Freq_def)
  snp_tar_def$Chr <- as.factor(snp_tar_def$Chr)


  print("Plotting the Difference of SNPs between two versions of map...")
  png(filename = "convert_map/SNP_difference_by_chromosome.png", height = height_d, width = width_d, units = "cm", res = 350)
  snp_diff <- ggplot(snp_tar_def, aes(x = Chr)) +
    geom_bar(aes(y = SNP_Freq_tar, fill = tarMap_title), stat = "identity", position = "identity") +
    scale_fill_manual(values = col_3) +
    geom_point(aes(y = SNP_Freq_def, shape = defMap_title), stat = "identity", position = "identity", size = 4, color = col_5) +
    geom_text(aes(y = SNP_Freq_tar + 500, label = Difference)) +
    theme_classic() +
    labs(x = "Chromosome", y = "Number of SNP", shape = NULL, fill = NULL, caption = "**The number represents the number of different SNPS between two maps") +
    theme(legend.position = c(0.9, 0.9)) # legend position - 0 is left/bottom, 1 is top/right
  print(snp_diff)
  dev.off()

  if(file.exists("convert_map/SNP_difference_by_chromosome.png")){
    print("SNP difference plot was saved in working directory.")
  }

  # if(species == "Bovine"){
  #   chr_length_tar <- data.frame("Chr" = c(30:1), "Length" = c( 139.009144, 51.098607, 45.94015, 45.612108, 51.992305,
  #                                                                   42.350435, 62.317253, 52.498615, 60.773035, 69.862954,
  #                                                                   71.974595, 63.449741, 65.820629, 73.167244, 81.013979,
  #                                                                   85.00778, 82.403003, 83.472345, 87.216183, 106.982474,
  #                                                                   103.308737, 105.454467, 113.31977, 110.682743, 117.80634,
  #                                                                   120.089316, 120.000601, 121.005158, 136.231102, 158.53411))
  #
  #
  #   chr_length_def <- data.frame("Chr" = c(30:1), "Length" = c(148.823899, 51.505224, 46.312546, 45.407902, 51.681464, 42.90417, 62.71493, 52.530062,
  #                                                                  61.435874, 71.599096, 72.042655, 64.057457,
  #                                                                  66.004023,75.158596, 81.724687, 85.296676, 84.64839, 84.24035, 91.163125,
  #                                                                  107.310763, 104.305016, 105.70825, 113.384836, 112.638659, 119.458736,
  #                                                                  121.1914245, 120.829699, 121.430405, 137.060424, 158.337067))
  # } else if(!(species == "Bovine")){
    #construct the border of each chromosome from imput data
    chr_length_tar <- tar_map %>%
                      group_by(Chr_tar) %>%
                      summarise(Length_tar = max(Position_tar, na.rm = TRUE) / 1000000) %>%
                      filter(Chr_tar >= 1 & Chr_tar <= chr_set) %>%
                      arrange(Chr_tar) %>%
                      rename(Chr = Chr_tar)


    chr_length_def <- def_map %>%
                      group_by(Chr_def) %>%
                      filter(Chr_def >= 1 & Chr_def <= chr_set) %>%
                      summarise(Length_def = max(Position_def, na.rm = TRUE) / 1000000) %>%
                      arrange(Chr_def) %>%
                      rename(Chr = Chr_def)
  # }


  snp_tar_def_density <- merge(snp_tar_def, chr_length_tar, by = "Chr", all.x = TRUE)
  snp_tar_def_density <- merge(snp_tar_def_density, chr_length_def, by = "Chr", all.x = TRUE)
  #calculate SNP density by chromosome, here need to skip some cell which contain missing value
  snp_tar_def_density <- snp_tar_def_density %>%
                         mutate(tar_Density = if_else(is.na(Length_tar), true = 0, false = round(SNP_Freq_tar / Length_tar, digits = 0)),
                                def_Density = if_else(is.na(Length_def), true = 0, false = round(SNP_Freq_def / Length_def, digits = 0)))

  snp_tar_def_density$tar_Density <- na_if(snp_tar_def_density$tar_Density, y = 0) #replace SNP density with 0 to NA
  snp_tar_def_density$def_Density <- na_if(snp_tar_def_density$def_Density, y = 0)

  fwrite(snp_tar_def_density, file = "convert_map/snp_density_chromosome.txt", sep ="\t", quote = FALSE)



  second_y_value <- max(snp_tar_def_density$SNP_Freq_def) / max(snp_tar_def_density$tar_Density, na.rm = TRUE)
  #add manual color for bar
  #bar_name = as.factor(tarMap_title)
  #color_bar <- c(bar_name = col_3)
  png(filename = "convert_map/SNP_difference_density.png", height = height_d, width = width_d, units = "cm", res = 350)
  diff_density <- ggplot(snp_tar_def_density, aes(x = Chr)) +
    geom_bar(aes(y = SNP_Freq_tar, fill = as.factor(tarMap_title)), stat = "identity", position = "identity") +
    scale_fill_manual(values = col_3) +
    geom_point(aes(y = SNP_Freq_def, shape = defMap_title), stat = "identity", position = "identity", size = 4, color = col_5) +
    geom_text(aes(y = SNP_Freq_tar + 500, label = Difference), size = 3, color = "red") +
    geom_line(aes(y = def_Density * second_y_value, group = 1), lwd = 1.2, color = col_6) +
    geom_line(aes(y = tar_Density * second_y_value, group = 1), lwd = 1.2, color = col_4) +
    geom_text(aes(y = tar_Density * second_y_value, label = tar_Density)) +
    geom_text(aes(y = def_Density * second_y_value, label = def_Density)) +
    scale_y_continuous(name = "Number of SNP", sec.axis = sec_axis(~./ second_y_value, name = "SNP Density")) +
    theme_classic() +
    labs(x = "Chromosome", y = "Number of SNP", shape = NULL, fill = NULL, caption = "**Red number indicates the number of different SNPS between two maps") +
    theme(legend.position = legend_position,
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0)) # legend position - 0 is left/bottom, 1 is top/right
  print(diff_density)
  dev.off()

}



