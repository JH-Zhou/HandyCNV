# Prepare map file for PennCNV and Plink with UMD_3.1 and ARS 1.2 assembly

#' Title Prepare SNP Map
#'
#' @param map_from_GenomeStudio
#' @param convert_map
#'
#' @return
#' @export
#'
#' @examples
convert_map <- function(umd_map, standard_ARS_map){
  umd_map <- fread(umd_map)
  names(umd_map) <- c("Chr_UMD", "Name","Morgan_UMD", "Position_UMD")

  ars_map <- fread(standard_ARS_map)
  names(ars_map) <- c("Chr_ARS","Name", "Mb_ARS", "Position_ARS")

  umd_ars_map <- merge(umd_map, ars_map, by = "Name", all.x = TRUE, sort = FALSE)
  #repalce X, Y, MT to 30, 31, 33
  umd_ars_map$Chr_UMD <- sub("X", "30", umd_ars_map$Chr_UMD)
  umd_ars_map$Chr_UMD <- sub("Y", "31", umd_ars_map$Chr_UMD)
  umd_ars_map$Chr_UMD <- sub("MT", "33", umd_ars_map$Chr_UMD)
  umd_ars_map$Chr_UMD <- as.numeric(umd_ars_map$Chr_UMD) # convert chr to numeric

  #compare snp difference between two vertions of map
  #summarize how many SNP on each Chromosome
  snp_umd <- data.frame(table(umd_ars_map$Chr_UMD))
  names(snp_umd) <- c("Chr", "SNP_Freq_UMD")
  #print("The number of SNP on UMD Chromosome:")
  #print(snp_umd)


  snp_ars <- data.frame(table(umd_ars_map$Chr_ARS))
  names(snp_ars) <- c("Chr", "SNP_Freq_ARS")
  #print("The number of SNP on UMD Chromosome:")
  #print(snp_ars)

  snp_ars_umd <- merge(snp_ars, snp_umd, all = TRUE)
  snp_ars_umd[is.na(snp_ars_umd)] <- 0
  snp_ars_umd$Difference <- abs(snp_ars_umd$SNP_Freq_ARS - snp_ars_umd$SNP_Freq_UMD)

  print("Ploting the SNP difference plot between two versions of map......")
  #snp_ars_umd$Chr <- as.numeric(snp_ars_umd$Chr)
  png("SNP_difference_by_chromosome.png", res = 300, width = 3500, height = 2000)
  snp_diff <- ggplot(snp_ars_umd, aes(x = Chr)) +
    geom_bar(aes(y = SNP_Freq_ARS, fill = "ARS-UCD 1.2"), stat = "identity", position = "identity") +
    geom_point(aes(y = SNP_Freq_UMD, shape = "UMD 3.1"), stat = "identity", position = "identity", size = 4, color = "turquoise3") +
    geom_text(aes(y = SNP_Freq_ARS + 500, label = Difference)) +
    theme_classic() +
    labs(x = "Chromosome", y = "Number of SNP", shape = NULL, fill = NULL, caption = "**The number represents the number of different SNPS in the two versions") +
    theme(legend.position = c(0.9, 0.9)) # legend position - 0 is left/bottom, 1 is top/right
  print(snp_diff)
  dev.off()
  if(file.exists("SNP_difference_by_chromosome.png")){
    print("SNP difference plot was saved in your working directory.")
  }

  #check if SNP on the same Chr between tow version of MAP
  umd_ars_map$Match <- umd_ars_map$Chr_UMD == umd_ars_map$Chr_ARS
  print("Ploting SNP Posotion comparison results.....")
  png(filename = "comparison_SNP_Position_before_QC.png", res = 300, height = 3500, width = 4000)
  comparison_bf_qc <- ggplot(umd_ars_map, aes(x = Position_UMD, y = Position_ARS, color = Match)) +
    geom_point() +
    theme_bw() +
    scale_y_continuous(labels = unit_format(unit = "Mb", scale = 1e-6)) +
    scale_x_continuous(labels = unit_format(unit = "Mb", scale = 1e-6)) +
    facet_wrap(~Chr_ARS) +
    labs(x = "Position on UMD 3.1", y = "Position on ARS-UCD 1.2")
  print(comparison_bf_qc)
  dev.off()
  if (file.exists("comparison_SNP_Position_before_QC.png")) {
    print("Comparison of SNP Position plot was saved in your working directory.")
  }

  #check how many negtive postion in new assembly and replace it to positive
  print("Checking if negtive position exsits in new converted map....")
  negtive_position <- table(umd_ars_map$Position_ARS < 0)
  print(negtive_position)
  print(paste0("There are ", negtive_position[2], " SNPs with negtive postion in new assembly." ))

  #replace the negtive position to positive
  umd_ars_map$Position_ARS <- sub("-", "", umd_ars_map$Position_ARS)
  umd_ars_map$Position_ARS <- as.numeric(umd_ars_map$Position_ARS)

  umd_ars_map$Match <- umd_ars_map$Chr_UMD == umd_ars_map$Chr_ARS

  print("Ploting SNP Position comparison results after correctting the negtive position.....")
  png(filename = "comparison_SNP_Position_after_QC.png", res = 300, height = 3500, width = 4000)
  comparison_after_qc <- ggplot(umd_ars_map, aes(x = Position_UMD, y = Position_ARS, color = Match)) +
    geom_point() +
    theme_bw() +
    scale_y_continuous(labels = unit_format(unit = "Mb", scale = 1e-6)) +
    scale_x_continuous(labels = unit_format(unit = "Mb", scale = 1e-6)) +
    facet_wrap(~Chr_ARS) +
    labs(x = "Position on UMD 3.1", y = "Position on ARS-UCD 1.2")
  print(comparison_after_qc)
  dev.off()
  if (file.exists("comparison_SNP_Position_after_QC.png")) {
    print("Comparison of SNP Position plot after correction was saved in your working directory.")
  }

  check_match <- umd_ars_map %>% group_by(Match) %>% count(Match, name = "Count")
  print("Comparasion of Chrosomes matching coditions as following:")
  print(check_match)
  check_chr <- umd_ars_map %>% group_by(Chr_ARS) %>% count(Chr_UMD, name = "Count_ARS")

  #how many type of SNP between two version
  #1 same chr
  #2 same unkonw chr
  #3 diferenct chr
  #4. UMD unkown position
  #5. ARS unkown position
  #6. wrong chr
  #7. miss value

  print("Checking how many SNPs are on the different chromosome between two version of map....")
  check_chr$Type <- ""
  for (i in 1:nrow(check_chr)) {
    if (is.na(check_chr$Chr_ARS[i])) {
      check_chr$Type[i] <- "Missing SNP on ARS"
    }
    else if (!(check_chr$Chr_UMD[i] == 0) & !(check_chr$Chr_ARS[i] == 0) & (check_chr$Chr_UMD[i] == check_chr$Chr_ARS[i])) {
      check_chr$Type[i] <- "Same Chromosome"
    }
    else if (check_chr$Chr_ARS[i] == 0 & check_chr$Chr_UMD[i] == 0) {
      check_chr$Type[i] <- "Both Unknown Position"
    }
    else if (check_chr$Chr_ARS[i] == 0 & !(check_chr$Chr_UMD[i] == 0)) {
      check_chr$Type[i] <- "Unknown Position in ARS"
    }
    else if (!(check_chr$Chr_ARS[i] == 0) & check_chr$Chr_UMD[i] == 0) {
      check_chr$Type[i] <- "Unknown Position in UMD"
    }
    else {
      check_chr$Type[i] <- "Different Chromosome"
    }
  }

  comparasion_snp <- check_chr %>% group_by(Type) %>% summarise("Count" = sum(Count_ARS))
  print("The comparsion details of SNP between two Assembly as following: ")
  print(comparasion_snp)

  #SNP density on each map


  print("Preparing the map files for PennCNV and Plink....")
  # prepare map file for pennCNV and plink
  penncnv_map_umd <- umd_map[, c("Name", "Chr_UMD", "Position_UMD")]
  names(penncnv_map_umd) <- c("Name", "Chr", "Position")

  penncnv_map_ars <- umd_ars_map[, c("Name", "Chr_ARS", "Position_ARS")]
  names(penncnv_map_ars) <- c("Name", "Chr", "Position")

  plink_map_ars <- umd_ars_map[, c("Chr_ARS", "Name","Mb_ARS", "Position_ARS")]
  plink_map_ars[is.na(plink_map_ars)] <- 0
  plink_map_ars$Mb_ARS <- 0 #plink map format the third column is Morgans position, if not sure better replace to 0

  fwrite(umd_ars_map, file = "umd_ars_map.map", sep ="\t", quote = FALSE)
  fwrite(snp_ars_umd, file = "number_of_snp_on_chromosome.txt", sep ="\t", quote = FALSE)
  fwrite(penncnv_map_ars, file = "ars.map", sep = "\t", quote = FALSE, col.names = TRUE)
  fwrite(penncnv_map_umd, file = "umd.map", sep = "\t", quote = FALSE, col.names = TRUE)
  fwrite(plink_map_ars, file = "ars_plink.map", sep = "\t", quote = FALSE, col.names = FALSE)
  fwrite(comparasion_snp, file = "diffirence_in_two_map.txt", sep ="\t", quote = FALSE, col.names = TRUE)

  if (file.exists("ars.map") & file.exists("umd.map") & file.exists("ars_plink.map")) {
    print("Target map was converted to PennCNV map file, file name as 'umd.map'.")
    print("New converted map for PennCNV was saved in working directory, file name as 'ars.map'.")
    print("New converted map for Plink was saved in working directory, file name as 'ars_plink.map'.")
  }

  chr_length_ars <- data.frame("Chr" = c(30:1), "ARS_Length" = c( 139.009144, 51.098607, 45.94015, 45.612108, 51.992305,
                                                                  42.350435, 62.317253, 52.498615, 60.773035, 69.862954,
                                                                  71.974595, 63.449741, 65.820629, 73.167244, 81.013979,
                                                                  85.00778, 82.403003, 83.472345, 87.216183, 106.982474,
                                                                  103.308737, 105.454467, 113.31977, 110.682743, 117.80634,
                                                                  120.089316, 120.000601, 121.005158, 136.231102, 158.53411))
  snp_ars_umd_density <- merge(snp_ars_umd, chr_length_ars, all.x = TRUE)

  chr_length_umd <- data.frame("Chr" = c(30:1), "UMD_Length" = c(148.823899, 51.505224, 46.312546, 45.407902, 51.681464, 42.90417, 62.71493, 52.530062,
                                                                 61.435874, 71.599096, 72.042655, 64.057457,
                                                                 66.004023,75.158596, 81.724687, 85.296676, 84.64839, 84.24035, 91.163125,
                                                                 107.310763, 104.305016, 105.70825, 113.384836, 112.638659, 119.458736,
                                                                 121.1914245, 120.829699, 121.430405, 137.060424, 158.337067))
  snp_ars_umd_density <- merge(snp_ars_umd_density, chr_length_umd, all.x = TRUE)
  snp_ars_umd_density$ARS_Density <- round(snp_ars_umd_density$SNP_Freq_ARS / snp_ars_umd_density$ARS_Length, digits = 0)
  snp_ars_umd_density$UMD_Density <- round(snp_ars_umd_density$SNP_Freq_UMD / snp_ars_umd_density$UMD_Length, digits = 0)

  second_y_value <- max(snp_ars_umd_density$SNP_Freq_UMD) / max(snp_ars_umd_density$ARS_Density,na.rm = TRUE)
  png("SNP_difference_density.png", res = 300, width = 3500, height = 2000)
  diff_density <- ggplot(snp_ars_umd_density, aes(x = Chr)) +
    geom_bar(aes(y = SNP_Freq_ARS, fill = "SNP Number on ARS"), stat = "identity", position = "identity") +
    geom_point(aes(y = SNP_Freq_UMD, shape = "SNP Number on UMD"), stat = "identity", position = "identity", size = 4, color = "turquoise3") +
    geom_text(aes(y = SNP_Freq_ARS + 500, label = Difference)) +
    geom_line(aes(y = UMD_Density * second_y_value, group = 1), lwd = 1.2, color = "turquoise3") +
    geom_line(aes(y = ARS_Density * second_y_value, group = 1), lwd = 1.2, color = "red") +
    geom_text(aes(y = ARS_Density * second_y_value, label = ARS_Density)) +
    geom_text(aes(y = UMD_Density * second_y_value, label = UMD_Density)) +
    scale_y_continuous(name = "Number of SNP", sec.axis = sec_axis(~./ second_y_value, name = "SNP Density")) +
    theme_classic() +
    labs(x = "Chromosome", y = "Number of SNP", shape = NULL, fill = NULL, caption = "**The number represents the number of different SNPS in the two versions") +
    theme(legend.position = c(0.9, 0.9)) # legend position - 0 is left/bottom, 1 is top/right
  print(diff_density)
  dev.off()

}



