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
  umd_ars_map$Position_ARS <- sub("-", "", umd_ars_map$Position_ARS)
  umd_ars_map$Match <- umd_ars_map$Chr_UMD == umd_ars_map$Chr_ARS
  check_match <- umd_ars_map %>% group_by(Match) %>% count(Match, name = "Count")
  print("Comparasion of Chrosomes matching coditions as following:")
  print(check_match)
  check_chr <- umd_ars_map %>% group_by(Chr_ARS) %>% count(Chr_UMD, name = "Count_ARS")

  #png("check_SNP_match.png", res = 300)
  #plot(umd_ars_map$Position_UMD, umd_ars_map$Position_ARS)
  #ggplot(umd_ars_map, aes(x = Position_UMD, y = Position_ARS, color = Match)) +
  #  geom_point() +
  #  theme_bw() +
  #  facet_wrap(~Chr_ARS)
  #dev.off()

  penncnv_map_umd <- umd_map[, c("Name", "Chr_UMD", "Position_UMD")]
  names(penncnv_map_umd) <- c("Name", "Chr", "Position")

  penncnv_map_ars <- umd_ars_map[, c("Name", "Chr_ARS", "Position_ARS")]
  names(penncnv_map_ars) <- c("Name", "Chr", "Position")

  plink_map_ars <- umd_ars_map[, c("Chr_ARS", "Name","Mb_ARS", "Position_ARS")]
  plink_map_ars[is.na(plink_map_ars)] <- 0
  plink_map_ars$Mb_ARS <- 0 #plink map format the third column is Morgans position, if not sure better replace to 0

  fwrite(penncnv_map_ars, file = "ars.map", sep = "\t", quote = FALSE, col.names = TRUE)
  fwrite(penncnv_map_umd, file = "umd.map", sep = "\t", quote = FALSE, col.names = TRUE)
  fwrite(plink_map_ars, file = "ars_plink.map", sep = "\t", quote = FALSE, col.names = FALSE)
  fwrite(check_chr, file = "diffirence_in_two_map.txt", sep ="\t", quote = FALSE, col.names = TRUE)
}

