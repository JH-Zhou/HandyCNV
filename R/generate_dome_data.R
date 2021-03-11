#' Title get_demo
#' The purpose of this function is to prepare demo data for 'HandyCNV'. The main reason is because of the SNP Signal Intensity and Genotype files always too large,
#' when we making plots can select the data inside all CNVR regions instead of using the whole file. The output will extract all Signal Intensity, genotype and map
#' information of given CNVRs. The new generated ped and map files keeps the same format with Plink Software, so we can --make-bed or --recode in Plink as well.
#'
#' @param intensity Four columns by order with
#' @param ped Export from GenomeStudio
#' @param map Export from GenomeStudio, requires four columns without header.
#' @param cnvr Generate by 'call_cnvr', requires at least have four columns: Chr, Start, End, Frequency
#' @param freq_threshold This used to filter the CNVRs, only prepare the demo data for CNVRs which passed this threshold
#' @param folder Set name of new folder to save results
#' @import dplyr tidyr
#' @importFrom data.table fread fwrite
#'
#' @return Demo data of Intensity, Map, Ped file which appears only in CNVRs
#' @export get_demo
#'
get_demo <- function(intensity = NULL, ped = NULL, map = NULL, cnvr = NULL, freq_threshold = 0, folder = "demo_data"){
  if (!(file.exists(folder))){
    dir.create(folder)
    print(paste0("New folder ", folder, " was created in working directory."))
  }
  print("Skip the first 9 lines to read the Signal Intensity file...")
  penn_intensity <- fread( file = intensity, skip = 9, header = TRUE)
  intensity_demo <- penn_intensity

  print("Reading the ped file...")
  ped <- fread(file = ped)

  print("Reading the map file...")
  map_snp <- fread(file = map)

  print("Reading the Interval list...")
  cnvr <- fread(file = cnvr)
  cnvr <- cnvr %>%
          filter(Frequent > freq_threshold)

  print(paste0(nrow(cnvr), " CNVRs were select to prepare Demo Data..."))

  #extract row names from map file
  get_map <- function(map_snp, chr, start_pos, end_pos){

    map_target <- map_snp %>%
      mutate(row_num = row.names(.)) %>%
      filter(V1 == chr & V4 > start_pos & V4 < end_pos) %>%
      arrange(row_num)

    return(map_target)
  }

  #creat an empty table to merge data in for loop. Note: usiang rbind() function requirs each columns in two table has the same type!
  lapply(map_snp, class)
  map_target <- data.frame("V1" = "",
                           "V2" = "",
                           "V3" = 1,
                           "V4" = 1,
                           "row_num" = "",
                           stringsAsFactors = FALSE)
  #start looping for each CNVR and merge all SNP togther
  for (i in 1:nrow(cnvr)){
    map_1 <- get_map(map_snp = map_snp, chr = cnvr$Chr[i], start_pos = cnvr$Start[i], end_pos = cnvr$End[i])
    map_target <- base::rbind(map_target, map_1)
    print(paste0("Processing ", round(i / nrow(cnvr), 4) * 100, " %"))
  }

  #remind to remove the fisrt row we generated before
  map_target <- map_target[-1, ]
  lapply(map_target, class)
  map_target$row_num <- as.integer(map_target$row_num) #we need to keep the map and ped file has the same order
  map_target <- map_target %>%
                arrange(row_num)
  #prepare the column ID according to row_num in map_target file  to extract genotype from ped file
  row_num <- sort(as.vector(map_target$row_num))

  #prepare columns ID for Ped file
  col_id <- paste0("V", sort(c(seq(1, 6, 1),  #first 6 information columns
                               row_num * 2 + 5,   #column name of first allele
                               row_num * 2 + 6))) #column name of second allele

  #check if any duplicated casued by mistake
  which(duplicated(col_id))

  print(paste0(nrow(map_target), " SNPs in total were selected from the original files."))

  #pick up the genotype for all SNPs
  ped_tartget <- ped %>%
                 select(c(col_id))

  intensity_target <- penn_intensity %>%
                      filter(`SNP Name` %in% map_target$V2)

  fwrite(intensity_target, file = paste0(folder, "/demo_intensity.txt"), sep = "\t", quote = FALSE, col.names = TRUE)
  fwrite(ped_tartget, file = paste0(folder, "/demo.ped"), sep = "\t", quote = FALSE, col.names = FALSE)
  fwrite(map_target[, -5], file = paste0(folder, "/demo.map"), sep = "\t", quote = FALSE, col.names = FALSE)

  if(file.exists(paste0(folder, "/demo_intensity.txt")) &
     file.exists(paste0(folder, "/demo.ped")) &
     file.exists(paste0(folder, "/demo.map"))){
    print("Task done.")
  } else {
    print("No output file detected, please check the input files carefully!")
  }
}

