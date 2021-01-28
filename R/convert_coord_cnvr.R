
#only need Chr, Start and End
convert_coord <- function(input_ars =NULL, input_umd = NULL, map){
  map <- fread(map)

  #convert from ARS to UMD
  if (is.null(input_umd)) {
    cnvr <- fread(input_ars)
    #first step, matching the start position
    cnvr_start <- dplyr::left_join(x = cnvr, y = map, by =c("Chr" = "Chr_ARS", "Start" = "Position_ARS")) %>%
      rename(Start_UMD = Position_UMD, Start_Match = Match) %>%
      unique()

    #sencond step, matching the end position
    cnvr_end <- dplyr::left_join(x = cnvr_start, y = map, by = c("Chr" = "Chr_ARS", "End" = "Position_ARS")) %>%
      rename(End_UMD = Position_UMD, End_Match = Match)

    #thord step, extract all information we need
    cnvr_new_coord <- cnvr_end %>%
      select(c("CNVR_ID", "Chr", "Start", "End", "Start_UMD", "End_UMD", "Start_Match", "End_Match")) %>%
      unique()

    #summary the quality of convertion
    print("The quality of converion of Start Position as below: FALSE means not match between two version.")
    print(table(cnvr_new_coord$Start_Match))
    print("The quality of converion of End Position as below: FALSE means not match between two version.")
    print(table(cnvr_new_coord$End_Match))
  } else{
    cnvr = fread(input_umd)
    #convert from UMD to ARS
    cnvr_start <- dplyr::left_join(x = cnvr, y = map, by =c("Chr" = "Chr_UMD", "Start" = "Position_UMD")) %>%
      rename(Start_ARS = Position_ARS, Start_Match = Match)

    #sencond step, matching the end position
    cnvr_end <- dplyr::left_join(x = cnvr_start, y = map, by = c("Chr" = "Chr_UMD", "End" = "Position_UMD")) %>%
      rename(End_ARS = Position_ARS, End_Match = Match)

    #thord step, extract all information we need
    cnvr_new_coord <- cnvr_end %>%
      select(c("CNVR_ID", "Chr", "Start", "End", "Start_ARS", "End_ARS", "Start_Match", "End_Match")) %>%
      unique()

    #summary the quality of convertion
    print("The quality of converion of Start Position as below: FALSE means not match between two version.")
    print(table(cnvr_new_coord$Start_Match))
    print("The quality of converion of End Position as below: FALSE means not match between two version.")
    print(table(cnvr_new_coord$End_Match))
  }
  cnvr_right <- cnvr_new_coord %>%
    filter(Start_Match == "TRUE" & End_Match == "TRUE") %>%
    filter(Start_UMD < End_UMD) %>%
    select(c("CNVR_ID", "Chr", "Start_UMD", "End_UMD")) %>%
    rename(Start = "Start_UMD", End = "End_UMD") %>%
    setDT()
  fwrite(cnvr_new_coord, "cnvr_convert_coord.txt", sep ="\t", quote = FALSE)
  fwrite(cnvr_right, "cnvr_convert_coord.correct", sep = "\t", quote = FALSE)

}


#compare cnvr, each input file should contain Chr, Start and End columns,
#the Chr column should be the number only, for example: 1 not chr1
compare_cnvrs <- function(cnvr_1, cnvr_2, map = NULL){
  cnvr_1 <- fread(cnvr_1)
  cnvr_1$Chr <- as.factor(cnvr_1$Chr)
  cnvr_2 <- fread(cnvr_2)
  cnvr_2$Chr <- as.factor(cnvr_2$Chr)
  #set key for cnvr_1
  setkey(cnvr_1, Chr, Start, End)
  names(cnvr_2) <- paste0(colnames(cnvr_2),"_2")

  #find overlap
  overlap <- data.table::foverlaps(x = cnvr_2, y = cnvr_1, by.x = c("Chr_2", "Start_2", "End_2"), type = "any")
  overlap <- overlap %>%
    filter(Start > 0) %>%
    mutate(overlap_len = pmin(End_2, End) - pmax(Start, Start_2) + 1)
  fwrite(overlap, file = "overlap.txt", sep ="\t", quote = FALSE)
  print("task done")
  return(overlap)
}
