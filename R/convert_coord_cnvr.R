
#' Title convert_coord
#' only need Chr, Start and End
#' @param input_ars
#' @param input_umd
#' @param map
#'
#' @import data.table dplyr
#'
#' @return
#' @export convert_coord
#'
#' @examples
convert_coord <- function(input_ars =NULL, input_umd = NULL, map){
  map <- fread(map)

  #convert from ARS to UMD
  if (is.null(input_umd)) {
    cnvr <- fread(input_ars)
    #first step, matching the start position
    cnvr_start <- dplyr::left_join(x = cnvr, y = map, by =c("Chr" = "Chr_tar", "Start" = "Position_tar")) %>%
      rename(Start_UMD = Position_def, Start_Match = Match) %>%
      unique()

    #second step, matching the end position
    cnvr_end <- dplyr::left_join(x = cnvr_start, y = map, by = c("Chr" = "Chr_tar", "End" = "Position_tar")) %>%
      rename(End_UMD = Position_def, End_Match = Match)

    #third step, extract all information we need
    cnvr_new_coord <- cnvr_end %>%
      select(c("CNVR_ID", "Chr", "Start", "End", "Start_UMD", "End_UMD", "Start_Match", "End_Match")) %>%
      unique()

    #summary the quality of convertion
    print("The quality of converion of Start Position as below: FALSE means not match between two version.")
    print(table(cnvr_new_coord$Start_Match))
    print("The quality of converion of End Position as below: FALSE means not match between two version.")
    print(table(cnvr_new_coord$End_Match))

    cnvr_right <- cnvr_new_coord %>%
      filter(Start_Match == "Match" & End_Match == "Match") %>%
      filter(Start_UMD < End_UMD) %>%
      select(c("CNVR_ID", "Chr", "Start_UMD", "End_UMD")) %>%
      rename(Start = "Start_UMD", End = "End_UMD") %>%
      setDT()
    fwrite(cnvr_new_coord, "cnvr_convert_coord.txt", sep ="\t", quote = FALSE)
    fwrite(cnvr_right, "cnvr_convert_coord.correct", sep = "\t", quote = FALSE)
  } else {
    cnvr = fread(input_umd)
    #convert from UMD to ARS
    cnvr_start <- dplyr::left_join(x = cnvr, y = map, by =c("Chr" = "Chr_def", "Start" = "Position_def")) %>%
      rename(Start_ARS = Position_tar, Start_Match = Match)

    #sencond step, matching the end position
    cnvr_end <- dplyr::left_join(x = cnvr_start, y = map, by = c("Chr" = "Chr_def", "End" = "Position_def")) %>%
      rename(End_ARS = Position_tar, End_Match = Match)

    #thord step, extract all information we need
    cnvr_new_coord <- cnvr_end %>%
      select(c("CNVR_ID", "Chr", "Start", "End", "Start_ARS", "End_ARS", "Start_Match", "End_Match")) %>%
      unique()

    #summary the quality of convertion
    print("The quality of converion of Start Position as below: FALSE means not match between two version.")
    print(table(cnvr_new_coord$Start_Match))
    print("The quality of converion of End Position as below: FALSE means not match between two version.")
    print(table(cnvr_new_coord$End_Match))

    cnvr_right <- cnvr_new_coord %>%
      filter(Start_Match == "Match" & End_Match == "Match") %>%
      filter(Start_ARS < End_ARS) %>%
      select(c("CNVR_ID", "Chr", "Start_ARS", "End_ARS")) %>%
      rename(Start = "Start_ARS", End = "End_ARS") %>%
      setDT()
    fwrite(cnvr_new_coord, "cnvr_convert_coord.txt", sep ="\t", quote = FALSE)
    fwrite(cnvr_right, "cnvr_convert_coord.correct", sep = "\t", quote = FALSE)
  }
}



#' Title compare_interval
#' compare cnvr, each input file should contain Chr, Start and End columns,
#' the Chr column should be the number only, for example: 1 not chr1
#' @param interval_1
#' @param interval_2
#' @param map
#'
#' @import data.table
#'
#' @return overlap
#' @export compare_interval
#'
#' @examples
compare_interval <- function(interval_1, interval_2, map = NULL){
  cnvr_1 <- fread(interval_1)
  cnvr_1$Chr <- as.factor(cnvr_1$Chr)
  cnvr_2 <- fread(interval_2)
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

