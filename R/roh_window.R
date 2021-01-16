#This function is used to summarise ROH by window size
#Contain three functions
#First step is generate a window size for all Chromosme
#Second step is find overlap between window size and roh by Chromosome, Start and End position
#Third step is summary results, count the roh frequency in each windowsize

#1.set module of windows for single chromosome

require(data.table, quietly = TRUE)
require(dplyr, quietly = TRUE)

window_module <- function(chr, size, chr_length){
  n_interval <- ceiling(chr_length/size)

  Chr <- rep(x = chr, times = n_interval)

  if (chr_length %% size == 0){
    Start <- seq(from = 0, to = chr_length - size, by = size)
    End <- seq(from = size, to = chr_length, by = size)
  } else{
    Start <- seq(from = 0, to = chr_length, by = size)
    End <- seq(from = size, to = chr_length + size, by = size)
  }

  window <- data.frame(Chr, Start, End)
  return(window)
}

#2.call all windows
set_window <- function(win_size){
  ars <- data.frame("chr" = c(29:1),
                    "length" = c( 51.098607, 45.94015, 45.612108, 51.992305,
                                  42.350435, 62.317253, 52.498615, 60.773035, 69.862954,
                                  71.974595, 63.449741, 65.820629, 73.167244, 81.013979,
                                  85.00778, 82.403003, 83.472345, 87.216183, 106.982474,
                                  103.308737, 105.454467, 113.31977, 110.682743, 117.80634,
                                  120.089316, 120.000601, 121.005158, 136.231102, 158.53411))

  windows <- data.frame()
  for(i in 1:nrow(ars)){
    window <- window_module(chr = ars$chr[i], chr_length = ars$length[i], size = win_size)
    windows <- rbind(windows, window)
  }

  windows <- dplyr::arrange(windows, Chr, Start) %>%
    mutate(window_id = paste0("Win_", seq(from = 1, to = nrow(windows), by = 1))) %>%
    mutate(Start_win = Start*1000000, End_win = End*1000000) %>%
    select(c("Chr", "Start_win", "End_win", "window_id"))
  windows <- data.table::setDT(windows)
  return(windows)
}


#mian function
roh_window <- function(roh, window_size = 5, length_autosomal = 2489.386) {
  roh <- fread(roh)

  #check if the input is a Plink results
  #convert to the stardards format if it is
  plink_roh_names <- c("FID", "IID", "PHE", "CHR", "SNP1", "SNP2", "POS1", "POS2",
                       "KB", "NSNP", "DENSITY", "PHOM", "PHET")
  if(all(colnames(roh) == plink_roh_names)){
    print("ROH with input file in PLINK format was detected, coverting to HandyCNV standard formats")
    colnames(roh) <- c("FID", "Sample_ID", "PHE", "Chr", "Start_SNP", "End_SNP",
                       "Start", "End", "Length", "Num_SNP", "DENSITY", "PHOM", "PHET")
    handycnv_name <- c("Sample_ID",	"Chr", "Start", "End", "Num_SNP",	"Length", "Start_SNP",	"End_SNP")
    roh <- roh %>% select(handycnv_name) %>% mutate(CNV_Value = 2)
    roh$Length <- roh$Length * 1000 #because the unit is kb in default plink results
    roh <- data.table::setDT(roh)
  }

  #summary of roh for individual
  indiv_roh <- roh %>%
    group_by(Sample_ID) %>%
    summarise(n = n_distinct(Start), total_length = sum(Length)/1000000, min_length = min(Length)/1000000, max_length = max(Length)/1000000) %>%
    mutate(means = total_length/n, F_roh = total_length/length_autosomal)
  print("Individual situations are shown below:")
  print(indiv_roh)
  fwrite(indiv_roh, file = "indiv_roh.txt", sep = "\t", quote = FALSE)

  #plot roh


  w <- set_window(win_size = window_size)

  setkey(w, Chr, Start_win, End_win)
  win_over <- data.table::foverlaps(x = roh, y = w, by.x = c("Chr", "Start", "End"), type = "any")

  roh_freq <- win_over %>%
    group_by(window_id) %>%
    summarise(number_roh = n_distinct(Sample_ID)) %>%
    merge(., y = w, by = "window_id") %>% #reassign the position for windows
    arrange(-number_roh)
  print("The top 10 ROHs by frequency are as follows:")
  print(roh_freq[1:10, ])
  fwrite(roh_freq, file = "roh_freq.txt", sep = "\t", quote = FALSE)
}
