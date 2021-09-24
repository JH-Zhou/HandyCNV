#' Prepare the phased genotype into the standard format
#'
#' This function used to prepare phased genotype file to standard format in 'HandyCNV' to get haplotype for any interested genomic region,
#' The input file should phased by Beagle 5.1.
#'
#' @param phased_geno phased genotype file from Beagle 5.1, the first nine columns contain basic information, the rest are genotype
#' @param convert_letter convert allele from 0, 1 to letter in REF and ALF
#'
#' @importFrom data.table fread
#' @importFrom purrr map
#' @import dplyr tidyr
#'
#' @return A list contain map and genotype information
#' @export prep_phased
#'
prep_phased <- function(phased_geno, convert_letter = FALSE){
  cat("Reading input file...\n")
  #1.read the phased chromosome VCF file from Beagle
  phased <- fread(file = phased_geno, header = T) #columns are sample ID, rows are SNP ID

  if(convert_letter == 'TRUE'){
    cat("Coverting haplotype from 0, 1 to Ref and ALT...\n")
    phased <- phased %>%
      gather(key, value, -c('#CHROM', POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)) %>%
      rowwise() %>%
      mutate(value = gsub("0", REF, gsub("1", ALT, value))) %>%
      spread(key, value)
  }


  #2.replace Sample name of datafame, because cannot read _ in name...
  name_sample <- colnames(phased) #original input name
  name_order <- paste0("V",c(1:length(name_sample))) #generate new name of samples
  name_index <- data.frame(name = name_sample, new_name = name_order) #build name index between 'name_sample' and 'name_order'
  colnames(phased) <- paste0("V",c(1:ncol(phased))) #replace columns name with new name in phased data

  phased_base <- phased[, 1:9]

  #3.load function to split columns
  split_vcf <- function(x){
    phased %>%
      select(all_of(x)) %>% #Use `all_of(x)` instead of `x` for the external variable.
      separate_(x, paste0(x, c(".1", ".2")))
  }

  #4.run split function for all columns
  cat("Spliting Sample columns...\n")
  target_hap <- names(phased)[10:ncol(phased)] %>% #first nine columns are default in Beagle
    map(split_vcf) %>%
    as.data.frame() #change to data.frame

  #5.convert SNP became to column and sample convert to row
  cat("Tansposing genotype...\n")
  target_hap_2 <- data.frame(t(target_hap))

  phase_list <- list("map" = phased_base, "geno" = target_hap_2, "name_index" = name_index)
  return(phase_list)
}


#' Pick out the names and positions of start and end SNPs of a genomic interval
#'
#' Get the closest SNPs within the given interval, working on the 'map' file generated from 'pre_phased function'.
#'
#' @param phased_input map file
#' @param chr chromosome ID
#' @param start start position, unite in bp'
#' @param end end position, unite in 'bp'
#'
#' @import dplyr
#' @importFrom tidyr replace_na
#'
#' @return SNP position
#' @export closer_snp
closer_snp <- function(phased_input, chr = 14, start, end){
  position <- phased_input[, 1:3] %>%
    filter(V1 == chr & V2 >= start*1000000 & V2 <= end*1000000) %>%
    summarise(min = min(V2), max = max(V2))
  print(position)
  return(position)
}


#' Get haplotype
#'
#' Generate haplotypes from Phased genotype file by given region.
#'
#' @param geno phase genotype and map file from 'prep_phased' function
#' @param pos start and end position from 'closer_snp' function
#'
#' @import dplyr
#' @importFrom tidyr unite replace_na
#'
#' @return A list contain Sample with recoded Haplotype and Index of Haplotype and Recode ID
#' @export get_haplotype
#'
get_haplotype <- function(geno, pos){

  #get the order of row of target region in phased data
  start_tar <- which(geno$map$V2 == pos$min)
  print(start_tar)
  end_tar <- which(geno$map$V2 == pos$max)
  print(end_tar)

  #extract map for roh region
  roh_map <- geno$map %>%
             filter(V2 >= pos$min & V2 <= pos$max)
  cat(paste0("There are ", nrow(roh_map), " SNPs in this region...\n"))

  #construct haplotype
  hap_aim <- geno$geno %>%
    select(paste0("X", start_tar):paste0("X", end_tar)) %>%
    na_if(., y = "") %>% #replace missing value as NA
    unite(col = Haploid, sep = "", remove = TRUE)
  hap_aim$Haploid <- gsub(pattern = "NA", replacement = "N", x = hap_aim$Haploid) #replace all 'NA' as 'N'

  #get the sample index from column for haplotype
  hap_aim_name <- data.frame(id = row.names(hap_aim)) %>%
    separate(id, into = c("id_1", "suf"))
  #merge sample index and Haplotype
  hap_aim_f <- cbind(hap_aim_name, hap_aim)
  #get sample original ID to Haplotype
  hap_aim_f <- merge(hap_aim_f, geno$name_index, by.x = "id_1", by.y =  "new_name", all.x = T) %>%
    separate(col = name, into = c("A","Sample_ID"), sep = "_", extra = "merge")

  #prepare genotype of Haplotype
  cat("Recoding Haploid...\n")
  uniq_hap <- hap_aim %>%
    group_by(Haploid) %>%
    summarise(Freq = n()) %>%
    arrange(-Freq) %>%
    mutate(Re_code = seq(1, nrow(.), 1))

  cat("Generating Haplotype...\n")
  hap_type <- hap_aim_f %>%
    left_join(uniq_hap, by = "Haploid") %>%
    group_by(Sample_ID) %>%
    arrange(Re_code) %>% #make sure no mirror type generated in collapse
    mutate(Recode_Type = paste0(Re_code, collapse = "-")) %>%
    select(Sample_ID, Recode_Type) %>%
    unique()
  #data.frame(table(Recode_type = hap_type$Recode_Type))

  my_list <- list("hap_type" = hap_type, "hap_index" = uniq_hap, "roh_map" = roh_map)

  return(my_list)
}


# Visualize haplotype
#' The order of haplotype present in plot are start from high frequency to low from top to bottom
#'
#' @param haplotype a list of haplotype information that generate by get_haplotype function
#' @param filter_freq used to filter haplotype by frequency, only the selected haplotype will present in plot
#' @param letter_size integer number to adjust the letter size of GCTA etc
#' @param show_letter assign "True" to present base pairs in plot, "False" to hide text, is useful on large plot
#' @param xlab_text add additional text into X label, for example, plot haplotype for specific genes
#' @importFrom stringr str_split
#' @importFrom dplyr filter
#' @importFrom stats reorder
#' @import ggplot2
#' @return haplotype plot
#' @export haplo_visual
#'
haplo_visual <- function(haplotype, filter_freq = 0, letter_size = 2, show_letter = TRUE, xlab_text = ""){

  haplo <- haplotype$hap_index %>% filter(Freq >= filter_freq)
  map <- haplotype$roh_map

  hap_heat_data <- data.frame()
  for (i in 1:nrow(haplo)){
    hap_b <- data.frame("Bases" = unlist(stringr::str_split(string = haplo$Haploid[i], pattern = "")),
                        "Order" =  map$V2,
                        #"Order" = seq(1, nchar(haplo$Haploid[i])),
                        "HapID" = paste0("H", i))
    hap_heat_data <- rbind(hap_heat_data, hap_b)
  }

  xlab = paste0(xlab_text,nrow(map), " SNPs from chr", map$V1[1], ":",map$V2[1], "-", map$V2[nrow(map)],"bp")
  hap_plot <- ggplot(hap_heat_data, aes(as.factor(Order), reorder(HapID, desc(HapID)), fill= Bases)) +
    geom_tile() +
    {if(show_letter == "TRUE") geom_text(aes(label=Bases), size = letter_size)} +
    #geom_vline(xintercept = gene_pos$Start[1], linetype = "dashed") +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    labs(x = xlab, y = "Haplotype ID")
  return(hap_plot)
}

