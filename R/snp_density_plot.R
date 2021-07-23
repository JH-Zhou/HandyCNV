#' Title  SNP density plot
#'
#' Plot SNP density form SNP genotyping map file, the default map file in Plink format. The density will calculate by number of SNPs per Mb.
#'
#' @param map Four columns are mandatory in the order of Chr, SNP, CM, Position. No header requires
#' @param max_chr The maximum Chromosome present in SNP density plot. Default value = 29
#' @param top_density Threshold of higher density region, integer number, unit in SNPs/Mb. Default value = 72, suggest setting as the value of Mean + SD of SNP density.
#' @param low_density Threshold of lower density region, integer number, unit in SNPs/Mb. Default value = 33suggest setting as the value of Mean - SD of SNP density.
#' @param color_top Point color of higher density SNPs
#' @param color_low Point color of lower density SNPs
#' @param color_mid Point color of middle regions density SNPS
#' @param legend_position Vectors, default values are c(0.9, 0.1), indicate left bottom corner
#' @param x_label Default label is "Physical distance (Mb)"
#' @param y_label Default label is "Number of SNPs per Mb"
#' @param ncol_1 Customize the number of column in final plot, default value is 5.
#'
#' @import dplyr ggplot2
#' @importFrom tidyr separate
#' @importFrom data.table fread
#' @importFrom stats sd
#'
#' @return SNP density plot
#' @export plot_snp_density
#'
plot_snp_density <- function(map, max_chr = 29, top_density = 72, low_density = 33, color_top = "yellow", color_low = "blue", color_mid = "red", legend_position = c(0.9, 0.1), x_label = "Physical distance (Mb)", y_label = "Number of SNPs per Mb", ncol_1 = 5){
  if(typeof(map) == "character"){
    snp_map <- fread(file = map)
  } else {
   snp_map = map
  }

  if(!ncol(snp_map) == 4){
    stop("The input map file is incorrectly formatted\n",
         "Four columns are mandatory in the order of Chr, SNP, CM, Position.")
  }

  names(snp_map) <- c("Chr", "SNP", "CM", "Pos")

  #count how many SNP on autosome
  num_snp_autochr <- snp_map %>% filter(Chr %in% c(1:max_chr)) %>% nrow()
  cat(paste0("The total number of SNPs on selected chromosomes are ", num_snp_autochr, "\n"))

  chr_bin <- snp_map %>%
    group_by(Chr) %>%
    mutate(bins = cut(Pos, breaks = seq(from = 0, to = max(snp_map$Pos), by = 1000000), include.lowest = T)) %>% #add 1Mb interval group to each SNP
    separate(col = bins, into = c("min_edge", "max_edge"), sep = ",", remove = F)

  #chr_bin$bins <- as.character(chr_bin$bins) #To avoid warnings: Factor `bins` contains implicit NA, consider using `forcats::fct_explicit_na`

  chr_bin$max_edge <- round(as.integer(sub(x = chr_bin$max_edge, pattern = ']', replacement = "")) / 1000000, 0) #used in final plot

  #plot the SNP density
  top <- paste0(">=", top_density)
  mid <- paste0(low_density, "-", top_density)
  low <- paste0("<=", low_density)

  plot_data <- chr_bin %>%
    group_by(Chr, max_edge) %>%
    summarise(num = n()) %>%
    #mutate(group_color = if_else(num >= 60, ">=60", "<60")) %>%
    mutate(group_color = case_when(num >= top_density ~ top,
                                   num >= low_density & num < top_density ~ mid,
                                   num <= low_density ~ low)) %>%
    filter(Chr %in% c(1:max_chr)) #only plot 1 to 29 chr

  cat("The number of regions on different SNP density groups as following:\n")
  print(table(plot_data$group_color))

  cat("The Mean of SNPs density is ", round(mean(plot_data$num, na.rm = T), digits = 0), " SNPs/Mb\n")
  cat("The Mean + SD of SNP density is ", round(mean(plot_data$num, na.rm = T)  + sd(plot_data$num, na.rm = T), digits = 0), " SNPs/Mb\n")
  cat("The Mean - SD of one SD is ", round(mean(plot_data$num, na.rm = T) - sd(plot_data$num, na.rm = T), digits = 0), " SNPs/Mb\n")

  plot_data$group_color <- factor(plot_data$group_color, levels = c(low, mid, top))

  snp_density <- ggplot(data = plot_data, aes(x = max_edge, y = num)) +
    geom_line(aes(group=1), linetype = "dashed") +
    geom_point(aes(color = as.factor(group_color))) +
    scale_color_manual(name = "Density (SNP/Mb)",
                       labels = c(low, mid, top), # the order of color are sorted by alphabetical
                       values = c(color_low, color_mid, color_top)) +
    facet_wrap(~as.numeric(Chr), scales = "free", ncol = ncol_1) +
    theme_classic() +
    theme(legend.position = legend_position) +
    labs(y = y_label, x = x_label)

  return(snp_density)
}
