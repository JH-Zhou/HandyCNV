plot_cnv_source <- function(cnvr, clean_cnv, pedigree, Frequent_threshold) {
  cnvr <- fread(cnvr)
  high_freq <- filter(cnvr, Frequent >= Frequent_threshold)

  for (i in 1:nrow(high_freq)){
    cnv_visual(clean_cnv = clean_cnv, chr_id = high_freq$Chr[i],
               start_position = high_freq$Start[i]/1000000, end_position = high_freq$End[i]/1000000,  report_id = 1,
               pedigree = pedigree)
  }
}
