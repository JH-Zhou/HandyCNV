gene_compare <- function(gene_freq_1, gene_freq_2, gene_freq_3, common_gene_threshold = 3){
  #Comparation of gene frequancy between two CNV annoation results
  ars_gene <- fread(gene_freq_1)
  names(ars_gene)[2] <- "ARS_Freq"
  umd_gene <- fread(gene_freq_2)
  names(umd_gene)[2] <- "UMD_Freq"
  ars_umd_gene <- merge(ars_gene, umd_gene, by = "name2", all = TRUE) #combine all gene together between two version of results
  ars_umd_gene <- ars_umd_gene[-1, ] #delete the row with empty gene name
  ars_umd_gene[is.na(ars_umd_gene)] <- 0 #replece all missing value as 0

  #adding new column to indicate the diffirence between two Gene frequent within CNVs
  ars_umd_gene$Common_Gene <- ""
  for (i in 1:nrow(ars_umd_gene)) {
    if (ars_umd_gene$ARS_Freq[i] < common_gene_threshold & ars_umd_gene$UMD_Freq[i] < common_gene_threshold) {
      ars_umd_gene$Common_Gene[i] = "Common_Low_Freq"
    } else if (ars_umd_gene$ARS_Freq[i] < common_gene_threshold & ars_umd_gene$UMD_Freq[i] >= common_gene_threshold) {
      ars_umd_gene$Common_Gene[i] = "UMD_High_Freq"
    } else if(ars_umd_gene$ARS_Freq[i] >= common_gene_threshold & ars_umd_gene$UMD_Freq[i] < common_gene_threshold){
      ars_umd_gene$Common_Gene[i] = "ARS_High_Freq"
    } else {
      ars_umd_gene$Common_Gene[i] = "Common_High_Freq"
    }
  }

  gene_summary <- ars_umd_gene %>% group_by(Common_Gene) %>% count(Common_Gene)
  fwrite(gene_summary, file = "gene_comparation_summary.txt", sep = "\t", quote = FALSE)
  fwrite(ars_umd_gene, file = "gene_comparation.txt", sep = "\t", quote = FALSE)

  png(filename = "gene_comparation.png", res = 300, width = 3000, height = 2000)
  compare_plot <- ggplot(data = ars_umd_gene, aes(x = UMD_Freq, y = ARS_Freq, color = Common_Gene)) +
    geom_point(pch = 19, size = 3, position = position_jitter(height = 0.15, width = 0.15)) +
    theme_bw() +
    labs( x = "UMD Gene Frequency", y = "ARS Gene Frequency")
  print(compare_plot)
  dev.off()

  if (!is.null(gene_freq_3)) {
    part_gene <- fread(gene_freq_3)
    names(part_gene)[2] <- "Part_gene"
    three_gene <- merge(ars_umd_gene, part_gene, all = TRUE)
    three_gene[is.na(three_gene)] <- 0
    three_gene <- three_gene[-1,]

    for (i in 1:nrow(three_gene)) {
      if (three_gene$ARS_Freq[i] < common_gene_threshold & three_gene$UMD_Freq[i] < common_gene_threshold & three_gene$Part_gene[i] < common_gene_threshold) {
        three_gene$Common[i] = "1"
      }
      else if (three_gene$ARS_Freq[i] >= common_gene_threshold & three_gene$UMD_Freq[i] >= common_gene_threshold & three_gene$Part_gene[i] >= common_gene_threshold) {
        three_gene$Common[i] = "2"
      }
      else {
        three_gene$Common[i] = "3"
      }
    }

    three_gene_summary <- three_gene %>% group_by(Common) %>% count(Common)

    png("3d_gene_comparation.png", res = 300, height = 2000, width = 2000, bg = "transparent")
    three_gene <- plot3d(z = three_gene$ARS_Freq, x = three_gene$UMD_Freq, y = three_gene$Part_gene,
           zlab = "ARS_Freq", xlab = "UMD_Freq", ylab = "Part_Freq",
           col = three_gene$Common, type = "s", radius = 0.5)
    print(three_gene)
    dev.off()

    fwrite(three_gene, file = "three_gene_comparation.txt", sep = "\t", quote = FALSE)
    fwrite(three_gene_summary, file = "three_gene_comparation_summary.txt", sep = "\t", quote = FALSE)

  }

}
