compare_gene <- function(gene_freq_1, gene_freq_2, gene_freq_3 = NULL, common_gene_threshold = 3){
  #Comparation of gene frequancy between two CNV annoation results
  ars_gene <- fread(gene_freq_1)
  names(ars_gene)[2] <- "gene_list_1"
  umd_gene <- fread(gene_freq_2)
  names(umd_gene)[2] <- "gene_list_2"
  ars_umd_gene <- merge(ars_gene, umd_gene, by = "name2", all = TRUE) #combine all gene together between two version of results
  ars_umd_gene <- ars_umd_gene[-1, ] #delete the row with empty gene name
  ars_umd_gene[is.na(ars_umd_gene)] <- 0 #replece all missing value as 0

  #adding new column to indicate the diffirence between two Gene frequent within CNVs
  ars_umd_gene$Common_Gene <- ""
  for (i in 1:nrow(ars_umd_gene)) {
    if (ars_umd_gene$gene_list_2[i] < common_gene_threshold & ars_umd_gene$gene_list_1[i] < common_gene_threshold) {
      ars_umd_gene$Common_Gene[i] = "Common_Low_Freq"
    } else if (ars_umd_gene$gene_list_2[i] < common_gene_threshold & ars_umd_gene$gene_list_1[i] >= common_gene_threshold) {
      ars_umd_gene$Common_Gene[i] = "High_Freq_list_1"
    } else if(ars_umd_gene$gene_list_2[i] >= common_gene_threshold & ars_umd_gene$gene_list_1[i] < common_gene_threshold){
      ars_umd_gene$Common_Gene[i] = "High_Freq_liste_2"
    } else {
      ars_umd_gene$Common_Gene[i] = "Common_High_Freq"
    }
  }

  gene_summary <- ars_umd_gene %>% group_by(Common_Gene) %>% count(Common_Gene)
  fwrite(gene_summary, file = "two_genelists_comparation_summary.txt", sep = "\t", quote = FALSE)
  fwrite(ars_umd_gene, file = "two_genelists_comparation.txt", sep = "\t", quote = FALSE)

  png(filename = "two_genelists_comparation.png", res = 300, width = 3000, height = 2000)
  compare_plot <- ggplot(data = ars_umd_gene, aes(x = gene_list_1, y = gene_list_2, color = Common_Gene)) +
    geom_point(pch = 19, size = 3, position = position_jitter(height = 0.15, width = 0.15)) +
    theme_bw() +
    labs( x = "Gene Frequency in First List", y = "Gene Frequency in Second List")
  print(compare_plot)
  dev.off()

  if (!is.null(gene_freq_3)) {
    part_gene <- fread(gene_freq_3)
    names(part_gene)[2] <- "gene_list_3"
    three_gene <- merge(ars_umd_gene, part_gene, all = TRUE)
    three_gene[is.na(three_gene)] <- 0
    three_gene <- three_gene[-1,]

    for (i in 1:nrow(three_gene)) {
      if (three_gene$gene_list_2[i] < common_gene_threshold & three_gene$gene_list_1[i] < common_gene_threshold & three_gene$gene_list_3[i] < common_gene_threshold) {
        three_gene$Common[i] = "1"
      }
      else if (three_gene$gene_list_2[i] >= common_gene_threshold & three_gene$gene_list_1[i] >= common_gene_threshold & three_gene$gene_list_3[i] >= common_gene_threshold) {
        three_gene$Common[i] = "2"
      }
      else {
        three_gene$Common[i] = "3"
      }
    }

    three_gene_summary <- three_gene %>% group_by(Common) %>% count(Common)

    plot3d(z = three_gene$gene_list_2, x = three_gene$gene_list_1, y = three_gene$gene_list_3,
           zlab = "gene_list_2", xlab = "gene_list_1", ylab = "gene_list_3",
           col = three_gene$Common, type = "s", radius = 0.5)
    rgl.postscript("three_genelists_comparison.pdf", "pdf")


    fwrite(three_gene, file = "three_genelists_comparation.txt", sep = "\t", quote = FALSE)
    fwrite(three_gene_summary, file = "three_genelists_comparation_summary.txt", sep = "\t", quote = FALSE)

  }

}
