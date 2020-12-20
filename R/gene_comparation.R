require(rgl)
require(data.table)
require(ggplot2)
#' Title Gene Comparison
#' There are two application scenarios:Compare two gene lists or three gene lists
#'
#' @param gene_freq_1
#' @param gene_freq_2
#' @param gene_freq_3
#' @param common_gene_threshold
#'
#' @return
#' @export
#'
#' @examples
compare_gene <- function(gene_freq_1, gene_freq_2, gene_freq_3 = NULL, common_gene_threshold = 3, title_1 = "genelist_1", title_2 = "genelist_2", title_3 = "genelist_3"){
  if (is.null(gene_freq_3)){
    #Comparation of gene frequancy between two CNV annoation results
    ars_gene <- fread(gene_freq_1)
    names(ars_gene)[2] <- "gene_list_1"
    umd_gene <- fread(gene_freq_2)
    names(umd_gene)[2] <- "gene_list_2"
    ars_umd_gene <- merge(ars_gene, umd_gene, by = "name2", all = TRUE) #combine all gene together between two version of results
    ars_umd_gene <- ars_umd_gene[-1, ] #delete the row with empty gene name
    ars_umd_gene[is.na(ars_umd_gene)] <- 0 #replece all missing value as 0

    #adding new column to indicate the diffirence between two Gene frequent within CNVs
    print("Clustering common gene between the two input files...")
    ars_umd_gene$Common_Gene <- ""
    for (i in 1:nrow(ars_umd_gene)) {
      if (ars_umd_gene$gene_list_2[i] < common_gene_threshold & ars_umd_gene$gene_list_1[i] < common_gene_threshold) {
        ars_umd_gene$Common_Gene[i] = "Common_Low"
      } else if (ars_umd_gene$gene_list_2[i] < common_gene_threshold & ars_umd_gene$gene_list_1[i] >= common_gene_threshold) {
        ars_umd_gene$Common_Gene[i] = "High_Freq_list_1"
      } else if(ars_umd_gene$gene_list_2[i] >= common_gene_threshold & ars_umd_gene$gene_list_1[i] < common_gene_threshold){
        ars_umd_gene$Common_Gene[i] = "High_Freq_liste_2"
      } else {
        ars_umd_gene$Common_Gene[i] = "Common_High"
      }
    }

    gene_summary <- ars_umd_gene %>% group_by(Common_Gene) %>% count(Common_Gene)
    print("Brief summary of comparison as below:")
    print(gene_summary)
    fwrite(gene_summary, file = "two_genelists_comparison_summary.txt", sep = "\t", quote = FALSE)
    fwrite(ars_umd_gene, file = "two_genelists_comparison.txt", sep = "\t", quote = FALSE)

    png(filename = "two_genelists_comparison.png", res = 300, width = 3000, height = 2000)
    compare_plot <- ggplot(data = ars_umd_gene, aes(x = gene_list_1, y = gene_list_2, color = Common_Gene)) +
      geom_point(pch = 19, size = 3, position = position_jitter(height = 0.15, width = 0.15)) +
      geom_text_repel(data = subset(ars_umd_gene, Common_Gene == "Common_High"), aes(label = name2), nudge_x = 0.3) +
      theme_classic() +
      theme(legend.position = c(0.9, 0.9)) +
      labs( x = title_1, y = title_2)
    print(compare_plot)
    dev.off()

    if(file.exists("two_genelists_comparison_summary.txt") & file.exists("two_genelists_comparison.txt") & file.exists("two_genelists_comparison.txt")) {
      print("Gene comparison list, brife summary and comparison plot was saved in working directory.")
    } else{
      print("Checking outputs files was faild, please check the input files and resetting a working directory!")
    }
  }

  else {

    #read all gene lists
    ars_gene <- fread(gene_freq_1)
    names(ars_gene)[2] <- "gene_list_1"

    umd_gene <- fread(gene_freq_2)
    names(umd_gene)[2] <- "gene_list_2"

    part_gene <- fread(gene_freq_3)
    names(part_gene)[2] <- "gene_list_3"

    #merge gene lists
    print("Starting merging gene lists...")
    two_gene <- merge(ars_gene, umd_gene, by = "name2", all = TRUE)
    three_gene <- merge(two_gene, part_gene, by = "name2", all = TRUE)
    three_gene[is.na(three_gene)] <- 0

    print("Clustering common gene among three gene lists...")
    # 2 = common_low, 1 = common_high, 3 = High_and_low
    for (i in 1:nrow(three_gene)) {
      if (three_gene$gene_list_2[i] < common_gene_threshold & three_gene$gene_list_1[i] < common_gene_threshold & three_gene$gene_list_3[i] < common_gene_threshold) {
        three_gene$Common_check[i] = "Common_low"
      }
      else if (three_gene$gene_list_2[i] >= common_gene_threshold & three_gene$gene_list_1[i] >= common_gene_threshold & three_gene$gene_list_3[i] >= common_gene_threshold) {
        three_gene$Common_check[i] = "Common_high"
      }
      else {
        three_gene$Common_check[i] = "High_and_Low"
      }
    }

    three_gene_summary <- three_gene %>% group_by(Common_check) %>% count(Common_check, name = "Number_Gene")
    print(paste0(length(which(three_gene$Common_check == "Common_high")), " genes with higher frequency in all three gene lists, they are: "))
    print(three_gene[three_gene$Common_check == "Common_high", ])
    print("Brief summary:")
    print(three_gene_summary)

    pdf("gene_comparison.pdf", width = 7, height = 6, onefile = T)
    #par(mfrow = c(2,2))
    three_gene$Common_check <- as.factor(three_gene$Common_check) #set color for factors
    color_custom = c("red", "orange", "blue")[three_gene$Common_check] #set color for factors
    gene_plot <- scatterplot3d::scatterplot3d(z = three_gene$gene_list_3, x = three_gene$gene_list_1, y = three_gene$gene_list_2,
                                 zlab = title_2, xlab = title_1, ylab = title_3,
                                 color = color_custom, pch = 16)
    common_high <- subset(three_gene, Common_check == "Common_high")#add common high gene name in plot
    #common_high <- subset(three_gene, !Common_check == "Common_low")#add common high gene name in plot
    #common_high_coord <- gene_plot$xyz.convert(x = common_high[,2], y = common_high[, 4], z = common_high[, 3]) #convert 3d coordinate to 2d
    common_high_coord <- gene_plot$xyz.convert(common_high[, 2:4])
    text(x = common_high_coord$x, y = common_high_coord$y,labels = common_high$name2, cex = 0.5, pos = 4) #add text
    legend("top", legend = c("Common High", "Common Low", "High and Low"), col =  c("red", "orange", "blue"), pch = 16, inset = -0.1, xpd = TRUE, horiz = TRUE, bty = "n")
    dev.off()

    fwrite(three_gene, file = "three_genelists_comparison.txt", sep = "\t", quote = FALSE)
    fwrite(three_gene_summary, file = "three_genelists_comparison_summary.txt", sep = "\t", quote = FALSE)
    if(file.exists("three_genelists_comparison.txt") & file.exists("three_genelists_comparison_summary.txt") & file.exists("gene_comparison.pdf")){
      print("Gene comparison list, brife summary and comparison plot was saved in working directory.")
    } else{
      print("Checking outputs files was faild, please check the input files and resetting a working directory!")
    }
  }

}
