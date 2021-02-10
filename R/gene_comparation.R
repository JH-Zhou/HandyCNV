#' Title compare_gene
#' Do comparison of genes between two annotated gene frequency list. The default genelist file was generate by call_gene function.
#' If the users coustomize the genelists, it requires at least two columns in each gene list, including fixed column names name2(Gene name) and Frequency(frequency of gene)
#'
#' @param gene_freq_1 first genelist, the default file was generated from call_gene function
#' @param gene_freq_2 second genelist, the default file was generated from call_gene function
#' @param gene_freq_3 third genelist, the default file was generated from call_gene function
#' @param gene_freq_4 fourth genelist, the default file was generated from call_gene function
#' @param common_gene_threshold integer input.could calculate by multiplying The Number of Sample by 0.05 or 0.1, et al
#' @param title_1 custom title(two and three lists) or label(four lists) in figture
#' @param title_2 custom title(two and three lists) or label(four lists) in figture
#' @param title_3 custom title(two and three lists) or label(four lists) in figture
#' @param title_4 custom title(two and three lists) or label(four lists) in figture
#' @param height_1 custom the height to save figure
#' @param width_1 custom the width to save figure
#'
#' @import dplyr ggplot2 scatterplot3d
#' @importFrom data.table fread fwrite setkey foverlaps setDT
#'
#' @return
#' @export compare_gene
#'
#' @examples
compare_gene <- function(gene_freq_1, gene_freq_2, gene_freq_3 = NULL, gene_freq_4 = NULL, common_gene_threshold = 3, title_1 = "list_1", title_2 = "list_2", title_3 = "list_3", title_4 = "list_4", height_1 = 10, width_1 =14){
  #check and creatfolder
  if(!file.exists("compare_gene")){
     dir.create("compare_gene")
    }

  if (is.null(gene_freq_3) & is.null(gene_freq_4)){
    #Comparation of gene frequancy between two CNV annoation results
    list_1 <- fread(gene_freq_1)
    colnames(list_1) <- paste(colnames(list_1), "1", sep = "_")

    list_2 <- fread(gene_freq_2)
    colnames(list_2) <- paste(colnames(list_2), "2", sep = "_")

    two_list <- merge(list_1, list_2, by.x = "name2_1", by.y = "name2_2", all = TRUE) #combine all gene together between two version of results
    two_list[is.na(two_list)] <- 0 #replece all missing value as 0

    #adding new column to indicate the diffirence between two Gene frequent within CNVs
    print("Clustering common gene between the two input files...")
    two_list <- two_list %>%
                mutate(Common_Gene = case_when(Frequency_1 < common_gene_threshold & Frequency_2 < common_gene_threshold ~ "Common_Low",
                                               Frequency_1 >= common_gene_threshold & Frequency_2 < common_gene_threshold ~ "High_Freq_list_1",
                                               Frequency_1 < common_gene_threshold & Frequency_2 >= common_gene_threshold ~ "High_Freq_list_2",
                                               TRUE ~ "Common_High"))

    gene_summary <- two_list %>% group_by(Common_Gene) %>% count(Common_Gene)
    print("Brief summary of comparison as below:")
    print(gene_summary)
    fwrite(gene_summary, file = "compare_gene/two_lists_comparison_summary.txt", sep = "\t", quote = FALSE)
    fwrite(two_list, file = "compare_gene/two_lists_comparison.txt", sep = "\t", quote = FALSE)

    #plot comparison
    compare_plot <- ggplot(data = two_list, aes(x = Frequency_1, y = Frequency_2, color = Common_Gene)) +
      geom_point(pch = 19, size = 3, position = position_jitter(height = 0.15, width = 0.15)) +
      geom_text_repel(data = subset(two_list, Common_Gene == "Common_High"), aes(label = name2_1), nudge_x = 0.3) +
      theme_classic() +
      theme(legend.position = c(0.85, 0.85)) +
      labs( x = title_1, y = title_2, color = NULL)
    ggsave(filename = "compare_gene/two_lists_comparison.png", plot = compare_plot, width = width_1, height = height_1, units = "cm", dpi = 300)

    if(file.exists("compare_gene/two_lists_comparison_summary.txt") & file.exists("compare_gene/two_lists_comparison.txt") & file.exists("compare_gene/two_lists_comparison.png")) {
      print("Gene comparison list, brife summary and comparison plot was saved in working directory.")
    } else{
      print("Checking outputs files was faild, please check the input files and resetting a working directory!")
    }
  }

  else if(is.null(gene_freq_4)){

    #read all gene lists
    list_1 <- fread(gene_freq_1)
    colnames(list_1) <- paste(colnames(list_1), "1", sep = "_")

    list_2 <- fread(gene_freq_2)
    colnames(list_2) <- paste(colnames(list_2), "2", sep = "_")

    list_3 <- fread(gene_freq_3)
    colnames(list_3) <- paste(colnames(list_3), "3", sep = "_")

    #merge gene lists
    print("Starting merging gene lists...")
    two_list <- merge(list_1, list_2, by.x = "name2_1", by.y = "name2_2", all = TRUE)
    three_list <- merge(two_list, list_3, by.x = "name2_1", by.y = "name2_3", all = TRUE)
    three_list[is.na(three_list)] <- 0

    print("Clustering common gene among three gene lists...")
    three_list <- three_list %>%
                  mutate(Common_Gene = case_when(Frequency_1 < common_gene_threshold & Frequency_2 < common_gene_threshold  & Frequency_3 < common_gene_threshold ~ "Common_Low",
                                                 Frequency_1 >= common_gene_threshold & Frequency_2 >= common_gene_threshold & Frequency_3 >= common_gene_threshold ~ "Common_High",
                                                 TRUE ~ "High_and_Low"))

    three_list_summary <- three_list %>% group_by(Common_Gene) %>% count(Common_Gene, name = "Number_Gene")
    print(paste0(length(which(three_list$Common_Gene == "Common_High")), " genes with higher frequency in all three gene lists, they are: "))
    print(subset(three_list[three_list$Common_Gene == "Common_High", ],
                 select = c("name2_1", "Frequency_1","Frequency_2", "Frequency_3", "Common_Gene"))) #here only plot selected columns
    print("Brief summary:")
    print(three_gene_summary)

    pdf("compare_gene/three_list_comparison.pdf", width = 7, height = 6, onefile = T)
    #par(mfrow = c(2,2))
    three_list$Common_Gene <- as.factor(three_list$Common_Gene) #set color for factors
    color_custom = c("red", "orange", "blue")[three_list$Common_Gene] #set color for factors
    gene_plot <- scatterplot3d::scatterplot3d(z = three_list$Frequency_3, x = three_list$Frequency_1, y = three_list$Frequency_2,
                                 zlab = title_3, xlab = title_1, ylab = title_2,
                                 color = color_custom, pch = 16)
    common_high <- subset(three_list, Common_Gene == "Common_High", select = c("name2_1", "Frequency_1","Frequency_2", "Frequency_3", "Common_Gene"))#add common high gene name in plot
    #common_high <- subset(three_gene, !Common_check == "Common_low")#add common high gene name in plot
    #common_high_coord <- gene_plot$xyz.convert(x = common_high[,2], y = common_high[, 4], z = common_high[, 3]) #convert 3d coordinate to 2d
    common_high_coord <- gene_plot$xyz.convert(common_high[, 2:4]) #convert 3d coordinates to 2d
    text(x = common_high_coord$x, y = common_high_coord$y,labels = common_high$name2_1, cex = 0.5, pos = 4) #add text
    legend("top", legend = c("Common High", "Common Low", "High and Low"), col =  c("red", "orange", "blue"), pch = 16, inset = -0.1, xpd = TRUE, horiz = TRUE, bty = "n")
    dev.off()

    fwrite(three_list, file = "compare_gene/three_genelists_comparison.txt", sep = "\t", quote = FALSE)
    fwrite(three_list_summary, file = "compare_gene/three_genelists_comparison_summary.txt", sep = "\t", quote = FALSE)
    if(file.exists("compare_gene/three_genelists_comparison.txt") & file.exists("compare_gene/three_genelists_comparison_summary.txt") & file.exists("compare_gene/gene_comparison.pdf")){
      print("Gene comparison list, brife summary and comparison plot was saved in working directory.")
    } else{
      print("Checking outputs files was faild, please check the input files and resetting a working directory!")
    }
  }
  else{
    #read all gene lists
    list_1 <- fread(gene_freq_1)
    colnames(list_1) <- paste(colnames(list_1), "1", sep = "_")

    list_2 <- fread(gene_freq_2)
    colnames(list_2) <- paste(colnames(list_2), "2", sep = "_")

    list_3 <- fread(gene_freq_3)
    colnames(list_3) <- paste(colnames(list_3), "3", sep = "_")

    list_4 <- fread(gene_freq_4)
    colnames(list_4) <- paste(colnames(list_4), "4", sep = "_")


    #merge gene lists
    print("Starting merging gene lists...")
    two_gene <- merge(list_1, list_2, by.x = "name2_1", by.y = "name2_2", all = TRUE)
    three_gene <- merge(two_gene, list_3, by.x = "name2_1", by.y = "name2_3", all = TRUE)
    four_gene <- merge(three_gene, list_4, by.x = "name2_1", by.y = "name2_4", all = TRUE)
    four_gene[is.na(four_gene)] <- 0

    print("Clustering common gene among four gene lists...")
    # 2 = common_low, 1 = common_high, 3 = High_and_low
    four_gene <- four_gene %>%
                 mutate(Common_check = case_when(Frequency_1 >= common_gene_threshold & Frequency_2 >= common_gene_threshold & Frequency_3 >= common_gene_threshold & Frequency_4 >= common_gene_threshold ~ "Common_high",
                                                Frequency_1 < common_gene_threshold & Frequency_2 < common_gene_threshold & Frequency_3 < common_gene_threshold & Frequency_4 < common_gene_threshold ~ "Common_low",
                                                TRUE ~ "High_and_Low"))

    four_gene_summary <- four_gene %>% group_by(Common_check) %>% count(Common_check, name = "Number_Gene")
    print(paste0(length(which(four_gene$Common_check == "Common_high")), " genes with higher frequency in all three gene lists, they are: "))
    print(subset(four_gene[four_gene$Common_check == "Common_high", ],
                 select = c("name2_1", "Frequency_1","Frequency_2", "Frequency_3", "Frequency_4", "Common_check"))) #here only plot selected columns
    print("Brief summary:")
    print(four_gene_summary)

    plot_data <- four_gene %>%
                 select(name2_1, Frequency_1, Frequency_2, Frequency_3, Frequency_4) %>%
                 filter(Frequency_1 >= common_gene_threshold | Frequency_2 >= common_gene_threshold | Frequency_3 >= common_gene_threshold | Frequency_4 >= common_gene_threshold)
    plot_data <- reshape2::melt(plot_data, id.vars = "name2_1")
    ggplot(plot_data, aes(x = variable , y = name2_1, fill = value)) +
      geom_tile() +
      scale_fill_gradientn(colours = c("yellow", "red"), na.value = "black") +
      theme_classic() +
      scale_x_discrete(labels = c(title_1, title_2, title_3, title_4)) +
      labs(x = NULL, y = "Gene", fill = "Quantity")
    ggsave(filename = "compare_gene/four_gene_heatmap.png", dpi = 300, height = height_1, width = width_1, units = "cm")

    fwrite(four_gene, file = "compare_gene/four_gene_comparison.txt", sep = "\t", quote = FALSE)
    fwrite(four_gene_summary, file = "compare_gene/four_gene_comparison_summary.txt", sep = "\t", quote = FALSE)
    if(file.exists("compare_gene/four_gene_comparison.txt") & file.exists("compare_gene/four_gene_comparison_summary.txt") & file.exists("compare_gene/four_gene_heatmap.png")){
      print("Gene comparison list, brife summary and comparison plot was saved in working directory.")
    } else{
      print("Checking outputs files was faild, please check the input files and resetting a working directory!")
    }
  }

}
