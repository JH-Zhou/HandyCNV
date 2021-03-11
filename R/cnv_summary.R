#' Title summarize CNVs
#' Make graphs for CNVs
#'
#' @param clean_cnv the clean cnv file was generate by clean_cnv function
#' @param plot_sum_1 first type of combination of summary plot
#' @param plot_sum_2 second type of combination of summary plot
#' @param length_group set group of vectors to divide CNV length, unit is Mb. such as vector of ‘c(0.05, 0.3,  0.6, 1)’, means divide the CNV length into four group: '<0.05Mb', '0.05 - 0.3Mb', '0.3-0.6Mb' and '>1Mb', maximum accept 11 values
#' @param folder set name of new create folder which used to save the results
#' @param col_0 set color for 0 copy of CNV
#' @param col_1 set color for 1 copy of CNV
#' @param col_2 set color for 2 copy of CNV (which might be ROH)
#' @param col_3 set color for 3 copy of CNV
#' @param col_4 set color for 4 copy of CNV
#' @param height_sum1 set height of Summary Plot 1
#' @param width_sum1 set width of Summary Plot 1
#' @param height_sum2 set height of Summary Plot 2
#' @param width_sum2 set width of Summary Plot 2
#'
#' @import ggplot2 dplyr reshape2 tidyr
#' @importFrom data.table fread fwrite
#' @importFrom scales unit_format
#'
#' @return Some summary pictures
#' @export cnv_summary_plot
#'
cnv_summary_plot <- function(clean_cnv, length_group = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1, 2, 5), plot_sum_1 = FALSE, plot_sum_2 = FALSE, folder = 'cnv_summary_plot', col_0 = "hotpink",  col_1 = "turquoise", col_2 = "gray", col_3 = "tomato", col_4= "deepskyblue", height_sum1 = 26, width_sum1 = 20,height_sum2 = 20, width_sum2 = 27) {
  if (!file.exists(folder)){
    dir.create(folder)
  }
  cnv_input <- fread(file = clean_cnv)

  #check input data completeness
  if(!("Length" %in% colnames(cnv_input))){
    print("Input data have no column of Length, it will be generated automatically")
    cnv_input <- cnv_input %>%
                 mutate(Length = End - Start + 1)
  }

  #divide the length group of CNV
  len_int <- as.vector(length_group)
  len_int_bp <- len_int * 1000000 #convert to Mb
  cnv_input <- cnv_input %>%
               mutate(group = case_when(Length <= len_int_bp[1] ~ paste0("<", len_int[1], "Mb"),
                                        Length > len_int_bp[1] & Length <= len_int_bp[2] ~ paste0(len_int[1], "-", len_int[2], "Mb"),
                                        Length > len_int_bp[2] & Length <= len_int_bp[3] ~ paste0(len_int[2], "-", len_int[3], "Mb"),
                                        Length > len_int_bp[3] & Length <= len_int_bp[4] ~ paste0(len_int[3], "-", len_int[4], "Mb"),
                                        Length > len_int_bp[4] & Length <= len_int_bp[5] ~ paste0(len_int[4], "-", len_int[5], "Mb"),
                                        Length > len_int_bp[5] & Length <= len_int_bp[6] ~ paste0(len_int[5], "-", len_int[6], "Mb"),
                                        Length > len_int_bp[6] & Length <= len_int_bp[7] ~ paste0(len_int[6], "-", len_int[7], "Mb"),
                                        Length > len_int_bp[7] & Length <= len_int_bp[8] ~ paste0(len_int[7], "-", len_int[8], "Mb"),
                                        Length > len_int_bp[8] & Length <= len_int_bp[9] ~ paste0(len_int[8], "-", len_int[9], "Mb"),
                                        Length > len_int_bp[9] & Length <= len_int_bp[10] ~ paste0(len_int[9], "-", len_int[10], "Mb"),
                                        Length > len_int_bp[10] & Length <= len_int_bp[11] ~ paste0(len_int[10], "-", len_int[11], "Mb"),
                                        Length > len_int_bp[11] ~ paste0(">", len_int[11], "Mb")))


  # cnv_input$group <- NA  #add a new column to make group of length
  # cnv_input$group[cnv_input$Length <= 50000] <- "1-50kb"
  # cnv_input$group[cnv_input$Length > 50000 & cnv_input$Length <= 100000] <- "50-100kb"
  # cnv_input$group[cnv_input$Length > 100000 & cnv_input$Length <= 200000] <- "100-200kb"
  # cnv_input$group[cnv_input$Length > 200000 & cnv_input$Length <= 300000] <- "200-300kb"
  # cnv_input$group[cnv_input$Length > 300000 & cnv_input$Length <= 400000] <- "300-400kb"
  # cnv_input$group[cnv_input$Length > 400000 & cnv_input$Length <= 500000] <- "400-500kb"
  # cnv_input$group[cnv_input$Length > 500000 & cnv_input$Length <= 600000] <- "500-600kb"
  # cnv_input$group[cnv_input$Length > 600000 & cnv_input$Length <= 700000] <- "600-700kb"
  # cnv_input$group[cnv_input$Length > 700000 & cnv_input$Length <= 1000000] <- "700-1000kb"
  # cnv_input$group[cnv_input$Length > 1000000 & cnv_input$Length <= 2000000] <- "1-2Mbp"
  # cnv_input$group[cnv_input$Length > 2000000 & cnv_input$Length <= 5000000] <- "2-5Mbp"
  # cnv_input$group[cnv_input$Length > 5000000] <- ">5Mbp"

  #customize color
  color_copy <- c("0" = col_0,
                  "1" = col_1,
                  "2" = col_2,
                  "3" = col_3,
                  "4" = col_4)
  #plot CNV distribution
  chr_freq <- cnv_input %>%
              group_by(Chr) %>%
              count(CNV_Value, name = "Freq")
  png(res = 350, filename = paste0(folder, "/f1_cnv_distribution.png"), height = 3500, width = 4000, bg = "transparent")
  adj_y <- max(cnv_input$Length)/max(chr_freq$Freq)
  cnv_distri <- ggplot(cnv_input) +
    geom_boxplot(aes(x = as.factor(Chr), y = Length/adj_y,  fill = as.factor(CNV_Value)), outlier.shape = NA) +
    scale_fill_manual(values = color_copy) +
    geom_line(aes(group = 1, x = as.factor(Chr), color = as.factor(CNV_Value)), stat = 'count') +
    scale_color_manual(values = color_copy) +
    geom_text(aes(x = as.factor(Chr), label = ..count..), stat = 'count') +
    scale_y_continuous(name = "Frequency of CNV (N)", sec.axis = sec_axis(~.*adj_y / 1000, name = "Length of CNV (Kb)")) +
    theme_classic() +
    theme(legend.position = "top",
          strip.text.x = element_blank()) +
    labs(x = "Chromosome", fill = "CNV Value", color = "Frequency of CNV") +
    facet_wrap(~CNV_Value, ncol = 1, scales = "free")

  print(cnv_distri)
  dev.off()

  if(file.exists(paste0(folder, "/f1_cnv_distribution.png"))){
    print("CNV distribution plot was saved in working directory.")
  } else {
    print("Task faild, please check your input data format and paramter used!")
  }


  #fig.1a Length group Vs CNV Count
  cnv_count <- cnv_input %>%
               group_by(CNV_Value) %>%
               count(group) #summarize the number of CNVs on each state each length group
  #cnv_count$group <- factor(cnv_count$group, levels = c("1-50kb", "50-100kb", "100-200kb", "200-300kb",
  #                                                      "300-400kb", "400-500kb", "500-600kb",
  #                                                      "600-700kb", "700-1000kb", "1-2Mbp", "2-5Mbp"))

  png(filename = paste0(folder, "/f2_length_group.png"), res = 350, width = 3500, height = 2000)
  f1a_lengthgroup <- ggplot(cnv_count, aes(fill = as.factor(CNV_Value), x = group, y = n)) +
    scale_fill_manual(values = color_copy) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 20, vjust = 0.8),
          axis.title.x = element_blank(),
          legend.key.size  = unit(0.5,"line"),
          legend.position = c(0.8, 0.9),
          legend.title = element_text(size = 6),
          legend.text  = element_text(size = 6),
          legend.direction = "horizontal") +
    labs(y = "Number of CNV", fill = "CNV Value") # col change title of legend
  print(f1a_lengthgroup)
  dev.off()

  if(file.exists(paste0(folder, "/f2_length_group.png"))){
    print("CNV length group plot was saved in working directory.")
  } else {
    print("Task faild, please check your input data format and paramter used!")
  }

  png(filename = paste0(folder, "/f3_lengthbox.png"), res = 350, width = 15, height = 14, units = "cm")
  #fig.1b CNV Vs Length
  f1b_lengthbox <- ggplot(cnv_input, aes(x = as.factor(CNV_Value), y = Length, color = as.factor(CNV_Value))) +
    geom_boxplot() +
    scale_color_manual(values = color_copy) +
    scale_y_continuous(labels = unit_format(unit = "" ,scale = 1e-3)) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9),
          legend.key.size  = unit(0.5,"line"),
          legend.title = element_text(size = 6),
          legend.text  = element_text(size = 6)) +
    labs(y = "Length (kb)", x = "Type of copy", col = "CNV Value")
  print(f1b_lengthbox)
  dev.off()
  if(file.exists(paste0(folder, "/f3_lengthbox.png"))){
    print("Box plot of CNV length was saved in working directory.")
  } else {
    print("Task faild, please check your input data format and paramter used!")
  }

  #fig.1c the count of CNV number distribute on each Chr
  cnv_chr <- cnv_input %>% group_by(CNV_Value) %>% count(Chr)
  cnv_chr$Chr <- factor(cnv_chr$Chr, levels = c(1:max(cnv_chr$Chr)))
  png(filename = paste0(folder, "/f4_chr.png"), res = 350, width = 15, height = 14, units = "cm")
  f1c_chr<- ggplot(cnv_chr, aes(x = Chr, y = n, group = CNV_Value)) +
    geom_line(linetype = "dashed", aes(color = as.factor(CNV_Value))) +
    geom_point(aes(color = as.factor(CNV_Value))) +
    scale_color_manual(values = color_copy) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9),
          legend.key.size  = unit(0.5,"line"),
          legend.title = element_text(size = 6),
          legend.text  = element_text(size = 6)) +
    labs(y = "Number of CNV", x = "Chr", col = "CNV Value") # col change title of legend
  print(f1c_chr)
  dev.off()
  if(file.exists(paste0(folder, "/f4_chr.png"))){
    print("Distribution of CNV type on chromosome plot was saved in working directory.")
  } else {
    print("Task faild, please check your input data format and paramter used!")
  }

  #plot individual CNV number
  png(filename = paste0(folder, "/f5_individual.png"), res = 350, width = 15, height = 12, units = "cm")
  cnv_indiv <- cnv_input %>%
               group_by(Sample_ID) %>%
               count(Sample_ID, name = "n_CNV") %>%
               arrange(-n_CNV)
  f1d_indiv <- ggplot(cnv_indiv, aes(x = as.factor(n_CNV))) +
    geom_bar(stat = "count",fill = col_0, color = "black") +
    #geom_text(stat = "count", aes(label = ..count..), vjust = -0.12) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "Number of CNVs per Sample", y = "Number of Samples")
  print(f1d_indiv)
  dev.off()
  fwrite(cnv_indiv, file = paste0(folder, "/individual_cnv_count.txt"), sep = '\t', quote = FALSE)

  if(file.exists(paste0(folder, "/f5_individual.png"))){
    print("The number of CNV in Individual plot was saved in working directory.")
  } else {
    print("Task faild, please check your input data format and paramter used!")
  }


  if(!missing(plot_sum_1)) {
    a_d <- plot_grid(f1d_indiv, f1a_lengthgroup, labels = c("a", "b"))

    png(res = 350, filename = paste0(folder, "/cnv_summary_plot_1.png"), height = height_sum1, width = width_sum1, units = "cm", bg = "transparent")
    ad_e <- plot_grid(a_d, cnv_distri, ncol = 1, rel_heights = c(1,3.3), labels = c("", "c"))
    print(ad_e)
    dev.off()
    if(file.exists(paste0(folder, "/cnv_summary_plot_1.png"))){
      print("cnv_summary_plot_1 was saved in working directory.")
    } else {
      print("Task faild, please check your input data format and paramter used!")
    }
  }

  if(!missing(plot_sum_2)){
      # merge four figure in one
    png(filename = paste0(folder, "/cnv_summaray_plot_2.png"), res = 350, height = height_sum2, width = width_sum2, units = "cm", bg = "transparent")

    all <-  plot_grid(f1a_lengthgroup, f1b_lengthbox, f1c_chr, f1d_indiv, labels = c("a", "b","c","d"))
    print(all)
    dev.off()
      if(file.exists(paste0(folder, "/cnv_summaray_plot_2.png"))){
        print("cnv_summaray_polt_2 was saved in working directory.")
        print("Task done.")
      } else {
        print("Task faild, please check your input data format and paramter used!")
      }
  }

  else {
    print("Task done. If you want to combine all these plot together, try arguments in function with plot_sum_1 = 1 and plot_sum_2 = 1")
  }
}




