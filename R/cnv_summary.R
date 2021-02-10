#' Title summarise CNVs
#' Make graphes for CNVs
#'
#' @param clean_cnv
#' @param plot_sum_1
#' @param plot_sum_2
#'
#' @import ggplot2 dplyr scales reshape2 tidyr
#' @importFrom data.table fread fwrite
#'
#' @return
#' @export cnv_summary_plot
#'
#' @examples
cnv_summary_plot <- function(clean_cnv, plot_sum_1 = NULL, plot_sum_2 = NULL) {
  if (!file.exists("cnv_summary_plot")){
    dir.create("cnv_summary_plot")
  }
  cnv_input <- fread(file = clean_cnv)
  cnv_input$group <- NA  #add a new column to make group of length
  cnv_input$group[cnv_input$Length <= 50000] <- "1-50kb"
  cnv_input$group[cnv_input$Length > 50000 & cnv_input$Length <= 100000] <- "50-100kb"
  cnv_input$group[cnv_input$Length > 100000 & cnv_input$Length <= 200000] <- "100-200kb"
  cnv_input$group[cnv_input$Length > 200000 & cnv_input$Length <= 300000] <- "200-300kb"
  cnv_input$group[cnv_input$Length > 300000 & cnv_input$Length <= 400000] <- "300-400kb"
  cnv_input$group[cnv_input$Length > 400000 & cnv_input$Length <= 500000] <- "400-500kb"
  cnv_input$group[cnv_input$Length > 500000 & cnv_input$Length <= 600000] <- "500-600kb"
  cnv_input$group[cnv_input$Length > 600000 & cnv_input$Length <= 700000] <- "600-700kb"
  cnv_input$group[cnv_input$Length > 700000 & cnv_input$Length <= 1000000] <- "700-1000kb"
  cnv_input$group[cnv_input$Length > 1000000 & cnv_input$Length <= 2000000] <- "1-2Mbp"
  cnv_input$group[cnv_input$Length > 2000000 & cnv_input$Length <= 5000000] <- "2-5Mbp"
  cnv_input$group[cnv_input$Length > 5000000] <- ">5Mbp"

  #plot CNV distribution
  chr_freq <- cnv_input %>% group_by(Chr) %>% count(CNV_Value, name = "Freq")
  png(res = 300, filename = "cnv_summary_plot/cnv_distribution.png", height = 3500, width = 4000, bg = "transparent")
  adj_y <- max(cnv_input$Length)/max(chr_freq$Freq)
  cnv_distri <- ggplot(cnv_input) +
    geom_boxplot(aes(x = as.factor(Chr), y = Length/adj_y,  fill = as.factor(CNV_Value)), outlier.shape = NA) +
    geom_line(aes(group = 1, x = as.factor(Chr), color = as.factor(CNV_Value)), stat = 'count') +
    geom_text(aes(x = as.factor(Chr), label = ..count..), stat = 'count') +
    scale_y_continuous(name = "Frequency of CNV (N)", sec.axis = sec_axis(~.*adj_y / 1000, name = "Length of CNV (Kb)")) +
    theme_classic() +
    theme(legend.position = "top", strip.text.x = element_blank()) +
    labs(x = "Chromosome", fill = "CNV Value", color = "Frequency of CNV") +
    facet_wrap(~CNV_Value, ncol = 1, scales = "free")

  print(cnv_distri)
  dev.off()

  if(file.exists("cnv_summary_plot/cnv_distribution.png")){
    print("CNV distribution plot was saved in working directory.")
  } else {
    print("Task faild, please check your input data format and paramter used!")
  }


  #fig.1a Length group Vs CNV Count
  cnv_count <- cnv_input %>% group_by(CNV_Value) %>% count(group) #summarise the number of CNVs on each state each length group
  cnv_count$group <- factor(cnv_count$group, levels = c("1-50kb", "50-100kb", "100-200kb", "200-300kb",
                                                        "300-400kb", "400-500kb", "500-600kb",
                                                        "600-700kb", "700-1000kb", "1-2Mbp", "2-5Mbp"))

  png(filename = "cnv_summary_plot/f1a_lengthgroup.png", res = 300, width = 3500, height = 2000)
  f1a_lengthgroup <- ggplot(cnv_count, aes(fill = as.factor(CNV_Value), x = group, y = n)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 20, vjust = 0.8), axis.title.x = element_blank(), legend.position = c(0.9, 0.9)) +
    labs(y = "Number of CNV", fill = "CNV Value") # col change title of legend
  print(f1a_lengthgroup)
  dev.off()

  if(file.exists("cnv_summary_plot/f1a_lengthgroup.png")){
    print("CNV length group plot was saved in working directory.")
  } else {
    print("Task faild, please check your input data format and paramter used!")
  }

  png(filename = "cnv_summary_plot/f1b_lengthbox.png", res = 300, width = 3500, height = 2000)
  #fig.1b CNV Vs Length
  f1b_lengthbox <- ggplot(cnv_input, aes(x = as.factor(CNV_Value), y = Length, color = as.factor(CNV_Value))) +
    geom_boxplot() +
    scale_y_continuous(labels = unit_format(unit = "" ,scale = 1e-3)) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9)) +
    labs(y = "Length (kb)", x = "Type of copy", col = "CNV Value")
  print(f1b_lengthbox)
  dev.off()
  if(file.exists("cnv_summary_plot/f1b_lengthbox.png")){
    print("Box plot of CNV length was saved in working directory.")
  } else {
    print("Task faild, please check your input data format and paramter used!")
  }

  #fig.1c the count of CNV number distribute on each Chr
  cnv_chr <- cnv_input %>% group_by(CNV_Value) %>% count(Chr)
  cnv_chr$Chr <- factor(cnv_chr$Chr, levels = c(1:29))
  png(filename = "cnv_summary_plot/f1c_chr.png", res = 300, width = 3500, height = 2000)
  f1c_chr<- ggplot(cnv_chr, aes(x = Chr, y = n, group = CNV_Value)) +
    geom_line(linetype = "dashed", aes(color = as.factor(CNV_Value))) +
    geom_point(aes(color = as.factor(CNV_Value))) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9)) +
    labs(y = "Number of CNV", x = "Chr", col = "CNV Value") # col change title of legend
  print(f1c_chr)
  dev.off()
  if(file.exists("cnv_summary_plot/f1c_chr.png")){
    print("Distribution of CNV type on chromosome plot was saved in working directory.")
  } else {
    print("Task faild, please check your input data format and paramter used!")
  }

  #plot individual CNV number
  png(filename = "cnv_summary_plot/f1d_indiv.png", res = 300, width = 3500, height = 2000)
  cnv_indiv <- cnv_input %>% group_by(Sample_ID) %>% count(Sample_ID, name = "n_CNV")
  f1d_indiv <- ggplot(cnv_indiv, aes(x = as.factor(n_CNV))) +
    geom_bar(stat = "count",fill = "lightblue", color = "black") +
    #geom_text(stat = "count", aes(label = ..count..), vjust = -0.12) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "The Number of CNV Detected in An Individual", y = "Number of Individual")
  print(f1d_indiv)
  dev.off()

  if(file.exists("cnv_summary_plot/f1d_indiv.png")){
    print("The number of CNV in Individual plot was saved in working directory.")
  } else {
    print("Task faild, please check your input data format and paramter used!")
  }


  if(!is.null(plot_sum_1)) {
    a_d <- plot_grid(f1d_indiv, f1a_lengthgroup, labels = c("a", "b"))

    png(res = 300, filename = "cnv_summary_plot/cnv_summary_plot_1.png", height = 4000, width = 3000, bg = "transparent")
    ad_e <- plot_grid(a_d, cnv_distri, ncol = 1, rel_heights = c(1,3.3), labels = c("", "c"))
    print(ad_e)
    dev.off()
    if(file.exists("cnv_summary_plot/cnv_summary_plot_1.png")){
      print("cnv_summary_plot_1 was saved in working directory.")
    } else {
      print("Task faild, please check your input data format and paramter used!")
    }
  }

  if(!is.null(plot_sum_2)){
      # set a big graph merge four figture in one
    png(filename = "cnv_summary_plot/cnv_summaray_plot_2.png",
        res = 300, # 300ppi
        width = 4300, height = 3400,
        bg = "transparent")

    all <-  plot_grid(f1a_lengthgroup, f1b_lengthbox, f1c_chr, f1d_indiv, labels = c("a", "b","c","d"))
    print(all)
    dev.off()
      if(file.exists("cnv_summary_plot/cnv_summaray_plot_2.png")){
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



