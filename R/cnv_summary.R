require(data.table, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(dplyr, quietly = TRUE)
require(scales, quietly = TRUE)
require(reshape2, quietly = TRUE)
require(tidyr, quietly = TRUE)
#3.1.1 Make graphes for CNVs from CNVPartition--------------------------
cnv_summary <- function(clean_cnv, plot_length = NULL, plot_copynumber =NULL, plot_chr = NULL, plot_individual = NULL, plot_merge = NULL) {
  cnv_input <- fread( file = clean_cnv)
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

  if (!is.null(plot_length)){
    #fig.1a Length group Vs CNV Count
    cnv_count <- cnv_input %>% group_by(CNV_Value) %>% count(group) #summarise the number of CNVs on each state each length group
    cnv_count$group <- factor(cnv_count$group, levels = c("1-50kb", "50-100kb", "100-200kb", "200-300kb",
                                                          "300-400kb", "400-500kb", "500-600kb",
                                                          "600-700kb", "700-1000kb", "1-2Mbp", "2-5Mbp"))

    png(filename = "f1a_length.png", res = 300, width = 3500, height = 2000)
    f1a_length <- ggplot(cnv_count, aes(fill = as.factor(CNV_Value), x = group, y = n)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 20, vjust = 0.8)) +
      labs(y = "Number of CNV", x = "Group of Length", fill = "CNV Value") # col change title of legend
    print(f1a_length)
    dev.off()
    if(file.exists("f1a_length.png")){
      print("Task done, plot saved in working directory.")
    } else {
      print("Task faild, please check your input data format and paramter used!")
    }
  }

  else if (!is.null(plot_copynumber)) {
    png(filename = "f1b_cnvnumber.png", res = 300, width = 3500, height = 2000)
    #fig.1b CNV Vs Length
    f1b_cnvnumber <- ggplot(cnv_input, aes(x = as.factor(CNV_Value), y = Length, color = as.factor(CNV_Value))) +
      geom_boxplot() +
      scale_y_continuous(labels = unit_format(unit = "" ,scale = 1e-3)) +
      theme_classic() +
      labs(y = "Length (kb)", x = "Number of copy", col = "CNV Value")
    print(f1b_cnvnumber)
    dev.off()
    if(file.exists("f1b_cnvnumber.png")){
      print("Task done, plot saved in working directory.")
    } else {
      print("Task faild, please check your input data format and paramter used!")
    }
  }

  else if (!is.null(plot_chr)) {
    #fig.1c the count of CNV number distribute on each Chr
    cnv_chr <- cnv_input %>% group_by(CNV_Value) %>% count(Chr)
    cnv_chr$Chr <- factor(cnv_chr$Chr, levels = c(1:29))
    png(filename = "f1c_chr.png", res = 300, width = 3500, height = 2000)
    f1c_chr<- ggplot(cnv_chr, aes(x = Chr, y = n, group = CNV_Value)) +
      geom_line(linetype = "dashed", aes(color = as.factor(CNV_Value))) +
      geom_point(aes(color = as.factor(CNV_Value))) +
      theme_classic() +
      labs(y = "Number of CNV", x = "Chr", col = "CNV Value") # col change title of legend
    print(f1c_chr)
    dev.off()
    if(file.exists("f1c_chr.png")){
      print("Task done, plot saved in working directory.")
    } else {
      print("Task faild, please check your input data format and paramter used!")
    }
  }

  else if(!is.null(plot_individual)) {
      #fig.1d the count of CNVs on individual level
      cnv_indiv <- cnv_input %>% group_by(Sample_ID) %>% count(Sample_ID) %>% group_by(n) %>% count(n)
      # n is the number of CNVs for each individual, nn is number of individuals on each CNV frequent group
      png(filename = "f1d_indiv.png", res = 300, width = 3500, height = 2000)
      f1d_indiv <- ggplot(cnv_indiv, aes(x = n, y = nn)) +
        geom_bar(stat = "identity",fill = "grey", color = "black") +
        geom_text(aes(label = nn), vjust = -0.3) +
        theme_classic() +
        labs(x = "Number of CNV", y = "Number of Individual")
      print(f1d_indiv)
      dev.off()
      if(file.exists("f1d_indiv.png")){
        print("Task done, plot saved in working directory.")
      } else {
        print("Task faild, please check your input data format and paramter used!")
      }
    }

  else if(!is.null(plot_merge)){
      # set a big graph merge four figture in one
      png(
        filename = "fig_1_part.png",
        type = "cairo", #
        res = 300, # 300ppi
        width = 3500, height = 2100,
        bg = "transparent" #
      )
      plot_grid(f1a_part, f1b_part, f1c_part, f1d_part, labels = c("a", "b","c","d"))
      dev.off()
      if(file.exists(f1d_indiv)){
        print("Task done, plot saved in working directory.")
      } else {
        print("Task faild, please check your input data format and paramter used!")
      }
    }

  else {
    print("Warning: Insufficient input parameters, please read help file for this fuction!")
  }

}



