#' Title
#' Visualizing CNV on both population level and individual level
#'
#' @param clean_cnv
#' @param max_chr_length
#' @param chr_id
#' @param chr_length
#' @param start_position
#' @param end_position
#' @param individual_id
#' @param plot_gene
#' @param plot_title
#'
#' @import
#' data.table,
#' dplyr,
#' ggplot2,
#' tidyr
#'
#' @return
#' CNV distribution plot
#'
#' @export
#'
#' @examples
#'
require(ggplot2, quietly = TRUE)
require(data.table, quietly = TRUE)
require(ggrepel, quietly = TRUE)
cnv_visual <- function(clean_cnv, max_chr_length = NULL, chr_id = NULL, chr_length = NULL, start_position = NULL, end_position = NULL, individual_id = NULL, plot_gene = NULL, plot_title = NULL, report_id = NULL, pedigree = NULL) {
  #myAgr <- formals(cnv_visual)
  #prepare for population data
  cnv <- fread(file = clean_cnv)

  #check if the input is a Plink results
  #convert to the stardards format if it is
  plink_roh_names <- c("FID", "IID", "PHE", "CHR", "SNP1", "SNP2", "POS1", "POS2",
                       "KB", "NSNP", "DENSITY", "PHOM", "PHET")
  if(all(colnames(cnv) == plink_roh_names)){
    print("ROH with input file in PLINK format was detected, coverting to HandyCNV standard formats")
    colnames(cnv) <- c("FID", "Sample_ID", "PHE", "Chr", "Start_SNP", "End_SNP",
                       "Start", "End", "Length", "Num_SNP", "DENSITY", "PHOM", "PHET")
    handycnv_name <- c("Sample_ID",	"Chr", "Start", "End", "Num_SNP",	"Length", "Start_SNP",	"End_SNP")
    cnv <- cnv %>% select(handycnv_name) %>% mutate(CNV_Value = 2)
  }

  id_coord <- data.frame("Sample_ID" = sort(unique(cnv$Sample_ID))) #extract unique ID prepare coordinate
  id_coord$Order <- seq(1,nrow(id_coord),1)
  id_coord$x <- 160
  id_coord$y <- (id_coord$Order-1)*5 + 1
  cnv_coord <- merge(cnv, id_coord, all.x = TRUE, sort = FALSE) #prepare original data
  #set length of chr
  chr_length_ars <- data.frame("Chr" = c(29:1), "Length" = c( 51.098607, 45.94015, 45.612108, 51.992305,
                                                                42.350435, 62.317253, 52.498615, 60.773035, 69.862954,
                                                                71.974595, 63.449741, 65.820629, 73.167244, 81.013979,
                                                                85.00778, 82.403003, 83.472345, 87.216183, 106.982474,
                                                                103.308737, 105.454467, 113.31977, 110.682743, 117.80634,
                                                                120.089316, 120.000601, 121.005158, 136.231102, 158.53411))
  #names(chr_length_ars) <- c("Chr", "Length")
  chr_length_ars <- chr_length_ars[order(chr_length_ars$Chr),]

  if(is.null(max_chr_length) == "FALSE") {
  #1.plot all CNV on all chromosome on population level
  cnv_pop <- cnv_coord
  png(res = 300, filename = "1_chr_all.png", width = 3000, height = 20000)
  popu_plot <- ggplot(cnv_pop, aes(xmin = Start/1000000, xmax = End/1000000, ymin = (Order-1)*5, ymax = (Order-1)*5 + 3)) +
    geom_rect(aes(fill = as.factor(CNV_Value))) +
    geom_text(aes(x,y, label = Sample_ID), size = 1.5) +
    theme_classic() +
    scale_y_continuous(labels = NULL) +
    scale_x_continuous(breaks = seq(0, max_chr_length +10, by = 10)) +
    facet_wrap(~as.numeric(Chr), ncol = 1) +
    labs(x = "Physical Position", y = "Individual ID", title = "CNV Distribution on Population Level", fill = "CNV_Num")
  print(popu_plot)
  dev.off()
  print("Task done, plot was stored in working directory.")
  }

  else if(is.null(chr_id) == "FALSE" & is.null(start_position) & is.null(plot_gene)){
  #2.plot for specific chromosome
  cnv_chr <- cnv[cnv$Chr == chr_id, ]
  id_coord <- data.frame("Sample_ID" = sort(unique(cnv_chr$Sample_ID))) #extract unique ID prepare coordinate
  id_coord$Order <- seq(1,nrow(id_coord),1)
  id_coord$x <- chr_length_ars[chr_id, 2]
  id_coord$y <- (id_coord$Order-1)*5 + 1
  cnv_chr <- merge(cnv_chr, id_coord, all.x = TRUE, sort = FALSE) #prepare original data

  chr_name <- paste0("Chr", chr_id, ".png")
  id_number <- nrow(id_coord)
  chr_title <- paste0("CNV on Chromosome ", chr_id, " with ", id_number, " Indiciduals")
  png(res = 300, filename = chr_name, width = 3500, height = 2000)
  #png(res = 300, filename = "10_chr.png", width = 3500, height = 2000)
  chr_plot <- ggplot(cnv_chr, aes(xmin = Start/1000000, xmax = End/1000000, ymin = (Order-1)*5, ymax = (Order-1)*5 + 3)) +
    geom_rect(aes(fill = as.factor(CNV_Value))) +
    #geom_text(aes(x,y, label = Sample_ID), size = 1.5, check_overlap = TRUE) +
    theme_bw() +
    scale_y_continuous(labels = NULL) +
    scale_x_continuous(breaks = seq(0, chr_length_ars[chr_id, 2], by = 5), limits = c(0, chr_length_ars[chr_id, 2])) +
    labs(x = "Physical Position (Mb)", y ="Individual ID", title = chr_title, fill = "CNV_Num")
  print(chr_plot)
  dev.off()
  print("Task done, plot was stored in working directory.")
  }

  #else if(is.null(start_position & end_position) == "FALSE")
  else if(is.null(start_position) == "FALSE" & is.null(end_position) == "FALSE" & is.null(plot_gene))
    {
  #3.zoom into specific region
  cnv_chr <- cnv[cnv$Chr == chr_id, ]
  cnv_chr_zoom <- filter(cnv_chr, Start >= start_position * 1000000 -1 & End <= end_position * 1000000 + 1)
  id_coord <- data.frame("Sample_ID" = sort(unique(cnv_chr_zoom$Sample_ID))) #extract unique ID prepare coordinate
  id_coord$Order <- seq(1,nrow(id_coord),1)
  id_coord$x <- chr_length_ars[chr_id, 2]
  id_coord$y <- (id_coord$Order-1)*5 + 1
  cnv_chr_zoom <- merge(cnv_chr_zoom, id_coord, all.x = TRUE, sort = FALSE) #prepare original data
  cnv_chr_zoom$zoom_x <- end_position
  zoom_name <- paste0("Chr", chr_id, "_",start_position,"-",end_position, "Mb", ".png")
  id_number <- nrow(id_coord)
  zoom_title <- paste0("CNV on Chromosome ", chr_id, ": ",start_position," - ",end_position, " Mb", " with ", id_number," Individual" ," - ", plot_title)
  png(res = 300, filename = zoom_name, width = 3500, height = 2000)
  zoom_plot <- ggplot(cnv_chr_zoom, aes(xmin = Start/1000000, xmax = End/1000000, ymin = (Order-1)*5, ymax = (Order-1)*5 + 3)) +
    geom_rect(aes(fill = as.factor(CNV_Value))) +
    #geom_text(aes(zoom_x, y, label = Sample_ID), size = 2.5) +
    theme_bw() +
    scale_y_continuous(labels = NULL) +
    scale_x_continuous(breaks = seq(start_position, end_position)) +
    labs(x = "Physical Position (Mb)", y ="Individual ID", title = zoom_title, fill = "CNV_Num")
  print(zoom_plot)
  dev.off()
  print("Task done, plot was stored in working directory.")
  if(!is.null(report_id)) {
    indiv_id <- cnv_chr_zoom$Sample_ID
    print("Individual ID in this CNVRs as following: ")
    print(indiv_id)
    #assign("indiv_id", value = cnv_chr_zoom, envir = .GlobalEnv)
    indiv <- cnv_chr_zoom
    #cnv_visual(clean_cnv = "clean_cnv/penncnv_clean.cnv", chr_id = 5, start_position = 93.6, end_position = 93.7, report_id = 1)
    pedb <- fread(pedigree)
    if (exists("pedb")) {
      print("Pedigree was read in.")
    }
    indiv_id_ped <- merge(indiv, pedb, by.x = "Sample_ID", by.y = "Chip_ID", all.x = TRUE)
    #print(colnames(indiv_id_ped))

    if ("Herd" %in% colnames(indiv_id_ped)){
      herd_cnv <- ggplot(indiv_id_ped, aes(x = as.factor(Herd), fill = as.factor(CNV_Value))) +
        geom_bar() +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 45), legend.position = "top") +
        labs(x = "Name of Herd", y = "Number of CNV", fill = "Copy of CNV")
    }

    if ("Sire_Source" %in% colnames(indiv_id_ped)){
      source_cnv <- ggplot(indiv_id_ped, aes(x = Sire_Source, fill = as.factor(CNV_Value))) +
        geom_bar() +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 45), legend.position = "top") +
        labs(x = "Bull Source", y = "Number of CNV", fill = "Copy of CNV")
    }

    if ("Sire_ID" %in% colnames(indiv_id_ped)){
      sire_cnv <- ggplot(indiv_id_ped, aes(x = Sire_ID, fill = as.factor(CNV_Value))) +
        geom_bar()+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 45), legend.position = "top") +
        labs(x = "Sire ID", y = "Number of CNV", fill = "Copy of CNV")
    }

    png(filename = paste0(zoom_name, "source.png"), res = 300, bg = "transparent", height = 2000, width = 2500)
    if (exists("herd_cnv") & exists("source_cnv") & exists("sire_cnv")){
      cnv_source <- plot_grid(herd_cnv, source_cnv, sire_cnv, ncol = 1)
      print(cnv_source)
      dev.off()
    } else if (exists("herd_cnv") & exists("source_cnv")) {
      cnv_source <- plot_grid(herd_cnv, source_cnv, sire_cnv, ncol = 1)
      print(cnv_source)
      dev.off()
    } else if (exists("herd_cnv") & exists("sire_cnv")) {
      cnv_source <- plot_grid(herd_cnv, sire_cnv, ncol = 1)
      print(cnv_source)
      dev.off()
    } else if (exists("source_cnv") & exists("sire_cnv")) {
      cnv_source <- plot_grid(source_cnv, sire_cnv, ncol = 1)
      print(cnv_source)
      dev.off()
    } else {
      print(sire_cnv)
      dev.off()
    }
  }
  }
  else if (is.null(start_position) == "FALSE" & is.null(start_position) == "FALSE" & is.null(end_position) == "FALSE" & is.null(plot_gene) == "FALSE") {
    cnv_chr <- cnv[cnv$Chr == chr_id, ]
    cnv_chr_zoom <- filter(cnv_chr, CNV_Start >= start_position * 1000000 -1 & CNV_End <= end_position * 1000000 + 1)
    id_coord <- data.frame("Sample_ID" = sort(unique(cnv_chr_zoom$Sample_ID))) #extract unique ID prepare coordinate
    try(id_coord$Order <- seq(1, nrow(id_coord),1), silent = FALSE)
    id_coord$x <- chr_length_ars[chr_id, 2]
    id_coord$y <- (id_coord$Order-1)*5 + 1
    cnv_chr_zoom <- merge(cnv_chr_zoom, id_coord, all.x = TRUE, sort = FALSE) #prepare original data
    cnv_chr_zoom$zoom_x <- end_position
    cnv_chr_zoom$gene_order <- max(cnv_chr_zoom$Order) + 3

    #gene_freq <- data.table(gene_freq) #convert it to data.table to set key
    #setkey(x = gene_freq, name2)
    gene_coord <- group_by(cnv_chr_zoom, name2) %>% slice(1) # generate gene data
    gene_coord$CNV_Start <- gene_coord$g_Start
    gene_coord$Order <- gene_coord$gene_order
    try(gene_coord$CNV_Value <- "5", silent = FALSE)

    gene_freq <- cnv %>% group_by(name2) %>% count(name2, name = "Frequent", sort = TRUE) #gene freq
    gene_coord_freq <- merge(gene_coord, gene_freq)
    gene_coord_freq <- gene_coord_freq[c(gene_coord_freq$Frequent >= 5), ]

    zoom_name <- paste0("Chr", chr_id, "_",start_position,"-",end_position, "Mb", ".png")
    id_number <- nrow(id_coord)
    zoom_title <- paste0("CNV on Chromosome ", chr_id, ": ",start_position," - ",end_position, " Mb", " with ", id_number," Individual" ," - ", plot_title)
    png(res = 300, filename = zoom_name, width = 3500, height = 2000)
    zoom_plot <- ggplot() +
      geom_rect(data = cnv_chr_zoom, aes(xmin = CNV_Start/1000000, xmax = CNV_End/1000000, ymin = (Order-1)*5, ymax = (Order-1)*5 + 3, fill = as.factor(CNV_Value))) +
      geom_rect(data = gene_coord, aes(xmin = g_Start/1000000, xmax = g_End/1000000, ymin = (Order-1)*5, ymax = (Order-1)*5 + 3), fill = "black") +
      geom_text_repel(data = gene_coord_freq, aes(x = g_Start/1000000, y = (Order-1)*5 + 10, label = name2)) +
      geom_hline(yintercept = (max(cnv_chr_zoom$Order) + 2)*5 - 2, linetype = "dashed") +
      #geom_text(aes(zoom_x, y, label = Sample_ID), size = 2.5) +
      #scale_color_manual(values = c("#F8766D", "#A3A500", "#00B0F6", "#E76BF3", "black")) +
      theme_bw() +
      scale_y_continuous(labels = NULL) +
      scale_x_continuous(breaks = seq(round(start_position,2), round(end_position,2), by = 0.2)) +
      labs(x = "Physical Position (Mb)", y ="Individual ID", title = zoom_title, fill = "CNV_Num")
    print(zoom_plot)
    dev.off()
    print("Task done, plot was stored in working directory.")
  }

  else if (is.null(individual_id) == "FALSE")
    {
  #plot on individual level
    #cnv_indiv <- cnv_p[which(cnv_p$Sample_ID == "204806050057_R01C01")]
  cnv_indiv <- cnv[which(cnv$Sample_ID == individual_id),]
  chr_coord <- data.frame("Chr" = seq(1,29,1))
  chr_coord$x <- 160 #adjust x for geom_text
  chr_coord$y <- (chr_coord$Chr-1)*5 + 1 #adjust y for geom_text
  cnv_indiv_coord <- merge(cnv_indiv, chr_coord, all.x = TRUE, sort = FALSE)

  #4.
  indiv_name <- paste0("CNV_of_", individual_id,".png")
  indiv_title <- paste0("CNV Distribution of Individual ", individual_id)
  png(res = 300, filename = indiv_name, width = 3500, height = 2000)
  indiv_plot <- ggplot(cnv_indiv_coord, aes(xmin = Start/1000000, xmax = End/1000000, ymin = (Chr-1)*5, ymax = (Chr-1)*5 + 3)) +
    geom_rect(aes(fill = as.factor(CNV_Value))) +
    geom_text(aes(x, y, label = Chr), size = 4) + theme_bw() +
    scale_y_continuous(breaks = seq(0, 160, by = 10),labels = NULL) +
    scale_x_continuous(breaks = seq(0, 160, by = 10)) +
    labs(x = "Physical Position (Mb)", y ="Autosome", title = indiv_title,  fill= "CNV_Num")
  print(indiv_plot)
  dev.off()
  print("Task done, plot was stored in working directory.")
  #return(cnv_coord)
  }

  else{
    print("Warning: Lack of input parameters!!!
           Warning: Please check input parameters carefully!!!")
  }
}
