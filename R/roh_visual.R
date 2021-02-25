#' Title plot_gene
#' The function used for plotting gene
#' @param refgene reference gene list.
#' @param chr_id chromosome ID, should be the integer only, like 1 or 29
#' @param start the start physical position used in the plot, the unit is Mb
#' @param end the end physical position used in plot, the unit is Mb
#' @param show_name default value is show_name = c(0, 160). accept the vectors only, unit is Mb. for example show_name = c(11.2, 12.4, 15.3, 18.4), means only plot the genes within the given interval
#' @param cnv is support for both CNV and ROH, set cnv argument as TRUE will plot gene for CNV, otherwise will plot gene for ROH
#' interval 11.2-12.4 Mb and 15.3-18.4 Mb, the maximum pairs of interval are three
#'
#' @import ggplot2 dplyr
#' @importFrom data.table fread fwrite setkey foverlaps setDT
#'
#' @return
#' @export plot_gene
#'
#' @examples
plot_gene <- function(refgene = NULL, chr_id, start, end, show_name = c(0,160), cnv = NULL){
  ###
  #read gene
  #if(refgene == "ARS_ens"){
  #  refgene = system.file("extdata", "Demo_data/gene_annotation/ensGene_ars_210202.txt", package = "HandyCNV")
  #  gene <- fread(file = refgene, header = TRUE)
  #} else if(refgene == "ARS_UCSC"){
  #  refgene = system.file("extdata", "Demo_data/gene_annotation/refGene_ars1.2.txt", package = "HandyCNV")
  #  gene <- fread(file = refgene, header = FALSE)
  #  names(gene) <- c("bin", "name", "Chr", "strand", "Start", "End",
  #                   "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
  #                   "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  #} else if(refgene == "UMD_UCSC"){
  #  refgene = system.file("extdata", "Demo_data/gene_annotation/refGene_umd3.1.txt", package = "HandyCNV")
  #  gene <- fread(file = refgene, header = FALSE)
  #  names(gene) <- c("bin", "name", "Chr", "strand", "Start", "End",
  #                   "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
  #                   "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  #}else{
  gene <- fread(file = refgene, header = TRUE)
   #gene <- refgene
   # names(gene) <- c("bin", "name", "Chr", "strand", "Start", "End",
    #                 "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
    #                 "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  #}

  #if(missing(gene)){
  #  gene <- fread(file = gene, header = T)
  #} else{
  #  gene <- fread(gene, header = F)
  #  names(gene) <- c("bin", "name", "Chr", "strand", "Start", "End",
  #                   "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
  #                   "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  #}

  gene$Chr <- sub("chr", "", gene$Chr)

  #extract gene list
  gene_sub <- gene %>%
    filter(Chr == chr_id & Start >= start * 1000000 & End <= end * 1000000) %>%
    arrange(Start) %>%
    mutate(y_min = case_when(length(Chr) <= 5 ~ rep(c(1, 1.1),length.out = length(Chr)),
                             length(Chr) > 5 & length(Chr) <= 15 ~ rep(c(1, 1.1, 1.2),length.out = length(Chr)),
                             length(Chr) > 15 & length(Chr) <= 25 ~ rep(c(1, 1.1, 1.2, 1.3),length.out = length(Chr)),
                             length(Chr) > 25 ~ rep(c(1, 1.1, 1.2, 1.3, 1.4),length.out = length(Chr))),
           y_max =  case_when(length(Chr) <= 5 ~ rep(c(1.1, 1.2),length.out = length(Chr)),
                              length(Chr) > 5 & length(Chr) <= 15 ~ rep(c(1.1, 1.2, 1.3),length.out = length(Chr)),
                              length(Chr) > 15 & length(Chr) <= 25 ~ rep(c(1.1, 1.2, 1.3, 1.4),length.out = length(Chr)),
                              length(Chr) > 25 ~ rep(c(1.1, 1.2, 1.3, 1.4, 1.5),length.out = length(Chr))))
  #mutate(y_min = case_when(row_number() %% 2 == 0 ~ "1",
  #                         row_number() %% 2 == 1 ~ "1.2"),
  #       y_max = case_when(row_number() %% 2 == 0 ~ "1.1",
  #                         row_number() %% 2 == 1 ~ "1.3")) #assign y axis by even and odd row number


  if(nrow(gene_sub) == 0){
    print("No gene in this region")
  } else{
    #plot the name of genes
    coord_name <- as.vector(show_name) * 1000000 #read the interval of gene want to present name
    if(length(coord_name) == 2){
      gene_present <- gene_sub%>%
        filter(Start > coord_name[1] & End < coord_name[2])
    } else if (length(coord_name) == 4){
      gene_present <- gene_sub %>%
        filter(Start > coord_name[1] & End < coord_name[2] | Start > coord_name[3] & End < coord_name[4])
    } else{
      gene_present <- gene_sub %>%
        filter(Start > coord_name[1] & End < coord_name[2] | Start > coord_name[3] & End < coord_name[4] | Start > coord_name[5] & End < coord_name[6])
    }

    #check if plot for roh or CNV, the diferrence is when gene figure combine to interval reduce the distance between the middle (from top to bottom) are reduced when set 'cnv'
    if(is.null(cnv)){
      ggplot() +
        geom_rect(data = gene_present, aes(xmin = Start/1000000, xmax = End/1000000, ymin = y_min, ymax = y_max, fill = as.character(name2)), show.legend = F) +
        #{if(nrow(gene_sub) < 50)geom_text_repel(aes(x = Start/1000000, y = y_max, label = name2))} +
        #{if(nrow(gene_present) < 50)geom_text_repel(data = gene_present, aes(x = Start/1000000, y = y_max, label = name2))} +
        geom_text_repel(data = gene_present, aes(x = Start/1000000, y = y_max, label = name2), size = 2.5) +
        scale_x_continuous(limits = c(start, end)) +
        #geom_vline(xintercept = coord_name/1000000, linetype = "dashed") +
        theme_bw() +
        theme(axis.text.y = element_blank(), axis.title.x = element_blank()) +
        labs(y = "Gene")
    } else {
      ggplot() +
        geom_rect(data = gene_present, aes(xmin = Start/1000000, xmax = End/1000000, ymin = y_min, ymax = y_max), fill = "gray50", show.legend = F) +
        #{if(nrow(gene_sub) < 50)geom_text_repel(aes(x = Start/1000000, y = y_max, label = name2))} +
        #{if(nrow(gene_present) < 50)geom_text_repel(data = gene_present, aes(x = Start/1000000, y = y_max, label = name2))} +
        geom_text_repel(data = gene_present, aes(x = Start/1000000, y = y_max, label = name2), size = 2.5) +
        scale_x_continuous(limits = c(start, end), labels = NULL) +
        #scale_x_continuous(limits = c(start, end)) +
        #geom_vline(xintercept = coord_name/1000000, linetype = "dashed") +
        theme_bw() +
        theme(axis.text.y = element_blank(), axis.title.x = element_blank(),
              plot.margin=unit(c(1,1,0,1), "cm"),
              axis.ticks.y = element_blank()) +
        labs(y = "Gene")
    }
  }
}


#' Title
#' Visualizing ROH on both population level and individual level
#'
#' @param chr_id the number of chromosome want to plot
#' @param chr_length set the length for the chromosome
#' @param start_position decimal digit, default unit is 'Mb'. such as 23.2112
#' @param end_position decimal digit, default unit is 'Mb'. such as 23.2112
#' @param individual_id reporting ID of individuals while plotting
#' @param plot_title add title in the plot
#' @param width_1 number to set the width of final plot size, unit is 'cm'
#' @param height_1 number to set the height of final plot size, unit is 'cm'
#' @param clean_roh support roh results from Plink or generated by cnv_clean form CNVPartition results
#' @param max_chr the maximum number of chromosomes to plot, it used for plot all chromosomes at once
#' @param gene if true, will plot gene above roh plot
#' @param report_id report the sample ID while plotting
#' @param pedigree pedigree list, require at least three columns, Sample_Id, Sire and Dam
#' @param show_name default value is show_name = c(0, 160). accept the vectors only, unit is Mb. for example show_name = c(11.2, 12.4, 15.3, 18.4), means only plot the genes within the given interval
#'
#' @import dplyr ggplot2 tidyr
#' @importFrom data.table fread fwrite setkey foverlaps setDT
#' @return
#' ROH distribution plot
#'
#' @export
#'
#' @examples
#'
roh_visual <- function(clean_roh, max_chr = NULL, chr_id = NULL, chr_length = NULL, start_position = NULL, end_position = NULL, individual_id = NULL, gene = NULL, plot_title = NULL, report_id = NULL, pedigree = NULL, show_name =c(0, 160), width_1 = 13, height_1 = 10) {

  #prepare for population data
  roh <- fread(file = clean_roh)

  #check if the input is a Plink results
  #convert to the standards format if it is
  plink_roh_names <- c("FID", "IID", "PHE", "CHR", "SNP1", "SNP2", "POS1", "POS2",
                       "KB", "NSNP", "DENSITY", "PHOM", "PHET")
  if(length(colnames(roh)) == length(plink_roh_names)){
    print("ROH with input file in PLINK format was detected, coverting to HandyCNV standard formats")
    colnames(roh) <- c("FID", "Sample_ID", "PHE", "Chr", "Start_SNP", "End_SNP",
                       "Start", "End", "Length", "Num_SNP", "DENSITY", "PHOM", "PHET")
    handycnv_name <- c("Sample_ID",	"Chr", "Start", "End", "Num_SNP",	"Length", "Start_SNP",	"End_SNP")
    roh <- roh %>% select(handycnv_name)
    roh$Chr <- as.numeric(roh$Chr)
    roh <- dplyr::filter(roh, roh$Chr %in% c(1:29))
  } else{
    print("Preparing plot data...")
    }

  id_coord <- data.frame("Sample_ID" = sort(unique(roh$Sample_ID))) #extract unique ID prepare coordinate
  id_coord$Order <- seq(1,nrow(id_coord),1)
  id_coord$x <- 160
  id_coord$y <- (id_coord$Order-1)*5 + 1
  cnv_coord <- merge(roh, id_coord, all.x = TRUE, sort = FALSE) #prepare original data
  #set length of chr
  chr_length_ars <- data.frame("Chr" = c(29:1), "Length" = c( 51.098607, 45.94015, 45.612108, 51.992305,
                                                              42.350435, 62.317253, 52.498615, 60.773035, 69.862954,
                                                              71.974595, 63.449741, 65.820629, 73.167244, 81.013979,
                                                              85.00778, 82.403003, 83.472345, 87.216183, 106.982474,
                                                              103.308737, 105.454467, 113.31977, 110.682743, 117.80634,
                                                              120.089316, 120.000601, 121.005158, 136.231102, 158.53411))
  #names(chr_length_ars) <- c("Chr", "Length")
  chr_length_ars <- chr_length_ars[order(chr_length_ars$Chr),]

  if(is.null(max_chr) == "FALSE") {
    #1.plot all CNV on all chromosome on population level
    cnv_pop <- cnv_coord
    cnv_pop <- dplyr::filter(cnv_pop, cnv_pop$Chr >=1 & cnv_pop$Chr <= max_chr)
    png(res = 300, filename = "1_chr_all_roh.png", width = 5000, height = 3000)
    popu_plot <- ggplot(cnv_pop, aes(xmin = Start/1000000, xmax = End/1000000, ymin = (Order-1)*5, ymax = (Order-1)*5 + 3)) +
      geom_rect(aes(fill = Length)) +
      scale_fill_gradientn(colours = c("red", "yellow", "blue")) +
      theme_classic() +
      scale_y_continuous(labels = NULL) +
      #scale_x_continuous(breaks = seq(0, max_chr_length +10, by = 10)) +
      scale_x_continuous(labels = NULL) +
      facet_wrap(~as.numeric(Chr), nrow = 1) +
      labs(title = "ROH Distribution on Population Level", fill = "Length")
    print(popu_plot)
    dev.off()
    print("Task done, plot was stored in working directory.")
  }

  else if(is.null(chr_id) == "FALSE" & is.null(start_position) & is.null(gene)){
    #2.plot for specific chromosome
    cnv_chr <- roh[roh$Chr == chr_id, ]
    id_coord <- data.frame("Sample_ID" = sort(unique(cnv_chr$Sample_ID))) #extract unique ID prepare coordinate
    id_coord$Order <- seq(1,nrow(id_coord),1)
    id_coord$x <- chr_length_ars[chr_id, 2]
    id_coord$y <- (id_coord$Order-1)*5 + 1
    cnv_chr <- merge(cnv_chr, id_coord, all.x = TRUE, sort = FALSE) #prepare original data

    chr_name <- paste0("Chr", chr_id, "_roh.png")
    id_number <- nrow(id_coord)
    chr_title <- paste0("ROH on Chromosome ", chr_id, " with ", id_number, " Individuals")
    png(res = 300, filename = chr_name, width = 3500, height = 2000)
    #png(res = 300, filename = "10_chr.png", width = 3500, height = 2000)
    chr_plot <- ggplot(cnv_chr, aes(xmin = Start/1000000, xmax = End/1000000, ymin = (Order-1)*5, ymax = (Order-1)*5 + 3)) +
      geom_rect(aes(fill = Length)) +
      scale_fill_gradientn(colours = c("red", "yellow", "blue")) +
      #geom_text(aes(x,y, label = Sample_ID), size = 1.5, check_overlap = TRUE) +
      theme_bw() +
      scale_y_continuous(labels = NULL) +
      scale_x_continuous(breaks = seq(0, chr_length_ars[chr_id, 2], by = 5), limits = c(0, chr_length_ars[chr_id, 2])) +
      labs(x = "Physical Position (Mb)", y ="Individual ID", title = chr_title, fill = "Length")
    print(chr_plot)
    dev.off()
    print("Task done, plot was stored in working directory.")
  }

  #else if(is.null(start_position & end_position) == "FALSE")
  else if(is.null(start_position) == "FALSE" & is.null(end_position) == "FALSE" & is.null(gene))
  {
    #3.zoom into specific region
    cnv_chr <- roh[roh$Chr == chr_id, ]
    cnv_chr_zoom <- filter(cnv_chr, Start >= start_position * 1000000 -1 & End <= end_position * 1000000 + 1)
    id_coord <- data.frame("Sample_ID" = sort(unique(cnv_chr_zoom$Sample_ID))) #extract unique ID prepare coordinate
    id_coord$Order <- seq(1,nrow(id_coord),1)
    id_coord$x <- chr_length_ars[chr_id, 2]
    id_coord$y <- (id_coord$Order-1)*5 + 1
    cnv_chr_zoom <- merge(cnv_chr_zoom, id_coord, all.x = TRUE, sort = FALSE) #prepare original data
    cnv_chr_zoom$zoom_x <- end_position
    zoom_name <- paste0("Chr", chr_id, "_",start_position,"-",end_position, "Mb", "_roh.png")
    id_number <- nrow(id_coord)
    zoom_title <- paste0("ROH on Chromosome ", chr_id, ": ",start_position," - ",end_position, " Mb", " with ", id_number," Individual" ," - ", plot_title)
    png(res = 300, filename = zoom_name, width = 3500, height = 2000)
    zoom_plot <- ggplot(cnv_chr_zoom, aes(xmin = Start/1000000, xmax = End/1000000, ymin = (Order-1)*5, ymax = (Order-1)*5 + 3)) +
      geom_rect(aes(fill = Length)) +
      scale_fill_gradientn(colours = c("red", "yellow", "blue")) +
      #geom_text(aes(zoom_x, y, label = Sample_ID), size = 2.5) +
      theme_bw() +
      scale_y_continuous(labels = NULL) +
      scale_x_continuous(breaks = seq(start_position, end_position)) +
      labs(x = "Physical Position (Mb)", y ="Individual ID", title = zoom_title, fill = "Length")
    print(zoom_plot)
    dev.off()
    print("Task done, plot was stored in working directory.")
    if(!is.null(report_id)) {
      indiv_id <- unique(cnv_chr_zoom$Sample_ID)
      print("Individual ID in this CNVRs as following: ")
      print(indiv_id)
      return(indiv_id)
      #assign("indiv_id", value = cnv_chr_zoom, envir = .GlobalEnv)
      indiv <- cnv_chr_zoom
      #cnv_visual(clean_cnv = "clean_cnv/penncnv_clean.cnv", chr_id = 5, start_position = 93.6, end_position = 93.7, report_id = 1)
      if(!(is.null(pedigree))){
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
  }
  else if (is.null(start_position) == "FALSE" & is.null(end_position) == "FALSE" & is.null(gene) == "FALSE") {
    #3.zoom into specific region
    coord_name <- as.vector(show_name) # add vline for plot
    cnv_chr <- roh[roh$Chr == chr_id, ]
    cnv_chr_zoom <- filter(cnv_chr, Start >= start_position * 1000000 -1 & End <= end_position * 1000000 + 1)
    id_coord <- data.frame("Sample_ID" = sort(unique(cnv_chr_zoom$Sample_ID))) #extract unique ID prepare coordinate
    id_coord$Order <- seq(1,nrow(id_coord),1)
    id_coord$x <- chr_length_ars[chr_id, 2]
    id_coord$y <- (id_coord$Order-1)*5 + 1
    cnv_chr_zoom <- merge(cnv_chr_zoom, id_coord, all.x = TRUE, sort = FALSE) #prepare original data
    cnv_chr_zoom$zoom_x <- end_position
    zoom_name <- paste0("Chr", chr_id, "_",start_position,"-",end_position, "Mb", "_roh.png")
    id_number <- nrow(id_coord)
    zoom_title <- paste0("Chr", chr_id, ": ",start_position," - ",end_position, " Mb", " with ", id_number," Individual" ," - ", plot_title)
    print(zoom_title)
    zoom_plot <- ggplot(cnv_chr_zoom, aes(xmin = Start/1000000, xmax = End/1000000, ymin = (Order-1)*5, ymax = (Order-1)*5 + 3)) +
      geom_rect(aes(fill = Length), show.legend = F) +
      scale_fill_gradientn(colours = c("red", "yellow", "blue")) +

      #geom_text(aes(zoom_x, y, label = Sample_ID), size = 2.5) +
      theme_bw() +
      theme(legend.position = "top") +
      scale_y_continuous(labels = NULL) +
      geom_vline(xintercept = coord_name, linetype = "dashed") +
      scale_x_continuous(limits = c(start_position, end_position)) +
      #labs(x = "Physical Position (Mb)", y ="Individual ID", title = zoom_title, fill = "Length")
      labs(x = "Physical Position (Mb)", y ="Individual ID")
    print("plotting gene....")
    gene_plot <- HandyCNV::plot_gene(chr_id = chr_id, start = start_position, end = end_position, show_name = show_name)
    #png(res = 300, filename = zoom_name, width = 3500, height = 2000)
    roh_gene <- plot_grid(gene_plot, zoom_plot, ncol = 1, rel_heights = c(1, 3))
    print(roh_gene)
    ggsave(filename = zoom_name, width = width_1, height = height_1, units = "cm", dpi = 300)
    #dev.off()
    print("Task done, plot was stored in working directory.")
  }

  else if (is.null(individual_id) == "FALSE")
  {
    #plot on individual level
    #cnv_indiv <- cnv_p[which(cnv_p$Sample_ID == "204806050057_R01C01")]
    cnv_indiv <- roh[which(roh$Sample_ID == individual_id),]
    chr_coord <- data.frame("Chr" = seq(1,29,1))
    chr_coord$x <- 160 #adjust x for geom_text
    chr_coord$y <- (chr_coord$Chr-1)*5 + 1 #adjust y for geom_text
    cnv_indiv_coord <- merge(cnv_indiv, chr_coord, all.x = TRUE, sort = FALSE)

    #4.
    indiv_name <- paste0("CNV_of_", individual_id,".png")
    indiv_title <- paste0("CNV Distribution of Individual ", individual_id)
    png(res = 300, filename = indiv_name, width = 3500, height = 2000)
    indiv_plot <- ggplot(cnv_indiv_coord, aes(xmin = Start/1000000, xmax = End/1000000, ymin = (Chr-1)*5, ymax = (Chr-1)*5 + 3)) +
      geom_rect(aes(fill = Length)) +
      scale_fill_gradient(low = "red", high = "blue") +
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
