cnv_visual <- function(clean_cnv, max_chr_length = NULL, chr_id = NULL, chr_length = NULL, start_position = NULL, end_position = NULL, individual_id = NULL, plot_title = NULL) {
  #myAgr <- formals(cnv_visual)
  #prepare for population data
  cnv <- fread(file = clean_cnv)
  id_coord <- data.frame("Sample_ID" = sort(unique(cnv$Sample_ID))) #extract unique ID prepare coordinate
  id_coord$Order <- seq(1,nrow(id_coord),1)
  id_coord$x <- 160
  id_coord$y <- (id_coord$Order-1)*5 + 1
  cnv_coord <- merge(cnv, id_coord, all.x = TRUE, sort = FALSE) #prepare original data

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
  print("Task done, plot was stored in setted directory.")
  }

  else if(is.null(chr_id) == "FALSE" & is.null(start_position))
   {

  #2.plot for specific chromosome
  cnv_chr <- cnv_coord[cnv_coord$Chr == chr_id, ]
  chr_name <- paste0("Chr", chr_id, ".png")
  chr_title <- paste0("CNV on Chromosome ", chr_id)
  png(res = 300, filename = chr_name, width = 3500, height = 2000)
  #png(res = 300, filename = "10_chr.png", width = 3500, height = 2000)
  chr_plot <- ggplot(cnv_chr, aes(xmin = Start/1000000, xmax = End/1000000, ymin = (Order-1)*5, ymax = (Order-1)*5 + 3)) +
    geom_rect(aes(fill = as.factor(CNV_Value))) +
    geom_text(aes(x,y, label = Sample_ID), size = 1.5, check_overlap = TRUE) +
    theme_bw() +
    scale_y_continuous(labels = NULL) +
    scale_x_continuous(breaks = seq(0, 160, by = 5), limits = c(0, 160)) +
    labs(x = "Physical Position (Mb)", y ="Individual ID", title = chr_title, fill = "CNV_Num")
  print(chr_plot)
  dev.off()
  print("Task done, plot was stored in setted directory.")
  }

  #else if(is.null(start_position & end_position) == "FALSE")
  else if(is.null(start_position) == "FALSE" & is.null(end_position) == "FALSE")
    {
  #3.zoom into specific region
  cnv_chr <- cnv_coord[cnv_coord$Chr == chr_id, ]
  cnv_chr_zoom <- filter(cnv_chr, Start >= start_position * 1000000 & End <= end_position * 1000000)
  cnv_chr_zoom$zoom_x <- end_position
  zoom_name <- paste0("Chr", chr_id, "_",start_position,"-",end_position, "Mb", ".png")
  zoom_title <- paste0("CNV on Chromosome ", chr_id, ": ",start_position," - ",end_position, " Mb", " - ", plot_title)
  png(res = 300, filename = zoom_name, width = 3500, height = 2000)
  zoom_plot <- ggplot(cnv_chr_zoom, aes(xmin = Start/1000000, xmax = End/1000000, ymin = (Order-1)*5, ymax = (Order-1)*5 + 3)) +
    geom_rect(aes(fill = as.factor(CNV_Value))) +
    geom_text(aes(zoom_x, y, label = Sample_ID), size = 2.5) + theme_bw() +
    scale_y_continuous(labels = NULL) +
    scale_x_continuous(breaks = seq(start_position, end_position, by = 0.25), limits = c(start_position, end_position + 1)) +
    labs(x = "Physical Position (Mb)", y ="Individual ID", title = zoom_title, fill = "CNV_Num")
  print(zoom_plot)
  dev.off()
  print("Task done, plot was stored in setted directory.")
  }

  else if(is.null(individual_id) == "FALSE")
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
  print("Task done, plot was stored in setted directory.")
  #return(cnv_coord)
  }

  else{
    print("Warning: Lack of input parameters!!!
           Warning: Please check input parameters carefully!!!")
  }
}

