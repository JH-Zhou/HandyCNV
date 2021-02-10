#' Title plot intensity of cnvr
#'
#' @param cnvr cnvr list, standard file was generated from call_cnvr function
#' @param cnv_annotation standard file was generaed from call_gene function
#' @param intensity the signal inyensity file contains LRR and BAF
#' @param map the map corresponding to the reference genome in cnvr file, standard file was generated from convert_map function
#' @param prefix_bed the prefix of bed, bim and fam file in plink format
#' @param sample_size the total number of unique samples in cnvr list
#' @param common_cnv_threshold two decimal places, combine with sample_size to set the common threshold
#' @param chr_id the number of chromosome, integer, such as chromosome 1 is: 1
#' @param start_position decimal digit, default unit is 'Mb'. such as 23.2112
#' @param end_position decimal digit, default unit is 'Mb'. such as 23.2112
#'
#' @return
#' @export
#'
#' @examples
  plot_intensity <- function(cnvr, cnv_annotation, intensity, map, prefix_bed, sample_size, common_cnv_threshold = 0.05, chr_id, start_position, end_position) {
  #1.Read the CNVR result----
  cnvr <- fread(file = cnvr)
  cnvr <- unite(cnvr, "title", names(cnvr[, c(2:4, 10)]), remove = FALSE) #generate a new columns names
  high_freq <- filter(cnvr, Frequent >= sample_size * common_cnv_threshold)
  #2.read cnv
  cnv <- fread(file = cnv_annotation)
  names(cnv)[names(cnv) == "CNV_Start"] <- "Start"
  names(cnv)[names(cnv) == "CNV_End"] <- "End"
  #merge cnv and cnvr
  setkey(cnv, Chr, Start, End)
  cnv_cnvr <- foverlaps(cnvr, cnv)
  #3. Read SNP intensity file----
  intensity <- fread(file = intensity, skip = 9, header = TRUE)  #
  #4.read map
  map <-fread(file = map, header = TRUE)
  names(map)[names(map) == "Name"] <- 'SNP Name'
  inten_chr <- merge(intensity, map, all.x =TRUE, sort =F) #add location information to intensity file

  #5. Read bed,bim and fam data from plink----
  x <- read.bed.matrix(basename = prefix_bed, rds = NULL) #path and prefix

  #zoom cnv, in order to match cnv value for each snp in Intensity
  chr_length_ars <- data.frame("chr" <- c(29:1), "length" <- c( 51.098607, 45.94015, 45.612108, 51.992305,
                                                                42.350435, 62.317253, 52.498615, 60.773035, 69.862954,
                                                                71.974595, 63.449741, 65.820629, 73.167244, 81.013979,
                                                                85.00778, 82.403003, 83.472345, 87.216183, 106.982474,
                                                                103.308737, 105.454467, 113.31977, 110.682743, 117.80634,
                                                                120.089316, 120.000601, 121.005158, 136.231102, 158.53411))
  names(chr_length_ars) <- c("Chr", "Length")
  chr_length_ars <- chr_length_ars[order(chr_length_ars$Chr),]


  cnv_chr <- cnv[cnv$Chr == chr_id, ]
  cnv_chr_zoom <- filter(cnv_chr, Start >= start_position - 1 & End <= end_position + 1)
  names(cnv_chr_zoom)[names(cnv_chr_zoom) == "Sample_ID"] <- "Sample ID"
  id_coord <- data.frame("Sample_ID" = sort(unique(cnv_chr_zoom$`Sample ID`))) #extract unique ID prepare coordinate
  try(id_coord$Order <- seq(1, nrow(id_coord),1), silent = FALSE)
  id_coord$x <- chr_length_ars[chr_id, 2]
  id_coord$y <- (id_coord$Order-1)*5 + 1
  names(id_coord)[names(id_coord) == "Sample_ID"] <- "Sample ID"
  cnv_chr_zoom <- merge(cnv_chr_zoom, id_coord, by = "Sample ID", all.x = TRUE, sort = FALSE) #prepare original data
  cnv_chr_zoom$zoom_x <- end_position
  cnv_chr_zoom$gene_order <- max(cnv_chr_zoom$Order) + 3
  gene_coord <- group_by(cnv_chr_zoom, name2) %>% slice(1) # generate gene data
  gene_coord$CNV_Start <- gene_coord$g_Start
  gene_coord$Order <- gene_coord$gene_order
  try(gene_coord$CNV_Value <- "5", silent = FALSE)

  #zoom_name <- paste0("Chr", chr_id, "_",start_position,"-",end_position, "Mb", ".png")
  #id_number <- nrow(id_coord)
  #zoom_title <- paste0("CNV on Chromosome ", chr_id, ": ",start_position," - ",end_position, " Mb", " with ", id_number," Individual")
  #png(res = 300, filename = zoom_name, width = 3500, height = 2000)

  #zoom_plot <- ggplot() +
  #  geom_rect(data = cnv_chr_zoom, aes(xmin = Start/1000000, xmax = End/1000000, ymin = (Order-1)*5, ymax = (Order-1)*5 + 3, fill = as.factor(CNV_Value))) +
  #  geom_rect(data = gene_coord, aes(xmin = g_Start/1000000, xmax = g_End/1000000, ymin = (Order-1)*5, ymax = (Order-1)*5 + 3), fill = "black") +
  #  geom_text_repel(data = gene_coord, aes(x = g_Start/1000000, y = (Order-1)*5 + 4, label = name2)) +
  #  geom_hline(yintercept = (max(cnv_chr_zoom$Order) + 2)*5 - 2, linetype = "dashed") +
    #geom_text(aes(zoom_x, y, label = Sample_ID), size = 2.5) +
    #scale_color_manual(values = c("#F8766D", "#A3A500", "#00B0F6", "#E76BF3", "black")) +
  #  theme_bw() +
  #  scale_y_continuous(labels = NULL) +
  #  scale_x_continuous(breaks = seq(round(start_position,2), round(end_position,2), by = 0.2)) +
  #  labs(x = "Physical Position (Mb)", y ="Individual ID", title = zoom_title, fill = "CNV_Num")

  #5.6 Plot BAF, LRR, LF MAF----
  #extract intensity for individual which with cnv in cnvr, so start and end position from cnvr
  sub_inten <- inten_chr[which(inten_chr$Chr == chr_id),]
  sub_inten <- sub_inten[order(sub_inten$Position), ]
  sub_inten <- filter(sub_inten, Position >= start_position & Position <= end_position)
  sub_inten <- sub_inten[which(sub_inten$`Sample ID` %in% cnv_chr_zoom$`Sample ID`), ]

  cnv_chr_zoom.byId = split(cnv_chr_zoom, cnv_chr_zoom$`Sample ID`)
  typeF = function(i) {
    id = sub_inten[i,'Sample ID']         ## get the ID in data_1 for row i
    pos = sub_inten[i,'Position']  ## get the Position in data_1 for row i
    tab = cnv_chr_zoom.byId[[id]]     ## get the subset of data_2 that matches this ID
    ## For each row in the data_2 subset, is the Position
    ## inside the interval [Start,End]
    ## idx is a Boolean (TRUE or FALSE) vector
    idx = pos >= tab$Start & pos <= tab$End
    ## Return the matching Type from the data_2 subset
    ## or return 2 if nothing matches
    ## any(idx): does any element of idx == TRUE?, i.e.,
    ## does the Position match any interval?
    ## Yes -> tab$Type[idx][1]: return the Type for the first match
    ## No  -> return 2
    ifelse( any(idx), tab$CNV_Value[idx][1], 2 )
  }

  CNV_Value = sapply( 1:nrow(sub_inten), typeF )
  sub_inten = cbind(sub_inten, CNV_Value)


  #id_all <- fread("1_CNV_Results/id_all.txt")
  #names(id_all)[1] <- "Sample ID"
  #sub_inten <- merge(sub_inten, id_all, all.x = TRUE, sort = FALSE)

  ld_x = LD(x, c(which(x@snps$id == sub_inten$`SNP Name`[1]), which(x@snps$id == sub_inten$`SNP Name`[nrow(sub_inten)])))
  ld_x[is.na(ld_x)] <- 0
  print(ld_x[1:2])
  snp_info <- x@snps[c(which(x@snps$id == sub_inten$`SNP Name`[1]):which(x@snps$id == sub_inten$`SNP Name`[nrow(sub_inten)])),]

  #show_col(hue_pal()(10))

  baf <- ggplot(sub_inten, aes(x = Position, y =`B Allele Freq`, color = as.factor(CNV_Value))) +
    scale_color_manual(values = c("#F8766D", "#A3A500","grey", "#00B0F6", "#E76BF3")) +
    theme_bw() + theme(legend.position = "top") +geom_point(shape = 1) +
    scale_x_continuous(breaks = seq(start_position, end_position, by = 250000), labels = unit_format(unit = "", scale = 1e-6)) +
    ylim(0,1) + labs(x = 'Position (Mb)', y = 'B Allele Freq', color = "Copy_Num")

  lrr <- ggplot(sub_inten, aes(x = Position, y = `Log R Ratio`, color = as.factor(CNV_Value))) +
    scale_color_manual(values = c("#F8766D", "#A3A500", "grey", "#00B0F6", "#E76BF3")) +
    theme_bw() + theme(legend.position = "top") + geom_point(shape = 1) +
    scale_x_continuous(breaks = seq(start_position, end_position, by = 250000), labels = unit_format(unit = "", scale = 1e-6)) +
    ylim(-2, 2) + labs(x = 'Position (Mb)', y = 'Log R Ratio', color = "Copy_Num")

  maf <- ggplot(data = snp_info) +
    geom_point(aes(x = pos, y = maf, color = "maf")) +
    geom_point(aes(x = pos, y = hz, color = "heterozygosity")) +
    geom_line(aes(x = pos, y = callrate, color = "callrate")) +
    scale_x_continuous(breaks = seq(start_position, end_position, by = 250000), labels = unit_format(unit = "", scale = 1e-6)) +
    ylim(0.0, 1.0) + theme_bw() + theme(legend.position = "top") +
    labs( x = "Position (Mb)", y = "Percentage") +
    scale_color_manual(values = c("maf" = "red", "heterozygosity" = "green", "callrate" = "purple"))

  #ld <- LD.plot(ld_x, snp.positions = snp_info$pos, draw.chr = FALSE)

  ld <- base2grob(~LD.plot(ld_x, snp.positions = snp_info$pos, draw.chr = FALSE))

  b <- cnvr$title[1]
  dir <- paste( b, ".png", sep = "")

  #plot_grid(zoom,zoom2, maf, baf, lrr, align = "v", axis = "t", ncol = 1, rel_heights = c(1,1,1,1,1))
  final_plot <- plot_grid(maf, baf, lrr, ld, align = "v", axis = "t", ncol = 1, rel_heights = c(1,1,1,1,1))
  png(dir, res = 300, width = 4000, height = 000, bg = "transparent")

  #save_plot(filename = dir, plot = final_plot)
}
