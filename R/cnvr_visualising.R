#' Title
#'
#' @param cnvr cnvr file generated from call_cnvr function
#' @param assembly which reference genome assenbly used in the data, each chromosome has different length between various assembly, now only support Bovine UMD3.1 and ARS1.2
#' @param cnv_annotation gene annotated cnv file, generated from call_gene function
#' @param sample_size integer, the total number of unique samples in the cnv result. combine with common_cnv_threshold to plot all CNVs which passed the threshold
#' @param common_cnv_threshold two decimal places, cobine with sample_size to plot all CNVs passed the common threshold
#' @param legend_x the x coordinate of legend, relative to the maximum length of Chromosome, unit is 'Mb'
#' @param legend_y the y coordinate of legend
#'
#' @import ggplot2 dplyr scales reshape2 tidyr
#'
#' @importFrom data.table fread fwrite
#'
#' @return A figture of CNVR distribution map and plot parameters.
#' If given cnv_annotation file, will plot all CNVRs which are passed common threhold.
#' @export cnvr_plot
#'
#' @examples
cnvr_plot <- function(cnvr, assembly = "ARS", legend_x = 127, legend_y = 30, cnv_annotation = NULL, sample_size = NULL, common_cnv_threshold = 0.05) {
  if (is.null(cnv_annotation)) {
    #Prepare parameters for CNVR distribution plot
    #prepare the Y axis for each chromosome
    cnvr_plot_part <- fread(file = cnvr)
    chr <- data.frame("Chr" = seq(1,max(cnvr_plot_part$Chr)))  #generate a chr order
    chr <- chr %>%
           mutate(chr_top = (max(cnvr_plot_part$Chr) + 1 - Chr)*3,
                  chr_bottom <- chr_top - 1.5)   #our plot depend on x and y coordinate axix, the toppest is chr 1 and bottom is chr29, interval of y axis is 3
                                                #bar plot depth is 1.5, this will using for polt rectangle for cnvr

    #chr$chr_top <- (max(cnvr_plot_part$Chr) + 1 -chr$Chr)*3 #our plot depend on x and y coordinate axix, the toppest is chr 1 and bottom is chr29, interval of y axis is 3
    #chr$chr_bottom <- chr$chr_top - 1.5 #bar plot deepth is 1.5, this will using for polt rectangle for cnvr
    #class(chr$Chr) #check the type
    chr$Chr = as.character(chr$Chr)#convert interger to charactor

    if (assembly == "UMD") {
      #UMD3.1 Chromsome Length
      chr_length <- c(51.505224, 46.312546, 45.407902, 51.681464, 42.90417, 62.71493, 52.530062,
                      61.435874, 71.599096, 72.042655, 64.057457,
                      66.004023,75.158596, 81.724687, 85.296676, 84.64839, 84.24035, 91.163125,
                      107.310763, 104.305016, 105.70825, 113.384836, 112.638659, 119.458736,
                      121.1914245, 120.829699, 121.430405, 137.060424, 158.337067)
      chr_name <- paste0("chr", seq(from = 29, to = 1, by = -1))
    } else if(assembly == "ARS"){
      # The length of X Chromsome is 139.009144 in ARS reference genome
      # ARS assembly
      chr_length <- c( 51.098607, 45.94015, 45.612108, 51.992305, 42.350435,
                       62.317253, 52.498615, 60.773035, 69.862954,
                           71.974595, 63.449741, 65.820629, 73.167244, 81.013979,
                           85.00778, 82.403003, 83.472345, 87.216183, 106.982474,
                           103.308737, 105.454467, 113.31977, 110.682743, 117.80634,
                           120.089316, 120.000601, 121.005158, 136.231102, 158.53411)
      chr_name <- paste0("chr", seq(from = 29, to = 1, by = -1))
    } else {
      #construct the border of each chromosome from imput data
      chr_length <- cnvr_plot_part %>%
                    group_by(Chr) %>%
                    summarise(length = max(End) / 1000000) %>%
                    arrange(-as.numeric(Chr)) %>%
                    select(length) %>%
                    as.matrix() %>%
                    t()
      chr_name <- paste0("chr", seq(from = max(cnvr_plot_part$Chr), to = 1, by = -1))
    }

    #4.6.1 prepare CNVPartition plot input data-----
    cnvr_plot_part$start_left <- cnvr_plot_part$Start/1000000 #convert bp to Mbp
    cnvr_plot_part$end_right <- cnvr_plot_part$End/1000000 #convert bp to Mbp
    cnvr_plot_part$Chr <- as.character(cnvr_plot_part$Chr)
    cnvr_plot_part <- merge(cnvr_plot_part, chr, all.x = TRUE)
    fwrite(cnvr_plot_part, file = "cnvr_plot.txt", sep ="\t", quote = FALSE)

    #this file need reread if you start from here
    #cnvr_plot_part <- fread("CNVRplot/cnvr_plot_part.txt", sep ="\t", quote = FALSE)
    gain_part <- cnvr_plot_part[cnvr_plot_part$Type == "Gain", ]
    loss_part <- cnvr_plot_part[cnvr_plot_part$Type == "Loss", ]
    mixed_part <- cnvr_plot_part[cnvr_plot_part$Type == "Mixed", ]

    #4.6.2 start make cnvr plot for CNVPartition
    png(filename = "cnvr_plot.png",
        res = 300,
        width = 3000, height = 2000,
        bg = "transparent"
    )

    #setting the graphic parameter
    par(lab=c(max(cnvr_plot_part$Chr),max(cnvr_plot_part$Chr),3),las=1,yaxt="s",xaxt="s",mar=c(7,7,4,5))

    #draw a bar plot and setup each bar name
    bar <- barplot(chr_length, horiz=TRUE,width=1.5,space=1,xlim=c(0,max(chr_length)),
                   names.arg= chr_name,
                   col=c("white"),xlab="Physical Position (Mbp)",cex.name=0.8)

    #read the gain type file, each file need prepare four columns as following order: start_left, end_right, chr_top and chr_bottom
    rect(xleft = gain_part$start_left, ybottom = gain_part$chr_bottom, xright = gain_part$end_right, ytop = gain_part$chr_top, col = "red", border = "red")

    #read the loss type file and draw a rectangle shape for each CNVR and color to green
    #loss <- cnvr_plot_part[cnvr_plot_part$Type == 'LOSS']
    rect(xleft = loss_part$start_left, ybottom = loss_part$chr_bottom, xright = loss_part$end_right, ytop = loss_part$chr_top, col = "green", border = "green")

    #read the mixed file and plot rectangle as red
    #mixed <- cnvr_plot_part[cnvr_plot_part$Type == 'mixed']
    rect(xleft = mixed_part$start_left, ybottom = mixed_part$chr_bottom, xright = mixed_part$end_right, ytop = mixed_part$chr_top, col = "blue", border = "blue")

    #prepare Y axis for CNVPartition overlap region, mark a under line on the overlap region
    #rect(xleft = overlap_part$start_left, ybottom = overlap_part$chr_bottom, xright = overlap_part$end_right, ytop = overlap_part$chr_top, col = "purple", border = "purple")


    #add legend
    legend(legend_x, legend_y,c("Gain","Loss","Mixed"), pch = c(15, 15, 15),
           col =c("red", "green", "blue"), bty = "n")

    #title(main = "CNVR Distribution on Population Level")

    dev.off()

    if (file.exists("cnvr_plot.png")) {
      print("Task done, CNVR Distribution Plot saved in working directory.")
    } else {
      print("Task faild, please check your CNVR input file carefully, the input cnvr file should be the CNVR results generated from call_cnvr function")
    }
  }



  else {
    cnvr <- fread(file = cnvr) # read cnvr result from call_cnvr function
    # high_freq <- cnvr[which(cnvr[, cnvr$Frequent >= sample_size * 0.05]), ] #call common CNVRs
    high_freq <- filter(cnvr, Frequent >= sample_size * common_cnv_threshold)
    print(paste0("There ", nrow(high_freq), " high frequent CNVR passed the customized threshold."))
    for (i in 1:nrow(high_freq)) {
      print(paste0("Ploting CNVR ", i, "..." ))
      cnv_visual(clean_cnv = cnv_annotation, chr_id = high_freq$Chr[i], start_position = high_freq$Start[i]/1000000, end_position = high_freq$End[i]/1000000, plot_gene = 1)
    }
  }
}

