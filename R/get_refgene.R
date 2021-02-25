#' Title clean_ensgene
#' clean refgene for ensembl gene.
#' Main works are add column names for table, add gene names to new table, replace missing gene names
#' with Ensembl ID then remove duplicated gene, then return the clean formatted gene list.
#' @param refgene ensembl reference gene download from UCSC website
#' @param ens_id ensembl id download from UCSC website
#'
#' @import dplyr
#'
#' @return formatted gene list.
#' @export clean_ensgene
#'
#' @examples
clean_ensgene <- function(refgene , ens_id){
  refgene = refgene
  ens_id = ens_id
  names(refgene) <- c("bin", "name", "Chr", "strand", "Start", "End",
                        "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
                        "score", "en_id", "cdsStartStat", "cdsEndStat", "exonFrames")
  names(ens_id) <- c("name", "name2")
  refgene <- merge(refgene, ens_id, by = "name", all.x = T)
  refgene_t <- refgene %>%
               mutate(gene = coalesce(name2, name)) %>% #fill the unknown gene names by Ensemble ID
               select(-name2) %>%
               rename(name2 = gene) %>%
               distinct(Chr, name2, .keep_all = T) #there are lots of duplicated genes with slightly different position, retain one only
  refgene_t$name2 <- sub("ENSBTAT000000", "", refgene_t$name2)
  return(refgene_t)
}

#' Title clean_ucsc
#' clean refgene for UCSC gene
#' @param refgene reference gene list with UCSC version which download from UCSC website
#'
#' @import dplyr
#'
#' @return formatted gene list.
#' @export clean_ucsc
#'
#' @examples
clean_ucsc <- function(refgene){
  refgene = refgene
  names(refgene) <- c("bin", "name", "Chr", "strand", "Start", "End",
                      "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
                       "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  refgene <- refgene %>%
              distinct(Chr, name2, .keep_all = T)
  return(refgene)
}


#' Title get_refgene
#' Download and make standard format of reference gene from UCSC website
#' @param gene_version The version of reference gene 'HandyCNV' support to download, if now don't have the version of the gene you want, please feel free to cotact the Maintainer to help you updating this function.
#'
#' @importFrom R.utils gunzip gunzip.default
#' @importFrom data.table fread
#'
#' @return Standard formatted reference gene list
#' @export get_refgene
#'
#' @examples
get_refgene <- function(gene_version = NULL){
  support_ver <- c("cow_ARS_UCSC", "cow_ARS_ENS", "cow_UMD_UCSC", "pig_susScr11_UCSC", "pig_susScr11_ENS", "Human_hg38")
  support_version <- data.frame(version = c("cow_ARS_UCSC",
                                            "cow_ARS_ENS",
                                            "cow_ARS_ENS_id",
                                            "cow_UMD_UCSC",
                                            "pig_susScr11_UCSC",
                                            "pig_susScr11_ENS",
                                            "pig_susScr11_ENS_id",
                                            "Human_hg38"),
                                URL = c("http://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/database/refGene.txt.gz",
                                        "http://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/database/ensGene.txt.gz",
                                        "http://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/database/ensemblToGeneName.txt.gz",
                                        "http://hgdownload.soe.ucsc.edu/goldenPath/bosTau8/database/refGene.txt.gz",
                                        "http://hgdownload.soe.ucsc.edu/goldenPath/susScr11/database/refGene.txt.gz",
                                        "http://hgdownload.soe.ucsc.edu/goldenPath/susScr11/database/ensGene.txt.gz",
                                        "http://hgdownload.soe.ucsc.edu/goldenPath/susScr11/database/ensemblToGeneName.txt.gz",
                                        "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
                                        ))

  #remind what's version we support to automatic download and formatting
  if(missing(gene_version)){
    cat("Now we support to download some of reference gene list of Bovine, Pig and Human from UCSC website. \nPlease select which version to download from following list:\n")
    print(support_ver)
  } else {
    if(!(file.exists("refgene"))){
      dir.create("refgene")
      print("New folder 'refgene' was created in working directory.")
    }

    #call clean function, check if it is ENS reference list?
   if(length(grep("ENS", gene_version)) == 1L){
      print("Converting formmat...")
      line <- which(support_version$version == gene_version)
      line_for_id <- which(support_version$version == paste0(gene_version, "_id"))
      refgene <- fread(input = as.character(support_version$URL[line]), header = FALSE)
      ens_id <- fread(input = as.character(support_version$URL[line_for_id]), header = FALSE)
      gene_version_clean <- clean_ensgene(refgene = refgene, ens_id = ens_id)
      fwrite(gene_version_clean, file = paste0("refgene/", gene_version, ".txt"), sep = "\t", quote = FALSE, col.names = TRUE)
    } else {
      print("Converting formmat...")
      line <- which(support_version$version == gene_version)
      refgene <- fread(input = as.character(support_version$URL[line]), header = FALSE)
      gene_version_clean <- clean_ucsc(refgene = refgene)
      fwrite(x = gene_version_clean, file = paste0("refgene/", gene_version, ".txt"), sep = "\t", quote = FALSE, col.names = TRUE)
    }

    if(file.exists(paste0("refgene/", gene_version, ".txt"))){
      print("Task done.")
    }
  }
}
