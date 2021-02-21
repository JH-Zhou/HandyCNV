#' Title clean_ensgene
#' clean refgene for ensembl gene.
#' Main works are add column names for table, add gene names to new table, replace missing gene names
#' with Ensembl ID then remove duplicated gene, then return the clean formatted gene list.
#' @param refgene ensembl reference gene download from UCSC website
#' @param ens_id enseml id download from UCSC website
#'
#' @return formated gene list.
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
    mutate(gene = coalesce(name2, name)) %>% #fill the unkown gene names by Ensemble ID
    select(-name2) %>%
    rename(name2 = gene) %>%
    distinct(Chr, name2, .keep_all = T) #there are lots of duplicated genes with slightly different position, retain one only
  return(refgene_t)
}

#' Title clean_ucsc
#' clean refgene for UCSC gene
#' @param refgene reference gene list with UCSC verison which download from UCSC website
#'
#' @return formated gene list.
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
#' @param gene_verison The verison of reference gene 'HandyCNV' support to download, if now don't have the version of the gene you want, please feel free to cotact the Maintainer to help you updating this function.
#'
#' @return Standard formatted reference gene list
#' @export get_refgene
#'
#' @examples
get_refgene <- function(gene_verison = NULL){
  support_ver <- c("cow_ARS_UCSC", "cow_ARS_ENS", "cow_UMD_UCSC", "pig_susScr11_UCSC", "pig_susScr11_ENS", "Human_hg38")
  support_verison <- data.frame(verison = c("cow_ARS_UCSC",
                                            "cow_ARS_ENS",
                                            "cow_ARS_Ens_id",
                                            "cow_UMD_UCSC",
                                            "pig_susScr11_UCSC",
                                            "pig_susScr11_ENS",
                                            "pig_susScr11_ENS_id",
                                            "Human_hg38"),
                                URL = c("http://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/database/refGene.txt.gz",
                                        "http://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/database/ensGene.txt.gz",
                                        "http://hgdownload.soe.ucsc.edu/goldenPath/bosTau9/database/ensemblToGeneName.txt",
                                        "http://hgdownload.soe.ucsc.edu/goldenPath/bosTau8/database/refGene.txt.gz",
                                        "http://hgdownload.soe.ucsc.edu/goldenPath/susScr11/database/refGene.txt.gz",
                                        "http://hgdownload.soe.ucsc.edu/goldenPath/susScr11/database/ensGene.txt.gz",
                                        "http://hgdownload.soe.ucsc.edu/goldenPath/susScr11/database/ensemblToGeneName.txt.gz",
                                        "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
                                        ))

  #remind what's version we support to automatic download and formatting
  if(missing(gene_verison)){
    cat("Now we support to download some of reference gene list of Bovine, Pig and Human from UCSC website. \nPlease select which version to download from following list:\n")
    print(support_ver)
  } else {
    if(!(file.exists("refgene"))){
      dir.create("refgene")
      print("New folder 'refgene' was created in working directory.")
    }

    #call clean function, check if it is ENS reference list?
   if(length(grep("ENS", gene_verison)) == 1L){
      line <- which(gene_verison %in% support_verison$verison)
      refgene <- fread(input = as.character(support_verison$URL[line]), header = FALSE)
      ens_id <- fread(input = paste0(verison, "_id"))
      gene_verison_clean <- clean_ensgene(refgene = ref, ens_id = ens_id)
      fwrite(gene_verison_clean, file = paste0("refgene/", gene_verison, ".txt"), sep = "\t", quote = FALSE, col.names = TRUE)
    } else {
      line <- which(gene_verison %in% support_verison$verison)
      refgene <- fread(input = as.character(support_verison$URL[line]), header = FALSE)
      gene_verison_clean <- clean_ucsc(refgene = refgene)
      fwrite(x = gene_verison_clean, file = paste0("refgene/", gene_verison, ".txt"), sep = "\t", quote = FALSE, col.names = TRUE)
    }
  }
}
