% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_haplotype.R
\name{get_haplotype}
\alias{get_haplotype}
\title{Get haplotype}
\usage{
get_haplotype(geno, pos)
}
\arguments{
\item{geno}{phase genotype and map file from 'prep_phased' function}

\item{pos}{start and end position from 'closer_snp' function}
}
\value{
A list contain Sample with recoded Haplotype and Index of Haplotype and Recode ID
}
\description{
Generate haplotypes from Phased genotype file by given region.
}
