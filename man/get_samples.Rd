% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_samples.R
\name{get_samples}
\alias{get_samples}
\title{Title get_samples}
\usage{
get_samples(annotated_cnvlist, gene_name)
}
\arguments{
\item{annotated_cnvlist}{CNV-gene annotated list, default format was generated by call_gene function.}

\item{gene_name}{Official gene name in text format.Standard parameter should set as gene_name = "GPB4".}
}
\value{
A datatable contain sample ID and CNV information
}
\description{
Used to get samples ID by searching interested gene from CNV annotation list.
}
