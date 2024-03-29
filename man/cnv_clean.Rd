% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnv_clean.R
\name{cnv_clean}
\alias{cnv_clean}
\title{Clean CNV}
\usage{
cnv_clean(
  cnvpartition = NULL,
  penncnv = NULL,
  standard_cnv = NULL,
  drop_length = 5,
  penn_id_sep = "cnv/",
  folder = "cnv_clean"
)
}
\arguments{
\item{cnvpartition}{load CNV results from CNVPartition}

\item{penncnv}{load CNV results from PennCNV}

\item{standard_cnv}{Load a user-generated CNV input file. The following columns must be present: Sample_ID, Chr, Start, End, CNV_Value}

\item{drop_length}{exclude CNVs longer than this threshold, unit is "Mb"}

\item{penn_id_sep}{the separator in the `Sample ID` column of PennCNV results. Useful if the ID is bound to the path}

\item{folder}{set the name of the output folder}
}
\value{
Formatted CNV results and brief summary files.
}
\description{
Import CNV call results produced by the software packages PennCNV and CNVPartition, and converts them into a standard format for use in other functions in the `HandyCNV` package.
It now support to accept a CNV list in the standard format, the standard format should have at least five columns with header: Sample_ID, Chr, Start, End, CNV_Value
}
