% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnv_summary.R
\name{cnv_summary_plot}
\alias{cnv_summary_plot}
\title{Summary plots for CNVs}
\usage{
cnv_summary_plot(
  clean_cnv,
  length_group = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1, 2, 5),
  plot_sum_1 = FALSE,
  plot_sum_2 = FALSE,
  folder = "cnv_summary_plot",
  col_0 = "hotpink",
  col_1 = "turquoise",
  col_2 = "gray",
  col_3 = "tomato",
  col_4 = "deepskyblue",
  height_sum1 = 26,
  width_sum1 = 20,
  height_sum2 = 20,
  width_sum2 = 27
)
}
\arguments{
\item{clean_cnv}{a clean CNV file, as generated by the `clean_cnv` function}

\item{length_group}{a vector of lengths (in Mb) to classify CNV lengths for plotting. For example, the vector `c(0.05, 0.3,  0.6, 1)`, instructs the function to divide the CNV lengths into four groups: '<0.05Mb', '0.05 - 0.3Mb', '0.3-0.6Mb' and '>1Mb'. Up to 11 values can be provided.}

\item{plot_sum_1}{plot the first type of summary plot combination}

\item{plot_sum_2}{plot the second type of summary plot combination}

\item{folder}{set the name of the output folder}

\item{col_0}{set colour for 0 copy of CNV}

\item{col_1}{set colour for 1 copy of CNV}

\item{col_2}{set colour for 2 copy of CNV (which might be ROH)}

\item{col_3}{set colour for 3 copy of CNV}

\item{col_4}{set colour for 4 copy of CNV}

\item{height_sum1}{height of Summary Plot 1 (in cm)}

\item{width_sum1}{width of Summary Plot 1 (in cm)}

\item{height_sum2}{height of Summary Plot 2 (in cm)}

\item{width_sum2}{width of Summary Plot 2 (in cm)}
}
\value{
Some summary pictures
}
\description{
Make summary graphs for CNVs. Types of plots produced include:
*  CNV counts by length and copy number
*  Length distribution by copy number
*  CNV counts by chromosome and copy number
*  CNV counts by individual
}
