% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fuzzy_filter.R
\name{closest_germline}
\alias{closest_germline}
\title{Determines distance to nearest germline event}
\usage{
closest_germline(bp = NULL, cores = 1, genome = NULL)
}
\arguments{
\item{bp}{\href{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}{Bedpe} from \link[InProgress]{svaba_vcf2bedpe} or \link[InProgress]{manta_vcf2bedpe}}

\item{cores}{Number of cores to run on, default is 1}

\item{genome}{run under hg19 or hg38}
}
\value{
\href{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}{Bedpe} with a column added for distance to nearest germline event
}
\description{
Uses \href{https://gnomad.broadinstitute.org/downloads#v2-structural-variants}{gnomAD} to annotate the nearest germline event to each structural variant.
For more information read \href{https://www.nature.com/articles/s41586-020-2287-8}{gnomAD blog}. Reference is in hg38.
}
