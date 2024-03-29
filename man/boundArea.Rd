% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boundArea.R
\name{boundArea}
\alias{boundArea}
\title{Annotations analysis around boundaries}
\usage{
boundArea(domain.gr, annot.gr, window.size = 50000, domain.boundary = "start")
}
\arguments{
\item{domain.gr}{\code{GRanges} with domains (TADs, compartments...).}

\item{annot.gr}{\code{GRanges} with genomic annotations (genes, repeat elements...).}

\item{window.size}{Window in base pair surrounding each side of the TAD boundaries.}

\item{domain.boundary}{Domain boundary to be analyzed: \code{"start"} or \code{"end"}. Default is the \code{"start"} of each domain to analyzed annotations around it. Note that it is also possible to take the center of the domains with \code{"center"}.}
}
\value{
\code{GRanges} object with all annotations surrounding all TAD boundaries. TAD boundaries are at the window.size position.
}
\description{
For each domain boundary (start or end), \code{boundArea()} return all annotations that are in the desired window of the boundary (i.e boundary +/- \code{window.size}).
The main use of this function is to be combined with \code{areaHist()} or \code{areaCov()} functions.
It is possible to use the \code{"start"}, \code{"end"} or \code{"center"} of the domains to analyze the distribution of annotations.
}
\examples{
# tad.gr with 2 TADs
tad.gr <- dataframes2grange(
  data.frame(chr = rep(1, 2), start = c(2e6, 4e6), end = c(3e6, 5e6)),
  data.frame(chr = "1", size = 6e6)
)


# annot.gr with 3 genes
annot.gr <- dataframes2grange(
  data.frame(
    chr = rep(1, 3),
    start = c(1.5e6, 3.5e6, 5.5e6),
    end = c(1.6e6, 3.6e6, 5.6e6),
    type = c("a", "b", "c")
  ),
  data.frame(chr = "1", size = 6e6),
  metadata.mcols = 4
)

boundArea(domain.gr = tad.gr, annot.gr = annot.gr, window.size = 2e6)

}
