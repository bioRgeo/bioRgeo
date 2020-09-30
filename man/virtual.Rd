\name{virtual}
\alias{virtual}
\title{Load virtual data}
\description{
Load a virtual dataset made of sites with species occurring in them. Each site has XY coordinates.
}
\usage{
data("virtual")
}
\arguments{

}
\details{

}
\format{
\code{virtual} is a long-format data.frame containing 40000 rows and nine columns. Each row describes a species occurrence within a site. This dataset has been constructed by simulating a virtual raster and virtual species with the \code{virtualspecies} R package. The columns are the following:
\describe{
\item{sp}{species occurring within a given site}
\item{x}{X coordinate of a site}
\item{y}{Y coordinate of a site}
\item{suitab}{suitability of the species in a site}
\item{pa}{presence/absence: is the species present (1) or absent (0) in a site.}
\item{opt}{Optimal environment for the species.}
\item{width}{Width of the environmental response curve of the species.}
\item{env}{Environmental value at a site.}
\item{site}{Sites where species occur.}
}

}

\section{Note}{
The construction of this dataset is detailed in the vignette of \code{bioRgeo}.
\code{virtualspecies} R package has been used to construct it.
}

\examples{
data("virtual")
}
