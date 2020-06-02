\name{run_infomap}
\alias{run_infomap}
\title{Performs the INFOMAP algorithm on a site-species dataframe.}
\usage{
run_infomap(dat, n_runs = 1, site, sp, ab = NULL, saving_directory, N = 10)
}
\arguments{
\item{dat}{a dataframe containing site column, species column, potentially
abundances of species.}

\item{site}{a character string corresponding to the column with sites.}

\item{sp}{a character string corresponding to the column with species.}

\item{ab}{a character string corresponding to the column with abundances of
species at sites.}

\item{saving_directory}{Path where the .rds output file containing
         the bioregions identified will be saved.}

\item{N}{Number of outer-most loops to run before picking the best solution.}

}
\description{
Run INFOMAP algorithm and save a .rds file with the identified bioregions in a
chosen directory.
See https://www.mapequation.org/infomap/ for more precisions.
}
\details{}
\section{Note}{}

\examples{
# Example dataset from bioRgeo:
data("virtual_sp")
sp_df <- virtual[[1]]

ex <- run_infomap(
  dat = sp_df[which(sp_df$pa > 0), ], site = "site", sp = "sp",
  ab = NULL, saving_directory = "local_path", N = 10)
}
