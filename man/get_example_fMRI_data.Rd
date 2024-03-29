% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_data.R
\name{get_example_fMRI_data}
\alias{get_example_fMRI_data}
\title{get_example_fMRI_data}
\usage{
get_example_fMRI_data(save_path = NULL, force = FALSE, subject = 1)
}
\arguments{
\item{save_path}{location where the data the example data is stored.}

\item{force}{force the download, even if the data already exists.}

\item{subject}{The example subject, must be 1 or 2.}
}
\value{
It returns an array of dimensions \code{[91, 109, 91, 310]}.
}
\description{
This function is used to download the example data used in the Vignettes.
}
\details{
The data for this example is related to an fMRI experiment where a sound stimulus is presented. 
That experiment is intended to offer a "voice localizer" scan, which allows rapid and reliable 
localization of the voice-sensitive "temporal voice areas" (TVA) of the human auditory cortex
\insertCite{pernet2015human}{BayesDLMfMRI}. The data of this "voice localizer" scan is freely  
available on the online platform OpenNEURO \insertCite{gorgolewski2017openneuro}{BayesDLMfMRI}.
}
\examples{
\dontrun{
# This example can take a long time to run.
fMRI.data  <- get_example_fMRI_data()
}
}
\references{
\insertRef{pernet2015human}{BayesDLMfMRI}

\insertRef{gorgolewski2017openneuro}{BayesDLMfMRI}
}
