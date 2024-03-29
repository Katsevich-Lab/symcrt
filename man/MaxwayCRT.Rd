% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{MaxwayCRT}
\alias{MaxwayCRT}
\title{Maxway CRT proposed in Li & Liu 2022.}
\usage{
MaxwayCRT(data, X_on_Z_reg, Y_on_Z_reg, test_hyperparams)
}
\arguments{
\item{data}{A named list with fields X, Y, Z.}

\item{X_on_Z_reg}{The regression method to apply for X|Z.}

\item{Y_on_Z_reg}{The regression method to apply for Y|Z.}

\item{test_hyperparams}{Additional test hyperparameters}
}
\value{
A data frame with columns "parameter," "target," "value" with p-value.
}
\description{
Maxway CRT proposed in Li & Liu 2022.
}
