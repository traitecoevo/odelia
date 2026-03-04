% Generated manually for odelia_load_dll
\name{odelia_load_dll}
\alias{odelia_load_dll}
\title{Load odelia shared library with configurable symbol visibility}
\usage{
odelia_load_dll(local = FALSE, now = TRUE)
}
\arguments{
  \item{local}{Passed to \code{dyn.load}. Default \code{FALSE} to expose symbols globally for downstream dynamic linking.}

  \item{now}{Passed to \code{dyn.load}. Default \code{TRUE}.}
}
\value{
Invisibly returns the resolved path to the \code{odelia} shared library.
}
\description{
Locates the installed \code{odelia} shared library and calls \code{dyn.load} on it.
Useful in \code{Rcpp::sourceCpp()} workflows where symbols from the package shared
object must be available to a temporary shared object.
}
