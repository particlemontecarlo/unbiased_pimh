% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_get_gss.R
\name{get_gss}
\alias{get_gss}
\title{Nonlinear growth model as in Gordon Salmond Smith 1993}
\usage{
get_gss()
}
\value{
A list
}
\description{
This function returns a list with objects such as
* rinit, rinit_rand to sample from the initial distribution
* rtransition, rtransition_rand to sample from the transition
* dtransition to evaluate the transition density
* dmeasurement to evaluate the measurement density
* dimension, which represents the dimension of the latent process
}
