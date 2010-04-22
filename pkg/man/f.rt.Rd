\name{f.rt}
\alias{f.rt}
\title{Function to compute the t statistic and p values for expression values corresponding to...}
\usage{f.rt(x, y, var.thresh=10^-4, num.thresh=10^(-6), denom.thresh=10^(-6))}
\description{Function to compute the t statistic and p values for expression values corresponding to 
two groups; the function applies a correction for small variances, t statistic numerators
and t statistic denominators. 
(this function is deprecated and can be replaced by the corresponding limma version)}
\value{matrix with four columns: numeratorTStatistic, denominatorTStatistic, tStatistic, 
and pValueTStatistic}
\arguments{\item{x}{matrix of expression values for group one}
\item{y}{matrix of expression values for group two}
\item{var.thresh}{threshold for the variance (a correction is applied if the variance is 
smaller than the number passed to var.thresh)}
\item{num.thresh}{threshold for the numerator of the t statistic (a correction is applied
if the numerator is smaller than the number passed to num.thresh)}
\item{denom.thresh}{threshold for the denominator of the t statistic (a correction is applied
if the denominator is smaller than the number passed to denom.thresh)}}
