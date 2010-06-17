\name{quantileCurves}
\alias{quantileCurves}
\title{Function used to calculate the quantile curves; the function is mainly used by ctpval1...}
\usage{quantileCurves(x, y, x0=x, y0=y, type=c("none", "dec", "inc"), m=20,
    lqi=0.05, hqi=0.95, sym=FALSE, plot=TRUE, flag=FALSE, df=15,
    logtran=FALSE)}
\description{Function used to calculate the quantile curves; the function is mainly used by ctpval1}
\value{matrix with the y0 values in the first column, the quantiles at x0 (as many columns as is specified in lqi and hqi) 
and the x0 values in the last column}
\arguments{\item{x}{x coordinates of the points used to generate the smoothed quantile curves}
\item{y}{y coordinates of the points used to generate the smoothed quantile curves}
\item{x0}{x coordinates for the observed data points; by default}
\item{y0}{y coordinates for the observed data points}
\item{type}{constraints to the smoothed quantile curve, one of "none" (default), "dec" (decreasing) or 
"inc" (increasing)}
\item{m}{block size for constructing the smoothed quantile curve (which defines the neighborhood
of the gene set); defaults to 20.}
\item{lqi}{vector of the percentiles at which the lower quantile curve is desired}
\item{hqi}{vector of the percentiles at which the upper quantile curve is desired}
\item{sym}{logical; should the curves be symmetric (TRUE) or not (FALSE); one does not expect it to
be symmetric, so the argument defaults to FALSE.}
\item{plot}{logical; should a plot be generated for the results ? Defaults to TRUE.}
\item{flag}{logical; if TRUE an alternative method is used to compute the quantiles which will be faster
but less accurate for large datasets than when using the apply function; defaults to FALSE.}
\item{df}{degrees of freedom for the smoothing spline}
\item{logtran}{whether whether or not to log-transform the x's and y's to calculate the smoothed
quantile curves}}
