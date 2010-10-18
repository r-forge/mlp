\name{smdecreasing1}
\alias{smdecreasing1}
\title{This function is used by quantileCurves when its type argument specifies monotonic non-increasing curves or...}
\usage{smdecreasing1(z, w, decreasing=TRUE)
}
\description{This function is used by quantileCurves when its type argument specifies monotonic non-increasing curves or
monotonic non-decreasing curves.}
\value{vector of y values with the monotonicity constraint imposed}
\arguments{\item{z}{output of smooth.spline with x component sorted; if z is a plain vector, w is used as the
y variable. The y component of z must be non-increasing, if not the function will make it non-increasing}
\item{w}{vector of y values in case z is a plain vector of x values}
\item{decreasing}{logical; if TRUE, specifies the curve to be non-increasing, if FALSE specifies the
curve to be non-decreasing}
}

