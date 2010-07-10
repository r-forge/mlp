#' Plot the Results of an MLP Run
#' @param x object of class 'MLP'
#' @param y argument added to comply with generic; not used, defaults to NULL
#' @param type character of length one; one of 'barplot', 'GOgraph' or 'quantileCurves'
#' @param ... further arguments for the plot functions for each type 
#' @return for type = "barplot", the midpoints of the barplot 
#' @S3method plot MLP
#' @export
plot.MLP <- function(x, y = NULL, type = c("barplot", "GOgraph", "quantileCurves") , ...){
  type <- match.arg(type)
  switch(type, 
      barplot = {
        mlpBarplot(object = x, ...)
      }, GOgraph = {
        cat("Not yet implemented\n")
      }, quantileCurves = {
        if (is.null(attr(x, "quantileCurveInformation"))){
          stop("A quantile curve plot can only be drawn for an MLP run for which smoothPValues = TRUE")
        } else {
           quantileCurveInformation <- attr(x, "quantileCurveInformation")
           plotQuantileCurves(x0 = quantileCurveInformation$x0, 
               y0 = quantileCurveInformation$y0,
               hqi = quantileCurveInformation$hqi, 
               xtp = quantileCurveInformation$xtp,
               qi = quantileCurveInformation$qi, 
               lqi = quantileCurveInformation$lqi,
               ...) # dots a.o. for sym argument
        }
      }) 
}
