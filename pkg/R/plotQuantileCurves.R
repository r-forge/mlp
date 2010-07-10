#' Plot Quantile Curves for an MLP Run 
#' @param x0 TODO
#' @param y0 TODO
#' @param hqi TODO
#' @param xtp TODO
#' @param qi TODO
#' @param lqi TODO
#' @param sym TODO defaults to TRUE
#' @return no return value; a quantile curve plot is drawn to the current device 
#' @export
plotQuantileCurves <- function(x0, y0, hqi, xtp, qi, lqi, sym = TRUE) {
  plot(x0, y0, xlab = "n", ylab = "MLP", axes = FALSE, col = "#08306B", pch = ".")
  axis(2, lwd = 1.5, col = "#08306B")
  atPositions <- axis(1, labels = FALSE)
  axis(1, lwd = 1.5, at = atPositions, labels = atPositions^2, las = 1, col = "#08306B")
  par(cex = 0.5)
  npc <- 16
  ccc <- "#2171B5"
  if (!is.null(hqi)){ 
    i2 <- y0 > xtp[,length(qi)]
    points(x0[i2], y0[i2], col = ccc, pch = npc) 
    if (length(i2[i2 == TRUE]) > 0){
      text(x = x0[i2], y = jitter(y0[i2], factor=4), labels = names(x0[i2]), cex = 1.25, col = "#08306B")
    }
  }
  if (!is.null(lqi)){ 
    i2 <- y0 < xtp[,1] 
    points(x0[i2], y0[i2], col = ccc, pch = npc) 
  }
  if (sym){ 
    i2 <- y0 < -xtp[,length(qi)] 
    points(x0[i2], y0[i2], col = ccc,pch = npc) 
  }
  par(cex = 1)
  matlines(x0[i2 <- sort.list(x0)],xtp[i2,],col= "#539ECC",lty = 1)
  if (sym) 
    matlines(x0[i2], -xtp[i2,], col= 3,lty = 1)
}
