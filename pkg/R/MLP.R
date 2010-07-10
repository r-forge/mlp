#' Compute column means of minus log10 of the data
#' @param x matrix
#' @return vector of column means of -log10(x) 
#' @export
computeMLP <- function(x){
  x <- as.matrix(x)
  if (ncol(x) > 1){
      retval <- c(nrow(x), colMeans(-log10(data.matrix(x))))
	  names(retval) <- c("n", paste("u", 1:ncol(x), sep = ""))
  } else {
     retval <- c(length(x), mean(-log10(x)))
	 names(retval) <- c("n", "u")
  }
  return(retval)
}

#' This function calculates the mean of the gene-statistic y[, 2] (the MLP statistic)
#'    for each unique gene set in geneSetPValues.
#' @param geneSetPValues is a list of gene sets; each component of the list (i.e. each gene
#'   set) contains a vector of p values for the genes in that gene set
#'   respectively
#' @return data frame with two columns; first column contains the number of genes
#'    in the gene set; the second column contains the MLP statistic (u) for the gene set.
#' @export
mlpStatistic <- function(geneSetPValues){
  retval <- sapply(geneSetPValues, computeMLP)
  return(retval)
}

#' Calculate all unique permutations for 1:n, with k = 2 in each group
#'   (e.g. wild-type vs. knock-out) ; used for column permutations only
#' @param n TODO
#' @param k TODO
permtwo = function (n, k) 
{
  ### From Javier 1/27/2006
  x <- c(0, 1)
  y <- NULL
  for (i in 2:n) {
    x <- rbind(cbind(x, 0), cbind(x, 1))
    j <- c(x %*% rep(1, i))
    jj <- j == k
    if ((ni <- sum(jj)) > 0) {
      xx <- matrix(0, ni, n)
      xx[, 1:i] <- x[jj, ]
      y <- rbind(y, xx)
    }
    x <- x[j < k & (j + n - i) >= k, ]
  }
  t(apply(y, 1, function(x) c((1:length(x))[x == 1], (1:length(x))[x == 
                        0])))
}


#' Calculate all permutations
#' (used for column permutations)
#' @param gr TODO
#' @export
permk <- function(gr) 
{
  ### From Javier 1/27/2006
  ff <- cumsum(gr)
  ttt <- permtwo(ff[2],ff[1])
  for(i in 3:length(gr)) 
  {
    tt <- ttt
    tt1 <- permtwo(ff[i],ff[i-1])
    ttt <- NULL
    for(i in ff[i-1]+1:gr[i]) tt <- cbind(tt,i)
    for(i in 1:nrow(tt1)) ttt <- rbind(ttt,cbind(tt,6,7)[,tt1[i,]])
  }
  ttt
}

###======================== MLP-CPERM===================================


#' Function that will replace all gene identifiers in a gene set with the p values associated to the gene identifiers
#' @param geneSet list of gene sets; each gene set (component of the list) is a character vector of gene identifiers
#' @param geneStatistic named vector; each name is a gene identifier and the corresponding element contains the p value
#'   corresponding to its gene identifier
#' @return list of gene sets; each component of the list (gene set) is a numeric vector of p values; the names
#'   of the numeric vector are the gene identifiers corresponding to the p values
#' @export 
mapGeneSetStatistic <- function(geneSet, geneStatistic){
	myfun <- function(x){
		rv <- geneStatistic[row.names(geneStatistic) %in% x, ]
		return(rv)
	}
	retval <- lapply(geneSet, myfun)
	return(retval)
}


#' Function to convert probeset datasets to gene datasets
#' (deprecated) 
#' TODO look for BioConductor based equivalent
#' @param xp input table of gene sets, probe.ids, and gene.ids, in that order.
#' @param yp input table of probe.ids and original and permuted probeset-statistics.
#' @return TODO fill out 
#' @export
pp2g <- function(xp, yp){
  
  x1  <- xp[, 1]
  x2  <- xp[, 2]
  x3  <- xp[, 3]
  
  pid <- yp[, 1]
  yp  <- yp[, -1]
  dimnames(yp)[[1]] <- pid # need for filtering for y downstream.
  
  y2xp <- mapGeneSetStatistic(xp[,2], pid) # map p-values in yp to probesets in xp.
  p1   <- data.frame(Probe.ID = pid[y2xp], P0 = yp[y2xp, 2])
  
  # Create gene-statistic table; take min p-value for all probesets for a given gene.
  x30 <- sort(unique(x3))
  x30 <- x30[x30 != ""]
  n30 <- length(x30)
  p30 <- vector("character", n30)
  
  for (j in 1:n30){ 
    # p30[j] <- (min(p1[x3==x30[j],2], na.rm = TRUE) ) # min p-value
    p30[j] <- p1[ p1[, 2] == min(p1[x3==x30[j], 2]), 1][1] # probe-set corr to min p-value
    # need checksum if genes have multiple mins.
  }
  y <- data.frame(Gene.ID = I(x30), yp[p30, ])
  
  # Create unique(GO.#, Gene.ID) combination table.
  y2x3<- mapGeneSetStatistic(x3,x30)
  x41 <- x1*10^10 + y2x3 # (GO number inflated  + gene number)
  n41 <- length(x41)
  x40 <- sort(unique(x41)) # multiple probe-sets corr. to gene can cause duplication in x41.
  n40 <- length(x40)  
  jj <- vector("numeric",n40)
  for (j in 1:n40){
    jj[j] <- (1:n41)[x41==x40[j]][1]
  }
  x <- data.frame(GO.Number=x1[jj],Gene.ID=I(x3[jj]))
  
  out <- list(x = x, y = y)
  return(out)
}


#' Function to define the cutoffs to be used for each individual gene set statistic;
#' @param observedGeneSetStats dataframe as returned by mlpStatistic 
#' @param permutationGeneSetStats output of MLP except for the gene set size column; 
#'   matrix of gene set statistics based on the permutation procedure
#' @return vector of p values for each gene set
#' @export
getIndividualPValues <- function(observedGeneSetStats, permutationGeneSetStats){
  # individual cut-offs
  w1 <- matrix(rep(observedGeneSetStats[, 2], ncol(permutationGeneSetStats)), nrow = nrow(observedGeneSetStats), ncol = ncol(permutationGeneSetStats))
  dw0w1 <-  ifelse(w1 - permutationGeneSetStats > 0, 1, 0)
  pw0 <- 1 - rowMeans(dw0w1) 
  return(pw0)
}

#' Function to define the smoothed monotonic cutoffs leveraging across gene sets of similar size;
#' it calls quantileCurves and ctpval1 to create smoothed quantile curves and calculate the corresponding
#' p values for the gene set.  
#' @param observedGeneSetStats dataframe as returned by mlpStatistic 
#' @param permutationGeneSetStats output of MLP except for the geneset size column; 
#'   matrix of gene set statistics based on the permutation procedure
#' @param q.cutoff vector of quantiles at which p values for each gene set are desired (to be
#'   specified by the users)
#' @param df degrees of freedom used for the smoothing splines used to calculate the smoothed quantile curves; defaults to 9.  
#' @return vector of p values for each gene set
#' @export
getSmoothedPValues <- function(observedGeneSetStats, permutationGeneSetStats, q.cutoff, df = 9){
  ### This imposes the decreasing criterion and has df =9
  w1 <- cbind(rep(observedGeneSetStats[,1], ncol(permutationGeneSetStats[,])), 
              as.vector(permutationGeneSetStats[,]))
  lqi <- NULL
  hqi <- q.cutoff

  work <- quantileCurves(x = sqrt(w1[, 1]), y = w1[, 2], x0 = sqrt(observedGeneSetStats[,1]), y0 = observedGeneSetStats[,2],
      type = "dec", m = 20, lqi = lqi, hqi = hqi, sym = FALSE, plot = TRUE, flag = FALSE, 
      df = df, logtran = FALSE)
  pval <- ctpval1(work, q.cutoff = q.cutoff)
  rm(work)
  return(pval)
}

#' Function used to calculate the quantile curves; the function is mainly used by ctpval1
#' @param x x coordinates of the points used to generate the smoothed quantile curves
#' @param y y coordinates of the points used to generate the smoothed quantile curves
#' @param x0 x coordinates for the observed data points; by default
#' @param y0 y coordinates for the observed data points
#' @param type constraints to the smoothed quantile curve, one of "none" (default), "dec" (decreasing) or 
#'   "inc" (increasing)
#' @param m block size for constructing the smoothed quantile curve (which defines the neighborhood
#'   of the gene set); defaults to 20.
#' @param lqi vector of the percentiles at which the lower quantile curve is desired
#' @param hqi vector of the percentiles at which the upper quantile curve is desired
#' @param sym logical; should the curves be symmetric (TRUE) or not (FALSE); one does not expect it to
#'   be symmetric, so the argument defaults to FALSE.
#' @param plot logical; should a plot be generated for the results ? Defaults to TRUE.
#' @param flag logical; if TRUE an alternative method is used to compute the quantiles which will be faster
#'   but less accurate for large datasets than when using the apply function; defaults to FALSE. 
#' @param df degrees of freedom for the smoothing spline
#' @param logtran whether whether or not to log-transform the x's and y's to calculate the smoothed
#'   quantile curves
#' @return matrix with the y0 values in the first column, the quantiles at x0 (as many columns as is specified in lqi and hqi) 
#'   and the x0 values in the last column
#' @export 
quantileCurves <- function(x, y, x0 = x, y0 = y, type = c("none", "dec", "inc"), m = 20, lqi = 0.05, hqi = 0.95,
    sym = FALSE, plot = TRUE, flag = FALSE, df = 15, logtran = FALSE) {
  
  type <- match.arg(type)
  if (is.null(lqi) & is.null(hqi))
	  stop("Please provide at least one of 'lqi' or 'hqi'.")
  
  qi <- if (!is.null(lqi)) c(sort(lqi), sort(hqi)) else sort(hqi)
  
  if (logtran){
    x <- log(x)
    x0 <- log(x0)
    y <- log(y)
    y0 <- log(y0)
  }
  
  if (sym) {
    qi <- unique(sort(pmax(qi, 1-qi)))
    lqi <- NULL
    y1 <- abs(y)
  } else {
    y1 <- y
  }
  n <- length(x)
  n1 <- n %/% m
  n0 <- m * n1
  i <- sort(sample(n, n0))    
  i <- sort.list(x)[i]
  xsp2 <- array(x[i], c(m, n1))
  xt2 <- array(y1[i], c(m, n1))
  if (sym) 
    xt2 <- abs(xt2) 
  zz <- csort(xt2)
  if (flag) 
    xq2 <- zz[round(qi * m), ]
  else xq2 <- apply(xt2, 2, quantile,qi)
  xp2 <- colMeans(xsp2)
  if(!sym) { 
    xmed <- zz[round(m/2), ]
    xq2  <- t(t(xq2) - xmed)
    ymed <- predict(smooth.spline(xp2,xmed,df=df), x)$y
    y1   <- y1-ymed
  } else { 
    xmed <- rep(0, ncol(zz))
    ymed <- rep(0, n) 
  } 
  
  tt <- function(xq, qq) {
    xxtp <- smooth.spline(xp2,xq,df=df)
    w <- predict(xxtp,x)$y
    kc <- quantile(y1/w, if(qq>0.5) qq else 1-qq)       
    w <- predict(xxtp, x)$y * kc + ymed
    as.matrix(predict(smooth.spline(x,w), x0)$y)
  }
  
  vv <- function(xq, qq) {
    xxtp <- smooth.spline(xp2,xq+xmed,df=df)
    w <- smdecreasing1(xxtp,x,decreasing=down)
    kc <- quantile((y1)/(w-ymed),if(qq>0.5) qq else 1-qq)       
    zz <- (w-ymed)*kc+ymed
    w <-  zz+ quantile(-zz+y1+ymed,qq)
    as.matrix(predict(smooth.spline(x,w),x0)$y) 
  }
  
  if (type == "none"){
    if (length(qi) == 1){
      xtp <- tt(xq2, qi)
    } else {
      xtp <- tt(c(xq2[1,]),qi[1])
      for(i in 2:length(qi)) 
        xtp <- cbind(xtp, tt(c(xq2[i,]), qi[i]))
    }    
  } else { 
    if (type == "inc") down <- FALSE
    else if(type=="dec") down <- TRUE
    if (length(qi) == 1) {
      xtp <- vv(xq2, qi)
    } else {
      xtp <- vv(c(xq2[1,]), qi[1])
      for (i in 2:length(qi))
        xtp <- cbind(xtp, vv(c(xq2[i,]), qi[i]))    
    }
  }
  
  if (logtran){
    xtp <- exp(xtp)
    x0 <-  exp(x0)
    y0 <-  exp(y0)
  }
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
  dimnames(xtp) <- list(dimnames(xtp)[[1]], paste("Curve", qi, sep=""))
  cbind(T = c(y0), xtp, Sp=x0)
}

# TODO: clarify the 'Extrapolations are linear' comment (and include in a Details section of the help page)

#' This function is used by quantileCurves when its type argument specifies monotonic non-increasing curves or
#' monotonic non-decreasing curves.
#'  
#' @param z output of smooth.spline with x component sorted; if z is a plain vector, w is used as the
#'   y variable. The y component of z must be non-increasing, if not the function will make it non-increasing  
#' @param w vector of y values in case z is a plain vector of x values 
#' @param decreasing logical; if TRUE, specifies the curve to be non-increasing, if FALSE specifies the
#'   curve to be non-decreasing
#' @return vector of y values with the monotonicity constraint imposed
#' @export
smdecreasing1 <- function(z, w, decreasing = TRUE) {
  
  x <- z$x
  if(decreasing) y <- z$y else y <- (-z$y)
  rx <- range(x)
  rw <- range(w)
  n <- length(x)
  while (any(diff(y) > 0)) 
    y <- pmax(y, y[c(2:n,n)])
  i1 <- (w < rx[1])
  i2 <- (w > rx[2])
  i <- !(i1|i2)
  w1 <- w
  if (any(i)) w1[i] <- approx(x, y, w[i])$y
  if (any(i1)) { 
    m1 <- lsfit(x[1:20],y[1:20])$coef[2]
    w1[i1] <- approx( c(rw[1],x[1]),c(y[1]-m1*(x[1]-rw[1]),y[1]),w[i1])$y
  }
  if (any(i2)) { 
    m2 <- lsfit(x[n-0:19],y[n-0:19])$coef[2]
    w1[i2] <- approx( c(x[n],rw[2]),c(y[n],y[n]-m2*(x[n]-rw[2])),w[i2])$y
  }
  if (decreasing) w1 else -w1
}

#' Function for column sorting of a matrix x
#' @param x matrix
#' @return sorted matrix 
#' @export
csort <- function(x) { 
  n  <- nrow(x)
  p  <- ncol(x)
  rx <- range(x)
  z  <- (x - rx[1])/(rx[2]-rx[1])*0.99
  k <- rep(1:p,rep(n,p))
  z  <- x[sort.list(z + k)]
  dim(z) <- dim(x)
  z
}

#' This function interpolates between the smoothed quantile curves to determine the p value for a given geneset
#' @param x.ct output of quantileCurves, i.e. the curves of the contours for the different quantiles 
#' @param q.cutoff critical values
#' @return actual p values for each gene set
#' @export
ctpval1 <- function(x.ct, q.cutoff) {
  p <- ncol(x.ct) 
  pr <- 2:(p-1)
  x <- log(x.ct[,pr])
  n <- nrow(x.ct)
  p2 <- p-2
  u <- q.cutoff # as.numeric(substring(dimnames(x.ct)[[2]][pr], nch))
  y <- log(-log(1-u))
  xm <- c(x %*% rep(1,p2)/p2)
  x2 <- (x-xm)^2 %*% rep(1,p2)
  xy <- (x-xm) %*% (y-mean(y))
  b <- xy / x2
  m <- mean(y) - b*xm
  exp(-exp(m + b*log(abs(x.ct[,1]))))
}

#' This function calculates p-values for each gene set based on row permutations
#' of the gene p values or column permutations of the expression matrix; the p values
#' can be obtained either as individual gene set p values or p values based on smoothing
#' across gene sets of similar size.
#' @param geneSet is the input list of gene sets (components) and gene IDs (character vectors). 
#' A gene set can, for example, be a GO category with for each category Entrez gene identifiers; 
#' The \link{getGeneSets} function can be used to construct the geneSet argument for different
#' pathway sources.
#' @param geneStatistic is either a named numeric vector (if rowPermutations is TRUE)
#'   or a numeric matrix of pvalues (if rowPermutations is FALSE). The names of the numeric vector
#'   or row names of the matrix should represent the gene IDs.
#' @param minGenes minimum number of genes in a gene set for it to be considered (lower threshold
#'    for gene set size)
#' @param maxGenes maximum number of genes in a gene set for it to be considered (upper threshold
#'    for gene set size)
#' @param rowPermutations logical indicating whether to use row permutations (TRUE; default) or column 
#'    permutations (FALSE) 
#' @param nPermutations is the number of simulations. By default 100 permutations are conducted.
#' @param smoothPValues logical indicating whether one wants to calculate smoothed cut-off thresholds (TRUE; default)
#'    or not (FALSE).
#' @param criticalValues vector of quantiles at which p values for each gene set are desired
#' @param df degrees of freedom for the smooth.spline function used in getSmoothedPValues 
#' @param addGeneSetDescription logical indicating whether a column with the gene set description be added to
#' the output data frame; defaults to TRUE.
#' @return data frame with four (or five) columns: totalGeneSetSize, testedGeneSetSize, geneSetStatistic and geneSetPValue
#' and (if addDescription is set to TRUE) geneSetDescription; the rows of the data frame are ordered by 
#' ascending geneSetPValue.
#'
#' @references Raghavan, Nandini et al. (2007). The high-level similarity of some disparate gene expression measures,
#' Bioinformatics, 23, 22, 3032-3038.
#' @export
MLP <- function (geneSet, geneStatistic, minGenes = 5, maxGenes = 100, 
        rowPermutations = TRUE, nPermutations = 100, smoothPValues = TRUE, 
        criticalValues = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999), 
        df = 9, addGeneSetDescription = TRUE) 
{
  if (!inherits(geneSet, "geneSetMLP")) {
    stop("The 'geneSet' should be an object of class 'geneSetMLP' as produced by getGeneSets")
  }
  species <- attr(geneSet, "species")
  if (is.vector(geneStatistic)) 
    if (!rowPermutations) 
      stop("If rowPermutations is set to FALSE, the 'geneStatistic' needs to be a matrix")
  geneStatistic <- as.matrix(geneStatistic)
  if (!is.numeric(geneStatistic) & !is.matrix(geneStatistic)) {
    warning("Argument 'geneStatistic' should be a numeric vector or matrix")
  }
  if (!rowPermutations) {
    nPermutations <- ncol(geneStatistic) - 2
  }
  
  totalGeneSetSize <- sapply(geneSet, length)
  geneSets <- sapply(geneSet, function(x) x[x %in% rownames(geneStatistic)])
  
  geneSetSize <- sapply(geneSets, function(x) length(x))
  geneSetIndices <- which(geneSetSize >= minGenes & geneSetSize <= 
          maxGenes)
  geneSets <- geneSets[geneSetIndices]
  geneStatistic <- na.omit(geneStatistic)
  geneStatistic <- geneStatistic[grep("^[[:digit:]]+$", rownames(geneStatistic)), 
      , drop = FALSE]
  mapResult <- mapGeneSetStatistic(geneSets, geneStatistic)
  pValueNAlist <- lapply(mapResult, function(x) which(!is.na(x)))
  filterFunction <- function(x, y) {
    x[y]
  }
  reducedMapResult <- mapply(filterFunction, mapResult)
  reducedGeneSet <- mapply(filterFunction, geneSets)
  n2 <- nrow(geneStatistic)
  w0 <- t(mlpStatistic(reducedMapResult))
  n11 <- length(geneSets)
  if (!rowPermutations) {
    w <- mlpStatistic(reducedMapResult)
  }
  else {
    p1 <- apply(matrix(1:n2, n2, nPermutations), 2, sample)
    y1 <- matrix(geneStatistic[p1, 1], n2, nPermutations)
    row.names(y1) <- row.names(geneStatistic)
    rm(p1)
    mapResultPermuted <- mapGeneSetStatistic(reducedGeneSet, 
        y1)
    w <- t(mlpStatistic(mapResultPermuted))
  }
  if (smoothPValues) {
    if (!rowPermutations) {
      warning("Cannot smooth p-values if ind.sim = FALSE")
    }
    pw0 <- getSmoothedPValues(w0, w[, -1], q.cutoff = criticalValues, 
        df = df)
  }
  else {
    pw0 <- getIndividualPValues(w0, w)
  }
  
  res <- data.frame(testedGeneSetSize = w0[, 1], geneSetStatistic = w0[, 
          2], geneSetPValue = pw0)
  res <- res[order(res$geneSetPValue), ]
  res <- data.frame(totalGeneSetSize = totalGeneSetSize[rownames(res)], res)
  class(res) <- c("MLP", class(res))
  if (addGeneSetDescription){
    pathwaySource <- attr(geneSet, "pathwaySource")
    res <- addGeneSetDescription(object = res, pathwaySource = pathwaySource)
  }
  if (is.null(attr(res, "species"))) 
    attr(res, "species") <- species
  return(res)
}

