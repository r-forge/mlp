###Note: 05-12-2009: Deleted GO Annotation functions. Instead use RFunctions_MLPSetup.txt 

# File: RFunctions_MLP.txt
# Author: Nandini Raghavan. Other Contributing Authors: Javier Cabrera, Dhammika Amaratunga.
# Last modified: March19, 2008

# Description: This code contains the necessary R functions to
#(i) create the necessary annotation data and expression data inputs 
#(ii) do an MLP analysis with options to use either gene or sample permutations. 
# Functions included here: 
# f.AffyP2GCleanup, f.GOTaxChip,f.GONames, f.AffyP2G,
# f.u, f.mlp0, f.mlp, f.permtwo, f.permk (incomplete),f.pperm
# f.rt, f.rtt, f.rsd,f.rsp,f.rspp, f.rvar, f.toarray, f.qnorm
# f.y2x, f.pp2g,
# f.cut0,f.cut1,f.cut2, f.ee1, f.smdecreasing1, f.ctpval1,f.csort
# f.GoFiSH
# f.GO2Name,f.GO2Gene, f.GO2G, f.GO2P,f.G2P,f.G2GO,f.P2G,f.P2GO




###====================== MLP STATISTIC ========================

f.u <- function(x){ 
  colMeans(-log10(x)) 
}

f.mlp0 <- function(x, y, y2x){
  ###This function calculates the mean of the gene-statistic y[,2] for each unique gene-set in x[,1].
  ###x is the input table of gene-sets and gene-names.
  ###y is the input table of gene-names and gene-statistics.
  ###y2x is the mapping of the gene-names in y to those in x.
  ###t1 returns the gene-sets and the mean statistic for the gene-set.
  
  n   <- table(x[,1])
  ###s1  <- split(y[y2x,-1],x[,1])
  ###u   <- t(sapply(s1,function(x){f.u(as.matrix(x))}))
  s1  <- split(y[y2x,-1],x[,1])
  u  <- t(sapply(s1, function(x){f.u(matrix(x, ncol=ncol(y)-1))}))
  ###dimnames(u)[[2]] <- dimnames(s1[[1]])[[2]]
  if (ncol(y) == 2){
    return(cbind(n, as.vector(u)))
  } else {
    return(cbind(n, u))
  }
}


f.mlp <- function(inputData, y, y2x){
  ###This function calculates the mean of the gene-statistic y[,2] for each unique gene-set in x[,1].
  ###x is the input table of gene-sets and gene-names.
  ###y is the input table of gene-names and gene-statistics.
  ###y2x is the mapping of the gene-names in y to those in x.
  ###t1 returns the gene-sets and the mean statistic for the gene-set.
  
  number   <- table(inputData[,1])
  ###s1  <- split(y[y2x,-1],x[,1])
  ###u   <- t(sapply(s1,function(x){f.u(as.matrix(x))}))
  s1  <- split(f.toarray(y[y2x,-1]), inputData[,1])
  u  <- t(sapply(s1,function(x){f.u(matrix(x,ncol=ncol(y)-1))}))
  ###dimnames(u)[[2]] <- dimnames(s1[[1]])[[2]]
  if (ncol(y) == 2){
    return(cbind(number, as.vector(u)))
  } else {
    return(cbind(number,u))
  }
}

###=============================== EXACT PERMUTATIONS===========================

f.permtwo = function (n, k) 
{
  ###From Javier 1/27/2006 to calculate all unique permutations for 1:n, with k in each group.
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


f.permk <- function(gr) 
{
  ###From Javier 1/27/2006 # To calculate all the permutations
  ff <- cumsum(gr)
  ttt <- f.permtwo(ff[2],ff[1])
  for(i in 3:length(gr)) 
  {
    tt <- ttt
    tt1 <- f.permtwo(ff[i],ff[i-1])
    ttt <- NULL
    for(i in ff[i-1]+1:gr[i]) tt <- cbind(tt,i)
    for(i in 1:nrow(tt1)) ttt <- rbind(ttt,cbind(tt,6,7)[,tt1[i,]])
  }
  ttt
}

###======================== MLP-CPERM===================================

f.pperm <- function(x, p0, p1, ind.exact = TRUE, np = 0){
  ### x is the expression data for which the original and column-permuted p-values are calculated.
  ### p0 and p1 are the original "treatment" and "control" columns.
  ### The original permutation is identified and removed from the simulated data.
  ### ind.exact: should exact permuations be used? logical defaulting to \code{TRUE}
  ### np is the number of random permutations, if ind.exact!=1.
  ### Example: ypp <- f.pperm(dma.norm, 1:4, 5:8, ind.exact = 1, 1))
  ### Example: ypp <- f.pperm(dma.norm, 1:4, 5:8, ind.exact = 0, 100))
  
  
  ### Calculating the p-value of the original t-statistic and replacing NA's with 1.
  n0 <- length(p0)
  n  <- ncol(x)
  nr <- nrow(x)
  
  ### Creating the random permutations:
  if (ind.exact){
    ### matrix of exact column permutations 
    pj <- f.permtwo(n,n0)
    p0p1String <- paste(c(p0, p1), collapse = "")
    pj <- pj[-pmatch(p0p1String,
            apply(pj, 1, function(x) paste(x, collapse = ""))),]
    np <- nrow(pj)
    pp <- matrix(NA,nr,(np+1)) 
    for (j in 1:np){
      j0 <- pj[j,][1:n0]
      j1 <- pj[j,][(n0+1):n]
      pp[,(j+1)] <- f.rt(x[,j0],x[,j1])[,4]
      
    }
  } else {  
    pp <- matrix(NA, nr, (np+1))
    for (j in 1:np){
      pj <- sample(n)
      
      p0p1String <- paste(c(p0, p1), collapse = "")
      pjString <- paste(pj, collapse = "")
      
      if (!is.na(pmatch(p0p1String, pjString))){
        pj <- sample(n)
      }
      j0 <- pj[1:n0]
      j1 <- pj[(n0+1):n]
      pp[, (j+1)] <- f.rt(x[,j0], x[,j1])[,4]
    }}
  
  pp[,1] <- f.rt(x[,p0],x[,p1])[,4]
  pp[is.na(pp)] <- 1
  
  work <- paste("P", 0:np, sep = "")
  dimnames(pp) <- list(I(dimnames(x)[[1]]), work)
  
  return(pp)
}


###==============================UTILITY FUNCTIONS=======================

f.rt <- function (x, y,var.thresh=10^-4,num.thresh=10^(-6),denom.thresh=10^(-6)) {
  nx <- (!is.na(x)) %*% rep(1, ncol(x))
  ny <- (!is.na(y)) %*% rep(1, ncol(y))
  t.numer <- (rowMeans(x) - rowMeans(y))
  t.denom <- f.rssp(x, y)
  
  i0 <- f.rvar(cbind(x, y)) < var.thresh
  i1 = abs(t.numer ) < num.thresh & t.denom < denom.thresh
  i <- i0 | i1
  t.numer[i] <-  0
  t.denom[i] <- 1
  t.stat <- t.numer/t.denom
  t.pval <- 2 * (1 - pt(abs(t.stat), nx + ny - 2))
  tt <- matrix(cbind(t.numer, t.denom, t.stat, t.pval), ncol = 4)
  dimnames(tt)[2] <- list(c("t.numer", "t.denom", "t.stat", 
          "t.pval"))
  return(tt)
}

f.rtt <- function(x1, x2, fu=0, var.thresh=10^-4, num.thresh=10^(-6), denom.thresh=10^(-6)) {
  i0 <- f.rvar(cbind(x1,x2)) < var.thresh
  num <- (rowMeans(x1) - rowMeans(x2))
  den <- (t(array(fu, c(nrow(x1), length(fu)))) + f.rssp(x1, x2))
  i1 <- abs(num) < num.thresh & den < denom.thresh
  i <- i0 | i1
  num[i] <- 0
  den[i] <- 1
  num / den
}

f.rsd <- function(x) sqrt(f.rvar(x))
f.rsp <- function(x1, x2) sqrt((f.rvar(x1) + f.rvar(x2))/2)

f.rvar <- function(x,n=ncol(x)) c((x-rowMeans(x))^2 %*%rep(1,n))/(n-1)
f.toarray <- function(x) array(unlist(x) , dim(x), dimnames(x))
f.rssp <- function(x1,x2) f.rsp(x1,x2)*sqrt(1/ncol(x1) +1/ncol(x2))


###==============================MAPPING FUNCTIONS=======================

###Pattern match (application: to see which row of permuted numbers corresponds to observed sequence)
# f.Vec2String <- function(x){
#  tmp0 <- NULL
#  for (j in 1:length(x)){
#    tmp0 <- paste(tmp0, x[j], sep = "")
#  }
#  return(tmp0)
#}
#### TV: not needed: plain paste(., collapse = "")

f.y2x <- function(x, y){
  ###Function to map the character vector y to character vector x.
  ###CAUTION: Only works correctly for unique entries in y. (07-28-2006)
  x <- as.character(x)
  y <- as.character(y)
  j <- seq_along(y)
  names(j) <- y
  y2x <- j[x]    
  return(y2x)
} 


f.pp2g <- function(xp,yp){
  ###Function to convert probeset datasets to gene datasets
  ###xp is the input table of gene-sets, probe.ids, and gene.ids, in that order.
  ###yp is the input table of probe.ids and original and permuted probeset-statistics.
  
  x1  <- xp[,1]
  x2  <- xp[,2]
  x3  <- xp[,3]
  
  pid <- yp[,1]
  yp  <- yp[,-1]
  dimnames(yp)[[1]] <- pid  ###need for filtering for y downstream.
  
  y2xp <- f.y2x(xp[,2],pid) ###map p-values in yp to probesets in xp.
  p1   <- data.frame(Probe.ID=pid[y2xp],P0=yp[y2xp,2])
  
  
  ###Create gene-statistic table; take min p-value for all probesets for a given gene.
  x30 <- sort(unique(x3))
  x30 <- x30[x30 !=""]
  n30 <- length(x30)
  p30 <- vector("character",n30)
  
  for (j in 1:n30){ 
    ###p30[j] <- (min(p1[x3==x30[j],2],na.rm=T) )     ### min p-value
    p30[j] <- p1[ p1[,2]==min(p1[x3==x30[j],2]),1][1] ###probe-set corr to min p-value
    ###need checksum if genes have multiple mins.
  }
  y <- data.frame(Gene.ID=I(x30),yp[p30,])
  
  
  ###Create unique(GO.#, Gene.ID) combination table.
  y2x3<- f.y2x(x3,x30)
  x41 <- x1*10^10 + y2x3 ###(GO number inflated  + gene number)
  n41 <- length(x41)
  x40 <- sort(unique(x41)) ###multiple probe-sets corr. to gene can cause duplication in x41.
  n40 <- length(x40)  
  jj <- vector("numeric",n40)
  for (j in 1:n40){
    jj[j] <- (1:n41)[x41==x40[j]][1]
  }
  x <- data.frame(GO.Number=x1[jj],Gene.ID=I(x3[jj]))
  
  out <- list(x=x, y=y)
  return(out)
}


###=======================================CUTOFFS========================

f.cut0 <- function(w0, w){
  ### individual cut-offs
  w1 <- matrix(rep(w0[,2], ncol(w)), nrow = nrow(w0), ncol = ncol(w))
  dw0w1 <-  ifelse(w1-w > 0, 1, 0)
  pw0   <- 1 - rowMeans(dw0w1) 
  return(pw0)
}

f.cut1 <- function(w0,w,q.cutoff){
  ### This imposes no criterion and has df = 15
  ### SMOOTHED CUTOFFS :Uses f.ee, f.ctpval and related functions from Javier's code. 
  w1 <- cbind(rep(w0[,1],ncol(w[,])),as.vector(w[,]))
  lqi = NULL; hqi= q.cutoff
  work <- f.ee1(x=sqrt(w1[,1]), y=w1[,2], x0 = sqrt(w0[,1]), y0 = w0[,2],type=c("dec", "inc", "none")[3], m = 20, lqi = lqi, hqi=hqi, sym = FALSE, plot = TRUE, flag = FALSE, dg = 15, logtran = FALSE)
  pval <- f.ctpval1(work,nch=6)
  rm(work)
  return(pval)
}


f.cut2 <- function(w0,w,q.cutoff){
  ###This imposes the decreasing criterion and has df =9
  ###SMOOTHED CUTOFFS :Uses f.ee, f.ctpval and related functions from Javier's code. 
  w1 <- cbind(rep(w0[,1],ncol(w[,])),as.vector(w[,]))
  lqi = NULL; hqi= q.cutoff
  ###work <- f.ee1(x=sqrt(w1[,1]), y=w1[,2], x0 = sqrt(w0[,1]), y0 = w0[,2],type=c("dec", "inc", "none")[1], m = 20, lqi = lqi, hqi=hqi, sym = F, plot = T, flag = F, dg = 15, logtran = F)
  work <- f.ee1(x=sqrt(w1[,1]), y=w1[,2], x0 = sqrt(w0[,1]), y0 = w0[,2],type=c("dec", "inc", "none")[1], m = 20, lqi = lqi, hqi=hqi, sym = FALSE, plot = TRUE, flag = FALSE, dg = 9, logtran = FALSE)
  pval <- f.ctpval1(work,nch=6)
  rm(work)
  return(pval)
}


f.ee1 <- function(x, y, x0 = x, y0 = y, type = c("none", "dec", "inc"), m = 20, lqi = 0.05, hqi = 0.95,
    sym = FALSE, plot = TRUE, flag = FALSE, dg = 15, logtran=FALSE) {
#  x, y are the coordinates of the points used to generate the
#       envelope.
#  x0, y0 are the points that are graphed with the envelop.
#        By default x0=x y0=y
#  type:  constraints to the envelop function, one of increasing ("inc"), decreasing ("dec") or none ("none", default). 
#  m: block size for constructing the quantile envelop
  
  type <- match.arg(type)
  
  qi <- c(sort(lqi), sort(hqi))
  if (logtran){
    x <- log(x)
    x0 <- log(x0)
    y <- log(y)
    y0 <- log(y0)
  }
  
  if (sym) {
    qi <- unique(sort(pmax(qi,1-qi)))
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
  zz <- f.csort(xt2)
  if (flag) 
    xq2 <- zz[round(qi * m), ]
  else xq2 <- apply(xt2, 2, quantile,qi)
  xp2 <- colMeans(xsp2)
  if(!sym) { 
    xmed <- zz[round(m/2), ]
    xq2  <- t(t(xq2) - xmed)
    ymed <- predict(smooth.spline(xp2,xmed,df=dg), x)$y
    y1   <- y1-ymed
  } else { 
    xmed <- rep(0, ncol(zz))
    ymed= rep(0, n) 
  } 
  
  tt <- function(xq, qq) {
    xxtp <- smooth.spline(xp2,xq,df=dg)
    w <- predict(xxtp,x)$y
    kc <- quantile(y1/w, if(qq>0.5) qq else 1-qq)       
    w <- predict(xxtp, x)$y * kc + ymed
    as.matrix(predict(smooth.spline(x,w), x0)$y)
  }
  
  vv <- function(xq, qq) {
    xxtp <- smooth.spline(xp2,xq+xmed,df=dg)
    w <- f.smdecreasing1(xxtp,x,decreasing=down)
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
      xtp <- vv(xq2,qi)
    } else {
      xtp <- vv(c(xq2[1,]), qi[1])
      for(i in 2:length(qi))
        xtp <- cbind(xtp, vv(c(xq2[i,]), qi[i]))    
    }
  }
  
  if (logtran){
    xtp <- exp(xtp)
    x0 <-  exp(x0)
    y0 <-  exp(y0)
  }
  plot(x0, y0, xlab="n", ylab="MLP", axes=FALSE, col="#08306B", pch=".")
  axis(2, lwd = 1.5, col="#08306B")
  atPositions <-(axis(1, labels = FALSE))
  axis(1,lwd = 1.5, at = atPositions, labels = atPositions^2, las = 1, col = "#08306B")
  par(cex = 0.5)
  npc <- 16
  ccc <- "#2171B5"
  if (!is.null(hqi)){ 
    i2 <- y0 > xtp[,length(qi)]
    print(length(i2[i2==TRUE]))
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

f.smdecreasing1 = function(z, w, decreasing=TRUE) {
# z ourput of smooth.spline with x component sorted.
# The y component of z must be non increasing, if not the 
#    function willmake it non-increasing
# Extrapolations are linear
  x <- z$x
  if(decreasing) y <- z$y else y <- (-z$y)
  rx <- range(x)
  rw <- range(w)
  n <- length(x)
  while (any(diff(y) > 0)) 
    y <- pmax(y,y[c(2:n,n)])
  i1 <- (w < rx[1])
  i2 <- (w > rx[2])
  i <- !(i1|i2)
  w1 <- w
  if(any(i)) w1[i] <- approx( x,y,w[i])$y
  if(any(i1)) { 
    m1 <- lsfit(x[1:20],y[1:20])$coef[2]
    w1[i1] <- approx( c(rw[1],x[1]),c(y[1]-m1*(x[1]-rw[1]),y[1]),w[i1])$y
  }
  if(any(i2)) { 
    m2 <- lsfit(x[n-0:19],y[n-0:19])$coef[2]
    w1[i2] <- approx( c(x[n],rw[2]),c(y[n],y[n]-m2*(x[n]-rw[2])),w[i2])$y
  }
  if (decreasing) w1 else -w1
}

f.csort = function(x) { 
  n  <- nrow(x)
  p  <- ncol(x)
  rx <- range(x)
  z  <- (x - rx[1])/(rx[2]-rx[1])*0.99
  k <- rep(1:p,rep(n,p))
  z  <- x[sort.list( z + k)]
  dim(z) <- dim(x)
  z
}

f.ctpval1 = function(x.ct,nch=6) {
  p <- ncol(x.ct) 
  pr <- 2:(p-1)
  x <- log(x.ct[,pr])
  n <- nrow(x.ct)
  p2 <- p-2
  u <- as.numeric(substring(dimnames(x.ct)[[2]][pr], nch))
  y <- log(-log(1-u))
  xm <- c(x%*%rep(1,p2)/p2)
  x2 <- (x-xm)^2 %*%rep(1,p2)
  xy <- (x-xm)%*%(y-mean(y))
  b <- xy / x2
  m <- mean(y) - b*xm
  exp(-exp(m + b*log(abs(x.ct[,1]))))
}


###==================================SIMULATION========================


f.GOFiSH <- function(x, y, nmin = 5, nmax = 100, ind.sim = TRUE, nsim = 10, ind.smooth = TRUE){
  
  if (!is.numeric(x)){
    warning("Input x should be numeric")
    x <- f.toarray(x)
  }
  if (!is.numeric(y)){
    warning("Input y should be numeric")
    y <- f.toarray(y)
  }
  
  ### x is the input (numeric)  table of gene-sets and gene-ids.
  ### y is the input (numeric) table of gene-ids and gene-pvalues (or other summary statistic).
  ### p is the user-defined threshold for significance. 
  ### nsim is the number of simulations. 
  ### ind.p is the indicator for calculating the p-value.
  ### ind.smooth calculates smoothed cut-off thresholds.
  
  if (!ind.sim){
    nsim <- ncol(y)-2
  }
  
  ### Remove missing values from x and y.
  y <- y[!is.na(y[,2]), ]
  y2x <- f.y2x(x[,2], y[,1])
  x <- x[!is.na(y2x), ]
  y2x <- f.y2x(x[,2], y[,1])
  
  n1 <- length(unique(x[,1]))
  n2 <- nrow(y)
  
  w0 <- f.mlp0(x,y[,1:2],y2x)
  i <- ((w0[,1] <= nmax) & (w0[,1] >= nmin))
  # Only relevant rows of x to reduce downstream computations.
  ix   <- i[match(x[,1],names(i))]
  y2xi <- f.y2x(x[ix,2],y[,1])
  n11 <- length(unique(x[ix,1]))
  
  ### Create simulation data:
  if (!ind.sim){
    ###Column Permutations
    w <- f.mlp(x[ix,], y[,-2], y2xi)[,-1]
  } else {
    ### Row Permutations
    p1 <- apply(matrix(1:n2,n2,nsim),2,sample) ##need "matrix" to resample each column separately.
    y1 <- data.frame(Gene=I(y[,1]),matrix(y[p1,2],n2,nsim))
    rm(p1)
    w <- f.mlp(x[ix,],y1,y2xi)[,-1]
  }
  
  ### Determine Cut-offs:
  if (ind.smooth){
    if (!ind.sim){
      warning("Cannot smooth p-values if ind.sim = FALSE")
    }
    pw0 <- f.cut2(w0[i,], w, q.cutoff = c(0.5,0.9,0.95,0.99,0.999,0.9999,0.99999))
  } else {
    pw0 <- f.cut0(w0[i,], w)
  }
  
  res <- data.frame(Geneset.Size = w0[i,1], MLP.Observed = w0[i,2], P.Value = pw0)
  return(res)
}


###===============================OUTPUT RETRIEVAL FUNCTIONS========================

f.GO2Name <- function(ANames, go.out, ind.p = TRUE, p.cutoff = 0.01){
  if (ind.p){
    ii <- row.names(go.out[go.out[,3]<= p.cutoff, ])
  } else{
    ii <- row.names(go.out)
  }
  res <- data.frame(go.out[ii,], Geneset.Description = I(ANames[ii]))
  return(res)
}

f.P2G <- function(p0, G){
  return(G[p0])
}

f.G2GO <- function(g0, A){
  lg0 <- length(g0)
  GO0 <- vector("list", lg0)
  names(GO0) <- g0
  for (j in 1:lg0){ 
    GO0[[j]] <- A[A[,2]==g0[j], 1]
    names(GO0[[j]]) <- NULL
  }
  return(GO0)
}

f.P2GO <- function(p0, G, A){
  g0 <- f.P2G(p0, G)
  g0 <- sort(unique(g0))
  GO0 <- f.G2GO(g0, A)
  return(GO0)
}

f.GO2G <- function(G0, A){
  lG0 <- length(G0)
  g0 <- vector("list", lG0)
  names(g0) <- G0
  for (j in 1:lG0){ 
    g0[[j]] <- A[A[,1]==G0[j],2]
    names(g0[[j]]) <- NULL
  }
  return(g0)
}

f.G2P <- function(g0, G){
  p0 <- NULL
  for (j in g0){
    p0 <- c(p0, names(G[G==j]))
  }
  return(p0)
}

f.GO2P <- function(G0, G, A){
  g0 <-  unlist(f.GO2G(G0, A))
  g0 <- sort(unique(g0))
  p0 <- f.G2P(g0, G)
  return(p0)
}


f.GO2Gene <- function(x, y, go.out, go.pthreshold = 0.01, gene.pthreshold = 0.05, ind.log = FALSE, fname = "foo.pdf"){

  ###Example:
  go  <- go.out[go.out[,3] < go.pthreshold, ]
  n <- go[,1]
  go <- as.numeric(row.names(go))
  
  i <- order(go)
  go <- go[i]
  n <- n[i]
   
  pdf(file = fname, bg = "white")
  
  par(mfrow=c(2, 3))
  
  xjj <- NULL
  for (j in go){
    i <- seq_along(go)[go==j]
    xj <- x[x[,1]==j, ]
    xj <- cbind(xj, PVal = y[match(xj[,2], y[,1]), 2])
    xj <- xj[order(xj[,3]),]
    xj <- xj[xj[,3] < gene.pthreshold, ]
    if (ind.log)
      xj[,3] <- -log10(xj[,3])
    barplot(xj[,3], names = xj[,2], las = 2, cex.names = 0.6)
    title(main=paste("GO :",j, 
            "n.Genes :", n[i], sep=" "), cex.main=0.75)
    xjj <- rbind(xjj, xj)
  }
  dev.off()
  dimnames(xjj) <- list(1:nrow(xjj), c("GO", "Gene", "PValue"))
  return(xjj)
}
