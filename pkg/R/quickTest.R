# TODO: Add comment
# 
# Author: Tobias Verbeke
###############################################################################


f.permtwo = function (n, k) 
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

f.permtwo(n = 5, k = 2)