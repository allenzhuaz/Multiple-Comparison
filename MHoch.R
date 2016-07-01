hyperpval <- function(n1,n0,x){
  x1 <- c(x:0)
  pval <- 1-phyper(x1-1,n1,n0,x)
  return(pval)
}

### Using summation of CDF to develop a modified Hochberg procedure
MHoch.discrete.adjusted <- function(p,p.set){
  o <- order(p); ro <- order(o); m <- length(p)
  sort.p <- p[o]; sort.p.set <- p.set[o]
  adjP <- numeric();pCDF <- matrix(NA,m,m)
  for(j in 1:m){
    for(i in m:1){
      pCDF[i,j] <- max(sort.p.set[[j]][sort.p.set[[j]] <= sort.p[i]],0)
      c <- sum(pCDF[i,i:m])
      adjP[i] <- ifelse(i==m,c,min(adjP[i+1],c))
    }
  }
  return(adjP[ro])
}
Roth.Rej <- function(p,p.set,alpha){
  m <- length(p); minP <- sapply(p.set,min)
  p_R1 <- p[minP < alpha]
  Q <-Inf
  for (j in seq_along(p_R1)){
    if(sort(p_R1,decreasing = T)[j] < alpha/j){
      Q <- j; break
    }
    else{next}
  }


  M <- c()
  for (j in seq_along(p)){
    M[j] <- sum(minP<alpha/j)
    if (M[j]<=j){
      K <- j; break
    }
    else {next}
  }
  p_RK <- if(M[K]<K) c(sort(p[minP < alpha/K],decreasing = T),rep(0,K-M[K])) else      sort(p[minP < alpha/K],decreasing = T)
  W <- Inf
  for (j in seq_along(p_RK)){
    # p_star[j] <- max(p_RK[j],p[minP<alpha/j & minP>=alpha/K])
    if(max(p_RK[j],p[minP<alpha/j & minP>=alpha/K]) < alpha/j){
      W <- j; break
    }
    else{next}
  }
  # p[p<alpha/W]
  return(c(Q,W))

}
