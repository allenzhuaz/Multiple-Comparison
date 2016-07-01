WHoch.thres <- function(p,p.set, alpha=0.05, t=0.5){
  m <- length(p); o <- order(p)
  sort.p <- p[o]; sort.p.set <- p.set[o]
  ## Two typical kinds of weights
  # first: w=F(t)
  w <- sapply(sort.p.set,function(x){max(x[x<=t],0)})

  #    w <- 1:m
  # second: w=I{F(t)>0}
  #    w <- sapply(sort.p.set,function(x){
  #          ifelse(max(x[x<=t],0)>0,1,0)
  #          })

  # third: w=I{F(t)>\lambda}
  #   lambda <- 0.1
  #    w <- sapply(sort.p.set,function(x){
  #        ifelse(max(x[x<=t],0)>lambda,1,0)
  #     })
  w0 <- min(w[w>0])
  weight.m <- sum(w)
  weight.r <- sapply(1:m, function(r){sum(w[1:r])})
  alpha.seq <- w0*alpha/(weight.m-weight.r+w0)
  R <- max(which(sort.p<=alpha.seq),0)
  return(max(alpha.seq[R],0))
}

WHoch.discrete.adjusted <- function(p,p.set,t=0.5,digits=4){
  m <- length(p)
  adjP <- c()
  for (i in 1:m){
    ##### setting the accuracy as 4 digitals
    init.alpha <- numeric(); eps <- 10^(-1*c(1:digits))

    init.alpha[1] <- 1
    while(p[i] <= WHoch.thres(p,p.set,init.alpha[1], t) & init.alpha[1] <= 1){
      init.alpha[1] <- init.alpha[1]-eps[1]
    }

    for (j in 2:digits){
      init.alpha[j] <- init.alpha[j-1]+eps[j-1]
      while(p[i] <= WHoch.thres(p,p.set,init.alpha[j], t) ){
        init.alpha[j] <- init.alpha[j]-eps[j]
      }
    }
    adjP[i] <- init.alpha[digits]+eps[digits]
  }
  return(adjP)
}


###########################################################################################
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

