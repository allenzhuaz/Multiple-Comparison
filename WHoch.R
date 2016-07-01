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
