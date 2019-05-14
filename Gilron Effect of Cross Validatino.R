suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(gap))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(Hotelling))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(RPushbullet))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plotly))
# suppressPackageStartupMessages(library(sparsediscrim))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(LiblineaR))
suppressPackageStartupMessages(library(data.table))
# source("https://bioconductor.org/biocLite.R")
# biocLite("globaltest")
suppressPackageStartupMessages(library('globaltest'))
suppressPackageStartupMessages(library(RobPer))
# suppressPackageStartupMessages(library(fungible))
suppressPackageStartupMessages(library(kernlab))
suppressPackageStartupMessages(library(energy))
suppressPackageStartupMessages(library(HDtest))
suppressPackageStartupMessages(library(emulator))



balance.log <- 'balance_log.txt'
file.remove(balance.log)
file.create(balance.log)
the.message <- paste(file.name, Sys.info()[["nodename"]], sep=' ')





tr <- function(A) sum(diag(A))
Frob <- function(A) sum(A^2)

Euclid <- function(x) sqrt(sum(x^2))

gcd <- function(x,y) {
  r <- x%%y;
  return(ifelse(r, gcd(y, r), y))
}


matrix2list <- function(A){
  split(A, rep(1:ncol(A), each = nrow(A)))
}


ar1_cov <- function(n, rho, sigma=1){
  x <- diag(n) 
  x <- sigma * rho^abs(row(x)-col(x))
  x
  return(x)
}
## Testing
# lattice::levelplot(ar1_cov(10,0.8))


ar1_cov2 <- function(n, rho, sigma=1){
  times <- 1:n
  H <- abs(outer(times, times, "-"))
  V <- sigma * rho^H
  p <- nrow(V)
  V[cbind(1:p, 1:p)] <- V[cbind(1:p, 1:p)] * sigma
  V
}
## Testing:
# lattice::levelplot(ar1_cov2(10,0.8))


ar1_preci <- function(p, rho){
  Q = matrix(0, p, p)
  diag(Q) = 1+rho^2
  for (i in 1:(p-1)) {
    Q[i, i+1] = -rho
    Q[i+1, i] = -rho
  }
  Q[1,1] = 1
  Q[p,p] = 1
  return(Q)
}
## Testing:
# lattice::levelplot(solve(ar1_preci(10,0.8)))





# The correlation implies by a brownian motion
browninan_cov <- function(n){
  Sigma <- matrix(NA, n, n)
  for(i in 1:n){
    for(j in 1:n){
      Sigma[i,j] <- min(i,j)
    }
  }
  D <- 1/sqrt(diag(Sigma))
  Corr <- diag(D) %*% Sigma %*% diag(D)
  return(Corr)
}
## Testing
# lattice::levelplot(browninan_cov(10,0.8))



browninan_cov2 <- function(n){
  Sigma <- matrix(NA, n, n)
  for(i in 1:n){
    for(j in 1:n){
      Sigma[i,j] <- min(i/n,j/n)
    }
  }
  return(Sigma)
}
## Testing
# lattice::levelplot(browninan_cov2(10))


seq_cov <- function(p){
  Sigma0 <- 1:p
  Sigma.0.trace <- sum(Sigma0)
  # Sigma <- diag(Sigma0/Sigma.0.trace)
  Sigma <- diag(Sigma0)
  return(Sigma)
}






balanced_folding <- function(labels, n.folds, balance){
  stopifnot(is.logical(balance))
  n <- length(labels)
  
  if(balance){
    folds.list <- createFolds(y=labels, n.folds)
    fold.inds <- rep(NA, n)
    for(v in seq_along(folds.list)) fold.inds[folds.list[[v]]] <- v 
  }
  else {
    fold.inds <- sample(rep(1:n.folds, length=n))
  }
  return(fold.inds)  
}


my.summary <- function(pvals){
  cat(
    sprintf('cost=%s, permutations=%d, replications=%d, n=%d, p=%d, train.prop=%s, effect=%s, n.folds=%d', 
            cost, n.permutations, n.replications, n, p, train.prop, effect, n.folds), 
    '\n',
    apply(pvals, 2, function(x) mean(x<0.05)))
}

my.ecdf <- function(x,t){
  # ecdf(-x)(-t) #+ 1/(length(x)+1)
  mean(x>=t)
}


randomizedTest <- function(alpha, pval, pval.strickt){
  stopifnot(pval>=pval.strickt)
  reject <- NULL
  if(pval.strickt> alpha) reject <- FALSE
  else if (pval<alpha) reject <- TRUE
  else {
    p.diff <- (alpha-pval.strickt)/(pval- pval.strickt)
    reject <- rbinom(1, 1, p.diff) %>% as.logical
  }
  return(reject)
}
## Testing:
# sum(replicate(1e3, randomizedTest(0.05, 0.06, 0.03)))


# statistic.levels <- c("Oracle", "Hotelling", "Hotelling.shrink", "Goeman", "sd", "MMD","dCOV", 
#                       "lda.CV.1", "lda.noCV.1", "svm.CV.1", "svm.CV.2", "svm.noCV.1", 
#                       "svm.noCV.2")


source('statistics.R')

makeSymmetric <- function(n){
  x <- matrix(rnorm(n*n), n) 
  ind <- lower.tri(x) 
  x[ind] <- t(x)[ind] 
  x
}
## Testing:
# makeSymmetric(3)



getIp <- function(){
  x <- system("ifconfig", intern=TRUE)
  x[grep("IPv4", x)]
}
## Testing
# getIp()

my.outer <- function(x){
  xx <- outer(x,x)
  xx[lower.tri(xx, diag = TRUE)]
}
##Testing:
# my.outer(1:3)

augmentDesign <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  x.aug <- cbind(x,t(apply(x,1,my.outer)))
  scale(x.aug)
}
## Testing:
# x <- matrix(1:9, ncol=3, byrow = TRUE)
# augmentDesign(x)

my.quad.form <- function(x,A){
  apply(x, 2, function (y) quad.form(A,y))
}
## Testing:


makeLogisticSetup <- function(p, beta00, B00, b00=-1){
  Sigma <<- diag(p)
  Sigma.inv <<- solve(Sigma)
  b0 <<- b00
  beta0 <<- rep(beta00,p)
  # B0 <<- matrix(B00, ncol = p, nrow = p)
  B0 <<- diag(B00, p)
}
## Testing


makeLogisticProbs <- function(beta0, B0, effect, noise){
  beta <- beta0*effect
  B <- B0*effect
  link0 <- noise%*%beta + colSums(t(noise) * (B %*% t(noise)))
  link <- link0-median(link0)
  probs <- plogis(link)
  return(probs)
}
## Testing


makeLogisticLabels <- function(probs,n, balance=FALSE){
  labels <- c(TRUE,FALSE)[rbinom(n,1, probs)+1]
  if(balance){
    while(sum(labels)<10 || sum(labels)>30){
      labels <- c(TRUE,FALSE)[rbinom(n,1, probs)+1] # group assignemt labels
    }
  }
  return(labels)
}