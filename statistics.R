
t_stat <- function(x,y){
  mean.x <- mean(x)
  mean.y <- mean(y)
  var.x <- var(x)
  var.y <- var(y)
  n.x <- length(x)
  n.x.1 <- n.x-1
  n.y <- length(y)
  n.y.1 <- n.y-1
  var.pooled <- (n.x.1 * var.x + n.y.1 * var.y)/(n.x.1+n.y.1)
  nom <- (mean.y-mean.x)^2
  denom <- var.pooled * (1/n.x + 1/n.y)
  t.stat <- nom / denom
  return(t.stat)
}
## Testing:
# t_stat(rnorm(10), rnorm(20,2))


# A cross validation wrapper.
# Constructed to allow validation within a permutation test (thus- "old.labels").
t_cv <- function(FUN, noise, new.labels, old.labels, fold.ids, ...){
  n.folds <- max(fold.ids)
  t.cv <- rep(NA, n.folds) # initialize output container
  
  for(v in 1:n.folds){
    # v <- 1
    test.ind <- fold.ids==v
    train.ind <- !test.ind
    
    train.noise <- noise[train.ind, ]
    train.labels <- new.labels[train.ind]
    test.noise <- noise[test.ind, ]
    test.labels <- old.labels[test.ind] 
    
    t.cv[v] <- FUN(train.noise, train.labels, test.noise, test.labels, ...)
  }
  result <- mean(t.cv)
  return(result) # average accuracy over folds
}







t_Oracle <- function(x,y,S.inv){
  x.bar <- colMeans(x)
  y.bar <- colMeans(y)
  delta <- x.bar - y.bar
  
  T2 <- delta %*% S.inv %*% delta
  return(as.numeric(T2))
}
### Testing:
# t_Oracle(x1,x2,Sigma.aug.inv)


t_Oracle_NP <- function(x,y, mu.x, mu.y, S){
  
  S.inv <- solve(S)
  x.bar <- colMeans(x)
  y.bar <- colMeans(y)
  
  l1.y <- (y.bar - mu.y) %*% S.inv %*% (y.bar - mu.y)
  l0.y <- (y.bar - mu.x) %*% S.inv %*% (y.bar - mu.x)
  
  T2 <- l1.y-l0.y
  return(T2)
}
### Testing:
# t_Oracle_NP(x = rmvnorm(1e2, rep(0,1e1)), 
            # y = rmvnorm(1e2, rep(1,1e1)),
            # mu.x=0, mu.y=1,
            # S=diag(1e1))




t_Hotelling <- function(x,y, shrinkage){
  T2 <- hotelling.stat(x,y, shrinkage=shrinkage)$statistic
  return(T2)
}
## Testing:
# t_Hotelling(x = rmvnorm(1e2, rep(0,1e1)), 
#             y = rmvnorm(1e2, rep(0,1e1)), 
#             shrinkage = TRUE)

# Goeman's high-dim test
t_goeman <- function(x,y){
  X <- rbind(x,y)
  dimnames(X) <- list(NULL,LETTERS[1:ncol(X)])
  y <- as.matrix(c(rep(FALSE, nrow(x)),rep(TRUE,nrow(y))))
  result <- globaltest::gt(y,X)
  result@result[,'Statistic'] %>% unname()
}
## Testing
# t_goeman(x1,x2)




# Linear SVM. Arguments self explanatory.
t_svm <- function(train.noise, train.labels, test.noise, test.labels, cost, type){
  svm.1 <- svm(x=train.noise, y=train.labels, type='C-classification', kernel='linear', cost=cost)
  accuracy <- mean(predict(svm.1, newdata=test.noise)==test.labels)
  
  if(type==1){
    statistic <- accuracy
  }
  else if(type==2){
    .p <- mean(test.labels)
    p <- max(.p,1-.p)
    statistic <- abs(accuracy-p)/sqrt(p*(1-p))
  }
  return(statistic)
}


t_svm_Gauss <- function(train.noise, train.labels, test.noise, test.labels, cost, type){
  svm.1 <- svm(x=train.noise, y=train.labels, type='C-classification', cost=cost)
  accuracy <- mean(predict(svm.1, newdata=test.noise)==test.labels)
  
  if(type==1){
    statistic <- accuracy
  }
  else if(type==2){
    .p <- mean(test.labels)
    p <- max(.p,1-.p)
    statistic <- abs(accuracy-p)/sqrt(p*(1-p))
  }
  return(statistic)
}


t_svm_cv_Gaus <- function(noise, new.labels, old.labels, fold.ids, cost, type){
  t_cv(FUN = t_svm, noise, new.labels, old.labels, fold.ids, cost, type)
}
## Testing:
# t_svm_cv_Gaus(noise, labels, labels, fold.ids, cost=100, type=1)
# svm(x=noise,y=labels, type='C-classification', cost=100)



# Compute cross validated test statistics
# noise: the predictors 
# labels: 
# fold.ids: the asignemt of observations to folds
# cost: the svm cost parameter. 
t_svm_cv <- function(noise, new.labels, old.labels, fold.ids, cost, type){
  t_cv(FUN = t_svm, noise, new.labels, old.labels, fold.ids, cost, type)
}
## Testing:
# t_svm_cv(noise, labels, noise, labels, cost=100, type=1)
# svm(x=noise,y=labels, type='C-classification', cost=100)





# Linear SVM. Arguments self explanatory.
t_svml2 <- function(train.noise, train.labels, test.noise, test.labels, cost, type){
  # svm.1 <- glmnet(x=train.noise, y=train.labels, family = 'binomial', alpha = 0)
  svm.1 <- LiblineaR(data=train.noise, target =train.labels, type=1, cost=cost)
  predict.labels <- predict(svm.1, newx=test.noise, type='class')$predictions
  accuracy <- mean(predict.labels==test.labels)
  
  if(type==1){
    statistic <- accuracy
  }
  else if(type==2){
    .p <- mean(test.labels)
    p <- max(.p,1-.p)
    statistic <- abs(accuracy-p)/sqrt(p*(1-p))
  }
  return(statistic)
}

# Compute cross validated test statistics
# noise: the predictors 
# labels: 
# fold.ids: the asignemt of observations to folds
# cost: the svm cost parameter. 
t_svml2_cv <- function(noise, new.labels, old.labels, fold.ids, cost, type){
  t_cv(FUN = t_svml2, noise, new.labels, old.labels, fold.ids, cost, type)
}



# Compute cross validated test statistics with cross validated regularization
# noise: the predictors 
# labels: 
# fold.ids: the asignemt of observations to folds
# cost: the svm cost parameter. 
t_svm_cvCV <- function(noise, new.labels, old.labels, fold.ids){
  t_cv(FUN = t_svm_cvLambda, noise, new.labels, old.labels, fold.ids)
}
## Testing:
# t_svm_cvCV(noise, labels, labels, fold.ids)



t_svm_cvLambda <- function(train.noise, train.labels, test.noise, test.labels, type=1){
  svm.1 <- best.svm(x=train.noise, y=train.labels, type='C-classification', kernel='linear', cost = 10^(-3:3))
  accuracy <- mean(predict(svm.1, newdata=test.noise)==test.labels)
  # capture.output(svm.1$cost, file='svmCV.txt', append=TRUE)
  
  if(type==1){
    statistic <- accuracy
  }
  else if(type==2){
    .p <- mean(test.labels)
    p <- max(.p,1-.p)
    statistic <- abs(accuracy-p)/sqrt(p*(1-p))
  }
  return(statistic)
}
## Testing:
# t_svm_cvLambda(train.noise, train.labels, test.noise, test.labels)





# Linear Discriminant analysis
t_lda <- function(train.noise, train.labels, test.noise, test.labels, type){
  lda.1 <- lda(x=train.noise, grouping =train.labels)
  predictions <- predict(lda.1, newdata=test.noise)$class
  accuracy <- mean(predictions==test.labels)
  
  if(type==1){
    statistic <- accuracy
  }
  else if(type==2){
    .p <- mean(test.labels)
    p <- max(.p,1-.p)
    statistic <- abs(accuracy-p)/sqrt(p*(1-p))
  }
  return(statistic)
}

t_lda_cv <- function(noise, new.labels, old.labels, fold.ids, type){
  t_cv(FUN = t_lda, noise, new.labels, old.labels, fold.ids, type)
}


t_dlda <- function(train.noise, train.labels, test.noise, test.labels, type){
  lda.1 <- dlda(x=train.noise, y =train.labels)
  predictions <- predict(lda.1, newdata=test.noise)$class
  accuracy <- mean(predictions==test.labels)
  
  if(type==1){
    statistic <- accuracy
  }
  else if(type==2){
    .p <- mean(test.labels)
    p <- max(.p,1-.p)
    statistic <- abs(accuracy-p)/sqrt(p*(1-p))
  }
  return(statistic)
}

t_dlda_cv <- function(noise, new.labels, old.labels, fold.ids, type){
  t_cv(FUN = t_dlda, noise, new.labels, old.labels, fold.ids, type)
}





t_hdrda <- function(train.noise, train.labels, test.noise, test.labels, type){
  lda.1 <- hdrda(x=train.noise, y =train.labels)
  predictions <- predict(lda.1, newdata=test.noise)$class
  accuracy <- mean(predictions==test.labels)
  
  if(type==1){
    statistic <- accuracy
  }
  else if(type==2){
    .p <- mean(test.labels)
    p <- max(.p,1-.p)
    statistic <- abs(accuracy-p)/sqrt(p*(1-p))
  }
  return(statistic)
}



t_hdrda_cv <- function(noise, new.labels, old.labels, fold.ids, type){
  t_cv(FUN = t_hdrda, noise, new.labels, old.labels, fold.ids, type)
}





t_sdlda <- function(train.noise, train.labels, test.noise, test.labels, type){
  lda.1 <- sdlda(x=train.noise, y =train.labels)
  predictions <- predict(lda.1, newdata=test.noise)$class
  accuracy <- mean(predictions==test.labels)
  
  if(type==1){
    statistic <- accuracy
  }
  else if(type==2){
    .p <- mean(test.labels)
    p <- max(.p,1-.p)
    statistic <- abs(accuracy-p)/sqrt(p*(1-p))
  }
  return(statistic)
}



t_sdlda_cv <- function(noise, new.labels, old.labels, fold.ids, type){
  t_cv(FUN = t_sdlda, noise, new.labels, old.labels, fold.ids, type)
}





t_SD <- function(x,y){
  x.bar <- colMeans(x)
  y.bar <- colMeans(y)
  delta <- x.bar - y.bar
  
  Sx <- cov(x)
  Sy <- cov(y)
  
  nx <- nrow(x)
  ny <- nrow(y)
  n <- nx + ny
  nx1 <- nx-1
  ny1 <- ny-1
  n1 <- n-2
  
  S <- (nx1*Sx + ny1*Sy)/n1
  
  D <- diag(S)
  D.inv <- if(length(D)>1) diag(1/D) else 1/D
  
  R <- sqrt(D.inv) %*% S %*% sqrt(D.inv)
  # R <- cor(rbind(x,y))
  p <- ncol(x)
  tr.R2 <- tr(R %*% R)
  d <- 1 + tr.R2/p^(3/2)
  
  nominator <- 1/(1/nx+1/ny) * delta %*% D.inv %*% delta - p
  denominator <- sqrt( 2 * d *(tr.R2 - p^2/n1))
  c(nominator/denominator)
  # list(nom=nominator, den=denominator, t=nominator/denominator)
}
## Testing:
# .p <- 1e0
# .nx <- 1e1
# .ny <- 1e1
# .x <- rmvnorm(.nx, rep(0,.p))
# .y <- rmvnorm(.ny, rep(0,.p))
# t2013(.x,.y)












# A Bootstrapping wrapper.
# Sketch:
## For B bootstrap samples:
## Compute statistic
## Average over samples
t_boot <- function(FUN, noise, labels, B, type2, ...){
  
  t.boot <- rep(NA, B) # initialize output container
  n.samples <- nrow(noise)
  
  for(b in 1:B){
    train.ind <- sample(1:n.samples, replace = TRUE)
    
    train.noise <- noise[train.ind, ]
    train.labels <- labels[train.ind]
    test.noise <- noise[-train.ind,]
    test.labels <- labels[-train.ind] 
    
    t.boot[b] <- FUN(train.noise, train.labels, test.noise, test.labels, ...)
  }
  error.boot <- mean(t.boot)
  
  if(type2==2) {
    result <- error.boot
  }
  if(type2==1) {
    error.resubstitute <- FUN(noise, labels, noise, labels, ...)
    result <- 0.368 * error.resubstitute + 0.632 * error.boot
  }
  
  return(result) 
}



t_svm_boot <- function(noise, labels, B, cost, type2, type){
  t_boot(FUN = t_svm, noise = noise, labels = labels, B=B, type2 = type2, cost, type)
}

t_lda_boot <- function(noise, labels, B, type2, type){
  t_boot(FUN = t_lda, noise = noise, labels = labels, B = B, type2, type)
}


t_sdlda_boot <- function(noise, labels, B, type2, type){
  t_boot(FUN = t_sdlda, noise = noise, labels = labels, B = B, type2=type2, type=type)
}




t_svm_highdim_boot <- function(noise, labels, B, cost, type2, type){
  t_boot(FUN = t_svml2, noise = noise, labels = labels, B=B, type2 = type2, cost=cost, type=type)
}


t_kmmd <- function(x,y,...){
  kmmd(x,y,...)@mmdstats[[1]]
}

t_dcov <- function(x,y){
  unname(eqdist.e(rbind(x,y), sizes = c(nrow(x), nrow(y))))
}
## Testing:
# t_dcov(rmvnorm(n = 20,rep(10,10)),rmvnorm(n = 20,rep(0,10)))



t_Simes <- function(x,y){
  stopifnot(ncol(y)==ncol(x))
  p <- ncol(x)
  group <- as.factor(c(rep('x', nrow(x)), rep('y',nrow(y))))
  my.t <- function(z) t.test(z~group)$p.value
  p.vals <- apply(rbind(x,y), 2, my.t) # Compute variable-wise pvalues
  p.Simes <- p * min(sort(p.vals)/seq_along(p.vals)) # Compute the Simes statistic
  return(as.numeric(-log(p.Simes)))
}
## Testing
# library(mvtnorm)
# t_simes(rmvnorm(n = 20,rep(0,10)),rmvnorm(n = 30,rep(0,10)))


t_CLX <- function(x,y,alpha=0.05,input='dummy'){
  CLX(X = t(x),Y = t(y),alpha = alpha,DNAME=input)$statistics %>% as.numeric()
}
## Testing:
# t_CLX(rmvnorm(n = 20,rep(0,10)),rmvnorm(n = 30,rep(10,10)))






#' Title Compute all statistics
#'
#' @param x1 Matrix of first group observations
#' @param x2 Matrix of second group observations
#' @param Sigma.inv Oracle precision matrix
#' @param noise All data
#' @param labels Group assignments
#' @param fold.ids 
#' @param shift.vec True shift. For NP Oracle Only.
#'
#' @return
#' @export
#'
#' @examples
statistics <- function(x1,x2,Sigma.inv,noise,labels,fold.ids,with.cv=FALSE,shift.vec){
  result <- c(
    ### Two Group Tests:
    Oracle=t_Oracle(x1, x2, Sigma.inv),
    # Oracle.Cov.Mu=t_Oracle_NP(x1, x2, mu.x=0, mu.y=shift.vec, Sigma),
    Hotelling=t_Hotelling(x1, x2, FALSE),
    Schafer=t_Hotelling(x1, x2, TRUE),
    Goeman=t_goeman(x1, x2),
    Srivastava=t_SD(x1, x2),
    Gretton=t_kmmd(x1, x2),
    dCOV=t_dcov(x1,x2),
    Simes=t_Simes(x1,x2),
    Cai=t_CLX(x1,x2),
    ### Accuracy Tests:
    lda.CV.1=t_lda_cv(noise, labels, labels, fold.ids, type=1),
    lda.noCV.1=t_lda(noise, labels, noise, labels, type=1),
    svm.CV.c100=t_svm_cv(noise, labels, labels, fold.ids, cost=100, type=1),
    svm.CV.c001=t_svm_cv(noise, labels, labels, fold.ids, cost=0.01, type=1),
    svm.noCV.c100=t_svm(noise, labels, noise, labels, cost=100, type=1),
    svm.noCV.c001=t_svm(noise, labels, noise, labels, cost=0.01, type=1)
  )
  
  
  if(with.cv) {
    result <- c(result, 
                svm.CV.cCV=t_svm_cvCV(noise, labels, labels, fold.ids)
                )
  }
  return(result)
}
## Testing:



#' Add a cross-validated accuracy test
statisticsWithCV <- function(x1,x2,Sigma,noise,labels,fold.ids,cost.1,cost.2){
  c(
    statistics(x1,x2,Sigma,noise,labels,fold.ids,cost.1,cost.2),
    
  )
}



statisticsBootstrap <- function(x1,x2,Sigma.inv,noise,labels,fold.ids){
  c(
    statistics(x1,x2,Sigma.inv,noise,labels,fold.ids),
    svm.Boot.c100.b10=t_svm_boot(noise, labels, B=10, type2=2, cost=100, type=1),
    svm.Boot.c001.b10=t_svm_boot(noise, labels, B=10, type2=2, cost=0.01, type=1),
    svm.Boot.c100.b50=t_svm_boot(noise, labels, B=50, type2=2, cost=100, type=1),
    svm.Boot.c001.b50=t_svm_boot(noise, labels, B=50, type2=2, cost=0.01,  type=1),
    lda.Boot.b10=t_lda_boot(noise, labels, B=10, type2=2, type=1)
  )
}


statisticsHighDim <- function(x1,x2,Sigma.inv,noise,labels,fold.ids){
  c(
    statistics(x1,x2,Sigma.inv,noise,labels,fold.ids),
    svm.Boot.c100.b50=t_svm_boot(noise, labels, B=50, type2=2, cost=10, type=1),
    svm.Boot.c001.b50=t_svm_boot(noise, labels, B=50, type2=2, cost=0.01,  type=1),
    lda.highdim.Dudoit.CV=t_dlda_cv(noise, labels, labels, fold.ids, type=1),
    lda.highdim.Ramey.CV=t_hdrda_cv(noise, labels, labels, fold.ids, type=1),
    lda.highdim.Pang.CV=t_sdlda_cv(noise, labels, labels, fold.ids, type=1),
    lda.highdim.Pang.b50=t_sdlda_boot(noise, labels, B=50, type2=2, type=1)
  )
}


statisticsLargeSample <- function(x1,x2,Sigma.inv,noise,labels,fold.ids){
  result <- c(
    ### Two Group Tests:
    Oracle.Cov=t_Oracle(x1, x2, Sigma.inv),
    Hotelling=t_Hotelling(x1, x2, FALSE),
    Schafer=t_Hotelling(x1, x2, TRUE),
    Goeman=t_goeman(x1, x2),
    Srivastava=t_SD(x1, x2),
    # Gretton=t_kmmd(x1, x2),
    dCOV=t_dcov(x1,x2),
    Simes=t_Simes(x1,x2),
    Cai=t_CLX(x1,x2),
    ### Accuracy Tests:
    lda.CV.1=t_lda_cv(noise, labels, labels, fold.ids, type=1),
    lda.noCV.1=t_lda(noise, labels, noise, labels, type=1),
    svm.CV.c100=t_svm_cv(noise, labels, labels, fold.ids, cost=100, type=1),
    svm.CV.c001=t_svm_cv(noise, labels, labels, fold.ids, cost=0.01, type=1),
    # svm.CV.cCV=t_svm_cvCV(noise, labels, labels, fold.ids),
    svm.noCV.c100=t_svm(noise, labels, noise, labels, cost=100, type=1),
    svm.noCV.c001=t_svm(noise, labels, noise, labels, cost=0.01, type=1)
  )
  return(result)
}
## Testing:


statisticsSparse <- function(){
  c(
    statistics(x1,x2,Sigma,noise,labels,fold.ids,cost.1,cost.2)
    # Add CART and other adaptive tests. 
  )
}

statisticsGausSVM <- function(x1,x2,Sigma,noise,labels,fold.ids,cost.1,cost.2){
  result <- c(
    statistics(x1,x2,Sigma,noise,labels,fold.ids,cost.1,cost.2),
    svmGAUS.CV.c100=t_svm_cv_Gaus(noise, labels, labels, fold.ids, cost=100, type=1),
    svmGAUS.CV.c001=t_svm_cv_Gaus(noise, labels, labels, fold.ids, cost=100, type=1)
  )
  return(result)
}

statisticsAugmented <- function(x1,x2,Sigma.inv,noise,labels,fold.ids){
  result <- c(
    ### Two Group Tests:
    Oracle.Cov=t_Oracle(x1, x2, Sigma.inv),
    Schafer=t_Hotelling(x1, x2, TRUE),
    Goeman=t_goeman(x1, x2),
    Srivastava=t_SD(x1, x2),
    Gretton=t_kmmd(x1, x2),
    dCOV=t_dcov(x1,x2),
    Simes=t_Simes(x1,x2),
    Cai=t_CLX(x1,x2),
    ### Accuracy Tests:
    svm.CV.c100=t_svm_cv(noise, labels, labels, fold.ids, cost=100, type=1),
    svm.CV.c001=t_svm_cv(noise, labels, labels, fold.ids, cost=0.01, type=1),
    # svm.CV.cCV=t_svm_cvCV(noise, labels, labels, fold.ids),
    svm.noCV.c100=t_svm(noise, labels, noise, labels, cost=100, type=1),
    svm.noCV.c001=t_svm(noise, labels, noise, labels, cost=0.01, type=1)
  )
  return(unlist(result))
}
## Testing
# fff <- statisticsAugmented(x1 = x1,
#                     x2 = x2,
#                     Sigma.inv = Sigma.aug.inv,
#                     noise = noise.augment,
#                     labels = labels,
#                     fold.ids = fold.ids)
