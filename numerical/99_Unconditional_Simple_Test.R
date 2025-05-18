library(MASS)
require(mixtools)
require(mclust)
library(ggplot2)
library(grid)
library(gridExtra)

# Computing the test statistic

#  Ri.compute <- function(datax = train.x, model1, model2, truncated = F){
#  class1pred = predict(model1, newdata = datax, what = "dens")
#  class2pred = predict(model2, newdata = datax, what = "dens")
#  Ri = log(class2pred) - log(class1pred)
#  return(Ri)
#}

fitting.mixture.model <- function(data = train.x, G = 1, modelNames, 
                                  starting.means = matrix(0, nrow = G, ncol = ncol(data)),
                                  starting.variance = diag(rep(1, ncol(data)))){
  iter = 0
  
  log.lik.old = -10^{10}
  density = numeric()
  for(i in 1:G){
    density = cbind(density, apply(data, 1, logdensitynormal, mu = t(starting.means[i,]),
                                   Sigma = starting.variance, diagonal = T))
  }
  log.lik.new = sum(apply(density, 1, max))
  
  while(iter < 100 && abs(log.lik.new - log.lik.old) > 10^(-2)){
    iter = iter + 1
    log.lik.old = log.lik.new
    which.class = apply(density, 1, which.max)
    
    variance.class = matrix(0, nrow = G, ncol = ncol(data))
    starting.means.new =  matrix(0, nrow = G, ncol = ncol(data))
    prop = rep(0, G)
    for(i in 1:G){
      prop[i] =  length(which(which.class == i))/nrow(data)
      print(prop[i])
      starting.means.new[i,] = colMeans(data[which(which.class == i),])
      variance.class[i,] = prop[i] * diag(var(data[which(which.class == i),]))
    }
    starting.variance.new = diag(colSums(variance.class))
    density = numeric()
    for(i in 1:G){
      density = cbind(density, apply(data, 1, logdensitynormal, mu = t(starting.means.new[i,]),
                                     Sigma = starting.variance.new, diagonal = T))
    }
    log.lik.new = sum(apply(density, 1, max))
  }
  if(iter > 100){print("Algorithm didn't converge!")}
  sigma = array(0,c(ncol(data),ncol(data),G))
  for(i in 1:G){
    sigma[,,i] = starting.variance.new
  }
  variance = list(sigma = sigma)
  parameters = list(mean = t(starting.means.new), variance = variance, pro = prop)
  return(list(parameters = parameters))
}


logdensitynormal <- function(x, mu, Sigma, diagonal = F){
  d = nrow(Sigma)
  if(diagonal == T){
    return(-0.5*d*log(2*pi) - 0.5*sum(log(diag(Sigma)))
           -0.5*(x - mu)%*% diag(1/diag(Sigma))%*% t(x - mu))
  }else{
    return(-0.5*d*log(2*pi) - 0.5*log(det(Sigma))
           -0.5*(x - mu)%*% solve(Sigma)%*% t(x - mu))
  }
  
}

logxpy <- function(lx,ly) max(lx,ly) + log1p(exp(-abs(lx-ly)))

Ri.compute <- function(datax = train.x, model1, model2, truncated = F, diagonal = F){
  
  class1pred2 = sapply(c(1:nrow(datax)), function(i){
    return( logdensitynormal(datax[i,],t(model1$parameters$mean), 
                             model1$parameters$variance$sigma[,,1], diagonal))
  })

  class2pred12 = sapply(c(1:nrow(datax)), function(i){
    return( logdensitynormal(datax[i,],t(model2$parameters$mean[,1]), 
                             model2$parameters$variance$sigma[,,1], diagonal))
  })
  class2pred12 = class2pred12 + log(model2$parameters$pro[1])
  
  class2pred22 = sapply(c(1:nrow(datax)), function(i){
    return( logdensitynormal(datax[i,],t(model2$parameters$mean[,2]), 
                             model2$parameters$variance$sigma[,,2], diagonal))
  })
  class2pred22 = class2pred22 + log(model2$parameters$pro[2])
  
  class2pred2 = sapply(c(1:length(class2pred22)), function(i){
    return(logxpy(class2pred12[i], class2pred22[i]))})
  Ri = class2pred2 - class1pred2
  return(Ri)
}

Ri.compute.mixtools <- function(datax = train.x, model1, model2, truncated = F, diagonal = F){
  
  class1pred2 = sapply(c(1:nrow(datax)), function(i){
    return( logdensitynormal(datax[i,],t(model1$parameters$mean), 
                             model1$parameters$variance$sigma[,,1], diagonal))
  })
  
  class2pred12 = sapply(c(1:nrow(datax)), function(i){
    return( logdensitynormal(datax[i,],t(model2$parameters$mean[,1]), 
                             model2$parameters$variance$sigma[,,1], diagonal))
  })
  class2pred12 = class2pred12 + log(model2$parameters$pro[1])
  
  class2pred22 = sapply(c(1:nrow(datax)), function(i){
    return( logdensitynormal(datax[i,],t(model2$parameters$mean[,2]), 
                             model2$parameters$variance$sigma[,,2], diagonal))
  })
  class2pred22 = class2pred22 + log(model2$parameters$pro[2])
  
  class2pred2 = sapply(c(1:length(class2pred22)), function(i){
    return(logxpy(class2pred12[i], class2pred22[i]))})
  Ri = class2pred2 - class1pred2
  return(Ri)
}


# Test 1 from a simple test

# Using mean

KL.Test1.mixture <- function(Ri, log.truncation.ratio = 0, alpha = 0.05){
  gamma.hat = mean(Ri) - log.truncation.ratio
  tau.hat = sd(Ri)
  reject = (gamma.hat > (-qnorm(alpha)*tau.hat)/sqrt(length(Ri)))
  test.stat = (gamma.hat*sqrt(length(Ri)))/tau.hat
  pvalue = 1 - pnorm(test.stat)
  return(list(reject = reject, pvalue = pvalue))
}

# Using median

KL.Test2.mixture.median <- function(Ri, log.truncation.ratio = 0, alpha = 0.05){
  Ri = Ri - log.truncation.ratio
  Tr = length(which(Ri > 0))
  pvalue = 1 - pbinom(Tr, length(Ri), 0.5)
  reject = (pvalue < alpha)
  return(list(reject = reject, pvalue = pvalue))
}

Simple.Test1 <- function(train.x, test.x, alpha = 0.05, use.median = FALSE,
                         log.truncation.ratio = 0){
  
  model1 = densityMclust(data = train.x, G = 1, model = "XXX")
  model2 = densityMclust(data = train.x, G = 2, model = "EEE")
  Ri = Ri.compute(test.x, model1, model2)
  if(use.median == F){
    return(KL.Test1.mixture(Ri, log.truncation.ratio))
  }else{
    return(KL.Test2.mixture.median(Ri, log.truncation.ratio))
  }
  
  #return(KL.Test2.mixture.exact)
}

# Test 2 from a simple test

gamma.null.compute <- function(Ri){
  ei = 2*rbinom(length(Ri), 1, 0.5) - 1
  gamma.null = mean(ei*Ri)
  return(gamma.null)
}

KL.Test2.mixture.exact <- function(Ri, alpha = 0.05){
  simul.gamma = replicate(1000,gamma.null.compute(Ri))
  gamma.hat = mean(Ri)
  pvalue = mean(gamma.hat <= simul.gamma)
  reject = (pvalue <= alpha)
  return(list(reject = reject, pvalue = pvalue))
}