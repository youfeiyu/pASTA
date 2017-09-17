## pASTA code for p-value estimation and subset identification

if(!require("smoothmest")) install.packages("smoothmest")
library(smoothmest)

cond.prob <- function(K, current.subset, study.size, cor){
  # K is the total number of traits
  # current.subset is a vector of length K, a row from all.combn
  # cor is the correlation matrix of each Z_k
  # return a conditional probability density function, which is a function of test statistic z
  
  current.size <- current.subset*study.size
  current.wt <- sqrt(current.size/sum(current.size)) # weights of studies in current subset
  
  cond.prob2 <- function(z)
  { # function: cond.prob
    tmp.prob <- 0
    
    for(k in 1:K)
    { # for1: go over all studies not included in current.subset
      if(k %in% which(current.subset==0))
      { # if1: if study k not included in current.subset
        e <- numeric(K)
        e[k] <- 1
        A <- rbind(e, current.wt)
        Sigma <- A%*%cor%*%t(A)
        cond.mean <- Sigma[1,2]/Sigma[2,2]*z
        cond.var <- Sigma[1,1] - Sigma[1,2]/Sigma[2,2]*Sigma[2,1]
        
        a <- sum(current.size)/(sum(current.size) + study.size[k])
        # conditional probability given Z_gamma=z
        tmp.prob <- tmp.prob + log(pnorm(z*(1-sqrt(a))/sqrt(1-a), mean=cond.mean, sd=sqrt(cond.var)))
      } # end if1
      else 
      { # else1: study k included in current.subset
        e <- numeric(K)
        e[k] <- 1
        A <- rbind(e, current.wt)
        Sigma <- A%*%cor%*%t(A)
        cond.mean <- Sigma[1,2]/Sigma[2,2]*z
        cond.var <- Sigma[1,1] - Sigma[1,2]/Sigma[2,2]*Sigma[2,1]
        b <- sum(current.size)/(sum(current.size) - study.size[k])
        if(sum(current.subset)>1)
        {
          tmp.prob <- tmp.prob + log(pnorm((z*(sqrt(b)-1))/sqrt(b-1), mean=cond.mean, sd=sqrt(cond.var), lower.tail=F))
        }
      } # end else1
    } # end for1
    return(exp(tmp.prob))
  }
  return(cond.prob2)
}

cond.prob.z <- function(K, current.subset, study.size, cor){ 
  # conditional probabilty multiply by standard normal distribution
  cond.prob.tmp = cond.prob(K, current.subset, study.size, cor)
  integrand.final <- function(z){
    current.size = current.subset*study.size
    current.wt = sqrt(current.size/sum(current.size)) # weights of studies in current subset
    Zvar=t(current.wt) %*% cor %*% current.wt
    cond.prob.tmp(z)*dnorm(z, sd=sqrt(Zvar))
  }
  return(integrand.final) # return a function of z
}

p.dlm <- function(K, test.stat, study.size, cor) {
  all.combn <- expand.grid(rep(list(0:1), K))
  all.combn <- all.combn[apply(all.combn, 1, sum) != 0,]
  p.val <- sum(apply(all.combn, 1, function(v){integrate(cond.prob.z(K=K, current.subset = v, study.size=study.size, cor=cor), lower = test.stat, upper = Inf)$value}))
  return(p.val)
}

test.statistic <- function(p.values, study.size)
{
  kk <- length(p.values) # total number of studies
  Z.stats <- (-1)*qnorm(p.values)
  all.combn <- expand.grid(rep(list(0:1), kk))
  all.combn <- all.combn[apply(all.combn, 1, sum) != 0,]
  Z.df <- cbind(all.combn, rep(NA, nrow(all.combn))) # data frame with the last column being the Z.meta of subset S
  colnames(Z.df) <- c(paste("Study", 1:kk, sep=""), "Z.S")
  rownames(Z.df) <- c(1:nrow(all.combn))
  
  for(r in 1:nrow(all.combn))
  { # for4: calculate Z(S) over all possible non-empty subsets
    current.size <- study.size*all.combn[r,] # sample size of studies in a subset
    current.wt <- sqrt(current.size/sum(current.size))
    Z.meta <- sum(current.wt*Z.stats)
    Z.df[r,(kk+1)] <- Z.meta
  } # end for4
  selected.subset=all.combn[which.max(Z.df$Z.S),]
  return(list(test.stat=max(Z.df$Z.S), Z.df=Z.df, selected.subset=selected.subset))
}
