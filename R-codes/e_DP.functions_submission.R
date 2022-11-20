e_fn.consistency.check <- function(parm)
{err <- 0
 
 
 if (max(sort(unique(parm$clust$s_e.v))) != parm$clust$K_e)
 {err <- 1
 }
 
 if (sum(parm$clust$n_e.v) != parm$n2)
 {err <- 2
 }
 
 if (length(parm$clust$n_e.v) != parm$clust$K_e) 
 {err <- 3
 }
 
 if (length(parm$clust$phi_e.v) != parm$clust$K_e) 
 {err <- 4
 }
 
 if (sum(is.na(parm$clust$phi_e.v)) > 0)
 {err <- 5
 }
 
 err
}


e_fn.swap.clusters <- function(parm, g1, g2)
{
  
  ####################################################
  # swap the cluster labels g1 and g2
  ####################################################
  
  ind1 <- parm$clust$s_e.v == g1
  ind2 <- parm$clust$s_e.v == g2
  parm$clust$s_e.v[ind1] <- g2
  parm$clust$s_e.v[ind2] <- g1
  
  buffer <- parm$clust$n_e.v[g1]
  parm$clust$n_e.v[g1] <- parm$clust$n_e.v[g2]
  parm$clust$n_e.v[g2] <- buffer
  
  buffer <- parm$clust$phi_e.v[g1]
  parm$clust$phi_e.v[g1] <- parm$clust$phi_e.v[g2]
  parm$clust$phi_e.v[g2] <- buffer
  
  
  parm
}

e_fn.clip.clusters <- function(parm, keep)
{
  parm$clust$n_e.v <- parm$clust$n_e.v[keep]
  
  parm$clust$phi_e.v <- parm$clust$phi_e.v[keep]
  
  
  parm
}



###########################################################


e_fn.gibbs <- function(k, parm)
{	
  old.s <- parm$clust$s_e.v[k]
  
  ###############
  
  parm$clust$n_e.v[old.s] <- parm$clust$n_e.v[old.s] - 1
  
  ### Has the current cluster just been emptied?
  empty.flag <- parm$clust$n_e.v[old.s]==0
  
  # actual number of clusters after removing kth case
  new.K <- parm$clust$K_e - as.numeric(empty.flag)
  
  #######################################################
  ### Likelihood for clusters (including possibly empty one)
  #######################################################
  
  Y.k.e<-mean(parm$Y.e[k,])
  
  L.v <- array(,parm$clust$K_e)
  
  for (je in 1:parm$clust$K_e)  ####caution!!!! different!!!!####9/24/2015
  {L.v[je] <- (dnorm(Y.k.e, mean=parm$clust$phi_e.v[je], sd=parm$tau/sqrt(parm$p), log=TRUE))
  }		
  
  # marginal likelihood of new cluster
  new.L <- (dnorm(Y.k.e, mean=parm$clust$mu_e, sd=sqrt(parm$tau^2/parm$p + parm$clust$tau_e^2), log=TRUE))
  
  L.v <- c(L.v, new.L)
  
  #######################################################
  ### Prior for clusters by Polya urn scheme
  ### Emptied cluster is "removed" by setting prob to -Inf
  #######################################################
  
  log.prior.v <- array(NA, (1+parm$clust$K_e))
  log.prior.v <- log(c(parm$clust$n_e.v, parm$clust$Me))
  
  #######################################################
  
  
  tmp2 <- log.prior.v + L.v 
  maxx <- max(tmp2)
  tmp2 <- tmp2 - maxx
  
  tmp2 <- exp(tmp2)
  
  parm$clust$post.k <- tmp2
  parm$clust$post.k <- parm$clust$post.k / sum(parm$clust$post.k)
  
  # sample s for kth case
  new.s <- sample(1:(parm$clust$K_e+1), size=1, replace=TRUE, prob=parm$clust$post.k)
  parm$clust$s_e.v[k] <- new.s
  new.flag <- new.s == (parm$clust$K_e+1)
  
  #######################
  
  if (!new.flag)
  {parm$clust$n_e.v[new.s] <- parm$clust$n_e.v[new.s] + 1			 
  }
  
  if (new.flag) 
  {parm$clust$K_e <- parm$clust$K_e + 1
   parm$clust$n_e.v <- c(parm$clust$n_e.v, 1)
   
   post.prec <- 1/parm$clust$tau_e^2 + parm$p*(1/parm$tau^2)
   post.sd <- sqrt(1/post.prec)
   post.mean <- post.sd^2 * (Y.k.e*parm$p/parm$tau^2 + parm$clust$mu_e/parm$clust$tau_e^2)
   new.phi <- rnorm(n=1, post.mean, post.sd)
   
   parm$clust$phi_e.v <- c(parm$clust$phi_e.v, new.phi)
   if (is.na(new.phi))
   {stop("Error in new.phi")
   }
  }
  
  parm
}

#########################################################

e_fn.drop <- function(parm)
{
  
  ##########################################
  ## Drop empty clusters:
  ## (i) Set parm$clust$K equal to number of non-empty clusters 
  ## (ii)  Move empty clusters to end by relabeling clusters
  ## (iii) Retain only clusters  1,...,parm$clust$K
  #########################################
  
  parm$clust$K_e <- sum(parm$clust$n_e.v>0)
  num.dropped <- sum(parm$clust$n_e.v==0)
  
  if (num.dropped > 0)
  {
    for (rr in 1:num.dropped)
    {
      old.label <-  min(which(parm$clust$n_e.v==0))
      new.label <- max(which(parm$clust$n_e.v>0))
      stopp <-  max(which(parm$clust$n_e.v>0)) == parm$clust$K_e
      if (stopp)
      {break
      }
      parm <- e_fn.swap.clusters(parm, g1=new.label, g2=old.label)
    }
  }
  
  keep <- 1:parm$clust$K_e
  parm <- e_fn.clip.clusters(parm, keep)
  
  parm
}





fn.epsilon.DP <- function(parm)
{		
  parm$Y.e<-parm$X.LANE-parm$theta.mt[parm$group,]
  for (k in 1:parm$n2)
  {parm  <- e_fn.gibbs(k, parm)
  }
  
  ##########################################
  ## Drop empty group clusters and relable clusters
  #########################################
  
  parm <- e_fn.drop(parm)
  
  ############################
  
  for (jje in 1:parm$clust$K_e)
  {indx.j <- parm$clust$s_e.v==jje
  n.j <- sum(indx.j)
  Y.e.j <- rowMeans(parm$Y.e[indx.j,,drop=F])
  
  post.prec <- 1/parm$clust$tau_e^2 + n.j*parm$p*(1/parm$tau^2)
  post.sd <- 1/sqrt(post.prec)
  post.mean <- post.sd^2 * (mean(Y.e.j)*n.j*parm$p/(parm$tau^2) + parm$clust$mu_e/parm$clust$tau_e^2)
  
  parm$clust$phi_e.v[jje] <- rnorm(n=1, post.mean, post.sd)
  }
  
  parm$clust$theta_e.v <- parm$clust$phi_e.v[parm$clust$s_e.v]
  
  parm$epsilon<-parm$clust$theta_e.v
  
  err <- e_fn.consistency.check(parm)
  if (err > 0)
  {stop(paste("Failed consistency check in epsilon: err=",err))
  }
  
  parm
}


