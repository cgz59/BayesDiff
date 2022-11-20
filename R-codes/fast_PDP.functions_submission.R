
PDP_fn.compact.consistency.check <- function(parm)
	{err <- 0

	if (max(sort(unique(parm$clust$c.v))) != parm$clust$G)
		{err <- 1
		}

	if (sum(colSums(parm$clust$C.m.mt)==0) > 0)
		{err <- 1.5
	}
	
	if (sum((parm$clust$C.m.vec)==0) > 0)
	{err <- 1.5
	}

	err <- PDP_fn.consistency.check(parm)


	err
	}

PDP_fn.check.nbhd <- function(parm)
{
  max_dist.v <- array(,length(parm$clust$col.nbhd.k))

  for (zz in 1:length(parm$clust$col.nbhd.k))
  {
    v1 <- parm$clust$col_post.prob.mt[,parm$clust$col.nbhd.k[zz]]
    m1 <- matrix(parm$clust$col_post.prob.mt[,parm$clust$col.nbhd[[zz]]],ncol=length(parm$clust$col.nbhd[[zz]]))
    max_dist.v[zz] <- max(2*(1-colSums(sqrt(v1*m1))))
  }

  parm$clust$nbhd_max_dist <- max(max_dist.v)

  parm

}

PDP_fn.consistency.check <- function(parm)
	{err <- 0


	if (sum(parm$clust$C.m.mt) != (parm$p))
		{err <- 4
	}
	
	if (sum(parm$clust$C.m.vec) != (parm$p))
	{err <- 4.4
	}
  
# 	if ((sum(parm$clust$n.vec)) != parm$tot.trt*parm$clust$G)
# 	{err <- 8
# 	}

	if (length(parm$clust$n.vec) != parm$clust$K)
		{err <- 9
		}

	if (length(parm$clust$phi.v) != parm$clust$K)
		{err <- 10
		}


	err
	}


PDP_fn.swap.clusters <- function(parm, g1, g2, computeMode, mt=T)
	{

	####################################################
	# swap the group labels g1 and g2
	####################################################

	ind1 <- parm$clust$c.v == g1
	ind2 <- parm$clust$c.v == g2
	parm$clust$c.v[ind1] <- g2
	parm$clust$c.v[ind2] <- g1

#	if (computeMode$computeR) {
	  buffer <- parm$clust$s.mt[,g1]
	  parm$clust$s.mt[,g1] <- parm$clust$s.mt[,g2] # HOT 3
	  parm$clust$s.mt[,g2] <- buffer
# 	} else {
# 	  .swapIntegerMatrix(parm$clust$s.mt, g1, g2, FALSE)
# 	}

if (mt){
  buffer <- parm$clust$C.m.mt[,g1]
  parm$clust$C.m.mt[,g1] <- parm$clust$C.m.mt[,g2]
  parm$clust$C.m.mt[,g2] <- buffer
}

  buffer <- parm$clust$C.m.vec[g1]
  parm$clust$C.m.vec[g1] <- parm$clust$C.m.vec[g2]
  parm$clust$C.m.vec[g2] <- buffer

  
  buffer <- parm$clust$gr.c.v[g1]
  parm$clust$gr.c.v[g1]<-parm$clust$gr.c.v[g2]
  parm$clust$gr.c.v[g2]<-buffer

	#####################

	buffer <- parm$clust$A.mt[,g1]
	parm$clust$A.mt[,g1] <- parm$clust$A.mt[,g2]
	parm$clust$A.mt[,g2] <- buffer

	parm
	}

PDP_fn.clip.clusters <- function(parm, keep, mt=T)
	{

	parm$clust$s.mt <- parm$clust$s.mt[,keep, drop = F]
  if (mt) {
    parm$clust$C.m.mt <- parm$clust$C.m.mt[,keep, drop = F]
  }
	parm$clust$C.m.vec <- parm$clust$C.m.vec[keep]
	parm$clust$gr.c.v <- parm$clust$gr.c.v[keep]
  #    parm$clust$small.indx <- parm$clust$small.indx[keep]
  #    parm$clust$order.v <- parm$clust$order.v[keep]

  parm$clust$A.mt <- parm$clust$A.mt[,keep, drop = F]

	parm
	}



###########################################################


PDP_fn.log.lik <- function(gg, x.mt, parm, colSums=TRUE)
	{
	#if (gg > 0)
		#{
      a2.v <- parm$clust$A.mt[,gg]
		#z.g.v <- parm$clust$s.mt[,gg] > 0
		log.lik.v <- rep(0,ncol(x.mt)) ## HOT

		#if (sum(z.g.v) > 0)
			{#a2.1.v <- parm$clust$A.mt[,gg]
			small.X.1 <- matrix(x.mt, ncol=ncol(x.mt)) # HOT

			if (colSums)
			{#if (!parm$flip.sign)
				  {tmp <- colSums(-.5*(parm$ni/parm$tau^2)*(small.X.1 - a2.v)^2)
				  }
# 			  if (parm$flip.sign)
# 			    {tmp <- colSums(-.5*(small.X.1 - a2.1.v)^2) + colSums(-.5*(-small.X.1 - a2.1.v)^2)
# 			  }
			} # end 	if (colSums)

			if (!colSums)
				{#if (!parm$flip.sign)
				  {tmp <- sum(-.5*(parm$ni/parm$tau^2)*(small.X.1 - a2.v)^2)
				  }
# 			  if (parm$flip.sign)
# 			    {tmp <- sum(-.5*(small.X.1 - a2.1.v)^2) + sum(-.5*(-small.X.1 - a2.1.v)^2)
# 			    }
				} # end if (!colSums)

			log.lik.v <- log.lik.v + tmp - parm$tot.trt*(.5*log(2*pi))-sum(log(parm$tau/sqrt(parm$ni)))
	    } # end if (sum(z.g.v) > 0)

# 		if (sum(1-z.g.v) > 0)
# 			{small.X.0 <- matrix(x.mt[!z.g.v,], ncol=ncol(x.mt))
# 			if (colSums)
# 				{if (!parm$flip.sign)
# 				  {tmp <- colSums(-.5*small.X.0^2)
# 				  }
# 			  if (parm$flip.sign)
# 			    {tmp <- 2*colSums(-.5*small.X.0^2)
# 			    }
# 				}
# 			if (!colSums)
# 				{if (!parm$flip.sign)
# 				  {tmp <- sum(-.5*small.X.0^2)
# 				  }
# 			  if (parm$flip.sign)
# 			    {tmp <- 2*sum(-.5*small.X.0^2)
# 			    }
# 				}
# 			log.lik.v <- log.lik.v + tmp/parm$tau_0^2 + sum(1-z.g.v)*(-.5*log(2*pi)-log(parm$tau_0))
# 			} # end if (sum(1-z.g.v) > 0)
		#} # end if (gg > 0)


	log.lik.v
	}



PDP_fn.nbhd <- function(relative_I, parm, max.col.nbhd.size)
{
  if (length(relative_I)>1)
  {relative_k <- sample(relative_I, size=1)
  }
  if (length(relative_I)==1)
  {relative_k <- relative_I
  }

  post.prob.mt <- parm$clust$col_subset_post.prob.mt

  tmp1.mt <- matrix(post.prob.mt[,relative_I], ncol = length(relative_I)) ## HOT
  tmp2.v <- post.prob.mt[,relative_k]
  tmp3.mt <- sqrt(tmp1.mt * tmp2.v) ## HOT
  H.v <-  2 * (1 - colSums(tmp3.mt))

  cutoff <- parm$col.delta
  flag.v <- which(H.v <= cutoff)
  relative_I.k <- relative_I[flag.v]

  if (length(relative_I.k) > max.col.nbhd.size) {
    relative_I.k <- relative_I[rank(H.v, ties="random") <= max.col.nbhd.size]
  }

  relative_I.k <- sort(relative_I.k)

  relative_I <- sort(setdiff(relative_I, relative_I.k))
  relative_I <- sort(relative_I)

  list(relative_k, relative_I.k, relative_I)

}


PDP_fn.post.prob.and.delta <- function(parm, max.col.nbhd.size, col.frac.probes, computeMode) {

  col.subset <- 1:parm$p

  if (computeMode$computeR) {

    ################################################
    ### Compute pmf of cluster variables w_1,...,w_p
    ###############################################

    prior.prob.v <- c(parm$clust$C.m.vec)
    small <- 1e-3 # compared to 1
    prior.prob.v[prior.prob.v < small] <- small

    subset_log.ss.mt <- array(, c((parm$clust$G), length(col.subset)))

    for (gg in 1:parm$clust$G)
    {subset_log.ss.mt[(gg),] <- PDP_fn.log.lik(gg, x.mt = parm$Z[,col.subset], parm)
    }

    subset_log.ss.mt <- subset_log.ss.mt + log(prior.prob.v)

    maxx.v <- apply(subset_log.ss.mt, 2, max)

    subset_log.ss.mt <- t(t(subset_log.ss.mt) - maxx.v)
    subset_ss.mt <- exp(subset_log.ss.mt)

    col.sums.v <- colSums(subset_ss.mt)
    subset_ss.mt <- t(t(subset_ss.mt)/col.sums.v)

    # replace zeros by "small"
    small2 <- 1e-5
    subset_ss.mt[subset_ss.mt < small2] <- small2

    # again normalize
    col.sums.v <- colSums(subset_ss.mt)
    subset_ss.mt <- t(t(subset_ss.mt)/col.sums.v)

    parm$clust$col_post.prob.mt <- array(,c((parm$clust$G), parm$p))
    parm$clust$col_post.prob.mt[,col.subset] <- subset_ss.mt

    dimnames(parm$clust$col_post.prob.mt) <- list(1:parm$clust$G, 1:parm$p)

    parm$clust$col_subset_post.prob.mt <- subset_ss.mt
    dimnames(parm$clust$col_subset_post.prob.mt) <- list(1:parm$clust$G, 1:length(col.subset))

    #########################################
    ### now compute the delta-neighborhoods
    #########################################

    if (computeMode$computeC) { # debugging
      savedSeed <- .GlobalEnv$.Random.seed # For debugging purposed only
    }

    parm$clust$col.nbhd <- NULL
    parm$clust$col.nbhd.k <- NULL
    relative_I <- 1:length(col.subset)

    while (length(relative_I) >= 1) {
      tmp <- PDP_fn.nbhd(relative_I, parm, max.col.nbhd.size) # HOT inside function
      relative_k <- tmp[[1]]
      relative_I.k <- tmp[[2]]
      relative_I <- tmp[[3]]
      #
      parm$clust$col.nbhd <- c(parm$clust$col.nbhd, list(col.subset[relative_I.k]))
      parm$clust$col.nbhd.k <- c(parm$clust$col.nbhd.k, col.subset[relative_k])
    }

    ### SG: how good are these nbhds?

    parm <- PDP_fn.check.nbhd(parm)

  }

  if (computeMode$computeC) {

    if (computeMode$computeR) { # debugging
      .GlobalEnv$.Random.seed <- savedSeed # Roll back PRNG
    }

    test <- .computeColumnPmfAndNeighborhoods(computeMode$device$engine,
                                        parm$clust$C.m0, parm$clust$C.m.vec, 1e-3, 1e-5,
                                        parm$clust$G, parm$n2,
                                        parm$Y, parm$X,
                                        parm$clust$A.mt, parm$clust$s.mt,
                                        col.subset,
                                        parm$clust$C.m.vec, parm$p,
                                        parm$clust$phi.v, parm$tau, parm$tau_0, parm$tau_int,
                                        max.col.nbhd.size, parm$col.delta, TRUE)

    if (computeMode$computeR) { # debugging
      assertEqual(test$index, parm$clust$col.nbhd.k)
      assertEqual(test$neighbor, unlist(parm$clust$col.nbhd))
      assertEqual(test$neighborhoodMax, parm$clust$nbhd_max_dist, computeMode$tolerance)
    }

    # Convert from simple flat format to list of int vectors
    end <- test$offset[-1] - 1
    begin <- test$offset
    length(begin) <- length(begin) - 1

    parm$clust$col.nbhd <- lapply(1:length(begin),FUN = function(x) {
                                    test$neighbor[begin[x]:end[x]]
                                  })
    parm$clust$col.nbhd.k <- test$index
    parm$clust$nbhd_max_dist <- test$neighborhoodMax

  } # computeMode

  ## END

  parm
}


###########################################################

PDP_fn.gibbs_order_1 <- function(k, parm, data, computeMode)
{	k <- parm$k
err <- PDP_fn.consistency.check(parm)
if (err > 0)
{stop(paste("GIBBS - 0: failed consistency check: err=",err))
}
# store current state
old.parm <- parm

empty.indx <- which(colSums(parm$clust$C.m.mt)==0)
#new.G <- parm$clust$G - length(empty.indx)
#new.G.gr<-rowSums(parm$clust$C.m.mt>0)

if (length(empty.indx) >0)
{
#  if (computeMode$computeR) {
    new.s.mt <- parm$clust$s.mt[,-empty.indx] # HOT 3
    new.gr.c.v <- parm$clust$gr.c.v[-empty.indx]
    new.s.mt_g1<-new.s.mt[,new.gr.c.v==1,drop=F]
    new.s.mt_g2<-new.s.mt[,new.gr.c.v==2,drop=F]
    if (sum(apply(new.s.mt_g1,2,max)-apply(new.s.mt_g1,2,min))>0) {stop(paste('problem in new.s.mt_g1'))}
    new.s.mt_g1_comp<-apply(new.s.mt_g1,2,unique)
    new.n.vec <- array(,parm$clust$K)
    for (pp in 1:parm$clust$K) {
      new.n.vec[pp] <- sum(new.s.mt_g2==pp)+sum(new.s.mt_g1_comp==pp) # HOT
    }
#  }
  
  # if (computeMode$computeC) {
  #   
  #   tab <- .fastTabulateExcludeEmptiedIndices(parm$clust$s.mt, emptied.indx, parm$clust$K, TRUE)
  #   
  #   if (computeMode$computeR) { # debugging
  #     assertEqual(tab[1], new.n0)
  #     assertEqual(tab[-1], new.n.vec)
  #   }
  #   
  #   #     	  new.n0 <- tab[1]
  #   new.n.vec <- tab[-1]
  #   
  # } # computeMode
  
  # emptied.s.indx <- which(new.n.vec == 0)
  # new.K <- parm$clust$K - length(emptied.s.indx)
  
}

if (length(empty.indx) ==0)
{
  new.s.mt <- parm$clust$s.mt
  new.gr.c.v <- parm$clust$gr.c.v
  new.s.mt_g1<-new.s.mt[,new.gr.c.v==1,drop=F]
  new.s.mt_g2<-new.s.mt[,new.gr.c.v==2,drop=F]
  if (sum(apply(new.s.mt_g1,2,max)-apply(new.s.mt_g1,2,min))>0) {stop(paste('problem in new.s.mt_g1'))}
  new.s.mt_g1_comp<-apply(new.s.mt_g1,2,unique)
  
  new.n.vec <- array(,parm$clust$K)
  for (pp in 1:parm$clust$K) {
    new.n.vec[pp] <- sum(new.s.mt_g2==pp)+sum(new.s.mt_g1_comp==pp) # HOT
  }
  # emptied.s.indx <- which(new.n.vec==0)
  # new.K <- parm$clust$K
}

## generate auxilliary P vector

tmp.M <- rep(parm$clust$M/parm$clust$K,parm$clust$K)
tmp.alpha <- tmp.M+new.n.vec
#tmp.alpha[emptied.s.indx] <- 0
P.aux <- rgamma(parm$clust$K,c(tmp.alpha),1)
P.aux <- P.aux/sum(P.aux)

c.k<-parm$clust$c.v[k]
g.k<-parm$clust$g.v[k]

# if (k>1){
   if (parm$clust$C.m.mt[g.k,c.k]<=0) {stop(paste("problem1"))}
  parm$clust$C.m.mt[g.k,c.k] <- parm$clust$C.m.mt[g.k,c.k] - 1
# }
#   if (k==1) {
#     parm$clust$C.m.mt[1,c.k] <- parm$clust$C.m.mt[1,c.k] - 1
#   }

  if (k<parm$p){
    c.k.right<-parm$clust$c.v[k+1]
    g.k.right<-parm$clust$g.v[k+1]
    if (parm$clust$C.m.mt[g.k.right,c.k.right]<=0) {stop(paste("problem2"))}
    parm$clust$C.m.mt[g.k.right,c.k.right]<-parm$clust$C.m.mt[g.k.right,c.k.right] - 1
  }
  
  
	###############
		
	x.mt <- matrix(parm$Z[,k], ncol=1)

	if (computeMode$computeR) {

	# intercept cluster or any existing cluster
	L.v <- sapply(1:parm$clust$G, PDP_fn.log.lik, x.mt, parm) # HOT

  }

	if (computeMode$computeC) {

    test <- .computePdpLogLikelihood(computeMode$device$engine, k, parm$X,
                                     parm$clust$A.mt, parm$clust$s.mt,
                                     parm$clust$G, parm$n2,
                                     parm$tau, parm$tau_0, parm$tau_int, FALSE)

    if (computeMode$computeR) { # debugging
      assertEqual(test$logLikelihood, L.v, computeMode$tolerance)
    }

    L.v <- test$logLikehood
  } # computeMode
  # NB: returned logLikelihood differ from those computed above by approx 1e-15.  I believe this is due to non-transitivity of FLOPs

	#######################################################
	### emptied clusters are gone forever under Gibbs sampling
	#######################################################
	new.G.gr<-rowSums(parm$clust$C.m.mt>0)

  ## marginal likelihood of new cluster

  if (computeMode$computeR) {
    
    #marginal likelihood of new cluster for clust$g.v=2
    marg.log.lik.v_g2 <- array(,length(x.mt))
    
    tmp.marg.v_g2.mt <- matrix(NA,length(P.aux),length(x.mt))
    
    for (tt in 1:length(x.mt))
    {
        tmp.log.lik.v_g2 <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=(parm$tau/sqrt(parm$ni[tt])),log=T) # HOT 3
        tmp.log.marg.v_g2 <- tmp.log.lik.v_g2+log(P.aux)
        marg.log.lik.v_g2[tt] <- log(sum(exp(tmp.log.marg.v_g2)))
        tmppost_g2<-tmp.log.marg.v_g2-max(tmp.log.marg.v_g2)
        tmppost_g2<-exp(tmppost_g2)
        tmp.marg.v_g2<-tmppost_g2/sum(tmppost_g2)
        tmp.marg.v_g2.mt[,tt]<-tmp.marg.v_g2
#       tmp.marg.v<-tmppost
#       tmp.marg.v[tmp.marg.v < 1e-300] <- 1e-300
#       tmp.marg.v[emptied.s.indx] <- 0
#       tmp.marg.v <- tmp.marg.v/sum(tmp.marg.v)
        
        #R2.new.star<-R2.new.star*tmp.marg.v[tmp.cand.s]
    }
    marg.log.lik_g2 <- sum(marg.log.lik.v_g2)
#     while (length(unique(cand.s.v.k_g2))<2) {
#       for (tt in 1:length(x.mt)) {
#       cand.s.v.k_g2[tt]<-sample(1:parm$clust$K, size=1, replace=TRUE, prob=tmp.marg.v_g2.mt[,tt])
#     }
#     }
    
    # cand.gr<-1
    # if (max(parm$clust$phi.v[cand.s.v.k])-min(parm$clust$phi.v[cand.s.v.k])>parm$delta_margin) {cand.gr<-2}
    
    #marginal likelihood of new cluster for clust$g.v=1
    tmp.log.lik.v_g1 <- array(,length(parm$clust$phi.v))
    for (pp in 1:length(parm$clust$phi.v)) {
      tmp.log.lik.v_g1[pp] <- sum(dnorm(x.mt,mean=parm$clust$phi.v[pp], sd=(parm$tau/sqrt(parm$ni)),log=T))
    }
    tmp.log.marg.v_g1 <- tmp.log.lik.v_g1+log(P.aux)
    marg.log.lik_g1 <- log(sum(exp(tmp.log.marg.v_g1)))
    tmppost_g1<-tmp.log.marg.v_g1-max(tmp.log.marg.v_g1)
    tmppost_g1<-exp(tmppost_g1)
    tmp.marg.v_g1<-tmppost_g1/sum(tmppost_g1)
    cand.s.v.k_g1<-sample(1:parm$clust$K, size=1, replace=TRUE, prob=tmp.marg.v_g1)
    cand.s.v.k_g1<-rep(cand.s.v.k_g1,length(x.mt))
    
  }
  

  if (computeMode$computeC) {

    test <- .computeMarginalLikelihood(computeMode$device$engine,
                                       x.mt,
                                       parm$clust$phi.v,
                                       P.aux,
                                       parm$tau, parm$tau_0,
                                       FALSE, # no sampling
                                       computeMode$exactBitStream)

    if (computeMode$computeR) { # debugging
      assertEqual(test$logMarginalLikelihood, marg.log.lik, computeMode$tolerance)
    }

    marg.log.lik <- test$logMarginalLikelihood

  } #end computemode

  L.v <- cbind(c(L.v, marg.log.lik_g1),c(L.v, marg.log.lik_g2))

  parm$marg.log.lik_open<-c(marg.log.lik_g1,marg.log.lik_g2,marg.log.lik_g2-marg.log.lik_g1)
  # allow -Inf's in Gibbs sampling log-prior (just emptied clusters)
  if (k>1){
    h.km1<-parm$clust$gr.v[k-1]#==h_{k-1}
    if (parm$eta>0) {
      tmp.prob.g.k<-parm$rho[h.km1]+(1-parm$rho[h.km1])*exp(-parm$distance[k-1]/parm$eta)
    }
    if (parm$eta==0) {
      tmp.prob.g.k<-parm$rho[h.km1]
    }
    # if (parm$eta==-1) {
    #   tmp.prob.g.k<-1
    # }
    prob.g.k.v<-rep(1-tmp.prob.g.k,2)
    prob.g.k.v[h.km1]<-tmp.prob.g.k
    tmp.prior <- prob.g.k.v*cbind((parm$clust$C.m.mt-c(0,parm$d)), (parm$b1+new.G.gr*c(0,parm$d)))
  }
  if (k==1) {
    tmp.prior <- parm$rho*cbind((parm$clust$C.m.mt-c(0,parm$d)), (parm$b1+new.G.gr*c(0,parm$d)))
    
  }
	tmp.prior[tmp.prior<0]<-0
	log.prior<-log(tmp.prior)

	tmp2 <- t(t(log.prior) + L.v)
	maxx <- max(tmp2)
	tmp2 <- tmp2 - maxx

	tmp2 <- exp(tmp2)

	prop.prob<-tmp2
	prop.prob<-prop.prob/sum(prop.prob)

	########################################################################

	########################################################################

	tmp.prop.c.g.k <- sample(c(1:(parm$clust$G+1),-1:-(parm$clust$G+1)), size=1, replace=TRUE, prob=t(prop.prob))
  prop.c.k <- abs(tmp.prop.c.g.k)
  prop.g.k <- 1+as.numeric(tmp.prop.c.g.k<0)
	new.flag <- prop.c.k == (parm$clust$G+1)
	
	match.s.idx<-which(cand.s.v.k_g1[1]==parm$clust$s.mt[1,parm$clust$gr.c.v==1])
	
	if (length(match.s.idx)>1) {stop(paste('problem in prop.g, same cluster exist before proposing new one'))}
	
	if (new.flag & prop.g.k==1 & length(match.s.idx)>0) {
	  tmp.idx_g1<-which(parm$clust$gr.c.v==1)
	  prop.c.k<-tmp.idx_g1[match.s.idx]
	  new.flag<-FALSE
	}
	
	if (new.flag) {
	  cand.s.v.k_g2 <- rep(-1,length(x.mt))
	  if (prop.g.k==2) {
	    tmp.count<-0
	    while (length(unique(cand.s.v.k_g2))<2 & tmp.count<100) {
	      for (tt in 1:length(x.mt)) {
	        cand.s.v.k_g2[tt]<-sample(1:parm$clust$K, size=1, replace=TRUE, prob=tmp.marg.v_g2.mt[,tt])
	        tmp.count<-tmp.count+1
	      }
	    }
	    if (length(unique(cand.s.v.k_g2))<2) {
	      prop.g.k<-1
	      if (length(match.s.idx)>0) {
	        tmp.idx_g1<-which(parm$clust$gr.c.v==1)
	        prop.c.k<-tmp.idx_g1[match.s.idx]
	        new.flag<-FALSE
	      }
	    }
	  }
	}
	
	#######################

	if (!new.flag)
		{gr.k<-parm$clust$gr.v[k]<-parm$clust$gr.c.v[prop.c.k]
		
#		if (k>1){
		  parm$clust$C.m.mt[prop.g.k,prop.c.k] <- parm$clust$C.m.mt[prop.g.k,prop.c.k] + 1
#		}
# 		if (k==1) {
# 		  parm$clust$C.m.mt[1,prop.c.k] <- parm$clust$C.m.mt[1,prop.c.k] + 1
# 		}
		#parm$clust$C.m.mt[parm$clust$gr.v[k],c.k.right]<-parm$clust$C.m.mt[parm$clust$gr.v[k],c.k.right] + 1
		#q1.star<-prop.prob[prop.c.k]
		}

	if (new.flag)
  {
	  cand.s.v.k<-cbind(cand.s.v.k_g1,cand.s.v.k_g2)
	  #q1.star<-prop.prob[prop.c.k]*(R2.new.star)
	  parm$clust$gr.v[k]<-prop.g.k
	  parm$clust$gr.c.v<-c(parm$clust$gr.c.v,prop.g.k)
    parm$cand$s.v.k <- cand.s.v.k[,prop.g.k]

#     parm$cand$n.vec.k <- array(,parm$clust$K)
#     if (prop.g.k==2) {
#       for (gg in 1:parm$clust$K)
#       {parm$cand$n.vec.k[gg] <- sum(cand.s.v.k[,prop.g.k]==gg)
#       }
#     }
# 
#     if (prop.g.k==1) {
#       parm$cand$n.vec.k<-rep(0,parm$clust$K)
#       parm$cand$n.vec.k[cand.s.v.k[,1]]<-1
#     }

  ##################
   parm$clust$G <- parm$clust$G + 1
   #new.G<-new.G+1
  #parm$clust$C.m.mt[cand.gr,c.k.right]<-parm$clust$C.m.mt[cand.gr,c.k.right] + 1
   tmp.C.m<-c(0,0)
#   if (k>1){
     tmp.C.m[prop.g.k]<-1
#   }
#    if (k==1) {
#      tmp.C.m[1]<-1
#    }
		parm$clust$C.m.mt <- cbind(parm$clust$C.m.mt,tmp.C.m)

#     parm$clust$s.v <- c(parm$clust$s.v, parm$cand$s.v.k)
# 
# 		parm$clust$s.mt <- matrix(parm$clust$s.v, nrow = parm$tot.trt)
		parm$clust$s.mt <- cbind(parm$clust$s.mt,parm$cand$s.v.k)
		

#   parm$clust$n.vec <- parm$clust$n.vec + parm$cand$n.vec.k

   tmp.a.v <- array(,parm$tot.trt)
   s.G.v <- parm$cand$s.v.k
		tmp.a.v <- parm$clust$phi.v[s.G.v]
		parm$clust$A.mt <- cbind(parm$clust$A.mt, tmp.a.v) # HOT

  } # end  if (new.flag)
	
	parm$clust$c.v[k] <- prop.c.k
	parm$clust$g.v[k] <- prop.g.k
	
#####Accounting for transitioning out, computing q2.star
	if (k<parm$p){
	x.mt.right <- matrix(parm$Z[,(k+1)], ncol=1)
	right.G.gr<-rowSums(parm$clust$C.m.mt>0)
	
	#if (computeMode$computeR) {
	  
	  # any existing cluster
	  L.v.right <- sapply(1:parm$clust$G, PDP_fn.log.lik, x.mt.right, parm) # HOT
	  
	#}
	  
	  # tmp.alpha.right <- tmp.M+parm$clust$n.vec
	  # tmp.alpha.right[emptied.s.indx] <- 0
	  # P.aux.right <- rgamma(parm$clust$K,c(tmp.alpha.right),1)
	  # P.aux.right <- P.aux.right/sum(P.aux.right)
	  
	  if (g.k.right==2) {
	  marg.log.lik.v.right_g2 <- array(,length(x.mt.right))
	  tmp.s.v.right<-parm$clust$s.mt[,c.k.right]

	  R2.new.star.right<-0
	  for (tt in 1:length(x.mt.right))
	  {
	    tmp.log.lik.v.right_g2 <- dnorm(x.mt.right[tt],mean=parm$clust$phi.v, sd=(parm$tau/sqrt(parm$ni[tt])),log=T) # HOT 3
	    tmp.log.marg.v.right_g2 <- tmp.log.lik.v.right_g2+log(P.aux)
	    marg.log.lik.v.right_g2[tt] <- log(sum(exp(tmp.log.marg.v.right_g2)))
	    tmppost.right_g2<-tmp.log.marg.v.right_g2-max(tmp.log.marg.v.right_g2)
	    tmppost.right_g2<-exp(tmppost.right_g2)
	    tmp.marg.v.right_g2<-tmppost.right_g2/sum(tmppost.right_g2)
	    #cand.s.v.k[tt]<-tmp.cand.s<-sample(1:parm$clust$K, size=1, replace=TRUE, prob=tmp.marg.v)
	    R2.new.star.right<-R2.new.star.right+log(tmp.marg.v.right_g2[tmp.s.v.right[tt]])
	  }
	  marg.log.lik.right <- sum(marg.log.lik.v.right_g2)
	  }

	  if (g.k.right==1) {
	    tmp.s.v.right<-unique(parm$clust$s.mt[,c.k.right])
	    if (length(tmp.s.v.right)>1) {stop(paste("problem in tmp.s.v.right when g.k.right=1"))}
	  #marginal likelihood of new cluster for clust$g.v=1
	  tmp.log.lik.v.right_g1 <- array(,length(parm$clust$phi.v))
	  for (pp in 1:length(parm$clust$phi.v)) {
	    tmp.log.lik.v.right_g1[pp] <- sum(dnorm(x.mt.right,mean=parm$clust$phi.v[pp], sd=(parm$tau/sqrt(parm$ni)),log=T))
	  }
	  tmp.log.marg.v.right_g1 <- tmp.log.lik.v.right_g1+log(P.aux)
	  marg.log.lik.right <- log(sum(exp(tmp.log.marg.v.right_g1)))
	  tmppost.right_g1<-tmp.log.marg.v.right_g1-max(tmp.log.marg.v.right_g1)
	  tmppost.right_g1<-exp(tmppost.right_g1)
	  tmp.marg.v.right_g1<-tmppost.right_g1/sum(tmppost.right_g1)
	  # cand.s.v.k_g1<-sample(1:parm$clust$K, size=1, replace=TRUE, prob=tmp.marg.v_g1)
	  # cand.s.v.k_g1<-rep(cand.s.v.k_g1,length(x.mt))
	  R2.new.star.right<-log(tmp.marg.v.right_g1[tmp.s.v.right])
	  }
	  
    L.v.right <- c(L.v.right, marg.log.lik.right)
  
	h.k<-parm$clust$gr.v[k]#==h_{k}
	if (parm$eta>0) {
	  tmp.prob.g.kp1<-parm$rho[h.k]+(1-parm$rho[h.k])*exp(-parm$distance[k]/parm$eta)
	}
	if (parm$eta==0) {
	  tmp.prob.g.kp1<-parm$rho[h.k]
	}
	# if (parm$eta==-1) {
	#   tmp.prob.g.kp1<-1
	# }
	prob.g.kp1.v<-rep(1-tmp.prob.g.kp1,2)
	prob.g.kp1.v[h.k]<-tmp.prob.g.kp1
	
	  tmp.prior.right <- (prob.g.kp1.v*cbind((parm$clust$C.m.mt-c(0,parm$d)), (parm$b1+right.G.gr*c(0,parm$d))))
	  tmp.prior.right[tmp.prior.right<0]<-0
	  log.prior.right<-log(tmp.prior.right)
	  
	  tmp2.right <- t(t(log.prior.right) + L.v.right)
	  maxx.right <- max(tmp2.right)
	  tmp2.right <- tmp2.right - maxx.right
	  
	  tmp2.right <- exp(tmp2.right)
	  
	  prop.prob.right<-tmp2.right
	  prop.prob.right<-prop.prob.right/sum(prop.prob.right)
	  
	  if (parm$clust$C.m.mt[g.k.right,c.k.right]==0){
	    q2.star<-log(prop.prob.right[g.k.right,parm$clust$G+1])+R2.new.star.right
	  }
	  if (parm$clust$C.m.mt[g.k.right,c.k.right]>0) {
	    q2.star<-log(prop.prob.right[g.k.right,c.k.right])
	  }
	  
	  parm$clust$C.m.mt[g.k.right,c.k.right]<-parm$clust$C.m.mt[g.k.right,c.k.right] + 1
	  #parm$G<-new.G
#####computing q2
	  if (!new.flag) {
	    old.L.v.right<-L.v.right
	  }
	  if (new.flag) {
	    old.L.v.right <- sapply(1:old.parm$clust$G, PDP_fn.log.lik, x.mt.right, old.parm)
	    old.L.v.right<-c(old.L.v.right,marg.log.lik.right)
	  }
	  
	  old.parm$clust$C.m.mt[g.k.right,c.k.right]<-old.parm$clust$C.m.mt[g.k.right,c.k.right] - 1
	  
	  old.G.gr<-rowSums(old.parm$clust$C.m.mt>0)

old.h.k<-old.parm$clust$gr.v[k]#==h_{k}
if (old.parm$eta>0) {
  old.tmp.prob.g.kp1<-old.parm$rho[old.h.k]+(1-old.parm$rho[old.h.k])*exp(-old.parm$distance[k]/parm$eta)
}
if (old.parm$eta==0) {
  old.tmp.prob.g.kp1<-old.parm$rho[old.h.k]
}
# if (old.parm$eta==-1) {
#   old.tmp.prob.g.kp1<-1
# }
old.prob.g.kp1.v<-rep(1-old.tmp.prob.g.kp1,2)
old.prob.g.kp1.v[old.h.k]<-old.tmp.prob.g.kp1

	  old.tmp.prior.right <- old.prob.g.kp1.v*cbind((old.parm$clust$C.m.mt-c(0,old.parm$d)), (old.parm$b1+old.G.gr*c(0,old.parm$d)))
	  old.tmp.prior.right[old.tmp.prior.right<0]<-0
	  old.log.prior.right<-log(old.tmp.prior.right)
	  
	  old.tmp2.right <- t(t(old.log.prior.right) + old.L.v.right)
	  old.maxx.right <- max(old.tmp2.right)
	  old.tmp2.right <- old.tmp2.right - old.maxx.right
	  
	  old.tmp2.right <- exp(old.tmp2.right)
	  
	  old.prop.prob.right<-old.tmp2.right
	  old.prop.prob.right<-old.prop.prob.right/sum(old.prop.prob.right)
	  
	  if (old.parm$clust$C.m.mt[g.k.right,c.k.right]==0){
	    q2<-log(old.prop.prob.right[g.k.right,old.parm$clust$G+1])+R2.new.star.right
	  }
	  if (old.parm$clust$C.m.mt[g.k.right,c.k.right]>0) {
	    q2<-log(old.prop.prob.right[g.k.right,c.k.right])
	  }
	  
	  old.parm$clust$C.m.mt[g.k.right,c.k.right]<-old.parm$clust$C.m.mt[g.k.right,c.k.right] + 1

####accept/reject proposal
  if (is.na(q2.star-q2)==F) {
	    mh_rate <- exp(min((q2.star-q2),0))
  }
  else {
    mh_rate <- 1
  }
	  
	  flip<-as.logical(rbinom(n=1,size=1,prob=mh_rate))
	  
	  if (!flip) {parm <- old.parm}
	}
	
	else {flip=T}
	
	parm$clust$C.m.vec<-colSums(parm$clust$C.m.mt)
	
	list(parm, new.flag, flip)
}

###########################################################
PDP_fn.gibbs <- function(k, parm, data, computeMode)
{	k <- parm$k
err <- PDP_fn.consistency.check(parm)
if (err > 0)
{stop(paste("GIBBS - 0: failed consistency check: err=",err))
}

## generate auxilliary P vector
tmp.M <- rep(parm$clust$M/parm$clust$K,parm$clust$K)
tmp.alpha <- tmp.M+parm$clust$n.vec
#tmp.alpha[emptied.s.indx] <- 0
P.aux <- rgamma(parm$clust$K,c(tmp.alpha),1)
P.aux <- P.aux/sum(P.aux)

old.c.k <- parm$clust$c.v[k]

###############

parm$clust$C.m.vec[old.c.k] <- parm$clust$C.m.vec[old.c.k] - 1

empty.indx <- which(parm$clust$C.m.vec==0)
#if (length(empty.indx)>0){stop(paste('drop problem inside, gibbs'))}

if (length(empty.indx) >0)
{
  #  if (computeMode$computeR) {
  new.s.mt <- parm$clust$s.mt[,-empty.indx] # HOT 3
  new.n.vec <- array(,parm$clust$K)
  for (pp in 1:parm$clust$K) {
    new.n.vec[pp] <- sum(new.s.mt==pp) # HOT
  }
  #new.n0 <- sum(new.s.mt == 0) # HOT 3
  #  }
  
  # if (computeMode$computeC) {
  #   
  #   tab <- .fastTabulateExcludeEmptiedIndices(parm$clust$s.mt, emptied.indx, parm$clust$K, TRUE)
  #   
  #   if (computeMode$computeR) { # debugging
  #     assertEqual(tab[1], new.n0)
  #     assertEqual(tab[-1], new.n.vec)
  #   }
  #   
  #   new.n0 <- tab[1]
  #   new.n.vec <- tab[-1]
  #   
  # } # computeMode
  
  # emptied.s.indx <- which(new.n.vec == 0)
  # new.K <- parm$clust$K - length(emptied.s.indx)
  
}

if (length(empty.indx) ==0)
{
  new.s.mt <- parm$clust$s.mt
  new.n.vec <- parm$clust$n.vec
  # emptied.s.indx <- which(new.n.vec==0)
  # new.K <- parm$clust$K
  #new.n0 <- parm$clust$n0
}
parm$clust$n.vec <- new.n.vec

x.mt <- matrix(parm$Z[,k], ncol=1)

if (computeMode$computeR) {
  
  # intercept cluster or any existing cluster
  L.v <- sapply(1:parm$clust$G, PDP_fn.log.lik, x.mt, parm) # HOT
  
}

if (computeMode$computeC) {
  
  test <- .computePdpLogLikelihood(computeMode$device$engine, k, parm$X,
                                   parm$clust$A.mt, parm$clust$s.mt,
                                   parm$clust$G, parm$n2,
                                   parm$tau, parm$tau_0, parm$tau_int, FALSE)
  
  if (computeMode$computeR) { # debugging
    assertEqual(test$logLikelihood, L.v, computeMode$tolerance)
  }
  
  L.v <- test$logLikehood
} # computeMode
# NB: returned logLikelihood differ from those computed above by approx 1e-15.  I believe this is due to non-transitivity of FLOPs

#######################################################
### emptied clusters are gone forever under Gibbs sampling
#######################################################

 emptied.indx <- which(parm$clust$C.m.vec==0)
 new.G <- parm$clust$G - length(emptied.indx)
# 
# if (length(emptied.indx) >0)
# {
#   if (computeMode$computeR) {
#     new.s.mt <- parm$clust$s.mt[,-emptied.indx] # HOT 3
#     new.n.vec <- array(,parm$clust$K)
#     for (pp in 1:parm$clust$K) {
#       new.n.vec[pp] <- sum(new.s.mt==pp) # HOT
#     }
#     #new.n0 <- sum(new.s.mt == 0) # HOT 3
#   }
#   
#   if (computeMode$computeC) {
#     
#     tab <- .fastTabulateExcludeEmptiedIndices(parm$clust$s.mt, emptied.indx, parm$clust$K, TRUE)
#     
#     if (computeMode$computeR) { # debugging
#       assertEqual(tab[1], new.n0)
#       assertEqual(tab[-1], new.n.vec)
#     }
#     
#     new.n0 <- tab[1]
#     new.n.vec <- tab[-1]
#     
#   } # computeMode
#   
#   emptied.s.indx <- which(new.n.vec == 0)
#   new.K <- parm$clust$K - length(emptied.s.indx)
#   
# }
# 
# if (length(emptied.indx) ==0)
# {
#   new.s.mt <- parm$clust$s.mt
#   new.n.vec <- parm$clust$n.vec
#   emptied.s.indx <- which(new.n.vec==0)
#   new.K <- parm$clust$K
#   #new.n0 <- parm$clust$n0
# }

## generate auxilliary P vector

# tmp.M <- rep(parm$clust$M/new.K,parm$clust$K)
# tmp.alpha <- tmp.M+new.n.vec
# tmp.alpha[emptied.s.indx] <- 0
# P.aux <- rgamma(parm$clust$K,c(tmp.alpha),1)
# P.aux <- P.aux/sum(P.aux)

## marginal likelihood of new cluster

if (computeMode$computeR) {
  
  marg.log.lik.v <- array(,length(x.mt))
  cand.s.v.k <- array(,length(x.mt))
  for (tt in 1:length(x.mt))
  {
    #if (!parm$flip.sign)
    {tmp.log.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=(parm$tau/sqrt(parm$ni[tt])),log=T) # HOT 3
    #tmp.lik.v <- c(dnorm(x.mt[tt],mean=0, sd=parm$tau_0),tmp.lik.v)
    }
    #       if (parm$flip.sign)
    #         {tmp.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=parm$tau) + dnorm(-x.mt[tt],mean=parm$clust$phi.v, sd=parm$tau)# HOT 3
    #         tmp.lik.v <- c((dnorm(x.mt[tt],mean=0, sd=parm$tau_0)+dnorm(-x.mt[tt],mean=0, sd=parm$tau_0)),tmp.lik.v)
    #         }
    
    #       marg.log.lik.v[tt] <- log(sum(tmp.lik.v*P.aux))
    #       tmp.prob.v <- tmp.lik.v*P.aux
    #       prob.gen.v <- tmp.prob.v/sum(tmp.prob.v)
    #       cand.s.v.k[tt]<-sample(1:parm$clust$K, size=1, replace=TRUE, prob=prob.gen.v) # HOT
    tmp.log.marg.v <- tmp.log.lik.v+log(P.aux)
    marg.log.lik.v[tt] <- log(sum(exp(tmp.log.marg.v)))
    tmppost<-tmp.log.marg.v-max(tmp.log.marg.v)
    tmppost<-exp(tmppost)
    tmp.marg.v<-tmppost/sum(tmppost)
    #       tmp.marg.v<-tmppost
    #       tmp.marg.v[tmp.marg.v < 1e-300] <- 1e-300
    #       tmp.marg.v[emptied.s.indx] <- 0
    #       tmp.marg.v <- tmp.marg.v/sum(tmp.marg.v)
    cand.s.v.k[tt]<-sample(1:parm$clust$K, size=1, replace=TRUE, prob=tmp.marg.v)
    
  }
  marg.log.lik <- sum(marg.log.lik.v)
  cand.gr<-1
  if (max(parm$clust$phi.v[cand.s.v.k])-min(parm$clust$phi.v[cand.s.v.k])>parm$delta_margin) {cand.gr<-2}
}

if (computeMode$computeC) {
  
  test <- .computeMarginalLikelihood(computeMode$device$engine,
                                     x.mt,
                                     parm$clust$phi.v,
                                     P.aux,
                                     parm$tau, parm$tau_0,
                                     FALSE, # no sampling
                                     computeMode$exactBitStream)
  
  if (computeMode$computeR) { # debugging
    assertEqual(test$logMarginalLikelihood, marg.log.lik, computeMode$tolerance)
  }
  
  marg.log.lik <- test$logMarginalLikelihood
  
}

L.v <- c(L.v, marg.log.lik)

# log.prior.v <- array(NA, (1+parm$clust$G))

# allow -Inf's in Gibbs sampling log-prior (just emptied clusters)
# if (length(emptied.indx) >0)
# {log.prior.v[-(emptied.indx)] <- log(c((parm$clust$C.m.vec[-emptied.indx]-parm$d), (parm$b1+new.G*parm$d)))
# log.prior.v[emptied.indx] <- -Inf
# }
# 
# if (length(emptied.indx) ==0)
# {log.prior.v <- log(c((parm$clust$C.m.vec-parm$d), (parm$b1+new.G*parm$d)))
# }
tmp.prior.v <- (c((parm$clust$C.m.vec-parm$d), (parm$b1+new.G*parm$d)))
tmp.prior.v[tmp.prior.v<0]<-0
log.prior.v<-log(tmp.prior.v)

tmp2 <- log.prior.v + L.v
maxx <- max(tmp2)
tmp2 <- tmp2 - maxx

tmp2 <- exp(tmp2)

parm$clust$post.k <- tmp2
parm$clust$post.k <- parm$clust$post.k / sum(parm$clust$post.k)

########################################################################
# store current state
old.parm <- parm

########################################################################

new.c.k <- sample(1:(parm$clust$G+1), size=1, replace=TRUE, prob=parm$clust$post.k)
parm$clust$c.v[k] <- new.c.k

new.flag <- new.c.k == (parm$clust$G+1)

#######################

# 	count.0 <- sum(new.c.k==0)
# 	parm$clust$C.m0 <- parm$clust$C.m0 + count.0

if ((!new.flag))
{parm$clust$C.m.vec[new.c.k] <- parm$clust$C.m.vec[new.c.k] + 1
parm$clust$gr.v[k]<-parm$clust$gr.c.v[new.c.k]
#parm$clust$n.vec<-new.n.vec
#parm$N <- sum(parm$clust$n.vec)
#if (parm$N!=parm$tot.trt * (new.G)){stop(paste("problem in n.vec"))}
}

if (new.flag)
{
  ###generate the latent vector first, condition on the single kth column
  #     cand.s.v.k <- array(,length(x.mt))
  # 
  #     for (tt in 1:length(x.mt))
  # 	  {
  #       if (!parm$flip.sign)
  #         {tmp.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=parm$tau)
  #       	tmp.lik.v <- c(dnorm(x.mt[tt],mean=0, sd=parm$tau_0),tmp.lik.v)
  #         }
  #       if (parm$flip.sign)
  #         {tmp.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=parm$tau) + dnorm(-x.mt[tt],mean=parm$clust$phi.v, sd=parm$tau)# HOT 3
  #         tmp.lik.v <- c((dnorm(x.mt[tt],mean=0, sd=parm$tau_0)+dnorm(-x.mt[tt],mean=0, sd=parm$tau_0)),tmp.lik.v)
  #         }
  # 
  #       tmp.prob.v <- tmp.lik.v*P.aux
  #       prob.gen.v <- tmp.prob.v/sum(tmp.prob.v)
  #       cand.s.v.k[tt]<-sample(0:parm$clust$K, size=1, replace=TRUE, prob=prob.gen.v) # HOT
  #     } # end for loop
  parm$clust$gr.v[k]<-cand.gr
  parm$clust$gr.c.v<-c(parm$clust$gr.c.v,cand.gr)
  parm$cand$s.v.k <- cand.s.v.k
  
  parm$cand$n.vec.k <- array(,parm$clust$K)
  for (gg in 1:parm$clust$K)
  {parm$cand$n.vec.k[gg] <- sum(cand.s.v.k==gg)
  }

  ##################
  parm$clust$G <- parm$clust$G + 1
  parm$clust$C.m.vec <- c(parm$clust$C.m.vec, 1)
  # 		parm$clust$beta.v <- c(parm$clust$beta.v, 0)
  # 		parm$clust$gamma.v <- c(parm$clust$gamma.v,0)
  
  parm$clust$s.v <- c(parm$clust$s.v, parm$cand$s.v.k)
  
  parm$clust$s.mt <- matrix(parm$clust$s.v, nrow = parm$tot.trt)
  #parm$clust$n.vec <- new.n.vec + parm$cand$n.vec.k
  parm$clust$n.vec <- parm$clust$n.vec + parm$cand$n.vec.k
  #parm$clust$n0 <- parm$clust$n0 + parm$cand$n0.k
  
  parm$N <- sum(parm$clust$n.vec) 
#  if (parm$N!=parm$tot.trt * sum(parm$clust$C.m.vec>0)){stop(paste("problem in n.vec"))}
  
  tmp.a.v <- array(,parm$tot.trt)
  s.G.v <- parm$cand$s.v.k
  #indxx <- s.G.v==0
  #tmp.a.v[indxx] <- 0
  tmp.a.v <- parm$clust$phi.v[s.G.v]
  #
  parm$clust$A.mt <- cbind(parm$clust$A.mt, tmp.a.v) # HOT
  #parm$clust$B.mt <- cbind(parm$clust$B.mt, tmp.a.v) # HOT
  
  # 		if (parm$tBB_flag)
  # 		{
  # 		if (computeMode$computeR) {
  # 		  parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt # HOT
  # 		} else {
  # 		  parm$clust$tBB.mt <- .fastXtX(parm$clust$B.mt)
  # 		}
  # 		}
  
} # end  if (new.flag)


list(parm, new.flag)
}

###########################################################

PDP_fn.fast_col <- function(cc, parm, data, computeMode)
{
  k <- parm$k <- parm$clust$col.nbhd.k[[cc]]
  I.k <- parm$clust$col.nbhd[[cc]]

   err <- PDP_fn.consistency.check(parm)
   if (err > 0)
   {stop(paste("FAST - 0: failed consistency check: err=",err))
   }
   # store so that we can revert to this state if MH propsal is rejected
   init.cc.parm <- parm
   

   ## generate auxilliary P vector
   tmp.M <- rep(parm$clust$M/parm$clust$K,parm$clust$K)
   tmp.alpha <- tmp.M+parm$clust$n.vec
   #tmp.alpha[emptied.s.indx] <- 0
   P.aux <- rgamma(parm$clust$K,c(tmp.alpha),1)
   P.aux <- P.aux/sum(P.aux)
   

  old.c.k <- parm$clust$c.v[I.k]

  # Note to MS: I've ignored the branching conditions
  #             in "elementwise_DP.functions.R" based on computeMode$useR

  if (computeMode$computeR) {

    parm$clust$C.m.vec.k <- array(,parm$clust$G)
    for (gg in 1:parm$clust$G) {
      parm$clust$C.m.vec.k[gg] <- sum(old.c.k == gg)
    }

    #parm$clust$C.m0.k <- sum(old.c.k == 0)
  }

  if (computeMode$computeC) {

    all.n.vec <- .fastTabulateVector(old.c.k, parm$clust$G, TRUE)

    if (computeMode$computeR) { # debugging
      assertEqual(all.n.vec[1], parm$clust$C.m0.k)
      assertEqual(all.n.vec[-1], parm$clust$C.m.vec.k)
    }

    parm$clust$C.m0.k <- all.n.vec[1] # test1
    parm$clust$C.m.vec.k <- all.n.vec[-1] # test2
  }

  parm$clust$C.m.vec.k.comp <- parm$clust$C.m.vec - parm$clust$C.m.vec.k
  #parm$clust$C.m0.k.comp <- parm$clust$C.m0 - parm$clust$C.m0.k  

  x.mt <- matrix(parm$Z[,k], ncol=1)

  if (computeMode$computeR) {

    # intercept cluster or any existing cluster
    L.v <- sapply(1:parm$clust$G, PDP_fn.log.lik, x.mt, parm) # HOT
  }

  if (computeMode$computeC) {
    test <- .computePdpLogLikelihood(computeMode$device$engine, k, parm$X,
                                     parm$clust$A.mt, parm$clust$s.mt,
                                     parm$clust$G, parm$n2,
                                     parm$tau, parm$tau_0, parm$tau_int, FALSE)

    if (computeMode$computeR) { # debugging
      assertEqual(test$logLikelihood, L.v, computeMode$tolerance)
    }

    L.v <- test$logLikehood
  } # computeMode
  # NB: returned logLikelihood differ from those computed above by approx 1e-15.  I believe this is due to non-transitivity of FLOPs

    emptied.indx <- which(parm$clust$C.m.vec.k.comp==0)
    new.G <- parm$clust$G - length(emptied.indx)
  
  if (length(emptied.indx) >0)
  {
    if (computeMode$computeR) {
      new.s.mt <- parm$clust$s.mt[,-emptied.indx] # HOT 3
      new.n.vec <- array(,parm$clust$K)
      for (pp in 1:parm$clust$K) {
        new.n.vec[pp] <- sum(new.s.mt==pp) # HOT
      }
      #new.n0 <- sum(new.s.mt == 0) # HOT 3
    }
  
  #   if (computeMode$computeC) { # computeMode
  #     tab <- .fastTabulateExcludeEmptiedIndices(parm$clust$s.mt, emptied.indx, parm$clust$K, TRUE)
  # 
  #     if (computeMode$computeR) { # debugging
  #       assertEqual(tab[1], new.n0)
  #       assertEqual(tab[-1], new.n.vec)
  #     }
  # 
  #     new.n0 <- tab[1]
  #     new.n.vec <- tab[-1]
  # 
  #   } # computeMode
  # 
  #   emptied.s.indx <- which(new.n.vec == 0)
  #   new.K <- parm$clust$K - length(emptied.s.indx)
  }
  # 
  if (length(emptied.indx) ==0)
  {
    new.s.mt <- parm$clust$s.mt
    new.n.vec <- parm$clust$n.vec
#    emptied.s.indx <- which(new.n.vec==0)
#    new.K <- parm$clust$K
    #new.n0 <- parm$clust$n0
  }

parm$clust$n.vec <- new.n.vec

  ## generate auxilliary P vector

  # tmp.M <- rep(parm$clust$M/new.K,parm$clust$K)
  # tmp.alpha <- tmp.M+new.n.vec
  # tmp.alpha[emptied.s.indx] <- 0
  # P.aux <- rgamma(parm$clust$K,c(tmp.alpha),1)
  # P.aux <- P.aux/sum(P.aux)

  # START

  if (computeMode$computeR) {

    if (computeMode$computeC) { # debugging
      savedSeed <- .GlobalEnv$.Random.seed # For debugging purposed only
    }

  ## marginal likelihood of new cluster
  marg.log.lik.v <- cand.s.v.k <- array(,length(x.mt))
  for (tt in 1:length(x.mt))
  {
    #if (!parm$flip.sign)
      #{
        tmp.log.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=(parm$tau/sqrt(parm$ni[tt])),log=T) # HOT 3
      #tmp.lik.v <- c(dnorm(x.mt[tt],mean=0, sd=parm$tau_0),tmp.lik.v)
      #}
#     if (parm$flip.sign)
#       {tmp.lik.v <- dnorm(x.mt[tt],mean=parm$clust$phi.v, sd=parm$tau) + dnorm(-x.mt[tt],mean=parm$clust$phi.v, sd=parm$tau)# HOT 3
#       tmp.lik.v <- c((dnorm(x.mt[tt],mean=0, sd=parm$tau_0)+dnorm(-x.mt[tt],mean=0, sd=parm$tau_0)),tmp.lik.v)
#       }

    tmp.log.marg.v <- tmp.log.lik.v+log(P.aux)
    marg.log.lik.v[tt] <- log(sum(exp(tmp.log.marg.v)))
    tmppost<-tmp.log.marg.v-max(tmp.log.marg.v)
    tmppost<-exp(tmppost)
    tmp.marg.v<-tmppost/sum(tmppost)
    #tmp.marg.v<-tmppost
# tmp.marg.v[tmp.marg.v < 1e-300] <- 1e-300
# tmp.marg.v[emptied.s.indx] <- 0
# tmp.marg.v <- tmp.marg.v/sum(tmp.marg.v)
    cand.s.v.k[tt]<-sample(1:parm$clust$K, size=1, replace=TRUE, prob=tmp.marg.v)
  }

  parm$cand$s.v.k <- cand.s.v.k
  marg.log.lik <- sum(marg.log.lik.v)
  cand.gr<-1
  if (max(parm$clust$phi.v[cand.s.v.k])-min(parm$clust$phi.v[cand.s.v.k])>parm$delta_margin) {cand.gr<-2}
  
  # SWTICH
  }

  if (computeMode$computeC) {

    if (computeMode$computeR) { # debugging
      .GlobalEnv$.Random.seed <- savedSeed # Roll back PRNG
    }

  test <- .computeMarginalLikelihood(computeMode$device$engine,
                                     x.mt,
                                     parm$clust$phi.v,
                                     P.aux,
                                     parm$tau, parm$tau_0,
                                     TRUE, # do sampling
                                     computeMode$exactBitStream)

    if (computeMode$computeR) { # debugging
      assertEqual(test$logMarginalLikelihood, marg.log.lik, computeMode$toleranace)
      assertEqual(test$sVk, cand.s.v.k)
    }

  parm$cand$s.v.k <- test$sVk
  marg.log.lik <- test$logMarginalLikelihood

  }

  # END

  L.v <- c(L.v, marg.log.lik)

 #log.prior.v <- array(NA, (1+parm$clust$G))

  spread.mass <- (parm$b1+new.G*parm$d)/(1+length(emptied.indx))

   # if (length(emptied.indx) >0)
   # {
   #  log.prior.v[-(emptied.indx)] <- log(c( (parm$clust$C.m.vec[-emptied.indx]-parm$d), spread.mass))
   #  log.prior.v[emptied.indx] <- log(spread.mass)
   # }
   # 
   # if (length(emptied.indx) ==0)
   # {log.prior.v <- log(c( (parm$clust$C.m.vec.k.comp-parm$d), spread.mass))
   # }
  tmp.prior.v <- (c( (parm$clust$C.m.vec.k.comp-parm$d), spread.mass))
  tmp.prior.v[tmp.prior.v<0]<-spread.mass
  log.prior.v<-log(tmp.prior.v)

   tmp2 <- log.prior.v + L.v
   maxx <- max(tmp2)
   tmp2 <- tmp2 - maxx

   tmp2 <- exp(tmp2)

   parm$clust$post.k <- tmp2
   parm$clust$post.k <- parm$clust$post.k / sum(parm$clust$post.k)

   ########################################################################
   # store current state
   old.parm <- parm

   ########################################################################

   new.c.k <- sample(1:(parm$clust$G+1), size=length(I.k), replace=TRUE, prob=parm$clust$post.k)

  exit <- (sum(new.c.k != old.c.k)==0)
  flip <- TRUE
  new.flag <- FALSE

  if (!exit) # CONTINUE W/O EXITING FUNCTION
  {

    parm$clust$c.v[I.k] <- new.c.k

      new.count <- sum(new.c.k == (parm$clust$G+1))

      #######################

      new.prop <- 0

      if (computeMode$computeR) {

#       count.0 <- sum(new.c.k==0)
#       parm$clust$C.m0 <- parm$clust$C.m0.k.comp + count.0
#       if (count.0 >0)
#         {new.prop <- new.prop + log(parm$clust$post.k[1])*count.0
#         }

      for (gg in 1:parm$clust$G)
      {count.gg <- sum(new.c.k==gg)
       parm$clust$C.m.vec[gg] <- parm$clust$C.m.vec.k.comp[gg] + count.gg

       if (count.gg > 0)
          {new.prop <- new.prop + log(parm$clust$post.k[gg])*count.gg
          }
      }
#        parm$clust$n.vec<-new.n.vec
        
#        parm$N <- sum(parm$clust$n.vec) 
        #if (parm$N!=parm$tot.trt * (new.G)){stop(paste("problem in n.vec"))}
      }

      if (computeMode$computeC) {

        all.count <- .fastTabulateVector(new.c.k, parm$clust$G + 1, TRUE)
        test1 <- parm$clust$C.m0.k.comp + all.count[1] # test1
        test2 <- parm$clust$C.m.vec.k.comp + all.count[2:(parm$clust$G + 1)] # test2
        test3 <- .fastSumSafeLog(parm$clust$post.k, all.count, parm$clust$G + 1) # test3

        if (computeMode$computeR) { # debugging
          assertEqual(test1, parm$clust$C.m0)
          assertEqual(test2, parm$clust$C.m.vec)
          assertEqual(test3, new.prop)
        }

        parm$clust$C.m0 <- test1
        parm$clust$C.m.vec <- test2
        new.prop <- test3


      }

    if (new.count > 0)
      {parm$clust$gr.c.v<-c(parm$clust$gr.c.v,cand.gr)
#      parm$clust$gr.v[I.k]<-parm$clust$gr.c.v[new.c.k]
      parm$clust$G <- parm$clust$G + 1

      parm$clust$C.m.vec <- c(parm$clust$C.m.vec, new.count)

        new.prop <- new.prop + log(parm$clust$post.k[parm$clust$G])*new.count

        parm$cand$n.vec.k <- array(,parm$clust$K)
        for (ss in 1:parm$clust$K)
        {parm$cand$n.vec.k[ss] <- sum(parm$cand$s.v.k==ss)
        }
        #parm$cand$n0.k <- sum(parm$cand.s.v.k==0)

        ##################


#         parm$clust$beta.v <- c(parm$clust$beta.v, 0)
#         parm$clust$gamma.v <- c(parm$clust$gamma.v,0)

        parm$clust$s.v <- c(parm$clust$s.v, parm$cand$s.v.k)

        parm$clust$s.mt <- matrix(parm$clust$s.v, nrow = parm$tot.trt)

        parm$clust$n.vec <- parm$clust$n.vec + parm$cand$n.vec.k
        #parm$clust$n0 <- parm$clust$n0 + parm$cand$n0.k

        parm$N <- sum(parm$clust$n.vec) #+ parm$clust$n0
#        if (parm$N!=parm$tot.trt * sum(parm$clust$C.m.vec>0)){stop(paste("problem in n.vec"))}
        
        tmp.a.v <- array(,parm$tot.trt)
        s.G.v <- parm$cand$s.v.k
        #indxx <- s.G.v==0
        #tmp.a.v[indxx] <- 0
        tmp.a.v <- parm$clust$phi.v[s.G.v]
        #
        parm$clust$A.mt <- cbind(parm$clust$A.mt, tmp.a.v) # HOT
        #parm$clust$B.mt <- cbind(parm$clust$B.mt, tmp.a.v)

#         if (parm$tBB_flag)
#         {
#         if (computeMode$computeR) {
#           parm$clust$tBB.mt <- t(parm$clust$B.mt) %*% parm$clust$B.mt # HOT
#         } else {
#           parm$clust$tBB.mt <- .fastXtX(parm$clust$B.mt)
#         }
#         }

      } # end  if (new.count > 0)
      parm$clust$gr.v[I.k]<-parm$clust$gr.c.v[new.c.k]
    ######################

    # sum(parm$clust$C.m.vec) + parm$clust$C.m0 == parm$p



    ##########################################
    ##########################################
    ##### Computing proposal prob of reverse move
    ##########################################
    ##########################################

    old.prop <- 0

    for (gg in 1:old.parm$clust$G)
    {flag.gg <- (old.c.k==gg)
     count.gg <- sum(flag.gg)

     if (count.gg > 0)
     {old.prop <- old.prop + log(old.parm$clust$post.k[gg])*count.gg
     }
    }

    rho.prop <- new.prop - old.prop

    #######################################################
    #######################################################
    ########## computing true log-ratio: 2015
    #######################################################
    #######################################################

    # formula on page 1 of 07/27/15 notes

    # need to ensure there are no empty clusters in either
    # init.cc.parm or parm, otherwise likelihood formula
    # doesn't work (lgamma of negative values)

    tmp.new.parm <- parm
    indx.new <- parm$clust$C.m.vec > 0
    tmp.new.parm$clust$G <- sum(indx.new)
    tmp.new.parm$clust$C.m.vec <- parm$clust$C.m.vec[indx.new]
    tmp.new.parm$clust$s.mt <- parm$clust$s.mt[,indx.new]
#    tmp.new.parm$clust$n.vec <- tmp.new.parm$clust$n.vec[tmp.new.parm$clust$n.vec>0]
    #if (new.count>0) {tmp.new.parm$clust$n.vec<-tmp.new.parm$clust$n.vec+parm$can$n.vec.k}
    #
    tmp.old.parm <- init.cc.parm
    indx.old <- init.cc.parm$clust$C.m.vec > 0
    tmp.old.parm$clust$G <- sum(indx.old)
    tmp.old.parm$clust$C.m.vec <- init.cc.parm$clust$C.m.vec[indx.old]
    tmp.old.parm$clust$s.mt <- init.cc.parm$clust$s.mt[,indx.old]
#    tmp.old.parm$clust$n.vec <- tmp.old.parm$clust$n.vec[tmp.old.parm$clust$n.vec>0]


    rho.tru <- fn.d(d=parm$d, tmp.new.parm)[[2]] - fn.d(d=parm$d, tmp.old.parm)[[2]]

    new.log.lik <- 0

    for (gg in new.c.k)
    {indx.gg <- new.c.k==gg
     x_gg.mt <- matrix(parm$Z[,I.k[indx.gg]], ncol=sum(indx.gg))
     new.log.lik <- new.log.lik + sum(PDP_fn.log.lik(gg, x.mt=x_gg.mt, parm))
    }

    old.log.lik <- 0
    for (gg in old.c.k)
    {indx.gg <- old.c.k==gg
     x_gg.mt <- matrix(parm$Z[,I.k[indx.gg]], ncol=sum(indx.gg))
     old.log.lik <- old.log.lik + sum(PDP_fn.log.lik(gg, x.mt=x_gg.mt, old.parm))
    }

    rho.tru <- rho.tru + new.log.lik - old.log.lik

    ########## toss a coin #################

    prob <- exp(min((rho.tru-rho.prop),0))

    flip<-as.logical(rbinom(n=1,size=1,prob=prob))

if (!flip) {parm <- init.cc.parm}


  } # end BIG if (!exit) loop


   list(parm, new.flag, exit, flip)
}


#####################################


PDP_fn.drop <- function(parm, computeMode, mt=T) 
 	{

	##########################################
	## Drop empty clusters:
	## (i)  Move empty clusters to end by relabeling clusters
	## (ii) Set parm$clust$G equal to number of non-empty clusters
	## (ii) Retain only clusters  1,...,parm$clust$G
	#########################################

  parm$clust$G <- sum(parm$clust$C.m.vec>0)
  num.dropped <- sum(parm$clust$C.m.vec==0)

  if (parm$clust$G > 0)
  {
	if (num.dropped > 0)
	{
	for (rr in 1:num.dropped)
		{
		old.label <-  min(which(parm$clust$C.m.vec==0))
		new.label <- max(which(parm$clust$C.m.vec>0))
		stopp <-  max(which(parm$clust$C.m.vec>0)) == parm$clust$G
		if (stopp)
			{break
			}
    if (mt){
      parm <- PDP_fn.swap.clusters(parm, g1 = new.label, g2 = old.label, computeMode)
    }
    if (!mt) {
      parm <- PDP_fn.swap.clusters(parm, g1 = new.label, g2 = old.label, computeMode, mt=F)
    }
      
 		}
	}

	##########

	keep <- 1:parm$clust$G
  if (mt){
    parm <- PDP_fn.clip.clusters(parm, keep)
  }
  if (!mt) {
    parm <- PDP_fn.clip.clusters(parm, keep, mt=F)
  }

	###########

	# parm$clust$K does not change (possibly some empty elementwise clusters)

#	if (parm$N != parm$tot.trt * parm$clust$G){stop(paste("problem in n.vec"))}
	tmp.s.mt <- parm$clust$s.mt
	tmp.gr.c.v <- parm$clust$gr.c.v
	tmp.s.mt_g1<-tmp.s.mt[,tmp.gr.c.v==1,drop=F]
#	tmp.s.mt_g2<-tmp.s.mt[,tmp.gr.c.v==2]
	if (sum(apply(tmp.s.mt_g1,2,max)-apply(tmp.s.mt_g1,2,min))>0) {stop(paste('problem in PDP_fn.drop-tmp.s.mt_g1'))}
#	tmp.s.mt_g1_comp<-apply(tmp.s.mt_g1,2,unique)
	
	parm$N <- parm$tot.trt * sum(parm$clust$gr.c.v==2) + sum(parm$clust$gr.c.v==1)
#	parm$clust$s.v <- as.vector(parm$clust$s.mt)
	parm$clust$s.stretch <- parm$clust$s.v<-NULL
	for (g in 1:parm$clust$G) {
	  if (parm$clust$gr.c.v[g]==1) {
	    parm$clust$s.v<-c(parm$clust$s.v,unique(parm$clust$s.mt[,g]))
	    parm$clust$s.stretch<-c(parm$clust$s.stretch,1)
	  }
	  if (parm$clust$gr.c.v[g]==2) {
	    parm$clust$s.v<-c(parm$clust$s.v,(parm$clust$s.mt[,g]))
	    parm$clust$s.stretch<-c(parm$clust$s.stretch,rep(0,parm$tot.trt))
	  }
	}
	if (sum(c(length(parm$clust$s.v)!=parm$N,length(parm$clust$s.stretch)!=parm$N)>0))
	{stop(paste('problem in fast_PDP-PDP_fn.drop'))}
	

	#parm$clust$n0 <- sum(parm$clust$s.v==0)
	parm$clust$n.vec <- array(,parm$clust$K)
	for (ss in 1:parm$clust$K)
		{parm$clust$n.vec[ss] <- sum(parm$clust$s.v==ss)  # HOT
		}

  }
	parm
	}

PDP_fn.orientation <- function(parm, cc_subset)
  {
  X.mt <- matrix(parm$X[,cc_subset], ncol=length(cc_subset))
  c.v <- parm$clust$c.v[cc_subset]
  orient.v <- parm$clust$orient.v[cc_subset]

  # manipulate PDP_fn.log.lik to get log-likelihoods separately for each sign
  tmp.parm <- parm
  tmp.parm$flip.sign=FALSE

  log_lik.mt <- array(, c(2, length(cc_subset)))

  for (gg in 0:parm$clust$G)
  {indx.gg <- which(c.v==gg)

    if (length(indx.gg)>0)
      {X_gg.mt <- matrix(X.mt[,indx.gg], ncol=length(indx.gg))
      log_lik.mt[1,indx.gg] <- PDP_fn.log.lik(gg, x.mt = X_gg.mt, tmp.parm)
      log_lik.mt[2,indx.gg] <- PDP_fn.log.lik(gg, x.mt = -X_gg.mt, tmp.parm)
      }
  }

  maxx.v <- apply(log_lik.mt, 2, max)

  log_lik.mt <- t(t(log_lik.mt) - maxx.v)
  lik.mt <- exp(log_lik.mt)


  # assuming equal prior prob to each sign

  for (tt in 1:length(cc_subset))
    {orient.v[tt] <- sample(c(-1,1),size=1,prob=lik.mt[,tt])
    }

  parm$clust$orient.v[cc_subset] <- orient.v

    parm
  }

PDP_fn.eta.log.lik <- function(eta,parm,prec=F) {
  if (length(parm$clust$gr.v)!=parm$p) {stop(paste('error in length of gr.v'))}
  h.v <- parm$clust$gr.v[-parm$p]
  rho.v<-parm$rho[h.v]
  if (eta>0) {
    tmp.r.v<-exp(-parm$distance/eta)
  }
  if (eta==0) {tmp.r.v<-0}
  #if (eta==-1) {tmp.r.v<-1}
  tmp.prob.v1<-rho.v+(rep(1,length(rho.v))-rho.v)*tmp.r.v
  tmp.prob.v2<-1-tmp.prob.v1
  tmp.prob.mt<-rbind(tmp.prob.v1,tmp.prob.v2)
  tmp.idx.v<-rep(2,parm$p-1)
  tmp.idx.v[parm$clust$g.v[-1]==h.v]<-1
  prob.v<-diag(tmp.prob.mt[tmp.idx.v,])
  log.lik<-sum(log(prob.v))
  
  prop.sd<-NA
  if (prec==T) {
  if (eta>0) {
    tmp.b<-parm$distance/eta-2+2*tmp.r.v
    tmp.p1<-(rho.v*tmp.b-2*tmp.r.v)/(rho.v+(1-rho.v)*tmp.r.v)
    tmp.p2<-tmp.b/(1-tmp.r.v)
    tmp.2deri.v<-(tmp.r.v*(1-rho.v)*(parm$distance)*eta^(-3))*(tmp.p1-tmp.p2)
    #  tmp.2deri.1<-(parm$distance^2)*(1-rho.v)*rho.v*tmp.r.v*(rho.v+(1-rho.v)*tmp.r.v)^(-2)
    #  tmp.2deri.2<--(parm$distance^2)*tmp.r.v*(1-tmp.r.v)^(-2)
    #  tmp.2deri.mt<-rbind(tmp.2deri.1,tmp.2deri.2)
    prop.prec<--sum(tmp.2deri.v)
    if (prop.prec<0) {stop(paste('proposal precision negative'))}
    prop.sd<-(prop.prec)^(-0.5)
  }

  }
  return(list(log.lik,prop.sd))
}

fn.eta <- function(parm){
  ##gamma prior for eta
#  a<-1 ##shape
  b<-1/(parm$p-1) ##rate
  
#  tmp.old<-PDP_fn.eta.log.lik(parm$eta_R,parm)
#  tmp.marg_monte<-rtruncnorm(1000,a=0,b=2,mean=parm$eta_R,sd=parm$mult*tmp.old[[2]])
  grid.loglik.v<-NULL
  grid.len<-2500
  
  if(parm$iter_num==0|(parm$iter_num<=50 & parm$iter_num%%2==0)|(parm$iter_num>50 & parm$iter_num%%10==0)){
    grid.v<-seq(0,0.1,length.out=grid.len+1)
    for (ttt in 1:length(grid.v)) {
      grid.loglik.v<-c(grid.loglik.v,PDP_fn.eta.log.lik(grid.v[ttt],parm,prec=F)[[1]])
    }
#    log.grid.prior.v<-dgamma(grid.v,shape=a,rate=b,log=T)-log(sum(dgamma(grid.v,shape=a,rate=b)))
    pre.prior.v<-b*(grid.v[-1]^(-2))*exp(-b/grid.v[-1])
    pre.prior.v<-pre.prior.v/sum(pre.prior.v)
    log.grid.prior.v<-log(0.5)+c(0,log(pre.prior.v))
    log.grid.post.v<-grid.loglik.v+log.grid.prior.v
    log.grid.post.v<-log.grid.post.v-max(log.grid.post.v)
    
    grid.post.v<-exp(log.grid.post.v)
    grid.post.v<-grid.post.v/sum(grid.post.v)
    
    parm$eta<-sample(grid.v,size=1,prob=grid.post.v)
    
    tmp.grid.post_no0<-grid.post.v[-1]/sum(grid.post.v[-1])
    tmp.grid.p<-cumsum(tmp.grid.post_no0)
    # if (sum(tmp.grid.post_no0)!=1) {stop(paste('error in cumsum'))}
    tmp.diff.mt<-abs(matrix(rep(seq(0.005,0.995,by=0.005),each=grid.len),grid.len,199)-tmp.grid.p)
    grid.filtered.idx<-apply(tmp.diff.mt,2,which.min)
    grid.filtered<-grid.v[c(1,(grid.filtered.idx+1))]
    parm$eta_prob<-grid.post.v[c(1,(grid.filtered.idx+1))]
    parm$eta_grid.filtered<-grid.filtered
    if(length(parm$eta_prob)!=200) {stop(paste('error in eta_prob'))}
  }
  
  else {
    for (ttt in 1:length(parm$eta_grid.filtered)) {
      grid.loglik.v<-c(grid.loglik.v,PDP_fn.eta.log.lik(parm$eta_grid.filtered[ttt],parm,prec=F)[[1]])
    }
#    log.grid.prior.v<-dgamma(parm$eta_grid.filtered,shape=a,rate=b,log=T)-log(sum(dgamma(parm$eta_grid.filtered,shape=a,rate=b)))
    pre.prior.v<-b*(parm$eta_grid.filtered[-1]^(-2))*exp(-b/parm$eta_grid.filtered[-1])
    pre.prior.v<-pre.prior.v/sum(pre.prior.v)
    log.grid.prior.v<-log(0.5)+c(0,log(pre.prior.v))
    log.grid.post.v<-grid.loglik.v+log.grid.prior.v
    log.grid.post.v<-log.grid.post.v-max(log.grid.post.v)
    
    grid.post.v<-exp(log.grid.post.v)
    grid.post.v<-grid.post.v/sum(grid.post.v)
    
    parm$eta<-sample(parm$eta_grid.filtered,size=1,prob=grid.post.v)
    parm$eta_prob<-grid.post.v
    
  }
  
  
  log.grid.post.v.2<-c(log.grid.post.v[1],log(sum(exp(log.grid.post.v[-1]))))
  
  
  parm$eta.0_log.BF<-log.grid.post.v.2[1]-log.grid.post.v.2[2]
  
  
  parm
}

fast_PDP_fn.main <- function(parm, data, col.frac.probes, prob.compute.col.nbhd, max.col.nbhd.size, computeMode)
{
  p <- parm$p
  
  parm$Z.f<-pre.Z<-parm$X.LANE-matrix(rep(parm$epsilon,parm$p),ncol=parm$p)
  #data.X.tmp<-tmp.Z
  tmp.Z<-matrix(,parm$tot.trt,parm$p)
  for (ii in 1:parm$tot.trt){
    tmp.Z[ii,]<-colMeans(pre.Z[parm$group==ii,,drop=F])
  }
  parm$Z<-tmp.Z
  
  new.flag.v <- NULL
  col.mh.flip.v <- col.mh.exit.v <- order_1.mh.flip.v<- NULL
  
#  if (parm$order==1){
    
    for (cc in 1:parm$p){
      parm$k<-cc
      
      tmp <- PDP_fn.gibbs_order_1(k=parm$k, parm, data, computeMode)
      parm <- tmp[[1]]
      new.flag.v <- c(new.flag.v, tmp[[2]])
      order_1.mh.flip.v<-c(order_1.mh.flip.v, tmp[[3]])
    }
    
#  }
  
#   if (parm$order==0){
#     ##########################
#     # compute delta-neighborhoods
#     #########################
#     
#     # SG: prob.compute.col.nbhd is the probability of computing
#     # the neighborhoods, and col.frac.probes is the fraction of neighborhoods updated
#     
#    col_flip <- as.logical(rbinom(n=1,size=1,prob=prob.compute.col.nbhd))
#    if (is.null(parm$clust$col.nbhd.k) | col_flip){
#      parm <- PDP_fn.post.prob.and.delta(parm, max.col.nbhd.size, col.frac.probes, computeMode)
#    }
# 
# 
#  	if (col.frac.probes < 1)
#  	{num_nbhds <- max(1,round(col.frac.probes*length(parm$clust$col.nbhd.k)))
#  	  parm$subset_nbhd.indx <- sort(sample(1:length(parm$clust$col.nbhd.k), size=num_nbhds))
#  	}
# 
#  	if (col.frac.probes == 1)
#  	{parm$subset_nbhd.indx <- 1:length(parm$clust$col.nbhd.k)
#  	}
#     
# #   if (parm$standardize.X)
# #     {parm <- fn.standardize_orient.X(parm)
# #     }
# 
#      for (cc in parm$subset_nbhd.indx)
# 			{
#       previous.parm <- parm
# 
# 			parm$k <- parm$clust$col.nbhd.k[[cc]]
# 			
#       if (length(parm$clust$col.nbhd[[cc]])==1)
#       {tmp <- PDP_fn.gibbs(k=parm$k, parm, data, computeMode)
#       parm <- tmp[[1]]
#       new.flag.v <- c(new.flag.v, tmp[[2]])
#       }
# 
# 			if (length(parm$clust$col.nbhd[[cc]])>1)
# 			{
# 			  tmp <- PDP_fn.fast_col(cc, parm, data, computeMode)
# 			  parm <- tmp[[1]]
# 			  new.flag.v <- c(new.flag.v, tmp[[2]])
# 			  col.mh.exit.v <- c(col.mh.exit.v, tmp[[3]])
# 			  col.mh.flip.v <- c(col.mh.flip.v, tmp[[4]])
# 			}
# 			
# 			} # end for loop
#    
#    parm$clust$C.m.mt<- matrix(NA,2,parm$clust$G)
#    for (g in 1:parm$clust$G)
#    {I.g <- which(parm$clust$c.v==g)
#    p.I.g<-setdiff((I.g-1),0)
#    parm$clust$C.m.mt[1,g]<-sum(parm$clust$gr.v[p.I.g]==1)
#    parm$clust$C.m.mt[2,g]<-sum(parm$clust$gr.v[p.I.g]==2)
#    }
#    ##adding the count from the first column
#    parm$clust$C.m.mt[1,parm$clust$c.v[1]]<-parm$clust$C.m.mt[1,parm$clust$c.v[1]]+1
# 
# } #end order==0

#   if (parm$flip.sign)
#   {
#     ############
#     ## Update signs for updated columns
#     ############
#     parm <- PDP_fn.orientation(parm, cc_subset=parm$subset_nbhd.indx)
# 
#     ############
#     # re-orient columns
#     ############
#     parm <- fn.standardize_orient.X(parm)
#   }

	err <- PDP_fn.consistency.check(parm)
	if (err > 0)
			{stop(paste("LOOP: failed consistency check: err=",err))
	}

  parm<-fn.eta(parm)

	parm$clust$col.new.flag <- mean(new.flag.v)
  if (!is.null(col.mh.flip.v))
    {parm$clust$col.mh.flip <- mean(col.mh.flip.v)
    parm$clust$col.mh.exit <- mean(col.mh.exit.v)
    }
  else
    {parm$clust$col.mh.flip <- 1
    parm$clust$col.mh.exit <- 1
    }
	if (!is.null(order_1.mh.flip.v))
	{parm$clust$order_1.mh.flip <- mean(order_1.mh.flip.v)
	}
	else
	{parm$clust$order_1.mh.flip <- NA
	}
	

	##########################################
	## Drop empty group clusters:
	## (i)  Move empty clusters to end by relabeling clusters
	## (ii) Set parm$clust$G equal to number of non-empty clusters
	## (ii) Retain only clusters  1,...,parm$clust$G
	#########################################

	parm <- PDP_fn.drop(parm, computeMode)

	# now drop empty elementwise clusters

	parm <- element_fn.drop(parm)

  parm$theta.mt<-parm$clust$A.mt[,parm$clust$c.v]

	err <- PDP_fn.compact.consistency.check(parm)
	if (err > 0)
		{stop(paste("MAIN FUNCTION END: failed consistency check: err=",err))
		}

	parm
	}



