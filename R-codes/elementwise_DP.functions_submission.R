

element_fn.compact_consistency.check <- function(parm, computeMode)
	{err <- 0

	if (max(unique(parm$clust$s.v)) != parm$clust$K)
		{err <- 1
		}

	err <- element_fn.consistency.check(parm, computeMode)

	err
	}


element_fn.consistency.check <- function(parm, computeMode)
	{err <- 0

	if (!is.null(parm$clust$row.nbhd)) {

	  if (computeMode$computeR) {
	    if (sum(unlist(lapply(parm$clust$row.nbhd, length))) != length(parm$row.subset.I)) {
	      err <- 2
	    }
	  }

	  if (computeMode$computeC) {
      if (parm$clust$row.nbhd.flat.length != length(parm$row.subset.I))
        err <- 2
	  }
	}

	if (sum(parm$clust$n.vec) != parm$N)
		{err <- 4
		}

	if (length(parm$clust$n.vec) != parm$clust$K)
		{err <- 5
		}

	if (length(parm$clust$phi.v) != parm$clust$K)
		{err <- 6
		}

	err
	}


element_fn.swap.clusters <- function(parm, g1, g2)
	{

	####################################################
	# swap the clusters with labels g1 and g2
	####################################################

	ind1 <- parm$clust$s.v == g1
	ind2 <- parm$clust$s.v == g2
	parm$clust$s.v[ind1] <- g2
	parm$clust$s.v[ind2] <- g1

	buffer <- parm$clust$phi.v[g1]
	parm$clust$phi.v[g1] <- parm$clust$phi.v[g2]
	parm$clust$phi.v[g2] <- buffer

	buffer <- parm$clust$n.vec[g1]
      parm$clust$n.vec[g1] <- parm$clust$n.vec[g2]
      parm$clust$n.vec[g2] <- buffer

	#####################

	parm
	}




########################################################


element_fn.log.lik <- function(mean, sd, num, Y, X.sd)  
	{log.lik <- -num/2*log(2*pi*sd^2) - .5*num*(Y-mean)^2/sd^2 - .5*(num-1)*X.sd^2/sd^2  ## HOT (15-08-28)

	log.lik
	}




element_fn.gen.phi <- function(parm, Y, lik.sd) 
{
	if (is.null(Y))
		{lik.sd <- Inf
		Y <- 0 #dummy
		}


	post.sd  <- 1/sqrt(1/parm$clust$tau2^2 + 1/lik.sd^2)

	post.mean <- post.sd^2 * (parm$clust$mu2/parm$clust$tau2^2 + Y/lik.sd^2)

	new.phi <- rnorm(n=1, mean=post.mean, sd=post.sd)

	new.phi
}



element_fn.nbhd <- function(I, parm, max.row.nbhd.size)
	{if (length(I)>1)
		{k <- sample(I, size=1)
		}
	if (length(I)==1)
		{k <- I
		}

   H.v <- rep(Inf,parm$N)

	 post.prob.mt <- parm$clust$post.prob.mt

	 tmp1.mt <- matrix(post.prob.mt[,I], ncol=length(I)) ## HOT
	 tmp2.v <- post.prob.mt[,k]
	 tmp3.mt <- sqrt(tmp1.mt * tmp2.v)
	 H.v[I] <-  2*(1-colSums(tmp3.mt))

	 cutoff <- parm$row.delta
	 I.k <- which(H.v <= cutoff)

	if (length(I.k) > max.row.nbhd.size)
		{cutoff <- quantile(H.v[I.k], probs=max.row.nbhd.size/length(I.k))
	 	 I.k <- which(H.v <= cutoff)
		}

	 I.k <- sort(I.k)

       I <- sort(setdiff(I, I.k))
       I <- sort(I)

       list(k, I.k, I)

	}



fast_element_fn.nbhd <- function(relative_I, parm, max.row.nbhd.size)
{
  if (length(relative_I)>1)
{relative_k <- sample(relative_I, size=1)
}
if (length(relative_I)==1)
{relative_k <- relative_I
}

post.prob.mt <- parm$clust$subset_post.prob.mt

tmp1.mt <- matrix(post.prob.mt[,relative_I], ncol=length(relative_I))
tmp2.v <- post.prob.mt[,relative_k]
tmp3.mt <- sqrt(tmp1.mt * tmp2.v)
H.v <-  2*(1-colSums(tmp3.mt))

cutoff <- parm$row.delta
flag.v <- which(H.v <= cutoff)
relative_I.k <- relative_I[flag.v]

if (length(relative_I.k) > max.row.nbhd.size)
{#cutoff <- quantile(H.v[flag.v], probs=max.row.nbhd.size/length(relative_I.k))
 #relative_I.k <- relative_I[which(H.v <= cutoff)]
 relative_I.k <- relative_I[rank(H.v, ties="random") <= max.row.nbhd.size]
}

relative_I.k <- sort(relative_I.k)

relative_I <- sort(setdiff(relative_I, relative_I.k))
relative_I <- sort(relative_I)

list(relative_k, relative_I.k, relative_I)

}



fast_element_fn.post.prob.and.delta <- function(parm, max.row.nbhd.size, computeMode)

	{

  if (computeMode$computeR) {

	################################################
	### Compute pmf of cluster variables w_1,...,w_n
	###############################################

	prior.prob.v <- c(parm$clust$n.vec)
	small <- 1e-3 # compared to 1
	prior.prob.v[prior.prob.v < small] <- small

	#log.ss.mt <- array(,c(parm$clust$K, parm$N))
  subset_log.ss.mt <- array(,c(parm$clust$K, length(parm$row.subset.I)))

	#
	Y.v <- parm$Y[parm$row.subset.I]
	Z.sd.v <- parm$Z.sd[parm$row.subset.I]
	## g.v also equals parm$g[parm$row.subset.I]
	#g.v <- (parm$row.subset.I-1) %/% parm$n2 + 1
	num.k.v <-  parm$num.Y[parm$row.subset.I]
  sd.v<-parm$tau


	for (ss in 1:parm$clust$K)
		{subset_log.ss.mt[ss,] <- element_fn.log.lik(mean=parm$clust$phi.v[ss], sd=sd.v, num=num.k.v, Y=Y.v, X.sd=Z.sd.v)
    #log.ss.mt[ss,parm$row.subset.I] <- element_fn.log.lik(mean=parm$clust$phi.v[ss], sd=parm$tau, num=num.k.v, Y=Y.v, X.sd=X.sd.v)
		}

	subset_log.ss.mt <- subset_log.ss.mt + log(prior.prob.v)

	#maxx.v <- apply(log.ss.mt, 2, max)    # HOT
	subset_maxx.v <- apply(subset_log.ss.mt, 2, max)

	subset_log.ss.mt <- t(t(subset_log.ss.mt) - subset_maxx.v)
	subset_ss.mt <- exp(subset_log.ss.mt)

	subset_col.sums.v <- colSums(subset_ss.mt)
	subset_ss.mt <- t(t(subset_ss.mt)/subset_col.sums.v)

	# replace zeros by "small"
	small2 <- 1e-5
	subset_ss.mt[subset_ss.mt < small2] <- small2

	# again normalize
	subset_col.sums.v <- colSums(subset_ss.mt)
	subset_ss.mt <- t(t(subset_ss.mt)/subset_col.sums.v)

	parm$clust$post.prob.mt <- array(,c((parm$clust$K), parm$N))
	parm$clust$post.prob.mt[,parm$row.subset.I] <- subset_ss.mt

	dimnames(parm$clust$post.prob.mt) <- list(1:parm$clust$K, 1:parm$N)

  parm$clust$subset_post.prob.mt <- subset_ss.mt
	dimnames(parm$clust$subset_post.prob.mt) <- list(1:parm$clust$K, 1:length(parm$row.subset.I))

	#########################################
	### now compute the delta-neighborhoods
	#########################################

	  if (computeMode$computeC) { # debugging
      savedSeed <- .GlobalEnv$.Random.seed # For debugging purposed only
	  }

	# I <- parm$row.subset.I
  relative_I <- 1:length(parm$row.subset.I)

  first <- TRUE
	while (length(relative_I)>=1)
		{tmp <- fast_element_fn.nbhd(relative_I, parm, max.row.nbhd.size)
		relative_k <- tmp[[1]]
		relative_I.k <- tmp[[2]]
		relative_I <- tmp[[3]]
		#
		if (first) {
		  parm$clust$row.nbhd <- list(parm$row.subset.I[relative_I.k])
		  parm$clust$row.nbhd.k <- parm$row.subset.I[relative_k]
		  first <- FALSE
		} else {
		  parm$clust$row.nbhd <- c(parm$clust$row.nbhd, list(parm$row.subset.I[relative_I.k]))
		  parm$clust$row.nbhd.k <- c(parm$clust$row.nbhd.k, parm$row.subset.I[relative_k])
		}

		}
  }

  if (computeMode$computeC) {

	  if (computeMode$computeR) { # debugging
	    .GlobalEnv$.Random.seed <- savedSeed # Roll back PRNG
	  }

	test <- .computePmfAndNeighborhoods(computeMode$device$engine,
	                            parm$clust$n0, parm$clust$n.vec, 1e-3, 1e-5,
	                            parm$clust$K, parm$N,
	                            parm$Y, parm$X.sd, parm$row.subset.I,
	                            parm$clust$C.m.vec, parm$n2,
	                            parm$clust$phi.v, parm$tau, parm$tau_0,
	                            max.row.nbhd.size, parm$row.delta)

	if (computeMode$computeR) { # debugging
	  assertEqual(parm$clust$row.nbhd.k, test$index)
	  assertEqual(unlist(parm$clust$row.nbhd), test$neighbor)
	}

	# Convert from simple flat format to list of int vectors
	end <- test$offset[-1] - 1
	begin <- test$offset
	length(begin) <- length(begin) - 1

	parm$clust$row.nbhd.flat.length <- length(test$neighbor)
	parm$clust$row.nbhd <- lapply(1:length(begin),
	                              FUN = function(x) {
	                                test$neighbor[begin[x]:end[x]]
	                              })
	parm$clust$row.nbhd.k <- test$index

  } # computeMode

	## END

	parm
	}



###########################################################

element_fn.row.gibbs.DP <- function(parm, computeMode)
{
     k <- parm$k
	I.k <- parm$I.k
	s.k <- parm$clust$s.v[k]

	###############

	err <- element_fn.consistency.check(parm, computeMode)
	if (err > 0)
		{stop(paste("GIBBS STEP - 0: failed consistency check: err=",err))
		}

	# store just in case
	init.cc.parm <- parm

	###############

	#if (s.k > 0)
		{parm$clust$n.vec[s.k] <- parm$clust$n.vec[s.k] - 1
		}

# 	if (s.k == 0)
# 		{parm$clust$n0 <- parm$clust$n0 - 1
# 		}

	# sum(parm$clust$n.vec) + parm$clust$n0 == (parm$N-1)

	#######################################################
	### emptied clusters are gone forever under Gibbs sampling move
	#######################################################

	Y.k <- parm$Y[k]
	Z.sd.k <- parm$Z.sd[k]
	#g.k <- (k-1) %/% parm$n2 + 1
  num.k <-  parm$num.Y[k]
  sd.k<-parm$tau


	if (computeMode$computeR) {
	  L.v <- sapply(parm$clust$phi.v, element_fn.log.lik, sd = sd.k, num = num.k, Y = Y.k, X.sd = Z.sd.k)
	}

	if (computeMode$computeC) {
	  L.v.C <- .vectorizedElementFnLogLik(computeMode$device$engine, parm$clust$phi.v, parm$tau, num.k, Y.k, X.sd.k)

	  if (computeMode$computeR) { # debugging
	    assertEqual(L.v.C, L.v, computeMode$tolerance)
	  }

	  L.v <- L.v.C
	}

	# add possibility that s=0 in front
	#L.v <- c(element_fn.log.lik(mean=0, sd=parm$tau_0, num=num.k, Y=Y.k, X.sd=X.sd.k), L.v)

	new.phi <- element_fn.gen.phi(parm, Y=Y.k, lik.sd=sd.k/sqrt(num.k))  ##check later 12/22/2015 12:16pm
	# last element is for new cluster
	parm$clust$n.vec <- c(parm$clust$n.vec, 0)
	parm$clust$phi.v <- c(parm$clust$phi.v, new.phi)
	parm$clust$K <- parm$clust$K + 1

	new.L <- element_fn.log.lik(mean=new.phi, sd=sd.k, num=num.k, Y=Y.k, X.sd=Z.sd.k) + dnorm(new.phi, mean=parm$clust$mu2, sd=parm$clust$tau2, log=TRUE)
	L.v <- c(L.v, new.L)

	############

	# forcing singletons to leave empty clusters
	log.prior.v <- log(c(parm$clust$n.vec[-parm$clust$K], parm$clust$M))

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

	new.s.k <- sample(1:parm$clust$K, size=length(I.k), replace=TRUE, prob=parm$clust$post.k)

	parm$clust$s.v[k] <- new.s.k
	#if (new.s.k > 0)
		{parm$clust$n.vec[new.s.k] <- parm$clust$n.vec[new.s.k] + 1
		}
# 	if (new.s.k == 0)
# 		{parm$clust$n0 <- parm$clust$n0 + 1
# 		}

	# sum(parm$clust$n.vec) + parm$clust$n0 == parm$N

	err <- element_fn.consistency.check(parm, computeMode)
	if (err > 0)
		{stop(paste("GIBBS STEP - 1: failed consistency check: err=",err))
		}

	parm
}


###########################################################

element_fn.fast.DP.iter <- function(parm, computeMode)
{
     k <- parm$k
	I.k <- parm$I.k

	###############

	err <- element_fn.consistency.check(parm, computeMode)
	if (err > 0)
		{stop(paste("FORWARD STEP - 0: failed consistency check: err=",err))
		}

	# store so that we can revert to this state if MH propsal is rejected
	init.cc.parm <- parm

	###############

	subset.s <- parm$clust$s.v[I.k]

	if (computeMode$computeR) {

	  parm$clust$n.vec.k <- array(,parm$clust$K)
	  for (gg in 1:parm$clust$K) {
	    parm$clust$n.vec.k[gg] <- sum(subset.s==gg)
		}
	  #parm$clust$n0.k <- sum(subset.s==0)
	}

	if (computeMode$computeC) {

	  # If this gets hot again, could write a fast subsettedTabulate
	  all.n.vec <- .fastTabulateVector(subset.s, parm$clust$K, TRUE)

	  if (computeMode$computeR) { # debugging
	    assertEqual(parm$clust$n0.k, all.n.vec[1])
	    assertEqual(parm$clust$n.vec.k, all.n.vec[-1])
	  }

	  parm$clust$n0.k <- all.n.vec[1]
	  parm$clust$n.vec.k <- all.n.vec[-1]
	}

	parm$clust$n.vec.k.comp <- parm$clust$n.vec - parm$clust$n.vec.k
	#parm$clust$n0.k.comp <- parm$clust$n0 - parm$clust$n0.k

	#

	# sum(parm$clust$n.vec.k) + parm$clust$n0.k == length(I.k)
	# sum(parm$clust$n.vec.k.comp) + parm$clust$n0.k.comp == (parm$N-length(I.k))

	Y.k <- parm$Y[k]
	Z.sd.k <- parm$Z.sd[k]
	#g.k <- (k-1) %/% parm$n2 + 1
	num.k <-  parm$num.Y[k]
	sd.k<-parm$tau
	

	if (computeMode$computeR) {
	  L.v <- sapply(parm$clust$phi.v, element_fn.log.lik, sd = sd.k, num = num.k, Y = Y.k, X.sd = Z.sd.k)
	}

	if (computeMode$computeC) {
	  L.v.C <- .vectorizedElementFnLogLik(computeMode$device$engine, parm$clust$phi.v, parm$tau, num.k, Y.k, X.sd.k)

	  if (computeMode$computeR) { # debugging
	    assertEqual(L.v.C, L.v, computeMode$tolerance)
	  }

	  L.v <- L.v.C
	}

# 	  	  if (any(abs(L.v - L.v2) > 1e-10)) {
# 	  	    stop("C++ error")
# 	  	  }

	# add possibility that s=0 in front
	#L.v <- c(element_fn.log.lik(mean=0, sd=parm$tau_0, num=num.k, Y=Y.k, X.sd=X.sd.k), L.v)

	############

	log.prior.v <- log(c(parm$clust$n.vec.k.comp))
	#######################################################
	# Buffering empty clusters with positive masses using Neal's method
	# NECESSARY under MH algorithm otherwise reverse proposals are impossible
	# and proposal always accepted (bad)
	#######################################################

	newly.empty.indx <- which(parm$clust$n.vec.k.comp==0)

	# possibly overwrite some 0's
	if (length(newly.empty.indx) > 0)
		{log.prior.v[newly.empty.indx] <- log(parm$clust$M / length(newly.empty.indx))
		}

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

	new.s.k <- sample(1:parm$clust$K, size = length(I.k), replace = TRUE, prob = parm$clust$post.k)
	old.s.k <- old.parm$clust$s.v[I.k]

	exit <- (sum(new.s.k != old.s.k) == 0)
	flip <- 1

  if (!exit) # CONTINUE W/O EXITING FUNCTION
  {

	if (computeMode$computeR) {

	  parm$clust$s.v[I.k] <- new.s.k

	  rho.prop <- 0
	  parm$clust$n.vec <- array(0, parm$clust$K)

	  for (gg in 1:parm$clust$K) {
	    flag.gg <- (new.s.k == gg)
	    count.gg <- sum(flag.gg)

	    #if (gg > 0) {
	      parm$clust$n.vec[gg] <- parm$clust$n.vec.k.comp[gg] + count.gg
	   # }

# 	    if (gg == 0) {
# 	      parm$clust$n0 <- parm$clust$n0.k.comp + count.gg
# 	    }

	    if (count.gg > 0) {
	      rho.prop <- rho.prop + log(parm$clust$post.k[gg]) * count.gg
	    }
	  }

	}

  if (computeMode$computeC) {

	  .fastIndexedSetNoCopy(parm$clust$s.v, I.k, new.s.k) # Same as: parm$clust$s.v[I.k] <- new.s.k
	  all.count <- .fastTabulateVector(new.s.k, parm$clust$K, TRUE)
	  rho.prop.C <- .fastSumSafeLog(parm$clust$post.k, all.count,
	                                length(parm$clust$post.k))

	  if (computeMode$computeR) { # debugging
	    assertEqual(parm$clust$n0.k.comp + all.count[1], parm$clust$n0)
	    assertEqual(parm$clust$n.vec.k.com + all.count[-1], parm$clust$n.vec)
	    assertEqual(rho.prop.C, rho.prop, computeMode$tolerance)
	  }

	  parm$clust$n0 <- parm$clust$n0.k.comp + all.count[1]
	  parm$clust$n.vec <- parm$clust$n.vec.k.com + all.count[-1]
    rho.prop <- rho.prop.C
	}
	# sum(parm$clust$n.vec) + parm$clust$n0 == parm$N


	##########################################
	##########################################
	##### Computing proposal prob of reverse move
	##########################################
	##########################################

  initial.rho.prop <- rho.prop

  if (computeMode$computeR) {

    for (gg in 1:parm$clust$K) {
      flag.gg <- (old.s.k == gg)
      count.gg <- sum(flag.gg)

      if (count.gg > 0) {
        rho.prop <- rho.prop - log(old.parm$clust$post.k[gg]) * count.gg
      }
    }

  }

  if (computeMode$computeC) {
    all.count <- .fastTabulateVector(old.s.k, parm$clust$K, TRUE)
    rho.prop.C <- initial.rho.prop - .fastSumSafeLog(old.parm$clust$post.k, all.count,
                                        length(old.parm$clust$post.k))

    if (computeMode$computeR) { # debugging
      assertEqual(rho.prop.C, rho.prop, computeMode$tolerance)
    }

    rho.prop <- rho.prop.C
  }

	#######################################################
	#######################################################
	########## computing true log-ratio
	#######################################################
	#######################################################

	new.indx <- which(parm$clust$n.vec > 0)
	old.indx <- which(old.parm$clust$n.vec > 0)

	new.K <- length(new.indx)
	old.K <- length(old.indx)

# 	new.is.s0 <- parm$clust$n0 > 0
# 	old.is.s0 <- old.parm$clust$n0 > 0

	rho.tru <- sum(lgamma(parm$clust$n.vec[new.indx])) - sum(lgamma(old.parm$clust$n.vec[old.indx])) + log(parm$clust$M) * (new.K - old.K)
	#rho.tru <- rho.tru + sum(lgamma(parm$clust$M0 + parm$clust$n0)) - sum(lgamma(parm$clust$M0 + old.parm$clust$n0)) + log(parm$clust$M0) * (new.is.s0 - old.is.s0)

	###################

	if (computeMode$computeR) { # START

	  Y.k.v <- parm$Y[I.k]
	  Z.sd.k.v <- parm$Z.sd[I.k]
	  ## g.k.v also equals parm$g[I.k]
	  #g.k.v <- (I.k - 1) %/% parm$n2 + 1
	  #num.k.v <-  parm$clust$C.m.vec[g.k.v]
	  num.k.v <-  parm$num.Y[I.k]
	  sd.k.v<-parm$tau
	  

	  new.log.lik <- old.log.lik <- 0

# 	  new.s.pos.indx <- new.s.k[new.s.k > 0]
# 	  old.s.pos.indx <- old.s.k[old.s.k > 0]

	  #if (length(new.s.pos.indx) > 0) ### TODO Can vectorize element_fn.log.lik below
	  #{
      new.log.lik <- new.log.lik + sum(element_fn.log.lik(mean = parm$clust$phi.v[new.s.k], sd = sd.k.v, num = num.k.v, Y = Y.k.v, X.sd = Z.sd.k.v))
	  #}
	  #if (length(old.s.pos.indx) > 0)
	  #{
      old.log.lik <- old.log.lik + sum(element_fn.log.lik(mean = old.parm$clust$phi.v[old.s.k], sd = sd.k.v, num = num.k.v, Y = Y.k.v, X.sd = Z.sd.k.v))
	  #}

# 	  if (sum(new.s.k == 0) > 0)
# 	  {new.log.lik <- new.log.lik + sum(element_fn.log.lik(mean = 0, sd = parm$tau_0, num = num.k.v[new.s.k == 0], Y = Y.k.v[new.s.k == 0], X.sd = X.sd.k.v[new.s.k == 0]))
# 	  }
# 	  if (sum(old.s.k == 0) > 0)
# 	  {old.log.lik <- old.log.lik + sum(element_fn.log.lik(mean = 0, sd = parm$tau_0, num = num.k.v[old.s.k == 0], Y = Y.k.v[old.s.k == 0], X.sd = X.sd.k.v[old.s.k == 0]))
# 	  }

	}

	if (computeMode$computeC) {

	  result <- .computeDPAcceptanceRatio(computeMode$device$engine,
	                                    parm$Y, parm$X.sd, I.k, parm$clust$C.m.vec, parm$clust$phi.v,
	                                    new.s.k, old.s.k,
	                                    parm$tau, parm$tau_0, parm$n2)

	  if (computeMode$computeR) { # debugging
	    assertEqual(result$new, new.log.lik, computeMode$tolerance)
	    assertEqual(result$old, old.log.lik, computeMode$tolerance)
	  }

	  new.log.lik <- result$new
	  old.log.lik <- result$old

	} # END

	rho.tru <- rho.tru + new.log.lik - old.log.lik

	########## toss a coin #################

	prob <- exp(min((rho.tru - rho.prop),0))
	flip <- as.logical(rbinom(n = 1, size = 1, prob = prob))
	if (!flip) {
	  parm <- init.cc.parm

	  if (computeMode$computeC) {
	    # C++ version does not make a deep copy (via copy-on-write), so values must be manually restored
	    .fastIndexedSetNoCopy(parm$clust$s.v, I.k, old.s.k)
	  }

  }

 } # end BIG if (!exit) loop


	list(parm, flip)
}

########################################

element_fn.gen.new <- function(parm, num.gen)
	{
	parm$clust$K.old <- parm$clust$K

	for (tt in 1:num.gen)
		{
		# generate from prior
		new.phi <- element_fn.gen.phi(parm, Y=NULL, lik.sd=NULL)
		#
		parm$clust$n.vec <- c(parm$clust$n.vec, 0)
		parm$clust$phi.v <- c(parm$clust$phi.v, new.phi)
		parm$clust$K <- parm$clust$K + 1
		#
		}

	parm
	}

#####################################


element_fn.drop <- function(parm)
 	{

	##########################################
	## Drop empty clusters:
	## (i)  Move empty clusters to end by relabeling clusters
	## (ii) Set parm$clust$K equal to number of non-empty clusters
	## (ii) Retain only clusters  1,...,parm$clust$K
	#########################################

	parm$clust$K <- sum(parm$clust$n.vec>0)
	num.dropped <- sum(parm$clust$n.vec==0)

	if (num.dropped > 0)
	{
	for (rr in 1:num.dropped)
		{
		old.label <-  min(which(parm$clust$n.vec==0))
		new.label <- max(which(parm$clust$n.vec>0))
		stopp <-  max(which(parm$clust$n.vec>0)) == parm$clust$K
		if (stopp)
			{break
			}
		parm <- element_fn.swap.clusters(parm, g1=new.label, g2=old.label)
 		}
	}

	##########

	keep <- 1:parm$clust$K

	parm$clust$phi.v <- parm$clust$phi.v[keep]
	parm$clust$n.vec <- parm$clust$n.vec[keep]
#	parm$clust$s.mt <- matrix(parm$clust$s.v, nrow = parm$tot.trt)
	tmp.s.mt<-NULL
	for (ss in 1:length(parm$clust$s.v)) {
	  if (parm$clust$s.stretch[ss]==1) {
	    tmp.s.mt<-c(tmp.s.mt,rep(parm$clust$s.v[ss],parm$tot.trt))
	  }
	  if (parm$clust$s.stretch[ss]==0) {
	    tmp.s.mt<-c(tmp.s.mt,parm$clust$s.v[ss])
	  }
	}
	if (length(tmp.s.mt)!=parm$tot.trt*parm$clust$G) {stop(paste('problem in elementwise_main-fn1.update'))}
	parm$clust$s.mt <- array(tmp.s.mt, c(parm$tot.trt, parm$clust$G))
	

	parm
	}



element_fn.fast.DP <- function(parm, max.row.nbhd.size, row.frac.probes, computeMode)
{

	if (row.frac.probes < 1)
		{parm$row.subset.I <- sort(sample(1:parm$N, ceiling(row.frac.probes*parm$N)))
		}

	if (row.frac.probes == 1)
		{parm$row.subset.I <- 1:parm$N
		}

	##########################
	# compute delta-neighborhoods
	#########################
	parm <- fast_element_fn.post.prob.and.delta(parm, max.row.nbhd.size, computeMode)

	err <- element_fn.consistency.check(parm, computeMode)
	if (err > 0)
		{stop(paste("ELEMENT_DP MAIN: failed consistency check: err=",err))
		}

	#################################
	# generate new empty clusters from prior
	#################################

	parm <- element_fn.gen.new(parm, num.gen=5)

	##########################
	# M-H move for random effects of delta-nieghborhoods
	#########################

	num.nbhds <- length(parm$clust$row.nbhd)
	flip.v <- NULL

      for (cc in 1:num.nbhds)
		{parm$k <- parm$clust$row.nbhd.k[cc]
		parm$I.k <- parm$clust$row.nbhd[[cc]]

		 # parm$I.k will contain at least parm$k

		if (length(parm$I.k) > 1)
			{tmp <- element_fn.fast.DP.iter(parm, computeMode)
		 	parm <- tmp[[1]]
		 	flip.v <- c(flip.v, tmp[[2]])
			}

		if (length(parm$I.k) == 1)
			{parm <- element_fn.row.gibbs.DP(parm, computeMode)
			}
		 #
		}


	######################################
if (length(flip.v)>0){
	parm$clust$row.flip <- mean(flip.v)
}

	##########################################
	## Drop empty clusters:
	## (i)  Move empty clusters to end by relabeling clusters
	## (ii) Set parm$clust$K equal to number of non-empty clusters
	## (ii) Retain only clusters  1,...,parm$clust$K
	#########################################

	parm <- element_fn.drop(parm)

	#################################
	# uprow.data parm$phi.v conditional on clusters
	#################################


	### implements data reduction formula on page 1 of 11/14/11 notes

	#indxx <- parm$clust$s.v>0
	num.v <- parm$num.Y
	Y_num.v <- parm$Y * num.v
	#s.pos <- parm$clust$s.v[indxx]
	Y_num.j.v <- tapply(Y_num.v, parm$clust$s.v, sum)
	num.j.v <- tapply(num.v,parm$clust$s.v, sum)
	Y.j.v <- Y_num.j.v/num.j.v


	for (j in 1:parm$clust$K)
		{parm$clust$phi.v[j] <- element_fn.gen.phi(parm, Y=Y.j.v[j], lik.sd=parm$tau/sqrt(num.j.v[j]))
		}

	err <- element_fn.compact_consistency.check(parm, computeMode)
	if (err > 0)
		{stop(paste("ELEMENT_DP MAIN END: failed consistency check: err=",err))
		}

	parm
	}



