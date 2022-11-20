createComputeMode <- function(language = "R",
                              exactBitStream = FALSE,
                              extraSort = TRUE,
                              completeTest = FALSE,
                              tolerance = 1E-10,
                              test1 = FALSE,
                              test2 = FALSE,
                              test3 = FALSE) {
  if (!(language %in% c("C","R"))) {
    stop("Invalid language")
  }
  
  useR <- (language == "R")
  device <- NULL
  if (!useR) {
    doSort <- (exactBitStream | extraSort)
    device <- .createEngine(doSort)
  }
  
  object <- list(
    computeR = (language == "R" | completeTest),
    computeC = (language == "C"),
    device = device,
    exactBitStream = exactBitStream,
    extraSort = extraSort,
    tolerance = tolerance,
    test1 = test1,
    test2 = test2,
    test3 = test3
  )
  class(object) <- "computeMode"
  return(object)
}

fn.dmvnorm <- function(x, mean, sigma, inv.sigma, log=TRUE)
{
  
  # Computes multivariate normal density function
  # a little faster than dmvnorm function of R!
  
  if (missing(inv.sigma))
  {inv.sigma <- solve(sigma)
  }
  
  logdet <- as.numeric(determinant(inv.sigma, logarithm=TRUE)$mod)
  r <- length(x)
  Q <- colSums(inv.sigma * (x-mean))
  Q <- sum(Q * (x-mean))
  
  val <- -r/2*log(2*pi) + logdet/2 - Q/2
  if (!log)
  {val <- exp(val)
  }
  
  val
}



fn.quality.check <- function(parm)
{err <- 0

if (!is.null(parm$clust$col.nbhd))
{#if (sum(unlist(lapply(parm$clust$col.nbhd, length))) != length(parm$col.subset.I))
  #{err <- 1
  #}
}

if (ncol(parm$clust$A.mt) != parm$clust$G)
{err <- 3
}

if ((sum(parm$clust$C.m.mt)) != parm$p)
{err <- 5
}

if ((sum(parm$clust$C.m.vec)) != parm$p)
{err <- 5.5
}

if (length(parm$clust$n.vec) != parm$clust$K)
{err <- 6
}

if (length(parm$clust$phi.v) != parm$clust$K)
{err <- 7
}

if ((sum(parm$clust$n.vec)) != parm$N)
{err <- 8
}

# if ((sum(parm$clust$n.vec)) != parm$tot.trt*parm$clust$G)
# {err <- 8
# }

err
}

######################

fn.init.clusters <- function(parm)
{
  
  X.mt <- parm$Z
  num.centers <- parm$G.new
  
  options(warn=0)
  tmp2 <- kmeans(t(X.mt), iter.max=1000, centers=num.centers, nstart=2)
  options(warn=2)
  
  parm$clust$c.v <- tmp2$cluster
  parm$clust$G <- length(tmp2$size)
  
  # start from PDP model with d=0 (i.e DP) 
  parm$d <- 0
  
  #introduce correlation range parameter eta
  parm$eta<-parm$eta_R <- .01
  
  parm
}


fn.eda <- function(parm, data, computeMode)
{
  parm$epsilon<-rowMeans(parm$X.LANE)-mean(parm$X.LANE)
  parm$Z.f<-pre.Z<-parm$X.LANE-matrix(rep(parm$epsilon,parm$p),ncol=parm$p)
  tmp.Z<-matrix(,parm$tot.trt,parm$p)
  for (ii in 1:parm$tot.trt){
    tmp.Z[ii,]<-colMeans(pre.Z[parm$group==ii,,drop=FALSE])
  }
  parm$Z<-tmp.Z
  
  parm <- fn.init.clusters(parm)
  parm$G.max <- min(parm$p/2, round(parm$clust$G*1.1))
  
  parm$Y.de_naive.eda<-array(,c(parm$n2,parm$clust$G))
  parm$Y.eda <- parm$clust$A.mt <- array(,c(parm$tot.trt,parm$clust$G))
  parm$clust$C.m.mt<- matrix(NA,2,parm$clust$G)
  
  ##compute parm$Y
  for (g in 1:parm$clust$G)
  {I.g <- (parm$clust$c.v==g)
  m.g <- sum(I.g)
  x.g.v <- parm$Z[,I.g]
  if (m.g > 1)
  {x.g.v <- rowMeans(x.g.v)
  }
  parm$Y.eda[,g] <- x.g.v
  }
  
  for (gd in 1:parm$clust$G)
  {I.gd <- (parm$clust$c.v==gd)
  m.gd <- sum(I.gd)
  x.gd.v <- parm$Z.f[,I.gd]
  if (m.gd > 1)
  {x.gd.v <- rowMeans(x.gd.v)
  }
  parm$Y.de_naive.eda[,gd] <- x.gd.v
  }
  
  if (sum(table(parm$group)==1)>0) {
    tmp.diff<-apply(parm$Y.de_naive.eda,2,max)-apply(parm$Y.de_naive.eda,2,min)
    eda.de.idx<-tmp.diff>quantile(tmp.diff,probs=0.9)
  }
  else {
    tmp.p_val.v<-NULL
    for (jj in 1:ncol(parm$Y.de_naive.eda)) {
      naive_fit<-lm(parm$Y.de_naive.eda[,jj]~factor(parm$group))
      tmp.p_val.v<-c(tmp.p_val.v,anova(naive_fit)$"Pr(>F)"[1])
    }
    
    eda.de.idx<-tmp.p_val.v<0.05
  }
  
  parm$clust$gr.c.v<-rep(1,parm$clust$G)
  parm$clust$gr.c.v[eda.de.idx]<-2
  
  parm$clust$gr.v<-parm$clust$gr.c.v[parm$clust$c.v]
  
  parm$clust$g.v<-parm$clust$gr.v
  
  parm$clust$M <- parm$a.R
  
  parm$clust$mu2 <- mean(as.vector(parm$Z))
  parm$clust$tau2 <- diff(range(as.vector(parm$Z)))/6
  
  #################################
  
  parm$g <- rep(1:parm$clust$G,each=parm$tot.trt)
  parm$trt.exp<-rep(1:parm$tot.trt,parm$clust$G)
  parm$N <- parm$tot.trt * sum(parm$clust$gr.c.v==2) + sum(parm$clust$gr.c.v==1)
  
  tmp.Y_g1<-colMeans(parm$Y.eda[,!eda.de.idx,drop=F])
  tmp.Y_g2<-as.vector(parm$Y.eda[,eda.de.idx])
  parm$Y_clipped<-c(tmp.Y_g2,tmp.Y_g1)
  parm$clust$K <- min(length(parm$Y_clipped),data$K.max)
  
  tmp <- tryCatch({
    iter.max <- ifelse((length(parm$Y_clipped) > 1000), 10, 1000)
    kmeans(parm$Y_clipped, iter.max = iter.max, centers = parm$clust$K, nstart = 10,
           algorithm = "Hartigan-Wong" 
    )}, error = function(e) {
      print("Kmeans did not converge ... using random assignment")
      tmp.choose_K<-sample(1:length(parm$Y_clipped),size=parm$clust$K,replace = F)
      cluster <- array(NA,length(parm$Y_clipped))
      cluster[tmp.choose_K]<-1:parm$clust$K
      cluster[-tmp.choose_K]<-sample(1:parm$clust$K, size=length(parm$Y_clipped)-parm$clust$K, replace = TRUE)
      centers <- sapply(1:parm$clust$K, FUN = function(x) {
        mean(parm$Y_clipped[which(cluster == x)])
      })
      size <- sapply(1:parm$clust$K, FUN = function(x) {
        sum(cluster == x)
      })
      list(cluster = cluster, centers = centers, size = size)
    })
  
  
  tmp.clust.s.v <- tmp$cluster
  parm$clust$s.mt <- array(NA, c(parm$tot.trt,parm$clust$G))
  parm$clust$s.mt[,eda.de.idx]<-matrix(tmp.clust.s.v[1:length(tmp.Y_g2)],nrow=parm$tot.trt)
  parm$clust$s.mt[,!eda.de.idx]<-matrix(rep(tmp.clust.s.v[-(1:length(tmp.Y_g2))],each=parm$tot.trt),nrow=parm$tot.trt)
  
  parm$clust$phi.v <- as.vector(tmp$centers)
  
  parm$clust$n.vec <- tmp$size
  
  for (g in 1:parm$clust$G)
  {parm$clust$A.mt[,g] <- parm$clust$phi.v[parm$clust$s.mt[,g]]
  }
  tmp.s.v_g1<-parm$clust$s.mt[1,parm$clust$gr.c.v==1]
  g1_idx<-which(parm$clust$gr.c.v==1)
  if (length(unique(tmp.s.v_g1))<length(tmp.s.v_g1)) {
    for (mm in 1:length(unique(tmp.s.v_g1))) {
      tmp.same.idx_g1<-which(tmp.s.v_g1==unique(tmp.s.v_g1)[mm])
      if (length(tmp.same.idx_g1)>1) {
        tmp.merge.idx<-g1_idx[tmp.same.idx_g1[-1]]
        tmp.keep.idx<-g1_idx[tmp.same.idx_g1[1]]
        parm$clust$c.v<-replace(parm$clust$c.v,parm$clust$c.v %in% tmp.merge.idx,tmp.keep.idx)
        parm$clust$n.vec[unique(tmp.s.v_g1)[mm]]<-parm$clust$n.vec[unique(tmp.s.v_g1)[mm]]-(length(tmp.same.idx_g1)-1)
      }
    }
  }
  
  for (g in 1:parm$clust$G)
  {I.g <- which(parm$clust$c.v==g)
  parm$clust$C.m.mt[1,g]<-sum(parm$clust$g.v[I.g]==1)
  parm$clust$C.m.mt[2,g]<-sum(parm$clust$g.v[I.g]==2)
  }
  
  parm$clust$C.m.vec<-colSums(parm$clust$C.m.mt)
  
  parm <- PDP_fn.drop(parm, computeMode)
  
  
  Z.g.mt <- parm$Z.f
  parm$theta.mt<-a.g.v <- parm$clust$A.mt[,parm$clust$c.v]
  resid.g.mt <- Z.g.mt - a.g.v[parm$group,]
  sum.resid.sq <- sum(resid.g.mt^2)
  parm$tau<-sqrt(sum.resid.sq/parm$n2/parm$p)
  
  
  parm
  
}



fn.gen.clust <- function(parm, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, computeMode)
{
  
  parm$G.new <- data$G.max
  parm <- fn.eda(parm, data, computeMode)
  
  parm
  
}

fn.init <- function(true, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, tBB_flag, standardize.X, flip.sign, computeMode = "R")
{
  
  parm <- NULL
  
  parm$n2 <- dim(data$X)[1] # TODO Check
  parm$p <- dim(data$X)[2]  # TODO Check
  
  # mass parameter of elementwise(s) groups
  # stored later in parm$clust$M
  parm$a.R <- true$M
  
  # mass paramater of columns
  parm$b1 <- true$b1
  
  # mass parameter of epsilon
  parm$a.R.e <- true$Me
  
  ############################
  # For delta neighborhoods
  ############################
  
  parm$col.delta <- .1
  
  # delta-neighborhood threshold for elements
  parm$row.delta <- .1
  
  #########################################
  # generating the R- and C- clusters
  ########################################
  
  parm$shift <- true$shift
  
  parm$mult<-true$mult
  
  parm$iter_num<-0
  
  parm$X <- data$X
  parm$tot.trt<-data$tot.trt
  parm$distance<-data$distance
  parm$rho<-c(1-data$rho,data$rho) #% of DE probes/genes
  parm$delta_margin<-data$delta_margin
  
  parm$LANE<-data$LANE
  parm$X.LANE<-parm$X-matrix(rep(parm$LANE,parm$p),ncol=parm$p)
  parm$group<-data$group
  parm$ni<-rep(NA,parm$tot.trt)
  for (jj in 1:parm$tot.trt){
    parm$ni[jj]<-sum(parm$group==jj)
  }
  
  parm <- fn.assign.priors(parm, data)
  
  parm <- fn.gen.clust(parm, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, computeMode)
  
  parm <- fn.element.DP(data, parm, max.row.nbhd.size, row.frac.probes=1, computeMode)
  
  tmp.s.v_g1_test<-parm$clust$s.mt[1,parm$clust$gr.c.v==1]
  if (length(unique(tmp.s.v_g1_test))<length(tmp.s.v_g1_test)) {
    stop(paste('init problem,element.DP no OK'))
  }
  
  parm <- fn.epsilon.init(parm, data)
  
  parm <- fn.hyperparameters(data, parm)
  
  parm
  
}

fn.epsilon.init<-function(parm, data) 
{
  parm$clust$Me <- parm$a.R.e  
  
  parm$clust$K_e <- round((parm$n2)/2)
  
  tmp <- kmeans(rowMeans(parm$X.LANE-parm$theta.mt[parm$group,]), centers=parm$clust$K_e)
  parm$clust$s_e.v <- tmp$cluster
  parm$clust$phi_e.v <- as.vector(tmp$centers)
  parm$clust$n_e.vec <- tmp$size
  
  parm$clust$mu_e <- mean(parm$clust$phi_e.v)

  parm$clust$tau_e <- sd(parm$clust$phi_e.v)

  parm$clust$theta_e.v <- parm$clust$phi_e.v[parm$clust$s_e.v]
  
  parm$epsilon<-parm$clust$theta_e.v
  
  parm
}



fn.gen.missing.X <- function(data, parm)
{
  # impute missing X values

  X.mt <- data$X

  if (parm$num.X.miss > 0)
  {
    for (cc in 1:parm$num.X.miss)
    {i.cc <- parm$X.missing.x[cc]
    j.cc <- parm$X.missing.y[cc]
    c.cc <- parm$clust$c.v[j.cc]
    if (c.cc != 0)
    {mean.cc <- parm$clust$A.mt[i.cc, c.cc]
    }
    if (c.cc == 0)
    {mean.cc <- 1
    }
    X.mt[i.cc, j.cc] <- rnorm(n=1, mean=mean.cc, sd=parm$tau)
    }
  }
  parm$X <- X.mt

  parm
}


fn.standardize_orient.X <- function(parm)
{

  ####
  ## STANDARDIZE X columns to unit variance and zero mean
  #####
  # Do only for columns with NA's
  # For other columns, it's just a one-time calculation at the beginning of MCMC

  if (parm$num.X.miss > 0)
    {tmp.X <- matrix(parm$X[,parm$X.missing.y],col=parm$num.X.miss)
    mean.v <- colMeans(tmp.X)
    sd.v <- apply(tmp.X, 2, sd)
    parm$X[,parm$X.missing.y] <- t((t(tmp.X) - mean.v)/sd.v)
   }

  ####
  ## ORIENT X
  ####

  parm$X <- t(t(parm$X) * parm$clust$orient.v)

  parm
}


fn.assign.priors <- function(parm, data)
{
  ######for tau###########
  parm$prior$tau <- NULL
  parm$prior$tau$alpha.tau <- 1/2
  parm$prior$tau$beta.tau <- 1/2
  
  #parm$prior$tau$max <- sqrt(.75)*sd(as.vector(data$X), na.rm=TRUE)
  parm$prior$tau$max <- 2*sd(as.vector(data$X), na.rm=TRUE)
  parm$prior$tau$min <- 1e-10
  parm$prior$tau.sq$max <- parm$prior$tau$max^2
  parm$prior$tau.sq$min <- parm$prior$tau$min^2
  parm$prior$inv.tau.sq$max <- 1/parm$prior$tau.sq$min
  parm$prior$inv.tau.sq$min <- 1/parm$prior$tau.sq$max
  
  ######for parameters related to thetas#########
  ###for tau2
  parm$prior$tau2 <- NULL
  parm$prior$tau2$alpha <- 1/2
  parm$prior$tau2$beta <- 1/2
  
  ###for mu2
  parm$prior$mu2 <- NULL
  parm$prior$mu2$mean <- 0
  parm$prior$mu2$sd <- 10*sd(as.vector(data$X))
  
  ######for parameters related to epsilons########
  ###for tau_e
  parm$prior$tau_e <- NULL
  parm$prior$tau_e$alpha <- 1/10
  parm$prior$tau_e$beta <- 1/10
  
  ###for mu_e
  parm$prior$mu_e <- NULL
  parm$prior$mu_e$mean <- 0
  parm$prior$mu_e$sd <- 200*sd(rowMeans(data$X))
  
  parm
}



########################################

fn.gen.tau  <- function(data, parm)
{
  ###################
  # update tau
  ###################
  
  # only covariates assigned to non-zero row and non-zero group clusters matter
  
  parm$Z.f<-parm$X.LANE-matrix(rep(parm$epsilon,parm$p),ncol=parm$p)
  
  sum.resid.sq <- 0
  count <- 0
  
  for (g in 1:parm$clust$G)
  {flag.v <- parm$clust$c.v == g
  #z.g.v <- parm$clust$s.mt[,g] > 0
  
  if ((sum(flag.v)>0))
  {X.g.mt <- parm$Z.f[,flag.v]
  a.g.v <- parm$clust$A.mt[parm$group,g]
  resid.g.mt <- X.g.mt - a.g.v
  sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
  #count <- count + sum(z.g.v)*sum(flag.v)
  }
  
  }
  
  shape <- parm$prior$tau$alpha + parm$n2*parm$p/2
  rate <- parm$prior$tau$beta + sum.resid.sq/2
  
  u.min <- pgamma(parm$prior$inv.tau.sq$min,shape=shape, rate=rate)
  u.max <- pgamma(parm$prior$inv.tau.sq$max,shape=shape, rate=rate)
  gen.u <- runif(n=1, min=u.min, max=u.max)
  
  parm$tau <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))
  
  # overwrite to avoid zeros and Inf
  if (round(u.min, digits = 5) == 1) # really close to 1
  {parm$tau<- 1/sqrt(parm$prior$inv.tau.sq$min)
  }
  if (round(u.max, digits = 5) == 0) # really close to 0
  {parm$tau<- 1/sqrt(parm$prior$inv.tau.sq$max)
  }
  
  parm
}


fn.gen.tau_0  <- function(data, parm)
{
  ###################
  # update tau_0
  ###################
  
  sum.resid.sq <- 0
  count <- 0
  
  for (g in 1:parm$clust$G)
  {flag.v <- parm$clust$c.v == g
  z.g.v <- parm$clust$s.mt[,g] > 0
  
  if ((sum(1-z.g.v) > 0) & (sum(flag.v)>0))
  {X.g.mt <- parm$X[!z.g.v,flag.v]
  resid.g.mt <- X.g.mt
  sum.resid.sq <- sum.resid.sq + sum(resid.g.mt^2)
  count <- count + sum(1-z.g.v)*sum(flag.v)
  }
  }
  
  shape <- 1 + count/2
  rate <- 1 + sum.resid.sq/2
  
  # minimum possible value of parm$tau_0 = 1.5 * maximum possible value of parm$tau
  # maximum possible value of parm$tau_0 = 3 * sd(as.vector(data$X))
  u.min <- pgamma(1/9 / var(as.vector(data$X),na.rm=TRUE),shape=shape, rate=rate)
  u.max <- pgamma(1/1.5^2/parm$prior$tau.sq$min,shape=shape, rate=rate)
  gen.u <- runif(n=1, min=u.min, max=u.max)
  
  parm$tau_0 <- 1/sqrt(qgamma(gen.u,shape=shape, rate=rate))
  
  parm
}


fn.hyperparameters <- function(data, parm)
	{

	# also updates update tau_int
	parm <- fn.gen.tau(data, parm)
  
	parm <- fn.mu2tau2(data, parm)
	parm <- fn.muetaue(data, parm)

	#parm <- fn.gen.tau_0(data, parm)

	parm

	}

fn.mu2tau2<-function(data, parm)
{
  #sample parm$clust$mu2
  post.mu2.prec<-parm$clust$K/parm$clust$tau2^2+1/parm$prior$mu2$sd^2
  post.mu2.sd<-1/sqrt(post.mu2.prec)
  post.mu2.mean<-post.mu2.sd^2*(mean(parm$clust$phi.v)*parm$clust$K/parm$clust$tau2^2+parm$prior$mu2$mean/parm$prior$mu2$sd^2)
  
  parm$clust$mu2<-rnorm(n=1,post.mu2.mean,post.mu2.sd)
  
  #sample parm$clust$tau2
  tau2.shape<-parm$prior$tau2$alpha+parm$clust$K/2
  tau2.rate<-parm$prior$tau2$beta+sum((parm$clust$phi.v-parm$clust$mu2)^2)/2
  
  parm$clust$tau2<-sqrt(1/rgamma(n=1,shape=tau2.shape,rate=tau2.rate))
  
  parm
}

fn.muetaue<-function(data, parm)
{
  #sample parm$clust$mu_e
  post.mu_e.prec<-parm$clust$K_e/parm$clust$tau_e^2+1/parm$prior$mu_e$sd^2
  post.mu_e.sd<-1/sqrt(post.mu_e.prec)
  post.mu_e.mean<-post.mu_e.sd^2*(mean(parm$clust$phi_e.v)*parm$clust$K_e/parm$clust$tau_e^2+parm$prior$mu_e$mean/parm$prior$mu_e$sd^2)
  
  parm$clust$mu_e<-rnorm(n=1,post.mu_e.mean,post.mu_e.sd)
  
  #sample parm$clust$tau_e
  tau_e.shape<-parm$prior$tau_e$alpha+parm$clust$K_e/2
  tau_e.rate<-parm$prior$tau_e$beta+sum((parm$clust$phi_e.v-parm$clust$mu_e)^2)/2
  
  parm$clust$tau_e<-sqrt(1/rgamma(n=1,shape=tau_e.shape,rate=tau_e.rate))
  
  parm
}

fn.funky <- function(s,t)
{# on log scale
  lgamma(s+t) - lgamma(s)
}

fn.d <- function(d, parm) 
{
  tmp0<-(parm$clust$n.vec)
  tmp1<-matrix(log(tmp0[parm$clust$s.mt]),ncol=ncol(parm$clust$s.mt))
  tmpF<-exp(colSums(tmp1)-parm$tot.trt*log(sum(tmp0)))
  tmp.s.mt_g1<-parm$clust$s.mt[,parm$clust$gr.c.v==1,drop=F]
  if (sum(apply(tmp.s.mt_g1,2,max)-apply(tmp.s.mt_g1,2,min))>0) {stop(paste('problem in iterations-fn.d'))}
  tmp.s.mt_g1_comp<-apply(tmp.s.mt_g1,2,unique)
  tmp1_g1<-tmp0[tmp.s.mt_g1_comp]
  tmpF[parm$clust$gr.c.v==1]<-exp(tmp1_g1-log(sum(tmp0)))
  # formula in review paper by Lijoi and Prunster - Models Beyond - p.23, eq.36
  tmp.n<-rowSums(parm$clust$C.m.mt)
  nempty.idx<-parm$clust$C.m.mt>0
  tmp.G.gr<-rowSums(nempty.idx)
  #	log.lik_1order <- sum(log(parm$b1 + (1:(tmp.G.gr[1]-1))*d)) - fn.funky((parm$b1+1), (tmp.n[1]-1)) + sum(fn.funky((1-d+tmpF[nempty.idx[1,]]*parm$b1), (parm$clust$C.m.mt[1,nempty.idx[1,]]-1)))
  #	log.lik_1order <- log.lik_1order+sum(log(parm$b1 + (1:(tmp.G.gr[2]-1))*d)) - fn.funky((parm$b1+1), (tmp.n[2]-1)) + sum(fn.funky((1-d+tmpF[nempty.idx[2,]]*parm$b1), (parm$clust$C.m.mt[2,nempty.idx[2,]]-1)))
  log.lik_1order <- sum(log(parm$b1 + (1:(tmp.G.gr[2]-1))*d)) - fn.funky((parm$b1+1), (tmp.n[2]-1)) + sum(fn.funky((1-d+tmpF[nempty.idx[2,]]*parm$b1), (parm$clust$C.m.mt[2,nempty.idx[2,]]-1)))
  
  log.lik_0order <- sum(log(parm$b1 + (1:(parm$clust$G-1))*d)) - fn.funky((parm$b1+1), (parm$p-1)) + sum(fn.funky((1-d+tmpF*parm$b1), ((parm$clust$C.m.vec)-1)))
  return(list(log.lik_1order,log.lik_0order,tmpF))
}


fn.poissonDP.hyperparm <- function(data, parm, w=.01, max.d) 
{
  
  ## update parm$d conditional on parm$b1
  ## 1/w must be an integer
  
  d.v <- seq(0,max.d,by=w)
  len <- length(d.v)
  d.v <- d.v[-len]
  len <- len-1
  
  # 	if (parm$order==0){
  # 	  log.lik.v <- unlist(sapply(d.v, fn.d, parm)[2,])
  # 	}
  #	if (parm$order==1) {
  log.lik.v <- unlist(sapply(d.v, fn.d, parm)[1,])
  #	}
  # putting 1/2 prior mass on 0 and remaining spread uniformly on positive points in d.v
  log.p.v <- log(.5) + c(0,  rep(-log(len-1),(len-1)))
  
  log.post.v <- log.lik.v + log.p.v
  log.post.v <- log.post.v - max(log.post.v)
  
  log.post.2 <- c(log.post.v[1], log(sum(exp(log.post.v[-1]))))
  # conditional log-BF in favor of PDP models
  parm$PDP_log.BF <- log.post.2[2]-log.post.2[1]
  
  ## end added
  
  post.v <- exp(log.post.v)
  post.v <- post.v/sum(post.v)
  
  # plot(d.v, post.v, type="l")
  
  prop.d <- sample(d.v, size=1, prob=post.v)
  
  if (prop.d > 0)
  {prop.d <- runif(n=1, min=(prop.d-w), max=(prop.d+w))
  }
  
  if (prop.d != parm$d)
  {
    # MH ratio for independent proposals and
    # prior same for all d (which is true if 0 wp .5 and \in (0,max.d) wp .5)
    
    # 	  if (parm$order==0) {
    # 	    log.ratio <- fn.d(d=prop.d, parm)[[2]] - fn.d(d=parm$d, parm)[[2]]
    # 	  }
    #	  if (parm$order==1) {
    log.ratio <- fn.d(d=prop.d, parm)[[1]] - fn.d(d=parm$d, parm)[[1]]
    #		}
    prob <- min(1, exp(log.ratio))
    flip <- rbinom(n=1, size=1, prob=prob)
    if (flip==1)
    {parm$d <- prop.d
    }
  }
  
  parm$debug_tmpF<-fn.d(parm$d,parm)[[3]]
  parm
  
}

fn.order<-function (parm){
  log.prior.v<-rep(log(0.5),2)
  log.lik.v<-c(fn.d(d=parm$d, parm)[[1]],fn.d(d=parm$d, parm)[[2]])
  log.post.v <- log.lik.v + log.prior.v
  log.post.v <- log.post.v - max(log.post.v)
  post.v <- exp(log.post.v)
  post.v <- post.v/sum(post.v)
  
  parm$order<-rbinom(n=1, size=1, prob=post.v[1])
  parm$order_log.BF<-fn.d(d=parm$d, parm)[[1]]-fn.d(d=parm$d, parm)[[2]]
  parm
}

########################################

fn.iter <- function(data, parm, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, computeMode)
{
  
  
  parm <- fast_PDP_fn.main(parm, data, col.frac.probes, prob.compute.col.nbhd, max.col.nbhd.size, computeMode)
  
  parm <- fn.element.DP(data, parm, max.row.nbhd.size, row.frac.probes, computeMode)
  
  parm <- fn.poissonDP.hyperparm(data, parm, w=.01, max.d=1)
  
  #	parm <- fn.order(parm)
  
  parm <- fn.epsilon.DP(parm)
  
  parm <- fn.hyperparameters(data, parm)
  
  parm$iter_num<-parm$iter_num+1
  
  # 	flip <- rbinom(n=1, size=1, prob=.1)
  # 	if (flip==1)
  # 	  {parm <- fn.gen.missing.X(data, parm)
  # 	}
  
  err <- fn.quality.check(parm)
  if (err > 0)
  {stop(paste("failed QC: err=",err))
  }
  
  parm
}


fn.mcmc <- function(text, true, data, n.burn, n.reps, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, true_parm, dahl.flag=FALSE,
                    standardize.X=FALSE, flip.sign=FALSE, tBB_flag=FALSE, computeMode = "R")
{
  
  # initialize
  parm <- fn.init(true, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, tBB_flag, standardize.X, flip.sign, computeMode)
  init.parm <- parm
  
  err <- fn.quality.check(parm)
  if (err > 0)
  {stop(paste("failed QC at fn.init: err=",err))
  }
  
  for (cc in 1:n.burn)
  {parm <- fn.iter(data, parm, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, computeMode)
  
  if (cc %% 10 == 0)
  {print(paste(text, "BURN = ",cc,date(),"***"))
  }
  }
  
  ##########################################
  ## first get an estimated G cluster
  ##########################################
  
  
  All.Stuff <- NULL
  #
  All.Stuff$eta.v <- All.Stuff$d.v <- All.Stuff$tau.v <- All.Stuff$G.v <- All.Stuff$K.v <- All.Stuff$K_e.v <- All.Stuff$mu2.v <- All.Stuff$tau2.v <- All.Stuff$mu_e.v <- All.Stuff$tau_e.v <- array(,n.reps)
  All.Stuff$row.flip.v  <- array(0,n.reps)
  All.Stuff$nbhd_max <- All.Stuff$col_new_clust.v  <- All.Stuff$col_exit.v <- All.Stuff$col_flip.v <- All.Stuff$order_1.mh.flip.v <- array(0,n.reps)
  All.Stuff$PDP_log.BF.v<- All.Stuff$eta.0_log.BF.v<-array(0,n.reps)
  
  All.Stuff$pi.mt <- array(0,c(parm$p,parm$p))
  All.Stuff$clust$g.mt<-array(0,c(parm$p,n.reps))
  All.Stuff$clust$c.mt<-array(0,c(parm$p,n.reps))
  All.Stuff$clust$gr.mt<-array(0,c(parm$p,n.reps))
  All.Stuff$epsilon.mt<-array(,c(length(data$group),n.reps))
  All.Stuff$mean.taxicab.v_g2 <- All.Stuff$mean.taxicab.v  <- array(0,n.reps)
  All.Stuff$theta.mt<-array(,c(data$tot.trt,parm$p,n.reps))
  All.Stuff$eta_prob<-array(,c(200,n.reps))
  
  for (cc in 1:n.reps)
  {parm <- fn.iter(data, parm, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, computeMode)
  
  All.Stuff$G.v[cc] <- parm$clust$G
  All.Stuff$K.v[cc] <- parm$clust$K
  All.Stuff$K_e.v[cc] <- parm$clust$K_e
  All.Stuff$tau.v[cc] <- parm$tau
  All.Stuff$mu2.v[cc] <- parm$clust$mu2
  All.Stuff$tau2.v[cc] <- parm$clust$tau2
  All.Stuff$mu_e.v[cc] <- parm$clust$mu_e
  All.Stuff$tau_e.v[cc] <- parm$clust$tau_e
  
  All.Stuff$theta.mt[,,cc]<-parm$theta.mt
  #All.Stuff$tau_0.v[cc] <- parm$tau_0
  #All.Stuff$tau_int.v[cc] <- parm$tau_int
  All.Stuff$clust$g.mt[,cc]<-parm$clust$g.v
  All.Stuff$clust$gr.mt[,cc]<-parm$clust$gr.v
  All.Stuff$clust$c.mt[,cc]<-parm$clust$c.v
  All.Stuff$epsilon.mt[,cc] <- parm$epsilon
  
  All.Stuff$d.v[cc] <- parm$d
  All.Stuff$eta.v[cc] <- parm$eta
  All.Stuff$PDP_log.BF.v[cc] <- parm$PDP_log.BF
  All.Stuff$eta.0_log.BF.v[cc] <- parm$eta.0_log.BF
  All.Stuff$eta_prob[,cc]<-parm$eta_prob
  
  # 		if (dahl.flag)
  # 		{All.Stuff$c.matrix[cc,] <- parm$clust$c.v
  # 		}
  
  
  # summarizing elementwise DP in "fn.groupwise.updates"
  
  All.Stuff$row.flip.v[cc]  <- parm$clust$row.flip
  
  All.Stuff$col_new_clust.v[cc]  <- parm$clust$col.new.flag
  All.Stuff$order_1.mh.flip.v[cc]  <- parm$clust$order_1.mh.flip
  #    All.Stuff$eta_mh.flip.v[cc] <- parm$eta_mh.flip
  All.Stuff$col_flip.v[cc]  <- parm$clust$col.mh.flip
  All.Stuff$col_exit.v[cc]  <- parm$clust$col.mh.exit
  
  tmp.mat <- array(0,c(parm$p,parm$p))
  
  for (jj in 1:parm$clust$G)
  {indx.jj <- which(parm$clust$c.v==jj)
  tmp.mat[indx.jj,indx.jj] <- 1
  }
  
  All.Stuff$pi.mt <- All.Stuff$pi.mt + tmp.mat
  
  #All.Stuff$mean.taxicab.v[cc] <- mean(true_parm$clust$nbhd.matrix != tmp.mat)*parm$p/(parm$p-1)
  All.Stuff$mean.taxicab.v[cc] <- mean(true_parm$clust$nbhd.matrix != tmp.mat)
  All.Stuff$mean.taxicab.v_g2[cc] <- mean(true_parm$clust$nbhd.matrix[true_parm$clust$gr.v==2,true_parm$clust$gr.v==2] != tmp.mat[true_parm$clust$gr.v==2,true_parm$clust$gr.v==2])
  
  #All.Stuff$nbhd_max[cc] <- round(parm$clust$nbhd_max_dist, digits = 2)
  
  if (cc %% 10 == 0)
  {print(paste(text, "REPS = ",cc,date(),"***"))
  }
  
  } # end for loop in cc
  
  
  All.Stuff$pi.mt <- All.Stuff$pi.mt/n.reps
  
  All.Stuff$parm <- parm
  All.Stuff$init.parm <- init.parm
  
  All.Stuff
}

fn.mcmc_notaxi <- function(text, true, data, n.burn, n.reps, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, dahl.flag=FALSE,
                           standardize.X=FALSE, flip.sign=FALSE, tBB_flag=FALSE, computeMode = "R")
{
  
  # initialize
  message("BayesDiff initializing...")
  parm <- fn.init(true, data, max.row.nbhd.size, row.frac.probes, col.frac.probes, tBB_flag, standardize.X, flip.sign, computeMode)
  init.parm <- parm
  
  err <- fn.quality.check(parm)
  if (err > 0)
  {stop(paste("failed QC at fn.init: err=",err))
  }
  
  for (cc in 1:n.burn)
  {parm <- fn.iter(data, parm, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, computeMode)
  
  if (cc %% 10 == 0)
  {message(paste(text, "BURN = ",cc,date(),"***"))
  }
  }
  
  ##########################################
  ## first get an estimated G cluster
  ##########################################
  
  
  All.Stuff <- NULL
  #
  All.Stuff$eta.v <- All.Stuff$d.v <- All.Stuff$tau.v <- All.Stuff$G.v <- All.Stuff$K.v <- All.Stuff$K_e.v <- All.Stuff$mu2.v <- All.Stuff$tau2.v <- All.Stuff$mu_e.v <- All.Stuff$tau_e.v <- array(,n.reps)
  All.Stuff$row.flip.v  <- array(0,n.reps)
  All.Stuff$nbhd_max <- All.Stuff$col_new_clust.v  <- All.Stuff$col_exit.v <- All.Stuff$col_flip.v <- All.Stuff$order_1.mh.flip.v <- array(0,n.reps)
  All.Stuff$PDP_log.BF.v<- All.Stuff$eta.0_log.BF.v<-array(0,n.reps)
  
  All.Stuff$pi.mt <- array(0,c(parm$p,parm$p))
  All.Stuff$clust$g.mt<-array(0,c(parm$p,n.reps))
  All.Stuff$clust$c.mt<-array(0,c(parm$p,n.reps))
  All.Stuff$clust$gr.mt<-array(0,c(parm$p,n.reps))
  All.Stuff$epsilon.mt<-array(,c(length(data$group),n.reps))
  # All.Stuff$mean.taxicab.v_g2 <- All.Stuff$mean.taxicab.v  <- array(0,n.reps)
  All.Stuff$theta.mt<-array(,c(data$tot.trt,parm$p,n.reps))
  All.Stuff$eta_prob<-array(,c(200,n.reps))
  
  # 	if (dahl.flag)
  # 	  {All.Stuff$c.matrix <- array(0,c(n.reps,parm$p))
  # 	  }
  
  
  for (cc in 1:n.reps)
  {parm <- fn.iter(data, parm, max.row.nbhd.size, max.col.nbhd.size, row.frac.probes, col.frac.probes, prob.compute.col.nbhd, computeMode)
  
  All.Stuff$G.v[cc] <- parm$clust$G
  All.Stuff$K.v[cc] <- parm$clust$K
  All.Stuff$K_e.v[cc] <- parm$clust$K_e
  All.Stuff$tau.v[cc] <- parm$tau
  All.Stuff$mu2.v[cc] <- parm$clust$mu2
  All.Stuff$tau2.v[cc] <- parm$clust$tau2
  All.Stuff$mu_e.v[cc] <- parm$clust$mu_e
  All.Stuff$tau_e.v[cc] <- parm$clust$tau_e
  
  All.Stuff$theta.mt[,,cc]<-parm$theta.mt
  #All.Stuff$tau_0.v[cc] <- parm$tau_0
  #All.Stuff$tau_int.v[cc] <- parm$tau_int
  All.Stuff$clust$g.mt[,cc]<-parm$clust$g.v
  All.Stuff$clust$gr.mt[,cc]<-parm$clust$gr.v
  All.Stuff$clust$c.mt[,cc]<-parm$clust$c.v
  All.Stuff$epsilon.mt[,cc] <- parm$epsilon
  
  All.Stuff$d.v[cc] <- parm$d
  All.Stuff$eta.v[cc] <- parm$eta
  All.Stuff$PDP_log.BF.v[cc] <- parm$PDP_log.BF
  All.Stuff$eta.0_log.BF.v[cc] <- parm$eta.0_log.BF
  All.Stuff$eta_prob[,cc]<-parm$eta_prob
  
  # summarizing elementwise DP in "fn.groupwise.updates"
  
  All.Stuff$row.flip.v[cc]  <- parm$clust$row.flip
  
  All.Stuff$col_new_clust.v[cc]  <- parm$clust$col.new.flag
  All.Stuff$order_1.mh.flip.v[cc]  <- parm$clust$order_1.mh.flip
  #    All.Stuff$eta_mh.flip.v[cc] <- parm$eta_mh.flip
  All.Stuff$col_flip.v[cc]  <- parm$clust$col.mh.flip
  All.Stuff$col_exit.v[cc]  <- parm$clust$col.mh.exit
  
  tmp.mat <- array(0,c(parm$p,parm$p))
  
  for (jj in 1:parm$clust$G)
  {indx.jj <- which(parm$clust$c.v==jj)
  tmp.mat[indx.jj,indx.jj] <- 1
  }
  
  All.Stuff$pi.mt <- All.Stuff$pi.mt + tmp.mat
  
  if (cc %% 10 == 0)
  {message(paste(text, "REPS = ",cc,date(),"***"))
  }
  
  } # end for loop in cc
  
  
  All.Stuff$pi.mt <- All.Stuff$pi.mt/n.reps
  
  All.Stuff$parm <- parm
  All.Stuff$init.parm <- init.parm
  
  All.Stuff
}

