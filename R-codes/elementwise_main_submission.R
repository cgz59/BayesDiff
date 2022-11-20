
fn1.update.element.objects <- function(parm, computeMode)
	{

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
	
	for (g in 1:parm$clust$G)
	{s.g.v <- parm$clust$s.mt[,g]
	#s.pos.indx <- s.g.v > 0
	#
	#if (sum(s.pos.indx) > 0)
	{parm$clust$A.mt[,g] <- parm$clust$phi.v[s.g.v]
	}
	# 		if ((parm$n2-sum(s.pos.indx)) > 0)
	# 			{parm$clust$A.mt[!s.pos.indx,g] <- 0
	# 			}
	# 		parm$clust$B.mt[,(g+1)] <- parm$clust$A.mt[,g]
	}
	
	tmp.s.v_g1<-parm$clust$s.v[parm$clust$s.stretch==1]
	g1_idx<-which(parm$clust$gr.c.v==1)
	if (length(unique(tmp.s.v_g1))<length(tmp.s.v_g1)) {
	  for (mm in 1:length(unique(tmp.s.v_g1))) {
	    tmp.same.idx_g1<-which(tmp.s.v_g1==unique(tmp.s.v_g1)[mm])
	    if (length(tmp.same.idx_g1)>1) {
	      tmp.merge.idx<-g1_idx[tmp.same.idx_g1[-1]]
	      tmp.keep.idx<-g1_idx[tmp.same.idx_g1[1]]
	      parm$clust$c.v<-replace(parm$clust$c.v,parm$clust$c.v %in% tmp.merge.idx,tmp.keep.idx)
	      parm$clust$C.m.mt[,tmp.keep.idx]<-rowSums(parm$clust$C.m.mt[,c(tmp.keep.idx,tmp.merge.idx)])
	      parm$clust$C.m.mt[,tmp.merge.idx]<-0
	      parm$clust$n.vec[unique(tmp.s.v_g1)[mm]]<-parm$clust$n.vec[unique(tmp.s.v_g1)[mm]]-(length(tmp.same.idx_g1)-1)
	    }
	  }
	}
	parm$clust$C.m.vec<-colSums(parm$clust$C.m.mt)
	
	parm<-PDP_fn.drop(parm, computeMode)

	parm$clust$theta.v <- as.vector(parm$clust$A.mt)
  parm$theta.mt<-parm$clust$A.mt[,parm$clust$c.v]

	parm
	}


fn2.update.element.objects <- function(parm, computeMode) 
	{
  parm$N <- parm$tot.trt * sum(parm$clust$gr.c.v==2) + sum(parm$clust$gr.c.v==1)
  
	parm$Y <- parm$Z.sd <- parm$num.Y <- parm$clust$s.stretch <- NULL

	# group covariate tells which parm$clust$rho.g
	# to use for likelihood calculation
	parm$g <- rep(1:parm$clust$G,each = parm$tot.trt)
	parm$trt.exp<-rep(1:parm$tot.trt,parm$clust$G)
  

	for (g in 1:parm$clust$G) {
		I.g <- (parm$clust$c.v==g)
		
		if (parm$clust$gr.c.v[g]==1) {
		  x.tmp<-as.vector(parm$Z.f[,I.g])
		  parm$Y<-c(parm$Y,mean(x.tmp))
		  parm$num.Y<-c(parm$num.Y,length(x.tmp))
		  parm$Z.sd<-c(parm$Z.sd,sd(x.tmp))
		  parm$clust$s.stretch<-c(parm$clust$s.stretch,1)
		}
		
		if (parm$clust$gr.c.v[g]==2) {
		  x.tmp.mt<-(parm$Z.f[,I.g,drop=FALSE])
		  for (ii in 1:parm$tot.trt) {
		    x.tmp<-as.vector(x.tmp.mt[parm$group==ii,])
		    parm$Y<-c(parm$Y,mean(x.tmp))
		    parm$num.Y<-c(parm$num.Y,length(x.tmp))
		    parm$Z.sd<-c(parm$Z.sd,sd(x.tmp))
		    parm$clust$s.stretch<-c(parm$clust$s.stretch,0)
		  }
		}
		}

	if (sum(c(length(parm$Y)!=parm$N,length(parm$num.Y)!=parm$N,length(parm$Z.sd)!=parm$N),length(parm$clust$s.stretch)!=parm$N)>0) 
	{stop(paste('problem in elementwise_main-fn2.update'))}

	parm
	}



fn.element.DP <- function(data, parm, max.row.nbhd.size, row.frac.probes,
                          computeMode)
{
  parm$Z.f<-pre.Z<-parm$X.LANE-matrix(rep(parm$epsilon,parm$p),ncol=parm$p)
  #data.X.tmp<-tmp.Z
  tmp.Z<-matrix(,parm$tot.trt,parm$p)
  for (ii in 1:parm$tot.trt){
    tmp.Z[ii,]<-colMeans(pre.Z[parm$group==ii,,drop=F])
  }
  parm$Z<-tmp.Z
  

#   if (parm$standardize.X)
#   {parm <- fn.standardize_orient.X(parm)
#   }

	# essentially, a Bush-Mac move: given groups, the parm$N=n2XG number of
	# invidividual elements (summaries of microarray elements) belonging to group g>0
	# are updated for s (phi) and z

	parm <- fn2.update.element.objects(parm, computeMode)

	parm <- element_fn.fast.DP(parm, max.row.nbhd.size, row.frac.probes, computeMode)

	#############################
	## Important: do not remove call to fn1.update.element.objects
	## updates A.mt, theta.v, B.mt, tBB.mt, s.v, s.mt, n.vec, n0
	#############################

	parm <- fn1.update.element.objects(parm, computeMode)

  	parm
}



