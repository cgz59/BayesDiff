
gen.X <- function(trt,group, p, true_parm, LANE) {
  
  data <- NULL
  data$X.raw <- array(,c(length(group),p))
  data$mean.theta <- data$mean.all <- array(,c(length(group),p))
  data$mean.chi <- numeric(p)
  base_p <- c(0.9, 0.5, 0.1, 0.5)
  
  for (xx in 1:p)
  {
    x.v <- array(,length(group))
    mean.theta<-true_parm$clust$A.mt[,true_parm$clust$c.v[xx]]
    p_current <- base_p[true_parm$clust$h.v[xx]]
    mean.chi <- log(p_current/(1-p_current))
    
    mean.all<-mean.theta[group]+mean.chi+LANE
    
    x.v <- rnorm(n=length(group),mean=mean.all,sd=true_parm$tau[xx])
    
    data$X.raw[,xx] <- x.v
    data$mean.all[,xx] <- mean.all
    data$mean.theta[,xx] <- mean.theta
    data$mean.chi[xx] <- mean.chi
  }
  
  data$G.max <- round(dim(data$X)[2]/2)
  data$K.max <- round(data$G.max*trt*.1*.5+data$G.max*.9*.5)
  
  return(data)
}


gen.clust <- function(trt,group, p,true) {
  
  true_parm <- NULL
  #####################PDP column clusters######################
  true_parm$d <- 1/3
  #  true_parm$d <- 0
  true_parm$tot.trt<-trt
  
  true_parm$clust <- NULL
  true_parm$clust$c.v <- array(,p)
  true_parm$clust$c.v[1] <- 1
  true_parm$clust$C.m.mt<-matrix(c(1,0),ncol=1)
  true_parm$clust$G <- 1
  true_parm$b1 <- true$b1
  
  true_parm$delta_margin<-true$delta_margin
  true_parm$rho<-true$rho
  #  true_parm$rho<-c(0.9,0.1)
  true_parm$eta<-true$eta
  true_parm$distance<-true$distance
  
  ##################################################################################
  ##generate h.v == methylation level status, from HMM, this decides different chi_j----
  ##################################################################################
  pi_base <- c(0.05,0.95,0.08,0.95)
  
  true_parm$clust$h.v <- array(NA,p)
  true_parm$clust$h.v[1] <- sample(1:4, size = 1)
  for (hh in 2:p) {
    pi_hh <- pi_base*exp(-1.895e-2*log(true_parm$distance[hh-1]))
    h_prev <- true_parm$clust$h.v[hh-1]
    prob_out <- pi_hh[h_prev]
    if(runif(1)<prob_out) {
      h_current <- h_prev + 1
      h_current <- ifelse(h_current == 5, 1, h_current)
    } else {
      h_current <- h_prev
    }
    true_parm$clust$h.v[hh] <- h_current
  }
  
  true_parm$clust$mu2 <- 0
  true_parm$clust$tau2 <- 1
  
  #  true_parm$clust$gr.c.v<-rbinom(1,1,0.1)+1
  true_parm$clust$gr.c.v<-1
  true_parm$clust$gr.v<-true_parm$clust$gr.c.v
  
  if (true_parm$clust$gr.v==1){
    s.v.tmp <- rep(1,trt)
    true_parm$clust$n.vec <- 1
    true_parm$clust$K <- 1
    true_parm$clust$M <- true$M
  } #end if (gr.v==1)
  
  if (true_parm$clust$gr.v==2) {
    s.v.tmp <- 1
    true_parm$clust$n.vec <- 1
    true_parm$clust$K <- 1
    true_parm$clust$M <- true$M
    
    for (yy in 2:trt)
    {prob.v <- (true_parm$clust$n.vec)
    prob.v <- c(prob.v, true_parm$clust$M)
    new.s <- sample(1:(true_parm$clust$K+1), size=1, prob=prob.v)
    
    s.v.tmp <- c(s.v.tmp,new.s)
    new.flag.s <- (new.s > true_parm$clust$K)
    #zero.flag <- new.s==0
    
    if (new.flag.s)
    {true_parm$clust$K <- true_parm$clust$K + 1
    true_parm$clust$n.vec <- c(true_parm$clust$n.vec, 1)
    }
    if (!new.flag.s)
    {true_parm$clust$n.vec[new.s] <- true_parm$clust$n.vec[new.s] + 1
    }
    }
  } #end if (gr.v==2)
  
  true_parm$clust$s.mt<-matrix(s.v.tmp,ncol=1)
  true_parm$clust$phi.v <- rnorm(n=true_parm$clust$K, mean=true_parm$clust$mu2, sd=true_parm$clust$tau2)
  true_parm$clust$A.mt<-true_parm$clust$phi.v[s.v.tmp]
  
  true_parm$clust$g.v<-true_parm$clust$gr.v
  
  for (xx in 2:p)
  {tmp.h.left<-true_parm$clust$gr.v[xx-1]
  if (true_parm$eta>0) {
    tmp.r<-exp(-true_parm$distance[xx-1]/true_parm$eta)
  }
  if (true_parm$eta==0) {
    tmp.r<-0
  }
  # if (true_parm$eta==-1) {
  #   tmp.r<-1
  # }
  tmp.prob.g.xx<-true_parm$rho[tmp.h.left]+(1-true_parm$rho[tmp.h.left])*tmp.r
  tmp.prob.g.xx.v<-rep(1-tmp.prob.g.xx,2)
  tmp.prob.g.xx.v[tmp.h.left]<-tmp.prob.g.xx
  tmp.g<-sample(c(1,2),size=1,prob=tmp.prob.g.xx.v)
  
  token<-rbinom(1,1,(tmp.g==1)*0+(tmp.g==2)*1) #transition in/out prob
  tmp.gr<-token+1
  
  true_parm$clust$g.v<-c(true_parm$clust$g.v,tmp.g)
  tmp.C.m.vec<-true_parm$clust$C.m.mt[tmp.g,]
  prob.v <- tmp.C.m.vec - true_parm$d*as.numeric(tmp.g==2)
  prob.v <- c(prob.v, (true_parm$b1 + sum(tmp.C.m.vec>0)*true_parm$d*as.numeric(tmp.g==2)))
  prob.v[prob.v<0]<-0
  new.c <- sample(1:(true_parm$clust$G+1), size=1, prob=prob.v)
  
  true_parm$clust$c.v[xx] <- new.c
  new.flag <- (new.c > true_parm$clust$G)
  
  
  if (new.flag) {
    true_parm$clust$G <- true_parm$clust$G + 1
    tmp.C.m<-c(0,0)
    tmp.C.m[tmp.g]<-1
    true_parm$clust$C.m.mt <- cbind(true_parm$clust$C.m.mt, tmp.C.m)
    
    if (tmp.gr==2) {
      s.v.tmp<-rep(-1,trt)
      old.K<-true_parm$clust$K
      while (length(unique(s.v.tmp))<2) {
        for (yy in 1:trt)
        {prob.v <- (true_parm$clust$n.vec)
        prob.v <- c(prob.v, true_parm$clust$M)
        new.s <- sample(1:(true_parm$clust$K+1), size=1, prob=prob.v)
        
        s.v.tmp[yy] <- new.s
        new.flag.s <- (new.s > true_parm$clust$K)
        #zero.flag <- new.s==0
        
        if (new.flag.s)
        {true_parm$clust$K <- true_parm$clust$K + 1
        true_parm$clust$n.vec <- c(true_parm$clust$n.vec, 1)
        }
        if (!new.flag.s)
        {true_parm$clust$n.vec[new.s] <- true_parm$clust$n.vec[new.s] + 1
        }
        } #end for (yy)
      }
    } #end if tmp.gr==2
    
    if (tmp.gr==1) {
      old.K<-true_parm$clust$K
      prob.v <- (true_parm$clust$n.vec)
      prob.v <- c(prob.v, true_parm$clust$M)
      # new.s<-true_parm$clust$s.mt[1,1]
      # while(sum(new.s==true_parm$clust$s.mt[1,true_parm$clust$gr.c.v==1])>0) {
      new.s <- sample(1:(true_parm$clust$K+1), size=1, prob=prob.v)
      # }
      s.v.tmp <- rep(new.s,trt)
      #       sum(colSums(apply(true_parm$clust$s.mt,2,'==',s.v.tmp))==trt)>0
      
      new.flag.s <- (new.s > true_parm$clust$K)
      if (new.flag.s)
      {true_parm$clust$K <- true_parm$clust$K + 1
      true_parm$clust$n.vec <- c(true_parm$clust$n.vec, 1)
      }
      if (!new.flag.s)
      {true_parm$clust$n.vec[new.s] <- true_parm$clust$n.vec[new.s] + 1
      }
      
    } #end tmp.gr==1
    true_parm$clust$s.mt<-cbind(true_parm$clust$s.mt,s.v.tmp)
    true_parm$clust$phi.v <- c(true_parm$clust$phi.v,rnorm(n=true_parm$clust$K-old.K, mean=true_parm$clust$mu2, sd=true_parm$clust$tau2))
    tmp.a.v<-true_parm$clust$phi.v[s.v.tmp]
    true_parm$clust$A.mt<-cbind(true_parm$clust$A.mt,tmp.a.v)
    
    true_parm$clust$gr.c.v<-c(true_parm$clust$gr.c.v,tmp.gr)
    true_parm$clust$gr.v<-true_parm$clust$gr.c.v[true_parm$clust$c.v]
    
  } #end if new.flag
  if (!new.flag)
  {true_parm$clust$C.m.mt[tmp.g,new.c] <- true_parm$clust$C.m.mt[tmp.g,new.c] + 1
  true_parm$clust$gr.v<-true_parm$clust$gr.c.v[true_parm$clust$c.v]
  }
  }
  
  
  
  ################DP clusters, latent elements#################
  true_parm$N <- trt * sum(true_parm$clust$gr.c.v==2) + sum(true_parm$clust$gr.c.v==1)
  true_parm$clust$s.stretch <- true_parm$clust$s.v<- NULL
  if (is.null(true_parm$clust$s.v)==F) {stop(paste('aaaa'))}
  for (g in 1:true_parm$clust$G) {
    if (true_parm$clust$gr.c.v[g]==1) {
      if (length(unique(true_parm$clust$s.mt[,g]))>1) {stop(paste('problem in gen.data with g1',g,unique(true_parm$clust$s.mt[,g])))}
      true_parm$clust$s.v<-c(true_parm$clust$s.v,unique(true_parm$clust$s.mt[,g]))
      true_parm$clust$s.stretch<-c(true_parm$clust$s.stretch,1)
    }
    if (true_parm$clust$gr.c.v[g]==2) {
      true_parm$clust$s.v<-c(true_parm$clust$s.v,(true_parm$clust$s.mt[,g]))
      true_parm$clust$s.stretch<-c(true_parm$clust$s.stretch,rep(0,trt))
    }
  }
  if (sum(c(length(true_parm$clust$s.v)!=true_parm$N,length(true_parm$clust$s.stretch)!=true_parm$N)>0))
  {stop(paste('problem in gen.data',true_parm$N,length(true_parm$clust$s.v),length(true_parm$clust$s.stretch)))}
  
  ######merge clusters#####
  tmp.s.v_g1<-true_parm$clust$s.v[true_parm$clust$s.stretch==1]
  g1_idx<-which(true_parm$clust$gr.c.v==1)
  if (length(unique(tmp.s.v_g1))<length(tmp.s.v_g1)) {
    for (mm in 1:length(unique(tmp.s.v_g1))) {
      tmp.same.idx_g1<-which(tmp.s.v_g1==unique(tmp.s.v_g1)[mm])
      if (length(tmp.same.idx_g1)>1) {
        tmp.merge.idx<-g1_idx[tmp.same.idx_g1[-1]]
        tmp.keep.idx<-g1_idx[tmp.same.idx_g1[1]]
        true_parm$clust$c.v<-replace(true_parm$clust$c.v,true_parm$clust$c.v %in% tmp.merge.idx,tmp.keep.idx)
        true_parm$clust$C.m.mt[,tmp.keep.idx]<-rowSums(true_parm$clust$C.m.mt[,c(tmp.keep.idx,tmp.merge.idx)])
        true_parm$clust$C.m.mt[,tmp.merge.idx]<-0
        true_parm$clust$n.vec[unique(tmp.s.v_g1)[mm]]<-true_parm$clust$n.vec[unique(tmp.s.v_g1)[mm]]-(length(tmp.same.idx_g1)-1)
      }
    }
  }
  true_parm$clust$C.m.vec<-colSums(true_parm$clust$C.m.mt)
  
  true_parm<-PDP_fn.drop(true_parm, computeMode)
  
  true_parm$num.dm<-sum(true_parm$clust$gr.v==2)
  
  ## neighborhood taxicab distances for column clusters
    tmp.mat <- array(0,c(p,p))
  
  for (jj in 1:true_parm$clust$G)
  {indx.jj <- which(true_parm$clust$c.v==jj)
  tmp.mat[indx.jj,indx.jj] <- 1
  }
  
  true_parm$clust$nbhd.matrix <- tmp.mat
  

  ###########################epsilon clusters#######################
  
  true_parm$n2<-length(group)
  true_parm$clust$s_e.v <- array(,true_parm$n2)
  true_parm$clust$s_e.v[1] <- 1
  true_parm$clust$n_e.v <- 1
  true_parm$clust$K_e <- 1
  true_parm$clust$Me <- true$Me
  
  
  for (xx in 2:true_parm$n2)
  {prob.v <- (true_parm$clust$n_e.v)
  prob.v <- c(prob.v, true_parm$clust$Me)
  new.s_e <- sample(1:(true_parm$clust$K_e+1), size=1, prob=prob.v)
  
  true_parm$clust$s_e.v[xx] <- new.s_e
  new.flag <- (new.s_e > true_parm$clust$K_e)
  #zero.flag <- new.s==0
  
  if (new.flag)
  {true_parm$clust$K_e <- true_parm$clust$K_e + 1
  true_parm$clust$n_e.v <- c(true_parm$clust$n_e.v, 1)
  }
  #if (zero.flag)
  #{true_parm$clust$n0 <- true_parm$clust$n0 + 1
  #}
  if (!new.flag)#&(!zero.flag)
  {true_parm$clust$n_e.v[new.s_e] <- true_parm$clust$n_e.v[new.s_e] + 1
  }
  }
  
  #sum(true_parm$clust$n.vec) + true_parm$clust$n0 == true_parm$N
  
  #true_parm$clust$s.mt <- matrix(true_parm$clust$s.v, nrow=T)
  
  ##########
  
  true_parm$clust$mu_e <- 0
  true_parm$clust$tau_e <- 0.35
  
  true_parm$clust$phi_e.v <- rnorm(n=true_parm$clust$K_e, mean=true_parm$clust$mu_e, sd=true_parm$clust$tau_e)
  
  ##################################################
  true_parm$tau<-true$tau
  
  #s.mt.c<-true_parm$clust$s.mt[,true_parm$clust$c.v]
  true_parm$theta.mt<-true_parm$clust$A.mt[,true_parm$clust$c.v]
  return (true_parm)
}

gen.data<-function(trt,group,p,true,LANE){
  true_parm<-gen.clust(trt,group, p,true)
  true_parm$p<-p
  #  true_parm$tot.trt<-trt
  #data<-NULL
  data<-gen.X(trt,group,p,true_parm,LANE)
  data$LANE<-LANE
  #data$H2<-H2
  data$group<-group
  data$tot.trt<-trt
  return(list(data=data,true_parm=true_parm))
}
