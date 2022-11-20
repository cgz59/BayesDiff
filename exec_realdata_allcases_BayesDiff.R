####preloading####
load('probe_gene_grp_info_50k.rda')
# load('meth_dataframe.rda')
# load('meth_datamt_STAD.rda')
# load('meth_datamt_LIHC.rda')
# load('meth_datamt_ESCA.rda')
# load('meth_datamt_PAAD.rda')
# library(SummarizedExperiment)
# library(biomaRt)
slurm_arrayid <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
message(slurm_arrayid)

# options(error=recover)
# options(warn=2)
# setwd("/home/cgz59/methTCGA/R_codes")
source("./R_codes/iterations_1order_stochastic_group_new_notaxi.R")
source("./R_codes/elementwise_main_1order_new.R")
source("./R_codes/elementwise_DP.functions_1order_new.R")
source("./R_codes/fast_PDP.functions_final.R")
source("./R_codes/e_DP.functions.R")
# source("gen.data_sigmaidentical_1order_stochastic_group_inverse_eta_control_chain_new.R")

true<-NULL
true$M<-20 #alpha_2,DP mass parameter(elementwise,phi's)
#true$a.R<-11
true$b1<-20 #alpha_1,PDP mass parameter
true$Me<-10 #alpha_3,epsilon_i DPmass parameter
true$shift <- 1e-4
true$FDR<-0.05 #desired FDR level

# mt_meth.all<-cbind(mt_meth.STAD, mt_meth.LIHC, mt_meth.ESCA, mt_meth.PAAD)
# rm(mt_meth.STAD, mt_meth.LIHC, mt_meth.ESCA, mt_meth.PAAD)
# gc()

n.burn <- 1000
n.reps <- 1000

gene_subset_idx<-which(global_realdata$num_probe>10)
# GeneListInfo$hgnc_symbol[gene_subset_idx]
# gene_subset_idx<-1:global_realdata$GG
# further analysis of top 16 genes
# load('real_data_gene_propdm_p5.rda')
# gene_subset_idx<-gene_propdm_p5

mean.discount.v<-NULL
discount.q1.v<-discount.q2.v<-NULL
eta.q1.v<-eta.q2.v<-NULL
eta.0.prop.v<-mean.eta.v<-median.eta.v<-NULL
lowerBound_PDP_BF.v<-se_LowerBound_PDP_BF.v<-NULL
eta.0_BF.v<-se_eta.0_BF.v<-NULL
# eta_prob_all<-array(NA,c(200,n.reps,length(gene_subset_idx)))
order_1.mh.rate<-NULL
de.gene.idx.ls<-de.post_prob.ls<-de.gene.FDR.idx.ls<-NULL
probe_level_ordered_diff_list.ls<-gene_level_max_grp_diff.ls<-NULL
prop_de.mt <- NULL
###########start#########
# library(foreach)
# library(doParallel)
# registerDoParallel(cores=2)
# results<-foreach(it = 1:2) %dopar% {
# for (it in (20*(slurm_arrayid-1)+1):min(length(gene_subset_idx),slurm_arrayid*20)) {
for (it in 1:209) {
  set.seed(124)
  tt<-gene_subset_idx[it]
  probe_idx_gene_tt<-unlist(global_realdata$probe_idx.allgenes[tt])
  load('meth_datamt_STAD.rda')
  data_touse.STAD <- mt_meth.STAD[probe_idx_gene_tt,]
  rm(mt_meth.STAD)
  gc()
  load('meth_datamt_LIHC.rda')
  data_touse.LIHC <- mt_meth.LIHC[probe_idx_gene_tt,]
  rm(mt_meth.LIHC)
  gc()
  load('meth_datamt_ESCA.rda')
  data_touse.ESCA <- mt_meth.ESCA[probe_idx_gene_tt,]
  rm(mt_meth.ESCA)
  gc()
  load('meth_datamt_PAAD.rda')
  data_touse.PAAD <- mt_meth.PAAD[probe_idx_gene_tt,]
  rm(mt_meth.PAAD)
  gc()
  data_touse<-cbind(data_touse.STAD,data_touse.LIHC,data_touse.ESCA,data_touse.PAAD)
  rm(data_touse.STAD,data_touse.LIHC,data_touse.ESCA,data_touse.PAAD)
  gc()
  # data_touse<-mt_meth.all[probe_idx_gene_tt,]
  data_touse<-replace(data_touse,is.na(data_touse),rowMeans(data_touse,na.rm=T)[which(is.na(data_touse)==T,arr.ind=T)[,1]])
  data_touse<-replace(data_touse,(data_touse==0),min(data_touse[data_touse>0])/2)
  data_touse<-replace(data_touse,(data_touse==1),mean(c(1,max(data_touse[data_touse<1]))))
  
  data<-NULL
  data$raw<-t(data_touse)
  data$X<-log(data$raw)-log(1-data$raw)
  #data$LANE<-log(calcNormFactors(t(data$raw),method="upperquartile"))
  data$LANE<-rep(0,length(global_realdata$type_grp))
  data$rho<-0.25
  data$distance<-unlist(global_realdata$distances.allgenes[tt])
  data$distance<-data$distance/sum(data$distance)
  data$group<-global_realdata$type_grp
  data$tot.trt<-length(unique(data$group))
  data$G.max <- round(dim(data$X)[2]/2)
  data$K.max <- round(data$G.max*data$tot.trt*.1*.5+data$G.max*.9*.5)
  
  all.stuff <- fn.mcmc(text="CLUST ANALYZE...",
                       true, data,
                       n.burn, n.reps, max.row.nbhd.size = nrow(data$X)*data$G.max/100, # should be small compared to n2*p^d (~ n2*G if d=.5)
                       max.col.nbhd.size = round(ncol(data$X)/5), # should be small compared to p
                       row.frac.probes = 0.25,
                       col.frac.probes = .5,
                       prob.compute.col.nbhd=.2,
                       true_parm,
                       dahl.flag=FALSE,
                       standardize.X=FALSE,
                       flip.sign=FALSE,
                       tBB_flag=FALSE,  computeMode= createComputeMode())
  
  mean.discount.v<-c(mean.discount.v,mean(all.stuff$d.v))
  discount.q1.v<-c(discount.q1.v,quantile(all.stuff$d.v,probs=c(0.05,0.95))[[1]])
  discount.q2.v<-c(discount.q2.v,quantile(all.stuff$d.v,probs=c(0.05,0.95))[[2]])
  order_1.mh.rate<-c(order_1.mh.rate,mean(all.stuff$order_1.mh.flip.v))
  
  #hist(all.stuff$eta.v,breaks=50)
  mean.eta.v<-c(mean.eta.v,mean(all.stuff$eta.v))
  median.eta.v<-c(median.eta.v,median(all.stuff$eta.v))
  eta.0.prop.v<-c(eta.0.prop.v,mean(all.stuff$eta.v==0))
  eta.q1.v<-c(eta.q1.v,quantile(all.stuff$eta.v,probs=c(0.05,0.95))[[1]])
  eta.q2.v<-c(eta.q2.v,quantile(all.stuff$eta.v,probs=c(0.05,0.95))[[2]])
  # eta_prob_all[,,it]<-all.stuff$eta_prob
  #all.stuff$eta_mh.flip.v
  #eta_mh.rate<-c(eta_mh.rate,mean(all.stuff$eta_mh.flip.v,na.rm=T))
  # estimated lower bound for log-BF in favor of PDP models
  # useful only if this value > log(10), i.e. strong evidence in favor of PDP
  lowerBound_PDP_BF.v<-c(lowerBound_PDP_BF.v,mean(all.stuff$PDP_log.BF.v))
  #quantile(all.stuff$PDP_log.BF.v,probs=c(0.05,0.95))
  # s.e. for estimate 
  se_LowerBound_PDP_BF.v<-c(se_LowerBound_PDP_BF.v,sd(all.stuff$PDP_log.BF.v)/sqrt(n.reps))
  
  eta.0_BF.v<-c(eta.0_BF.v,mean(all.stuff$eta.0_log.BF.v))
  se_eta.0_BF.v<-c(se_eta.0_BF.v,sd(all.stuff$eta.0_log.BF.v)/sqrt(n.reps))

  #de.gene.ind<-rowMeans(apply(try1,c(2,3),max)-apply(try1,c(2,3),min)>true$delta_margin)>0.5
  de.gene.idx<-which(rowMeans(all.stuff$clust$gr.mt)>1.5)
  # de.post_prob<-rowMeans(all.stuff$clust$gr.mt==2)
  de.gene.idx.ls[it]<-list(de.gene.idx)
  
  #get which sample group is largest/smallest for each mcmc sample 
  min_grp.mt <- apply(all.stuff$theta.mt, c(2,3), function(x) which(x == min(x)))
  max_grp.mt <- apply(all.stuff$theta.mt, c(2,3), function(x) which(x == max(x)))
  
  ##setting up some threshold for being dm in individual mcmc sample, get post_prob of dm and which pair gives largest difference
  threshold_list <- seq(from = 0.5, to = 1.5, by = 0.05)
  
  de.post_prob.mt <- NULL
  probe_level_ordered_diff.mt <- matrix(NA,ncol(data$X), length(threshold_list))
  pair_mt<-cbind(combn(1:data$tot.trt,2),combn(1:data$tot.trt,2)[2:1,])
  
  for(tt in 1:length(threshold_list)) {
    threshold <- threshold_list[tt]
    de.ind_post.mt <- apply(all.stuff$theta.mt, c(2,3), function(x) max(x)-min(x)>threshold) & all.stuff$clust$gr.mt==2
    de.post_prob.mt <- cbind(de.post_prob.mt,
                             rowMeans(de.ind_post.mt))
    
    for (pp in 1:nrow(de.ind_post.mt)) {
      
      min.pp_lst <- min_grp.mt[pp,de.ind_post.mt[pp,]==T]
      max.pp_lst <- max_grp.mt[pp,de.ind_post.mt[pp,]==T]
      
      if (length(min.pp_lst) != length(max.pp_lst)) {stop("min/max group length unequal")}
      
      if (length(min.pp_lst)>0 & length(max.pp_lst)>0) {
        ordered_diff_pair_count<-rep(0,12)
        for (iiit in 1:length(min.pp_lst)) {
          tmp_min.v<-unlist(min.pp_lst[iiit])
          tmp_max.v<-unlist(max.pp_lst[iiit])
          tmp_diff_pair<-expand.grid(tmp_max.v,tmp_min.v)
          if (nrow(tmp_diff_pair)>1) {
            for (pit in 1:nrow(tmp_diff_pair)) {
              tmp_pair_idx<-which(colSums(apply(pair_mt,2,'==',tmp_diff_pair[pit,]))==2)
              ordered_diff_pair_count[tmp_pair_idx]<-ordered_diff_pair_count[tmp_pair_idx]+1/nrow(tmp_diff_pair)
            }
          }
          else {
            tmp_pair_idx<-which(colSums(apply(pair_mt,2,'==',tmp_diff_pair))==2)
            ordered_diff_pair_count[tmp_pair_idx]<-ordered_diff_pair_count[tmp_pair_idx]+1/nrow(tmp_diff_pair)
          }
        } ##end looping through list of min and max
        
        probe_level_ordered_diff.mt[pp,tt] <- which.max(ordered_diff_pair_count)
      } else { ##end if there is any post burn mcmc sample indicate pp probe is dm
        probe_level_ordered_diff.mt[pp,tt] <- 0 ##if no difference, largest diff pair is arbitrarily set to "0"
      }
    } ##end looping through all probes
  } ##end looping through all threshold options
  
  gene_level_max_grp_diff.v <- apply(probe_level_ordered_diff.mt, 2, function(x) 
    names(which.max(table(x[x!=0]))))
  
  ##store results in master list
  probe_level_ordered_diff_list.ls[it]<-list(probe_level_ordered_diff.mt) #probe level max pairwise group difference for gene tt (vector of length num_probe[tt])
  gene_level_max_grp_diff.ls[it]<-list(gene_level_max_grp_diff.v) #gene level max pairwise group difference for gene tt
  # de.post_prob.ls[it]<-list(de.post_prob)
  de.post_prob.ls[it]<-list(de.post_prob.mt)

  de.post_prob.sorted<-apply(de.post_prob.mt, 2, sort, decreasing = T)
  de.post_prob.order<-apply(de.post_prob.mt, 2, order, decreasing = T)
  # de.post_prob.rank<-apply(de.post_prob.mt, 2, rank)
  de.post_prob.cummean<-apply(de.post_prob.sorted, 2, function(x) 
    cumsum(1-x)/(1:nrow(de.post_prob.mt)))
  # tmp <- de.post_prob.cummean<=true$FDR
  tmptmp.de.gene.FDR.idx <- matrix(F, nrow = nrow(de.post_prob.order), ncol = ncol(de.post_prob.order))
  tmp.de.gene.FDR.idx <-sapply(1:ncol(de.post_prob.order), function(x) {
    tmp.ind <- tmptmp.de.gene.FDR.idx[,x]
    tmp.ind[de.post_prob.order[(de.post_prob.cummean<=true$FDR)[,x],x]] <- T
    tmp.ind
  })
  
  can_include.idx<-de.post_prob.mt!=0
  de.gene.FDR.idx<-tmp.de.gene.FDR.idx & can_include.idx
  # de.gene.FDR.idx<-tmp.de.gene.FDR.idx
  de.gene.FDR.idx.ls[it]<-list(de.gene.FDR.idx)
  
  prop_de.mt <- cbind(prop_de.mt, colMeans(de.gene.FDR.idx))
  
  message(paste("cycle =",it,"finished @",Sys.time()))
}
# setwd("/home/cgz59/methTCGA")
dir.create("./Results_allCases_expandedThreshold")
save(mean.discount.v,discount.q1.v,discount.q2.v,eta.q1.v,eta.q2.v,eta.0.prop.v,mean.eta.v,median.eta.v,
     lowerBound_PDP_BF.v,se_LowerBound_PDP_BF.v,eta.0_BF.v,se_eta.0_BF.v,order_1.mh.rate,
     de.gene.idx.ls,de.post_prob.ls,de.gene.FDR.idx.ls, prop_de.mt,
     probe_level_ordered_diff_list.ls,gene_level_max_grp_diff.ls,file=paste('./Results_allCases_expandedThreshold/results_50kOut_pt',slurm_arrayid,'.rda',sep=''))

# save(mean.discount.v,discount.q1.v,discount.q2.v,eta.q1.v,eta.q2.v,eta.0.prop.v,mean.eta.v,median.eta.v,
#      lowerBound_PDP_BF.v,se_LowerBound_PDP_BF.v,eta.0_BF.v,se_eta.0_BF.v,order_1.mh.rate,
#      de.gene.idx.ls,de.post_prob.ls,de.gene.FDR.idx.ls, prop_de.mt,
#      probe_level_ordered_diff_list.ls,gene_level_max_grp_diff.ls,file=paste('./Results_newAllCases/results_50kOut_rerun359-360','.rda',sep=''))
