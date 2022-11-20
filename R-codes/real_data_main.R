source("R-codes/iterations_submission.R")
source("R-codes/elementwise_main_submission.R")
source("R-codes/elementwise_DP.functions_submission.R")
source("R-codes/fast_PDP.functions_submission.R")
source("R-codes/e_DP.functions_submission.R")

analyzeBayesdiff <- function(data_touse, n.burn = 1000, n.reps = 2000, FDR = 0.05) {
  true<-NULL
  true$M<-20 #alpha_2,DP mass parameter(elementwise,phi's)
  #true$a.R<-11
  true$b1<-20 #alpha_1,PDP mass parameter
  true$Me<-10 #alpha_3,epsilon_i DPmass parameter
  true$shift <- 1e-4
  true$FDR<-FDR #desired FDR level
  # n.burn <- 1000
  # n.reps <- 1000
  
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
  
  all.stuff <- fn.mcmc_notaxi(text="BayesDiff MCMC...",
                              true, data,
                              n.burn, n.reps, 
                              max.row.nbhd.size = nrow(data$X)*data$G.max/100, 
                              max.col.nbhd.size = round(ncol(data$X)/5),
                              row.frac.probes = 0.25,
                              col.frac.probes = .5,
                              prob.compute.col.nbhd=.2,
                              dahl.flag=FALSE,
                              standardize.X=FALSE,
                              flip.sign=FALSE,
                              tBB_flag=FALSE,  computeMode= createComputeMode())
  
  all_results <- list()
  
  all_results$model_parm_properties$mean.discount.v <- mean(all.stuff$d.v)
  all_results$model_parm_properties$discount.q1.v <- quantile(all.stuff$d.v,probs=c(0.05,0.95))[[1]]
  all_results$model_parm_properties$discount.q2.v <- quantile(all.stuff$d.v,probs=c(0.05,0.95))[[2]]
  all_results$model_parm_properties$order_1.mh.rate <- mean(all.stuff$order_1.mh.flip.v)
  all_results$model_parm_properties$mean.eta.v <- mean(all.stuff$eta.v)
  all_results$model_parm_properties$median.eta.v <- median(all.stuff$eta.v)
  all_results$model_parm_properties$eta.0.prop.v <- mean(all.stuff$eta.v==0)
  all_results$model_parm_properties$eta.q1.v <- quantile(all.stuff$eta.v,probs=c(0.05,0.95))[[1]]
  all_results$model_parm_properties$eta.q2.v <- quantile(all.stuff$eta.v,probs=c(0.05,0.95))[[2]]
  all_results$model_parm_properties$lowerBound_PDP_BF.v <- mean(all.stuff$PDP_log.BF.v)
  all_results$model_parm_properties$se_LowerBound_PDP_BF.v <- sd(all.stuff$PDP_log.BF.v)/sqrt(n.reps)
  all_results$model_parm_properties$eta.0_BF.v <- mean(all.stuff$eta.0_log.BF.v)
  all_results$model_parm_properties$se_eta.0_BF.v <- sd(all.stuff$eta.0_log.BF.v)/sqrt(n.reps)
  
  #get which sample group is largest/smallest for each mcmc sample 
  min_grp.mt <- apply(all.stuff$theta.mt, c(2,3), function(x) which(x == min(x)))
  max_grp.mt <- apply(all.stuff$theta.mt, c(2,3), function(x) which(x == max(x)))
  
  ##setting up some threshold for being dm in individual mcmc sample, get post_prob of dm and which pair gives largest difference
  threshold_list <- seq(from = 0.5, to = 1.5, by = 0.05)
  
  dm.post_prob.mt <- NULL
  probe_level_ordered_diff.mt <- matrix(NA,ncol(data$X), length(threshold_list))
  pair_mt<-cbind(combn(1:data$tot.trt,2),combn(1:data$tot.trt,2)[2:1,])
  
  for(tt in 1:length(threshold_list)) {
    threshold <- threshold_list[tt]
    dm.ind_post.mt <- apply(all.stuff$theta.mt, c(2,3), function(x) max(x)-min(x)>threshold) & all.stuff$clust$gr.mt==2
    dm.post_prob.mt <- cbind(dm.post_prob.mt,
                             rowMeans(dm.ind_post.mt))
    
    for (pp in 1:nrow(dm.ind_post.mt)) {
      
      min.pp_lst <- min_grp.mt[pp,dm.ind_post.mt[pp,]==T]
      max.pp_lst <- max_grp.mt[pp,dm.ind_post.mt[pp,]==T]
      
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
  all_results$probe_level_ordered_diff_list<-(probe_level_ordered_diff.mt[,5]) #probe level max pairwise group difference for gene tt (vector of length num_probe[tt])
  all_results$gene_level_max_grp_diff<-(gene_level_max_grp_diff.v[5]) #gene level max pairwise group difference for gene tt
  all_results$dm.post_prob<-(dm.post_prob.mt[,5])
  
  dm.post_prob.sorted<-apply(dm.post_prob.mt, 2, sort, decreasing = T)
  dm.post_prob.order<-apply(dm.post_prob.mt, 2, order, decreasing = T)
  # de.post_prob.rank<-apply(dm.post_prob.mt, 2, rank)
  dm.post_prob.cummean<-apply(dm.post_prob.sorted, 2, function(x) 
    cumsum(1-x)/(1:nrow(dm.post_prob.mt)))
  # tmp <- dm.post_prob.cummean<=true$FDR
  tmptmp.dm.probe.FDR.idx <- matrix(F, nrow = nrow(dm.post_prob.order), ncol = ncol(dm.post_prob.order))
  tmp.dm.probe.FDR.idx <-sapply(1:ncol(dm.post_prob.order), function(x) {
    tmp.ind <- tmptmp.dm.probe.FDR.idx[,x]
    tmp.ind[dm.post_prob.order[(dm.post_prob.cummean<=true$FDR)[,x],x]] <- T
    tmp.ind
  })
  
  can_include.idx<-dm.post_prob.mt!=0
  dm.probe.FDR.idx<-tmp.dm.probe.FDR.idx & can_include.idx
  # dm.probe.FDR.idx<-tmp.dm.probe.FDR.idx
  all_results$dm.probe.FDR.idx<-(dm.probe.FDR.idx[,5])
  
  all_results$prop_dm <- (colMeans(dm.probe.FDR.idx)[5])
  
  return(all_results)
}
