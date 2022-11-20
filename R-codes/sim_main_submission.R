library(MASS)
library(stats)
library(ggplot2)
source("https://bioconductor.org/biocLite.R")
# biocLite("BiSeq")
biocLite("bsseq")
library(bsseq)
if (!"BiocManager" %in% installed.packages()[,1]){
  install.packages("BiocManager")
}
BiocManager::install("BiSeq")
library(BiSeq)
library(abind)
library(lmtest)
library(betareg)

source("R-codes/iterations_submission.R")
source("R-codes/elementwise_main_submission.R")
source("R-codes/elementwise_DP.functions_submission.R")
source("R-codes/fast_PDP.functions_submission.R")
source("R-codes/e_DP.functions_submission.R")
source("R-codes/gen.data.submission.R")

load('probe_distances_allgenes.rda')

calcavgAUC<-function (FPR_roc.mt,TPR_roc.mt) {
  FPR_roc_avg<-rowMeans(FPR_roc.mt)
  TPR_roc_avg<-rowMeans(TPR_roc.mt)
  FPR_roc_avg<-sort(FPR_roc_avg, decreasing = T)
  TPR_roc_avg<-sort(TPR_roc_avg, decreasing = T)
  FPR_roc_avg_20<-FPR_roc_avg[FPR_roc_avg<=0.2]
  TPR_roc_avg_20<-TPR_roc_avg[FPR_roc_avg<=0.2]
  tmp_ratio_20<-(0.2-max(FPR_roc_avg[FPR_roc_avg<=0.2]))/(min(FPR_roc_avg[FPR_roc_avg>0.2])-max(FPR_roc_avg[FPR_roc_avg<=0.2]))
  FPR_roc_avg_10<-FPR_roc_avg[FPR_roc_avg<=0.1]
  TPR_roc_avg_10<-TPR_roc_avg[FPR_roc_avg<=0.1]
  tmp_ratio_10<-(0.1-max(FPR_roc_avg[FPR_roc_avg<=0.1]))/(min(FPR_roc_avg[FPR_roc_avg>0.1])-max(FPR_roc_avg[FPR_roc_avg<=0.1]))
  
  mbase<-rowMeans(cbind(TPR_roc_avg[-1],TPR_roc_avg[-length(TPR_roc_avg)]))
  AUC<-sum(-diff(FPR_roc_avg)*mbase) #AUC
  mbase_20<-rowMeans(cbind(TPR_roc_avg_20[-1],TPR_roc_avg_20[-length(TPR_roc_avg_20)]))
  AUC20<-5*(sum(-diff(FPR_roc_avg_20)*mbase_20)+(max(TPR_roc_avg[FPR_roc_avg<=0.2])+(min(TPR_roc_avg[FPR_roc_avg>0.2])-max(TPR_roc_avg[FPR_roc_avg<=0.2]))*tmp_ratio_20*0.5)*(0.2-max(FPR_roc_avg[FPR_roc_avg<=0.2]))) #AUC20*5
  mbase_10<-rowMeans(cbind(TPR_roc_avg_10[-1],TPR_roc_avg_10[-length(TPR_roc_avg_10)]))
  AUC10<-10*(sum(-diff(FPR_roc_avg_10)*mbase_10)+(max(TPR_roc_avg[FPR_roc_avg<=0.1])+(min(TPR_roc_avg[FPR_roc_avg>0.1])-max(TPR_roc_avg[FPR_roc_avg<=0.1]))*tmp_ratio_20*0.5)*(0.1-max(FPR_roc_avg[FPR_roc_avg<=0.1]))) #AUC10*10
  
  return(list(AUC=AUC,AUC20=AUC20,AUC10=AUC10))
}

simBayesdiff <- function(Sigma=1,eta=0.004,replic=20, seedAssigned=1107) {
  #Sigma=1.0 high noise -> SNR ~ 0.4
  #Sigma=0.6 low noise -> SNR ~ 0.7
  true<-NULL
  true$M<-20 
  true$b1<-20 
  true$Me<-10 
  true$shift <- 1e-4
  
  p<-500
  rho<-0.1
  true$tau <- rep(Sigma,p)
  true$rho<-c(1-rho,rho)
  true$eta<-eta
  
  group<-rep(1:5,each=4)
  trt<-length(unique(group))
  LANE<-rep(0,length(group))
  FDR_nominal <- seq(0.01,0.99,by=0.01)
  
  # For BiSeq
  metadata <- list(Sequencer = "Instrument", Year = "2017")
  colData <- DataFrame(group = factor(group),
                       row.names = 1:length(group))
  # End
  
  n.burn <- 2500
  n.reps <- 5000
  
  eta.0_BF.v<-se_eta.0_BF.v<-NULL
  data_X_array<-array(NA,c(length(group),p,replic))
  true_distance.mt<-true_dm.probe.ind.mt<-NULL
  dm.probe.post.prob.mt <- FDR_actual.mt <- NULL
  FPR_roc.mt<-TPR_roc.mt<-NULL
  FPR_roc_naive.mt<-TPR_roc_naive.mt<-FPR_roc_naive_ks.mt<-TPR_roc_naive_ks.mt<-NULL
  FPR_roc_BiSeq.mt<-TPR_roc_BiSeq.mt<-NULL
  FPR_roc_COHCAP.mt<-TPR_roc_COHCAP.mt<-NULL
  FPR_roc_RADMeth.mt<-TPR_roc_RADMeth.mt<-NULL
  FPR_roc_Methylkit.mt<-TPR_roc_Methylkit.mt<-NULL
  
  
  for (tt in 1:replic) {
    set.seed(seedAssigned+tt)
    print(paste("BayesDiff simulation, replication",tt,'start',date(),"******"))
    tmp.distance_random.gene<-unlist(distances.allgenes)
    tmp.distance<-sample(tmp.distance_random.gene,size=p-1,replace=F)
    true$distance<-tmp.distance/(sum(tmp.distance))
    true_distance.mt<-cbind(true_distance.mt,true$distance)
    num.dm<-0
    eta_lik.v<-0
    eta_tt<-seq(0,0.015,by=0.0001)
    tmpSNR <- 0
    while (num.dm>p*0.15|num.dm<p*0.05|eta_tt[which.max(eta_lik.v)]>(eta+0.001)*as.numeric(eta!=0)|eta_tt[which.max(eta_lik.v)]<(eta-0.001)){
      gen<-gen.data(trt,group,p,true,LANE)
      data<-gen$data
      true_parm<-gen$true_parm
      data$X<-data$X.raw
      data$rho<-rho
      data$delta_margin<-true$delta_margin
      data$distance<-true$distance
      num.dm<-true_parm$num.dm
      
      eta_lik.v<-NULL
      for (ttt in 1:length(eta_tt)){
        eta_lik.v<-c(eta_lik.v,PDP_fn.eta.log.lik(eta_tt[ttt],true_parm)[[1]])
      }
      eta_lik.v<-unlist(eta_lik.v)
      
      tmpSNR <- 1-(var(as.numeric(data$X.raw-data$mean.all))/var(as.numeric(t(data$X.raw)-(data$mean.chi))))
    }
    print(paste0("SNR=",tmpSNR))
    data_X_array[,,tt]<-data$X
    true_dm.probe.ind<-true_parm$clust$gr.v==2
    true_dm.probe.ind.mt<-cbind(true_dm.probe.ind.mt,true_dm.probe.ind)
    
    data$X <- t(t(data$X)-data$mean.chi) ##try this quick way first
    
    all.stuff <- fn.mcmc(text="BayesDiff...",
                         true, data,
                         n.burn, n.reps, max.row.nbhd.size = 20,
                         max.col.nbhd.size = 25,
                         row.frac.probes = 0.5,
                         col.frac.probes = .25,
                         prob.compute.col.nbhd=.2,
                         true_parm,
                         dahl.flag=FALSE,
                         standardize.X=FALSE,
                         flip.sign=FALSE,
                         tBB_flag=FALSE,  computeMode= createComputeMode())
    
    eta.0_BF.v<-c(eta.0_BF.v,mean(all.stuff$eta.0_log.BF.v))
    se_eta.0_BF.v<-c(se_eta.0_BF.v,sd(all.stuff$eta.0_log.BF.v)/sqrt(n.reps))
    
    dm.probe.post.prob <- rowMeans(all.stuff$clust$gr.mt==2)
    dm.probe.post.prob.mt <- cbind(dm.probe.post.prob.mt, dm.probe.post.prob)
    
    dm.post_prob.sorted<-sort(dm.probe.post.prob,decreasing = T)
    dm.post_prob.order<-order(dm.probe.post.prob,decreasing = T)
    dm.post_prob.cummean<-cumsum(1-dm.post_prob.sorted)/(1:length(dm.probe.post.prob))
    dm.probe.FDR.ind.mt <- sapply(FDR_nominal, function(xx) {
      if (sum(dm.post_prob.cummean<=xx)>0) {
        1:ncol(data$X) %in% dm.post_prob.order[1:max(which(dm.post_prob.cummean<=xx))]
      } else {
        rep(FALSE, ncol(data$X))
      }
    })
    
    FDR_actual <- apply(dm.probe.FDR.ind.mt, 2, function(yy) {
      if (sum(yy==T)>0) {
        sum(true_dm.probe.ind[yy==T]==F)/sum(yy==T)
      } else {
        0
      }
    })
    FDR_actual.mt <- cbind(FDR_actual.mt, FDR_actual)
    
    dm.probe.ind_roc<-sapply(seq(0.000,1.001,by=0.001),'<=',rowMeans(all.stuff$clust$gr.mt==2))
    FPR_roc<-colSums(dm.probe.ind_roc[true_dm.probe.ind==F,]==T)/sum(true_dm.probe.ind==F)
    TPR_roc<-colSums(dm.probe.ind_roc[true_dm.probe.ind==T,]==T)/sum(true_dm.probe.ind==T)
    FPR_roc.mt<-cbind(FPR_roc.mt,FPR_roc)
    TPR_roc.mt<-cbind(TPR_roc.mt,TPR_roc)
    
    print(paste("BayesDiff simulation, replication",tt,'complete',date(),"******"))
  }
  
  print(paste("simulation, alternative methods, start",date(),"******"))
  # Alternative methods
  for (tt in 1:replic) {
    tmp.data<-data_X_array[,,tt]
    data_adj<-tmp.data-(rowMeans(tmp.data)-(mean(tmp.data)))
    
    # Start ANOVA & Kruskal-Wallis analysis
    p_val.v<-NULL
    for (jj in 1:ncol(tmp.data)) {
      naive_fit<-lm(data_adj[,jj]~factor(group))
      p_val.v<-c(p_val.v,anova(naive_fit)$"Pr(>F)"[1])
    }
    p_val.v <- p.adjust(p_val.v, method = "BH")
    
    p_val.v_ks<-NULL
    for (jj in 1:ncol(tmp.data)) {
      p_val.v_ks<-c(p_val.v_ks,kruskal.test(data_adj[,jj]~factor(group))$p.value)
    }
    p_val.v_ks <- p.adjust(p_val.v_ks, method = "BH")
    
    tmp_roc<-sapply(seq(0.000,1,by=0.001),'>',p_val.v)
    FPR_roc_naive<-colSums((tmp_roc)[true_dm.probe.ind.mt[,tt]==F,]==T)/sum(true_dm.probe.ind.mt[,tt]==F) #FPR
    TPR_roc_naive<-colSums((tmp_roc)[true_dm.probe.ind.mt[,tt]==T,]==T)/sum(true_dm.probe.ind.mt[,tt]==T) #TPR
    FPR_roc_naive.mt<-cbind(FPR_roc_naive.mt,FPR_roc_naive)
    TPR_roc_naive.mt<-cbind(TPR_roc_naive.mt,TPR_roc_naive)
    
    tmp_roc_ks<-sapply(seq(0.000,1,by=0.001),'>',p_val.v_ks)
    FPR_roc_naive_ks<-colSums((tmp_roc_ks)[true_dm.probe.ind.mt[,tt]==F,]==T)/sum(true_dm.probe.ind.mt[,tt]==F) #FPR
    TPR_roc_naive_ks<-colSums((tmp_roc_ks)[true_dm.probe.ind.mt[,tt]==T,]==T)/sum(true_dm.probe.ind.mt[,tt]==T) #TPR
    FPR_roc_naive_ks.mt<-cbind(FPR_roc_naive_ks.mt,FPR_roc_naive_ks)
    TPR_roc_naive_ks.mt<-cbind(TPR_roc_naive_ks.mt,TPR_roc_naive_ks)
    # End ANOVA & Kruskal-Wallis analysis
    
    # Start BiSeq analysis
    rowRanges <- GRanges(seqnames = "chr1",
                         ranges = IRanges(start = cumsum(c(1,true_distance.mt[,tt]/min(true_distance.mt[,tt]))), end = cumsum(c(1,true_distance.mt[,tt]/min(true_distance.mt[,tt])))))
    methLevel <- t(exp(data_adj)/(1+exp(data_adj)))
    data_Obj<-BSrel(metadata = metadata,
                    rowRanges = rowRanges,
                    colData = colData,
                    methLevel = methLevel)
    pred.meth.part <- methLevel(data_Obj)
    p.val <- rep(NA,nrow(data_Obj))
    min.meth <- min(methLevel(data_Obj)[methLevel(data_Obj) > 0], na.rm=TRUE)
    max.meth <- max(methLevel(data_Obj)[methLevel(data_Obj) < 1], na.rm=TRUE)
    
    for(j in 1:nrow(pred.meth.part)){
      pred.meth <- pred.meth.part[j,]
      pred.meth[pred.meth == 0] <- min.meth
      pred.meth[pred.meth == 1] <- max.meth
      data_1probe <- cbind(pred.meth = pred.meth,
                           as.data.frame(colData(data_Obj)))
      options(show.error.messages = FALSE)
      suppressWarnings(
        lmodel <- try((betareg(formula = pred.meth ~ group,data=data_1probe, link = 'logit')), silent=TRUE)
      )
      suppressWarnings(
        lmodel_red <- try((betareg(formula = pred.meth ~ 1,data=data_1probe, link = 'logit')), silent=TRUE)
      )
      options(show.error.messages = TRUE)
      
      if((class(lmodel) == "try-error" | class(lmodel_red) == "try-error")){
        # p.val[j] <- NA
        stop(paste('BiSeq error in lmodel',j))
      } else{
        if(!lmodel$converged | !lmodel_red$converged){
          # p.val[j] <- NA
          stop(paste('error in lmodel$converged',j))
        } else{
          pVal<-waldtest(lmodel_red,lmodel)$`Pr(>Chisq)`[2]
          # Test auf 0: 2 * pnorm(-abs(Estimate / Std.Error))
          # Test auf min.diff > 0: 2 * min(0.5, pnorm( -(abs(Estimate)-min.diff)/Std.Error, mean=-min.diff))
          p.val[j] <- (pVal) # should not be zero
        }
      }
    }
    
    p_val.BiSeq<-p.adjust(p.val, method = "BH")
    tmp_roc<-sapply(seq(0.000,1,by=0.001),'>',p_val.BiSeq)
    FPR_roc_BiSeq<-colSums((tmp_roc)[true_dm.probe.ind.mt[,tt]==F,]==T)/sum(true_dm.probe.ind.mt[,tt]==F) #FPR
    TPR_roc_BiSeq<-colSums((tmp_roc)[true_dm.probe.ind.mt[,tt]==T,]==T)/sum(true_dm.probe.ind.mt[,tt]==T) #TPR
    FPR_roc_BiSeq.mt<-cbind(FPR_roc_BiSeq.mt,FPR_roc_BiSeq)
    TPR_roc_BiSeq.mt<-cbind(TPR_roc_BiSeq.mt,TPR_roc_BiSeq)
    # End BiSeq analysis
    
    # Start COHCAP analysis
    methLevel <- t(exp(data_adj)/(1+exp(data_adj)))
    p.val <- rep(NA,nrow(methLevel))
    for(j in 1:nrow(methLevel)){
      pred.meth <- methLevel[j,]
      tmp_fit<-lm(pred.meth~factor(group))
      p.val[j]<-anova(tmp_fit)$"Pr(>F)"[1]
      
    }
    p_val.COHCAP<-p.adjust(p.val, method = "BH")
    tmp_roc<-sapply(seq(0.000,1,by=0.001),'>',p_val.COHCAP)
    FPR_roc_COHCAP<-colSums((tmp_roc)[true_dm.probe.ind.mt[,tt]==F,]==T)/sum(true_dm.probe.ind.mt[,tt]==F) #FPR
    TPR_roc_COHCAP<-colSums((tmp_roc)[true_dm.probe.ind.mt[,tt]==T,]==T)/sum(true_dm.probe.ind.mt[,tt]==T) #TPR
    FPR_roc_COHCAP.mt<-cbind(FPR_roc_COHCAP.mt,FPR_roc_COHCAP)
    TPR_roc_COHCAP.mt<-cbind(TPR_roc_COHCAP.mt,TPR_roc_COHCAP)
    # End COHCAP analysis
    
    # Start RADMeth analysis
    methLevel <- t(exp(data_adj)/(1+exp(data_adj)))
    p.val <- rep(NA,nrow(methLevel))
    min.meth <- min(methLevel[methLevel> 0], na.rm=TRUE)
    max.meth <- max(methLevel[methLevel < 1], na.rm=TRUE)
    
    for(j in 1:nrow(methLevel)){
      pred.meth <- methLevel[j,]
      pred.meth[pred.meth == 0] <- min.meth
      pred.meth[pred.meth == 1] <- max.meth
      data_1probe <- cbind(pred.meth = pred.meth,
                           as.data.frame(colData))
      options(show.error.messages = FALSE)
      suppressWarnings(
        lmodel <- try((betareg(formula = pred.meth ~ group,
                               data=data_1probe, link = 'logit')), silent=TRUE)
      )
      suppressWarnings(
        lmodel_red <- try((betareg(formula = pred.meth ~ 1,
                                   data=data_1probe, link = 'logit')), silent=TRUE)
      )
      options(show.error.messages = TRUE)
      if((class(lmodel) == "try-error" | class(lmodel_red) == "try-error")){
        # p.val[j] <- NA
        stop(paste('RADMeth error in lmodel',j))
      } else{
        if(!lmodel$converged | !lmodel_red$converged){
          # p.val[j] <- NA
          stop(paste('error in lmodel$converged',j))
        } else{
          pVal<-lrtest(lmodel_red,lmodel)$`Pr(>Chisq)`[2]
          # Test auf 0: 2 * pnorm(-abs(Estimate / Std.Error))
          # Test auf min.diff > 0: 2 * min(0.5, pnorm( -(abs(Estimate)-min.diff)/Std.Error, mean=-min.diff))
          p.val[j] <- (pVal) # should not be zero
        }
      }
    }
    
    p_val.RADMeth<-p.adjust(p.val, method = "BH")
    tmp_roc<-sapply(seq(0.000,1,by=0.001),'>',p_val.RADMeth)
    FPR_roc_RADMeth<-colSums((tmp_roc)[true_dm.probe.ind.mt[,tt]==F,]==T)/sum(true_dm.probe.ind.mt[,tt]==F) #FPR
    TPR_roc_RADMeth<-colSums((tmp_roc)[true_dm.probe.ind.mt[,tt]==T,]==T)/sum(true_dm.probe.ind.mt[,tt]==T) #TPR
    FPR_roc_RADMeth.mt<-cbind(FPR_roc_RADMeth.mt,FPR_roc_RADMeth)
    TPR_roc_RADMeth.mt<-cbind(TPR_roc_RADMeth.mt,TPR_roc_RADMeth)
    # End RADMeth analysis
    
    # Start Methylkit analysis
    methLevel <- t(exp(data_adj)/(1+exp(data_adj)))
    
    p.val <- rep(NA,nrow(methLevel))
    
    for(j in 1:nrow(methLevel)){
      pred.meth <- methLevel[j,]
      data_1probe <- data.frame(pred.meth = pred.meth,group=factor(group))
      logimodel <- (glm(formula = pred.meth ~ group, #weights = rep(100,length(group)),
                        data=data_1probe, family = quasibinomial))
      phi<-summary(logimodel)$dispersion
      dev_diff<-logimodel$null.deviance-logimodel$deviance
      df_diff<-logimodel$df.null-logimodel$df.residual
      p.val[j]<-pchisq(dev_diff/phi, df_diff, lower.tail = FALSE)
      
    }
    
    p_val.Methylkit<-p.adjust(p.val, method = "BH")
    tmp_roc<-sapply(seq(0.000,1,by=0.001),'>',p_val.Methylkit)
    FPR_roc_Methylkit<-colSums((tmp_roc)[true_dm.probe.ind.mt[,tt]==F,]==T)/sum(true_dm.probe.ind.mt[,tt]==F) #FPR
    TPR_roc_Methylkit<-colSums((tmp_roc)[true_dm.probe.ind.mt[,tt]==T,]==T)/sum(true_dm.probe.ind.mt[,tt]==T) #TPR
    FPR_roc_Methylkit.mt<-cbind(FPR_roc_Methylkit.mt,FPR_roc_Methylkit)
    TPR_roc_Methylkit.mt<-cbind(TPR_roc_Methylkit.mt,TPR_roc_Methylkit)
    # End Methylkit analysis
    
  } # end alternative methods
  print(paste("simulation, alternative methods, end",date(),"******"))
  
  # Plot ROC curves
  pdf(file=paste('ROC_20_sigma',Sigma,'_eta',eta,"_",format(Sys.time(),"%Y%m%d-%H%M"),".pdf",sep='')
      ,width = 8, height = 6)
  par(mar=c(3,3,1.5,1))
  par(mgp=c(1.5,0.5,0))
  par(mfrow=c(4,5))
  
  for (ss in 1:replic) {
    plot(FPR_roc.mt[,ss],TPR_roc.mt[,ss],type='l',xlim=c(0,1),ylim=c(0,1),xlab='1-Specificity',ylab='Sensitivity',cex.lab=1.5,cex.axis=1.5,lwd=2)
    lines(FPR_roc_naive.mt[,ss],TPR_roc_naive.mt[,ss],lty=2,col='red')
    lines(FPR_roc_naive_ks.mt[,ss],TPR_roc_naive_ks.mt[,ss],lty=3,col='darkgreen')
    lines(FPR_roc_COHCAP.mt[,ss],TPR_roc_COHCAP.mt[,ss],lty=4,col='blue')
    lines(FPR_roc_Methylkit.mt[,ss],TPR_roc_Methylkit.mt[,ss],lty=5,col='darkred')
    lines(FPR_roc_BiSeq.mt[,ss],TPR_roc_BiSeq.mt[,ss],lty=6,col='darkmagenta')
    lines(FPR_roc_RADMeth.mt[,ss],TPR_roc_RADMeth.mt[,ss],lty=1,col='tan3')
  }
  dev.off()
  
  pdf(file=paste('ROC_avg_sigma',Sigma,'_eta',eta,"_",format(Sys.time(),"%Y%m%d-%H%M"),".pdf",sep='')
      ,width = 8, height = 6)
  par(mar=c(5.1, 4.6, 4.1, 2.1),mgp=c(3, 1, 0), las=0)
  
  plot(rowMeans(FPR_roc.mt),rowMeans(TPR_roc.mt),type='l',xlim=c(0,1),ylim=c(0,1),xlab='1-Specificity',ylab='Sensitivity',cex.lab=2,cex.axis=2,lwd=2)
  lines(rowMeans(FPR_roc_naive.mt),rowMeans(TPR_roc_naive.mt),lty=2,col='red',lwd=1.5)
  lines(rowMeans(FPR_roc_naive_ks.mt),rowMeans(TPR_roc_naive_ks.mt),lty=3,col='darkgreen',lwd=1.5)
  lines(rowMeans(FPR_roc_COHCAP.mt),rowMeans(TPR_roc_COHCAP.mt),lty=4,col='blue',lwd=1.5)
  lines(rowMeans(FPR_roc_Methylkit.mt),rowMeans(TPR_roc_Methylkit.mt),lty=5,col='darkred',lwd=1.5)
  lines(rowMeans(FPR_roc_BiSeq.mt),rowMeans(TPR_roc_BiSeq.mt),lty=6,col='darkmagenta',lwd=1.5)
  lines(rowMeans(FPR_roc_RADMeth.mt),rowMeans(TPR_roc_RADMeth.mt),lty=1,col='tan3',lwd=1.5)
  legend(0.55,0.65,c('BayesDiff','ANOVA','Kruskal-Wallis','COHCAP','Methylkit','BiSeq','RADMeth'),bty = "n",lty=c(1,2,3,4,5,6,1),col=c('black','red','darkgreen','blue','darkred','darkmagenta','tan3'),lwd = c(2,rep(1.5,6)),cex=1.7)
  dev.off()
  
  # Calculate average AUCs
  AUCs<-calcavgAUC(FPR_roc.mt,TPR_roc.mt)
  AUCs_naive<-calcavgAUC(FPR_roc_naive.mt,TPR_roc_naive.mt)
  AUCs_naive_ks<-calcavgAUC(FPR_roc_naive_ks.mt,TPR_roc_naive_ks.mt)
  AUCs_COHCAP<-calcavgAUC(FPR_roc_COHCAP.mt,TPR_roc_COHCAP.mt)
  AUCs_Methylkit<-calcavgAUC(FPR_roc_Methylkit.mt,TPR_roc_Methylkit.mt)
  AUCs_BiSeq<-calcavgAUC(FPR_roc_BiSeq.mt,TPR_roc_BiSeq.mt)
  AUCs_RADMeth<-calcavgAUC(FPR_roc_RADMeth.mt,TPR_roc_RADMeth.mt)
  
  return(list(FPR_roc.mt=FPR_roc.mt,TPR_roc.mt=TPR_roc.mt,FPR_roc_naive.mt=FPR_roc_naive.mt,TPR_roc_naive.mt=TPR_roc_naive.mt,
              FPR_roc_naive_ks.mt=FPR_roc_naive_ks.mt,TPR_roc_naive_ks.mt=TPR_roc_naive_ks.mt,FPR_roc_COHCAP.mt=FPR_roc_COHCAP.mt,TPR_roc_COHCAP.mt=TPR_roc_COHCAP.mt,
              FPR_roc_Methylkit.mt=FPR_roc_Methylkit.mt,TPR_roc_Methylkit.mt=TPR_roc_Methylkit.mt,FPR_roc_BiSeq.mt=FPR_roc_BiSeq.mt,TPR_roc_BiSeq.mt=TPR_roc_BiSeq.mt,
              FPR_roc_RADMeth.mt=FPR_roc_RADMeth.mt,TPR_roc_RADMeth.mt=TPR_roc_RADMeth.mt,eta.0_BF.v=eta.0_BF.v,se_eta.0_BF.v=se_eta.0_BF.v,
              AUCs=AUCs,AUCs_ANOVA=AUCs_naive,AUCs_KruskalWallis=AUCs_naive_ks,AUCs_COHCAP=AUCs_COHCAP,AUCs_Methylkit=AUCs_Methylkit,AUCs_BiSeq=AUCs_BiSeq,AUCs_RADMeth=AUCs_RADMeth
  ))
  
}

