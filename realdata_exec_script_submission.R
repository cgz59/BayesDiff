source("R-codes/real_data_main.R")

load('probe_gene_grp_info_50k.rda') ##meta data

gene_subset_idx<-which(global_realdata$num_probe>10)

results_all_genes <- list()
# for (it in 1:2) {
for (it in 1:length(gene_subset_idx)) {
  # set.seed(124)
  message("Loading data...")
  
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
  
  message("Data loaded.")
  
  results_all_genes[[it]] <- analyzeBayesdiff(data_touse, n.burn = 1000, n.reps = 2000, FDR = 0.05)
  
  message(paste("Analysis #",it,"finished @",Sys.time()))
}
