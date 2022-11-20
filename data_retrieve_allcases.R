source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")
biocLite("SummarizedExperiment")
biocLite("biomaRt")
# biocLite("Homo.sapiens")
# biocLite("GenomicRanges")
# library(stringi)
# library(digest)
library(TCGAbiolinks)
# library(dplyr)
# library(DT)
library(SummarizedExperiment)
library(biomaRt)
# library(Homo.sapiens)
# library(GenomicRanges)
options(error=recover)
options(warn=2)

##########start##########
query_meth.STAD <- GDCquery(project= "TCGA-STAD", 
                       data.category = "DNA methylation", 
                       platform = "Illumina Human Methylation 450", 
                       experimental.strategy = 'Methylation array',
                       #barcode = c("TCGA-HT-8111-01A-11D-2399-05","TCGA-HT-A5R5-01A-11D-A28N-05"), 
                       legacy = TRUE)
STAD.cases<-query_meth.STAD$results[[1]]$cases
query_meth.LIHC <- GDCquery(project= "TCGA-LIHC", 
                            data.category = "DNA methylation", 
                            platform = "Illumina Human Methylation 450", 
                            experimental.strategy = 'Methylation array',
                            #barcode = c("TCGA-HT-8111-01A-11D-2399-05","TCGA-HT-A5R5-01A-11D-A28N-05"), 
                            legacy = TRUE)
LIHC.cases<-query_meth.LIHC$results[[1]]$cases
query_meth.ESCA <- GDCquery(project= "TCGA-ESCA", 
                            data.category = "DNA methylation", 
                            platform = "Illumina Human Methylation 450", 
                            experimental.strategy = 'Methylation array',
                            #barcode = c("TCGA-HT-8111-01A-11D-2399-05","TCGA-HT-A5R5-01A-11D-A28N-05"), 
                            legacy = TRUE)
ESCA.cases<-query_meth.ESCA$results[[1]]$cases
query_meth.PAAD <- GDCquery(project= "TCGA-PAAD", 
                            data.category = "DNA methylation", 
                            platform = "Illumina Human Methylation 450", 
                            experimental.strategy = 'Methylation array',
                            #barcode = c("TCGA-HT-8111-01A-11D-2399-05","TCGA-HT-A5R5-01A-11D-A28N-05"), 
                            legacy = TRUE)
PAAD.cases<-query_meth.PAAD$results[[1]]$cases

####break up by cancer type due to capacity####
GDCdownload(query_meth.STAD)
data_meth.STAD<- GDCprepare(query_meth.STAD)
save(data_meth.STAD,file = 'meth_dataframe_STAD.rda')

GDCdownload(query_meth.LIHC)
data_meth.LIHC<- GDCprepare(query_meth.LIHC)
save(data_meth.LIHC,file = 'meth_dataframe_LIHC.rda')

GDCdownload(query_meth.ESCA)
data_meth.ESCA<- GDCprepare(query_meth.ESCA)
save(data_meth.ESCA,file = 'meth_dataframe_ESCA.rda')

GDCdownload(query_meth.PAAD)
data_meth.PAAD<- GDCprepare(query_meth.PAAD)
save(data_meth.PAAD,file = 'meth_dataframe_PAAD.rda')

# load("./meth_dataframe_STAD.rda")
mt_meth.STAD<-assay(data_meth.STAD)
save(mt_meth.STAD,file = 'meth_datamt_STAD.rda')
# load("./meth_datamt_STAD.rda")

# load("./meth_dataframe_LIHC.rda")
mt_meth.LIHC<-assay(data_meth.LIHC)
save(mt_meth.LIHC,file = 'meth_datamt_LIHC.rda')
# load("./meth_datamt_LIHC.rda")

# load("./meth_dataframe_ESCA.rda")
mt_meth.ESCA<-assay(data_meth.ESCA)
save(mt_meth.ESCA,file = 'meth_datamt_ESCA.rda')
# load("./meth_datamt_ESCA.rda")

# load("./meth_dataframe_PAAD.rda")
mt_meth.PAAD<-assay(data_meth.PAAD)
save(mt_meth.PAAD,file = 'meth_datamt_PAAD.rda')
# load("./meth_datamt_PAAD.rda")


###get methylation data matrix
# mt_meth.all<-cbind(mt_meth.STAD, mt_meth.LIHC, mt_meth.ESCA, mt_meth.PAAD)

