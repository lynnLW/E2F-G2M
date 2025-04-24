# E2F-G2M
`library(AutoML)`
load("../../sva/nocombat/list_data.Rdata")
train_data<-list_data$TCGA[[5]]
###geneset list
load("../../../geneset/glist3.Rdata")
###
source("R/feature_selection.R")
for (i in 1:3){
  genelist<-glist[[i]]
  g.name<-names(glist)[[i]]
  ####parameter setting
  InputMatrix=train_data
  seed = 123
  fold=5
  method="all"
  svm_method=F
  print("Setp 1 ----feature selection method 1")
  outdir=paste0("E2F_DFS/1_feature_selection/scale/",g.name)
  selected.feature<-feature_selection(InputMatrix,
                                       genelist,
                                       seed = 123,
                                       outdir=outdir,
                                       meta_time="none",
                                       fold=5,
                                       unicox_km=T,
                                       deg=T,
                                       up=F,
                                       method="all",
                                       svm_method=F)
}
`
