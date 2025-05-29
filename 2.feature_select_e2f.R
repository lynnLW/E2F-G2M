###feature selection
###training data
library(AutoML)
load("../../sva/nocombat/list_data.Rdata")
train_data<-list_data$TCGA[[5]]
###geneset list
load("../../../geneset/glist3.Rdata")
for (i in 1:3){
  genelist<-glist[[i]]
  g.name<-names(glist)[[i]]
  ####parameter setting
  selected.feature<-feature_selection(train_data,
                                      genelist=genelist,
                                      outdir=paste0("../E2F_DFS/1.feature_selection/",g.name))
}

