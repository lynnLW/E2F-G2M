library("AutoML")
## loading data list and training data
load("../dataset/set7/scale/list_vali_data.Rdata")
train_data<-list_vali_data[['TCGA']]
## geneset list
load("../../../geneset/glist3.Rdata")
### feature selection
for (i in 1:3){
  #####################
  print("Setp 2 ----final feature selection")
  ## feature selection
  g.name<-names(glist)[[i]]
  indir=paste0("../E2F_DFS/1.feature_selection/",g.name)
  load(paste0(indir,"/5_selected_feature.Rdata"))
  f<-top_feature_select(selected.feature = selected.feature,
                        nmethod = 7,
                        width=7.5,
                        height = 10,
                        outdir=indir)

  ## train models
  print("Setp 3 ----model training")
  candidate_genes<-f
  outdir=paste0("../E2F_DFS/2.train/",g.name)
  model.list<-ML.survival.model(train_data,
                                candidate_genes,
                                filter_OS_time=F,
                                meta_time="m",
                                cor=F,
                                cor_threshold=0.85,
                                fold=5,
                                rep=10,
                                p=0.75,
                                deep_method = F,
                                gbm_method = T, # memory consuming
                                outdir=outdir,
                                seed=5201314,
                                ncore=4)
}

## published signature
load("../../../geneset/ref_glist.Rdata")
### feature selection
for (i in 1:2){
  #####################
  print("Setp 2 ----final feature selection")
  ## feature selection
  g.name<-names(glist2)[[i]]
  genelist<-glist2[[i]]
  ## train models
  print("Setp 3 ----model training")
  candidate_genes<-genelist
  outdir=paste0("../E2F_DFS/2.train/",g.name)
  model.list<-ML.survival.model(train_data,
                                candidate_genes,
                                filter_OS_time=F,
                                meta_time="m",
                                cor=F,
                                cor_threshold=0.85,
                                fold=5,
                                rep=10,
                                p=0.75,
                                deep_method = F,
                                gbm_method = T, # memory consuming
                                outdir=outdir,
                                seed=5201314,
                                ncore=4)
}
