###geneset list
load("../../../geneset/glist3.Rdata")
library(AutoML)
for (i in 1:3){
    g.name<-names(glist)[i]
    ###model assess
    indir=paste0("../E2F_DFS/2.train/",g.name)
    load(paste0("../E2F_DFS/2.train/",g.name,"/10_5_model_list.RData"))

    out2=paste0("../E2F_DFS/3.figure/",g.name)
    ## evaluate cindex
    if(T){
      cindex_list<-lapply(model.list,function(x)x$metrics_list)
      save(cindex_list,file=paste0(indir,"/cindex_list.rdata"))
      cindex_rank2(cindex_list,order="valid",index="all",
                   outdir=out2,plot_type="boxplot")
      model_list<-lapply(model.list,function(x)x$final_model)
      cindex_rank2(model_list=model_list,outdir=out2)
      save(model_list,file=paste0(indir,"/model_list.rdata"))
    }

    ## evaluate performance in external cohorts
    if(T){
      ###features
      load(paste0(indir,"/train_features.Rdata"))
      load("../dataset/set7/scale/list_vali_data.Rdata")
      list_train_vali_Data<-list_vali_data
      out3=paste0("../E2F_DFS/4.test/",g.name)
      test.auc<-cal_vali_index(list_train_vali_Data,candidate_genes,
                               model_list,rep=1,outdir=out3)
      cindex_rank(vali_auc_list = test.auc,index="cindex",
                  plot_type="barplot",outdir=out3)
      cindex_rank(vali_auc_list = test.auc,index="km_auc_1",
                  plot_type="barplot",outdir=out3)
      cindex_rank(vali_auc_list = test.auc,index="km_auc_2",
                  plot_type="barplot",outdir=out3)
      cindex_rank(vali_auc_list = test.auc,index="km_auc_3",
                  plot_type="barplot",outdir=out3)
      cindex_rank(vali_auc_list = test.auc,index="km_auc_5",
                  plot_type="barplot",outdir=out3)
      cindex_rank(vali_auc_list = test.auc,index="km_auc_7",
                  plot_type="barplot",outdir=out3)
      roc_plot(vali_auc_list = test.auc,model="all",outdir=out3)
      surv_plot(vali_auc_list = test.auc,model="all",outdir=out3)
    }
}

