#########
load("../../../geneset/glist3.Rdata")
load("../dataset/set7/scale/list_vali_data.Rdata")
load("../dataset/set7/scale/list_vali_meta.Rdata")
for (i in 1:3){
    g.name=names(glist)[i]
    indir=paste0("../E2F_DFS/4.test/",g.name)
    out=paste0("../E2F_DFS/5.multicox/",g.name)
    load(paste0(indir,"/test_index.Rdata"))
    #model='SuperPC'
    model_name='RFRSF'
    rs_list=lapply(model_auc_list[[model_name]],function(x){x[[1]]$pred_df})
    outdir=paste0(out,"/",model_name,"/",g.name)
    ####
    cat_opti_summary<-c()
    for(n in 1:length(rs_list)){
      # risk table
      dataset_name<-names(rs_list)[n]
      print(dataset_name)
      rs<-rs_list[[n]]
      names(rs)[3]<-g.name
      # meta table
      meta<-list_vali_meta[[dataset_name]]
      # multicox
      features<-colnames(meta)
      features<-intersect(c("Age","PSA","GS","path_T"),features)
      meta<-meta[,features]
      rs<-rs[row.names(meta),]
      meta<-meta[row.names(rs),]
      combined_meta<-merge(meta,rs,by="row.names")
      ##multicox for continuous rs
      dir.create(paste0(outdir,"/best_cutoff"),recursive = T)
      sum.cox<-generate_multicox_analysis(
        data = combined_meta,
        features = features,
        gene = g.name,
        dataset_name = dataset_name,
        outdir = paste0(outdir,"/best_cutoff"),
        cut_type = NULL,
        plot_height = 4
      )
      cat_opti_summary<-rbind(cat_opti_summary,sum.cox)
  }
  write.table(cat_opti_summary,file=paste0(outdir,"/best_cutoff/cat_opti_summary.csv"),sep=",")
}
