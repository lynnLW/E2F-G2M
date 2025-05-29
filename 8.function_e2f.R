##loading validation data
load("../dataset/set7/scale/list_vali_data.Rdata")
load("../dataset/set7/scale/list_vali_meta.Rdata")
## separate samples into high and low groups
load("../E2F_DFS/4.test/all_test_list.rdata")
auc_list=all_test_list[['E2F_G2M']]$RFRSF
group_list<-list()
for (i in 1:length(auc_list)){
  dataset<-names(auc_list)[i]
  rs_df<-auc_list[[i]][[1]]$pred_df
  rs_df$group<-ifelse(rs_df$pred>median(rs_df$pred),"High","Low")
  group_list[[dataset]]<-rs_df
}
save(group_list,file="../E2F_DFS/8.function/group_list.rdata")
###function enrichment
library(dplyr)
pathways<-c("TILs","MPI","Immune_cycle","H","Metabolic","Cell_proliferation")
for (id in pathways){
  result_list<-list()
  for (i in 1:length(list_vali_data)){
    expr<-list_vali_data[[i]][,-c(1:3)] %>% t() %>% as.data.frame()
    dataset_name=names(list_vali_data)[i]
    gs<-geneset_cal(expr,category="TILs",prefix=paste0(dataset_name,"_",id),
                    output_dir = paste0("../E2F_DFS/8.function/",id,"/",dataset_name))
    result_list[[dataset_name]]<-gs

    ##calculate correlation
    rs_df<-group_list[[dataset_name]]
    rs_df<-rs_df[row.names(gs),]
    merge_df<-cbind(rs_df$pred,gs)
    names(merge_df)[1]<-"risk_score"
    analyze_drug_correlations(
      data = merge_df,
      risk_score_col = "risk_score",
      top_n = 10,
      ncol = 2,
      outdir=paste0("../E2F_DFS/8.function/",id,"/",dataset_name))

    ##calculate differences
    cal_diff(gs,rs_df$group,outdir=paste0("../E2F_DFS/8.function/",id,"/",dataset_name))
  }
  save(result_list,file=paste0("../E2F_DFS/8.function/",id,"/result_list.rdata"))
}


















