##calculate drug sensitivity
load("../../sva/nocombat/list_data.Rdata")
list_data_surv<-lapply(list_data,function(df){
  df[[1]]
})
###
library(dplyr)
library(limma)
dataset<-c("DKFZ","SU2C","TCGA","GSE46602","GSE70768","MSKCC")
for (i in 1:6){
  dataset_name=dataset[i]
  expr<-list_data_surv[[dataset_name]]
  expr<-normalizeQuantiles(expr)
  cal_drug_sensitive(expr,database ="CTRP2",
                     output_dir = "../E2F_DFS/9.drug/CTRP2/",
                     output_filename =paste0(dataset_name,"_drug_sensitivity"))

}
##calculate correlation and differences
library(dplyr)
load("../E2F_DFS/8.function/group_list.rdata")
drug_file<-list.files("../E2F_DFS/9.drug/CTRP2/",pattern="csv",full.names = T,recursive = T)
for (i in 1:length(group_list)){
  dataset<-names(group_list)[i]
  file<-drug_file[grep(dataset,drug_file)]
  result<-read.csv(file,row.names = 1)
  head(result[1:5,1:5])
  result<-scale(result) %>% as.data.frame()
  rs_df<-group_list[[dataset]]
  rs_df<-rs_df[row.names(result),]
  merge_df<-cbind(rs_df$pred,result)
  names(merge_df)[1]<-"risk_score"
  analyze_drug_correlations(
    data = merge_df,
    risk_score_col = "risk_score",
    r_threshold = 0,
    top_n = 10,
    ncol = 2,
    outdir=paste0("../E2F_DFS/9.drug/CTRP2/",dataset))
  # calculate differences
  cal_diff(result,rs_df$group,outdir=paste0("../E2F_DFS/9.drug/CTRP2/",dataset))
}
