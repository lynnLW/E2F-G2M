# loading external cohort
load("../dataset/set7/scale/list_vali_data.Rdata")
load("../dataset/set7/scale/list_vali_meta.Rdata")
# loading model list
load("../E2F_DFS/2.train/all_model_list.rdata")
# extract RFRSF of each signature from own model list
model_name="RFRSF"
own_model_list<-all_model_list[c(2,3,4)]
own_model_list<-lapply(own_model_list,function(x){x[[model_name]]})
# prepare meta information of training data
train_data<-list_vali_data[['TCGA']]
meta<-list_vali_meta[['TCGA']]
features<-colnames(meta)
features<-intersect(c("Age","PSA","GS"),features) # only continuous features
meta<-meta[,features]
meta<-meta[train_data$ID,]
meta$status<-train_data$OS_status
meta$time<-train_data$OS_time
# loading packages
library(randomForestSRC)
seed=5201314
param=own_model_list$`E2F-G2M`$best_param[[1]]
# export result
cli_model_list <- vector("list", length(features))
names(cli_model_list) <- features
# build model of each clinical features
for(k in 1:length(features)){
  feature=features[k]
  fm=as.formula(paste0("Surv(time, status) ~",feature))
  fit <-rfsrc(fm,
              data = meta,
              ntree = param$ntree,
              nodesize = param$nodesize,
              mtry = param$mtry,
              splitrule = 'logrank',
              importance = TRUE,
              proximity = TRUE,
              forest = TRUE,
              seed = seed)
  cli_model_list[[feature]] <- list(model = fit)
}
cli_model_list<-c(own_model_list,cli_model_list)
dir.create("../E2F_DFS/7.clinical/RFRSF/",recursive = T)
save(cli_model_list,file="../E2F_DFS/7.clinical/RFRSF/cli_model_list.rdata")
#######################################################
# calculate cindex list
load("../E2F_DFS/7.clinical/RFRSF/cli_model_list.rdata")
cli_auc_list<-vector("list", length(list_vali_data))
names(cli_auc_list) <- names(list_vali_data)
for(n in 1:length(list_vali_data)){
  dataset_name<-names(list_vali_meta)[n]
  print(dataset_name)
  dataset<-list_vali_data[[n]]
  names(dataset)[2]<-"time"
  names(dataset)[3]<-"status"
  #
  meta<-list_vali_meta[[n]]
  features<-colnames(meta)
  features<-intersect(c("Age","PSA","GS"),features)
  meta<-meta[,features]
  meta<-meta[dataset$ID,]
  meta$ID<-row.names(meta)
  combined_data<-merge(meta,dataset,by="ID")

  #calculate cindex
  for (i in 1:length(cli_model_list)){
    g.name=names(cli_model_list)[i]
    print(g.name)
    if(g.name %in% c(features,names(cli_model_list)[1:3])){
      model<-cli_model_list[[i]]
      common_feature<-colnames(model$model$xvar)
      common_feature<-c("time","status",common_feature)
      newdata<-combined_data[,common_feature]
      newdata<-na.omit(newdata)
      index<-cal_metrics(newdata,model$model,model_name="RFRSF")
      auc<-cal_multi_auc(newdata,model$model,model_name="RFRSF")
      cli_auc_list[[dataset_name]][[g.name]][[1]]<-c(index,auc)
    }
  }
}
save(cli_auc_list,file="../E2F_DFS/7.clinical/RFRSF/cli_auc_list.rdata")

# loading cindex list
load("../E2F_DFS/7.clinical/RFRSF/cli_auc_list.rdata")
# barplot of all features
out="../E2F_DFS/7.clinical/RFRSF/"
cindex_rank(vali_auc_list = cli_auc_list,
            index="cindex",plot_type="barplot",outdir=out)
cindex_rank(vali_auc_list = cli_auc_list,
            index="km_auc_1",plot_type="barplot",outdir=out)
cindex_rank(vali_auc_list = cli_auc_list,
            index="km_auc_2",plot_type="barplot",outdir=out)
cindex_rank(vali_auc_list = cli_auc_list,
            index="km_auc_3",plot_type="barplot",outdir=out)
cindex_rank(vali_auc_list = cli_auc_list,
            index="km_auc_5",plot_type="barplot",outdir=out)

# ROC curves
out="../E2F_DFS/7.clinical/RFRSF/"
roc_plot(vali_auc_list = cli_auc_list,model=NULL,
         cohort="all",auc_time=1,outdir=out)
roc_plot(vali_auc_list = cli_auc_list,model=NULL,
         cohort="all",auc_time=2,outdir=out)
roc_plot(vali_auc_list = cli_auc_list,model=NULL,
         cohort="all",auc_time=3,outdir=out)
roc_plot(vali_auc_list = cli_auc_list,model=NULL,
         cohort="all",auc_time=5,outdir=out)

