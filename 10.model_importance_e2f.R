load("../E2F_DFS//2.train/all_model_list.rdata")
##gene importance of each signature
for (i in 1:5){
  g.name<-names(all_model_list)[i]
  model_name='RFRSF'
  fit=all_model_list[[g.name]][[model_name]]$model
  generat_rfrsf_importance(fit,g.name,paste0("../E2F_DFS/10.importance/",model_name))
}
