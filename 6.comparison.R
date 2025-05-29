##megre model list and test index list with published gene list
all_test_list<-c()
indir=paste0("../E2F_DFS/4.test/")
files<-list.files(indir,recursive = T,pattern = "test_index.Rdata",full.names = T)
files
for (i in 1:length(files)){
  load(files[i])
  g.name=basename(dirname(files[i]))
  all_test_list[[i]]<-model_auc_list
  names(all_test_list)[[i]]<-g.name
}
save(all_test_list,file=paste0(indir,"/all_test_list.rdata"))

##comparison with published signatures
load("../E2F_DFS/4.test/all_test_list.rdata")
own_auc_list<-all_test_list[c(2,3,4)]
published_auc_list<-all_test_list[c(1,5)]
indexs=c("Cindex","AUC_1","AUC_2","AUC_3","AUC_5","AUC_7")
dataset<-names(own_auc_list[[1]][[1]])
outdir=paste0("../E2F_DFS/6.comparison/")
for (index in indexs){
  model_name="RFRSF"
  dir.create(paste0(outdir,"/",model_name),recursive = T)
  p<-index_comp(own_auc_list =own_auc_list,
                published_auc_list=published_auc_list,
                model_name=model_name,  # color value for cohort
                dataset=dataset, # input datasets name
                index=index)
  ggsave(p,file=paste0(outdir,"/",model_name,"/",index,".jpg"),dpi=600,height =6,width =30,units="cm")
}

# Calculate the significance of c-index
result_list<-indexC_comp(own_auc_list =all_auc_list,
                         published_auc_list=other_auc_list,
                         model_name=model_name,  # color value for cohort
                         dataset=dataset)
write.table(result_list,file=paste0(outdir,"/",model_name,"/",index,"_compare_sig.csv"),sep=",")
