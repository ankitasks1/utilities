#Make venn diagram for multiple vectors
v1 <- c("a","b","c","d")
v2 <- c("a","c")
v3 <- c("b", "a")
v4 <- c("d","a","b")
v5 <- c("d","c")
v6 <- c("c","a")
v7 <- c("d")
v8 <- c("a")

veclist <- list(v1=v1,v2=v2,v3=v3,v4=v4,v5=v5,v6=v6,v7=v7,v8=v8)

code_for_venn <- expand.grid(rep(list(0:1),length(names(veclist))))
colnames(code_for_venn) <- names(veclist)

#Prepare the max combination required
temp50 <- code_for_venn
for (names_i in names(code_for_venn)){
  print(names_i)
  temp50_names_i <- temp50[,names_i]
  #print(head(temp50_names_i))
  temp50[names_i] <- ifelse(temp50_names_i >= 1, names_i, 0)
}

temp50$code <- apply(temp50[ , colnames(temp50)] , 1 , paste , collapse = "_" )

#Prepare dataframe compatible for venn diagrame
# create a master list
master_list_all_degenes <- unique(as.character(unlist(veclist)))
master_sublist_all_degenes <- lapply(veclist, function(x) as.numeric(master_list_all_degenes%in%x))
#Bind the binary matrix
subDF_all_degenes <- do.call(cbind, master_sublist_all_degenes)
rownames(subDF_all_degenes) <- master_list_all_degenes
subDF_all_degenes <- data.frame(subDF_all_degenes)
temp51 <- subDF_all_degenes
for (names_i in colnames(subDF_all_degenes)){
  #print(names_i)
  temp51_names_i <- temp51[,names_i]
  #print(head(temp50_names_i))
  temp51[names_i] <- ifelse(temp51_names_i >= 1, names_i, 0)
}

temp51$code <- apply(temp51[ , colnames(temp51)] , 1 , paste , collapse = "_" )

#Extract combination/venn's
temp52gostlist <- list()
temp52list <- list()
for (code in temp50$code){
  print(code)
  temp52 <- temp51[which(temp51$code == code),]
  temp52list[[code]] <- temp52
}

assign("venn_outlist",temp52list)


