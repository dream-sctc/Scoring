### Define geometry
geo=data.frame(dm@geometry)
colnames(geo)=c("x","y","z")
rownames(geo)=paste0("L_",1:dim(geo)[1])

## Run the code for cv

for(sc in 1:3)
{
  for(cv in 1:10)
  {
    ### location of individual submission
    ### change this accordingly
    all_files <- c(paste0("submissions_cv/BCBU_bin/BCBU_",20*(4-sc),"_",cv,".csv"),
                   paste0("submissions_cv/challengers18/",20*(4-sc),"genes_",cv,".csv"),
                   paste0("submissions_cv/christoph_hafemeister/",20*(4-sc),"genes_",cv,".csv"),
                   paste0("submissions_cv/DeepCMC/",20*(4-sc),"genes_a_",cv,"_test.csv"),
                   paste0("submissions_cv/NAD_Result/",20*(4-sc),"genes_",cv,".csv"),
                   paste0("submissions_cv/OmicsEngineering_new/gene",20*(4-sc),"_testfold",cv,"_index.csv"),
                   paste0("submissions_cv/Team_MLB/",20*(4-sc),"genes_",cv,".csv"),
                   paste0("submissions_cv/thinnguyen_ten_fold_binarized/",20*(4-sc),"genes_",cv,".csv"),
                   paste0("submissions_cv/WhatATeam/",20*(4-sc),"genes_",cv,".csv"),
                   #"submissions_cv/WOC/",20*(4-sc),"_genes_",k,".csv",
                   paste0("submissions_cv/zho_team/",20*(4-sc),"genes_",cv,"_.csv"))
    
    ###### Get kmeans prediction
    d=find_aggregate_predictions_kmeans_cv(all_files,sc)
    
    ##### To have the same format let us read
    ##### a community prediction
    fname=paste0(all_files[5])
    ### woc_kmeans is the woc
    woc_kmeans=read.csv(fname,header=FALSE,stringsAsFactors = FALSE)
    lines_s=2*(4-sc)
    ###### get woc genes
    genes_sc=read.csv(paste0("cv_common_genes_sc_",sc,".csv"),header = FALSE,stringsAsFactors = FALSE)
    ###### replace woc genes
    ###### as well as locations of the community prediction
    woc_kmeans[1:lines_s,]=genes_sc
    woc_kmeans[(lines_s+1):(lines_s+dim(d)[1]),2:11]=d[paste0("C_",woc_kmeans[(lines_s+1):(lines_s+dim(d)[1]),1]),]
    ## Write woc prediction to file
    write.table(woc_kmeans,file=paste0(folder_cv,"/","WOC/",20*(4-sc),"_genes_",cv,".csv"),na=" ",
                col.names=FALSE,row.names = FALSE,sep=",")
    
    
    ################
    rm(fname,woc_kmeans,lines_s,all_files,d)
  }
}




############ functions used for the cv 
############
#############

find_aggregate_predictions_kmeans_cv <- function(all_files,sc)
{
  ### Total number of locations
  locations_all=3039*2
  ## just get number of cells and name of cells read  
  ## one of the community predictions
  xy=get_submission(all_files[3],sc)
  p_matrix=xy$prediction_matrix
  
  all_cells=rownames(p_matrix)
  
  ### define aggregate prediction matrix
  ### for each we have 10 woc location predictions
  
  aggregate_prediction_matrix=matrix(0,length(all_cells),10)
  rownames(aggregate_prediction_matrix)=all_cells
  
  for(cell in all_cells)
  {
    location_prediction_matrix=matrix(0,length(all_files),
                                      locations_all)
    rownames(location_prediction_matrix)=all_files
    ### Not confuse add L_ to locations
    colnames(location_prediction_matrix)=paste0("L_",1:locations_all)
    for(i in all_files)
    {
      xy=get_submission(i,sc)
      p_matrix=xy$prediction_matrix
      locations=p_matrix[cell,]
      #print(colnames(location_prediction_matrix))
      location_prediction_matrix[i,paste0("L_",locations)]=1
    }
    idx=which(rowSums(location_prediction_matrix)!=0)
    location_prediction_matrix=location_prediction_matrix[idx,]
    
    idx=which(colSums(location_prediction_matrix)!=0)
    location_prediction_matrix=location_prediction_matrix[,idx]
    
    d_pred=kmeans_location_arrangment(location_prediction_matrix)
    aggregate_prediction_matrix[cell,]=as.numeric(substr(d_pred,3,nchar(d_pred)))
  }
  
  return(aggregate_prediction_matrix)
}  

get_submission<-function(file_name, sub){
  if (!exists("dm")) initialize()
  
  #read all folds
  
  submission=read.csv(file_name,header=FALSE,stringsAsFactors = FALSE)
  #separate the gene names from the location predictions
  
  gene.lines <- (4-sub)*2
  genes <- submission %>% slice(1:gene.lines)
  locations <- submission %>% slice(-1:-gene.lines)
  
  #preprocess genes and locations, remove NAs, sort locations by cellid
  #mutations cause cancer, but anyway ...
  genes <- genes %>% select(-1) %>% unlist %>% as.character
  
  #fix incompatibility
  genes = gsub("-",".",genes,fixed = T)
  genes = gsub("(spl)",".spl.",genes,fixed = T)
  
  
  
  #### Find no na rows
  
  locations[locations==""]=NA
  idx_na=which(rowSums(is.na(locations))==0)
  
  
  locations <- locations[idx_na,]
  locations <- locations[order(locations[,1]),] 
  cell.ids <- locations[,1]
  loc_matrix=matrix(0,length(idx_na),10)
  
  locations <- locations %>% select(-1) %>% apply(2,as.numeric)
  
  loc_matrix=locations
  rownames(loc_matrix)=paste0("C_",cell.ids)
  
  
  
  # rownames(loc_matrix)=paste0("C_",locations[,1])
  #loc_matrix=loc_matrix[,-1]
  
  return(list(genes=genes,prediction_matrix=loc_matrix))
  
}

kmeans_woc <- function(data,
                        max_centers = 10, nstart = 5, iter.max = 100, 
                        epsilon=0.01,plot=FALSE) {
  
  
  library(ggplot2)
  wss = rep(NA,max_centers)
  set.seed(42)
  for (i in 1:max_centers) 
  {
    
    wss[i] <- sum(kmeans(x = data , iter.max = iter.max , centers = i , nstart = nstart)$withinss)
  }
  
  ## plot
  
  
  ## find best number of clusters 
  ## using elbow method
  wss_delta = rep(0,length(wss)-1)
  for (i in 1:length(wss)-1) 
  {
    wss_delta[i] = (wss[i+1] - wss[i])/wss[i]
  }
  best_n_cluters = min(which(wss_delta > epsilon)[1]+1,max_centers,
                       na.rm=TRUE)
  #print(best_n_cluters)
  kopt = kmeans(x = data , centers = best_n_cluters , nstart = nstart)
  
  if(plot==TRUE)
  {
    df=data.frame(x=1:max_centers,
                  y=wss,type="all",stringsAsFactors = FALSE)
    df[best_n_cluters,"type"]="Optimal"
    
    p=ggplot(df,aes(x,y))+geom_point(aes(color=type),size=2)+
      geom_line(size=0.75)
    plot(p)
  }
  return(list(max_centers = max_centers,
              best_n_cluters = best_n_cluters,
              wss = wss,kopt=kopt ))
}


kmeans_location_arrangment <-function(location_matrix)
{
  ## 35 is for number of teams
  locations_used=colnames(location_matrix)
  d_geo=geo[locations_used,]
  
  k_s=kmeans_woc(d_geo,max_centers = min(35,dim(location_matrix)[1]))
  k_opt=k_s$kopt$cluster
  clusters=unique(k_opt)
  
  cluster_scores=vector()
  for(i in clusters)
  {
    locations=names(k_opt)[which(k_opt==i)]
    if(length(locations)>1)
    {
      cluster_scores[i]=mean(colSums(location_matrix[,locations]))
    }
    else
    {
      cluster_scores[i]=sum(location_matrix[,locations])
    }
  }
  
  idx=order(cluster_scores,decreasing = TRUE)
  
  ## order locations
  ordered_locations=vector()
  for(i in idx)
  {
    locations_cluster=names(k_opt)[which(k_opt==i)]
    
    if(length(locations_cluster)>1)
    {
      scores=as.matrix(colSums(location_matrix[,locations_cluster]))
      idx_c=order(scores[,1],decreasing=TRUE)  
      ordered_locations=c(ordered_locations,rownames(scores)[idx_c])
    }
    else{
      ordered_locations=c(ordered_locations,locations_cluster[1])  
      
    }
  }
  
  return(ordered_locations[1:10])
}