# additional approach - First determine the number of subpopulations by subtypes the samples
#Clustering KNN rank selection
rm_outlier_replicate <- function(eset,druglist) {
  # out of 3 replicates, keep only the best 2 correlated
  data<-cor(exprs(eset))
  replicate_info<-pData(eset)
  time<-unique(eset$Time)
  concentration<-c("IC20","IC20_LOW")
  array_to_keep<-vector()
  for(i in 1:length(time)){
    for(j in 1:length(druglist)){
      for(z in 1:length(concentration)){
        arrays_to_compare<-rownames(subset(replicate_info , Time==time[i] & Drug==druglist[j] & Concentration==concentration[z]))
        if(length(arrays_to_compare)==1){
          array_to_keep<-append(array_to_keep,arrays_to_compare)
        } else{
          # Many replicate available, select best one
            cor_arrays<-data[arrays_to_compare,arrays_to_compare]
            diag(cor_arrays)<-0
            current<-rownames(cor_arrays[which(cor_arrays == max(cor_arrays), arr.ind = TRUE)[,1],])
            if(length(current)>2){
                remaining<-setdiff(arrays_to_compare,current)
                array_to_keep<-append(array_to_keep,names(which(cor_arrays[remaining,]==max(cor_arrays[remaining,]))))
            } else if(length(current)>1){
                array_to_keep<-append(array_to_keep,current[1])  
            } else {
                array_to_keep<-append(array_to_keep,current)
            }
        }
      }
    }
  }
  # Add dmso data
  for(i in 1:length(time)){
    arrays_to_compare<-rownames(subset(replicate_info , Time==time[i] & Drug=='DMSO'))
    if(length(arrays_to_compare)==1){
      array_to_keep<-append(array_to_keep,arrays_to_compare)
    } else{
      # Many replicate available, select best one
      cor_arrays<-data[arrays_to_compare,arrays_to_compare]
      diag(cor_arrays)<-0
      current<-rownames(cor_arrays[which(cor_arrays == max(cor_arrays), arr.ind = TRUE)[,1],])
      if(length(current)>2){
        remaining<-setdiff(arrays_to_compare,current)
        array_to_keep<-append(array_to_keep,names(which(cor_arrays[remaining,]==max(cor_arrays[remaining,]))))
      } else if(length(current)>1){
        array_to_keep<-append(array_to_keep,current[1])  
      } else {
        array_to_keep<-append(array_to_keep,current)
      }
    }
  }
  
  # Add media data
  for(i in 1:length(time)){
    arrays_to_compare<-rownames(subset(replicate_info , Time==time[i] & Drug=='Media'))
    if(length(arrays_to_compare)==1){
      array_to_keep<-append(array_to_keep,arrays_to_compare)
    } else{
      # Many replicate available, select best one
      cor_arrays<-data[arrays_to_compare,arrays_to_compare]
      diag(cor_arrays)<-0
      current<-rownames(cor_arrays[which(cor_arrays == max(cor_arrays), arr.ind = TRUE)[,1],])
      if(length(current)>2){
        remaining<-setdiff(arrays_to_compare,current)
        array_to_keep<-append(array_to_keep,names(which(cor_arrays[remaining,]==max(cor_arrays[remaining,]))))
      } else if(length(current)>1){
        array_to_keep<-append(array_to_keep,current[1])  
      } else {
        array_to_keep<-append(array_to_keep,current)
      }
    }
  }
  
  filtered_eset<-eset[, sampleNames(eset) %in% array_to_keep]
  return (filtered_eset)
}

knn_subpopulation<-function(eset){
  require(cluster)
  require(simpleaffy)
  subpopulation<-list()
  ICeset<-get.array.subset.exprset(eset,"Concentration", "IC20")  
  #Fct for KNN clustering based on silhouette
  OptK <- function(x, ns=3,nk=14,title,metric='manhattan',...) {
    asw<-numeric(nk-ns)
    for (k in ns:(nk-1)) {
      asw[k-ns+1] <- pam(x, k,metric=metric)$silinfo$avg.width
      k.best <- which.max(asw)+2
    }
    plot(ns:(nk-1), asw, type="s", main=title,
      xlab="K (number of clusters)", ylab = "mean silhouette width")
      
    axis(1, k.best, paste("best",k.best,sep="\n"), col="red", col.axis="red")
    return(pam(x, k.best, ...))
  }
    par(mfrow=c(1,3))
    time<-unique(as.vector(unlist(pData(ICeset)$Time)))
    for(t in 1:length(time)){
      subset<-get.array.subset.exprset(ICeset,"Time", time[t])
      clust <- OptK(t(exprs(subset)), ns=3,nk=14, title=paste("Clustering at ",time[t],"h",sep=""),metric='manhattan')$clustering
      subpopulation<-append(subpopulation,list(clust))
    }

  return (subpopulation)
}

# 1st step - Find differentially expressed genes
time_filter <- function(
eset, #From initialize_eset function
druglist,
ctrl='DMSO', # Either DMSO, Media or time (6hr used as control)
threshold=.95, # above which is interesting
time='all' # Either 1,2 or all
)
{
  output<-list()
  
	for(i in 1:length(druglist)){
		time_result <-loop_time_filter(eset,druglist[i],ctrl,threshold,time)
		output<-append(output,list(time_result))
    names(output)[i]<-druglist[i]
    print(paste(i,' - ',druglist[i],sep=''))
	}
	cat('End Time Filtering\n')
	return (output)
}

loop_time_filter <- function(eset,drug_test,ctrl,threshold,time) {
  require(simpleaffy)
  if(ctrl=='time'){
    eset <- get.array.subset.exprset(eset, "Concentration", c("IC20"))
  	eset <- get.array.subset.exprset(eset, "Drug", drug_test)
    dat <- exprs(eset)[,rank(pData(eset)$Time)]
    filter1<-apply(dat[,1:2],1,sd)
  	filter2<-apply(dat[,2:3],1,sd)
    dat1 <- dat[filter1>quantile(filter1,probs=threshold),]
  	dat2 <- dat[filter2>quantile(filter2,probs=threshold),]
  } else{
      subeset <- get.array.subset.exprset(eset, "Concentration", c("IC20"))
      subeset <- get.array.subset.exprset(subeset, "Drug", drug_test)
      sub_dat <- exprs(subeset)[,rank(pData(subeset)$Time)]
      
	    ctrleset <- get.array.subset.exprset(eset, "Concentration", c("NONE"))
	    ctrleset <- get.array.subset.exprset(ctrleset, "Drug", ctrl)
      ctrl_dat <- exprs(ctrleset)[,rank(pData(ctrleset)$Time)]
      
      if(!identical(match(rownames(sub_dat),rownames(ctrl_dat)),1:length(rownames(sub_dat)))){stop('error')}
      for(i in 1:nrow(sub_dat)){
        sub_dat[i,]<-sub_dat[i,]/ctrl_dat[i,]
      }
      filter1<-apply(sub_dat[,1:2],1,sd)
      filter2<-apply(sub_dat[,2:3],1,sd)
      dat1 <- sub_dat[filter1>quantile(filter1,probs=threshold),]
      dat2 <- sub_dat[filter2>quantile(filter2,probs=threshold),]
      }
  
  #Combining results according to parameters
  if(time=='all'){
    genes<-unique(append(rownames(dat1),rownames(dat2)))
  }
  else if(time=='1'){
    genes<-rownames(dat1)
  }
  else if(time=='2'){
    genes<-rownames(dat2)
  } else{
    stop('Selection parameters wrongs (1,2,all)')
  }
  return(genes)
}

#2nd step - Probeset-Drug concentration model fitting
model_fit <- function(
eset,druglist,relevant_probeset,
timepoint=6,
drug_flag=TRUE # Use drug data instead of assuming IC vector c(0,0.02,0.2)
)
{
  output<-list()
	for(i in 1:length(druglist)){
		fit_result<-loop_model_fit(eset,druglist[i],as.vector(unlist(relevant_probeset[druglist[i]])),
                               drug_flag,timepoint)
		output<-append(output,list(fit_result))
		names(output)[i]<-druglist[i]
		print(paste(i,' - ',druglist[i],sep=''))
	}
  cat('End Model Fitting\n')
	return(output)
}

loop_model_fit <- function(eset,drug,probeset,drug_flag,timepoint) {
	
	# Libraries and options
	require(IsoGene)
	require(drc)
  require(tcltk)			
	set.seed(1234)
  
  # Preparation Data
	eset<-get.array.subset.exprset(eset,"Drug", c("DMSO",drug))
	eset<-get.array.subset.exprset(eset,"Time", timepoint)
	Expression<-as.data.frame(exprs(eset[probeset,]))
  
	# Prepare dosage vectors with data instead of IC
	IC<-eset$Concentration
	if(drug_flag) {
		temp_data<-read.table(file='inputs/drug_input.txt',sep='\t',header=TRUE)
		drug_data<-as.vector(temp_data[,2])
		names(drug_data)<-temp_data[,1]
		experimental_ic<-as.numeric(drug_data[drug])
		for(z in 1:length(IC)) {
		  if(IC[z]=='NONE') {IC[z]=0}
		  if(IC[z]=='IC20') {IC[z]=experimental_ic}	
		  if(IC[z]=='IC20_LOW') {IC[z]=experimental_ic*.1}
		}
	} else {
	  for(z in 1:length(IC)) {
	    if(IC[z]=='NONE') {IC[z]=0}
	    if(IC[z]=='IC20') {IC[z]=0.2}	
	    if(IC[z]=='IC20_LOW') {IC[z]=0.02}
	  }
	}
	class(IC)<-"numeric"
	stat<-IsoGenemSAM(IC,Expression,fudge.factor=rep(1,length(rownames(Expression))))
	relevant_probeset<-names(stat$E2[stat$E2>0.05])
  slope_probeset<-as.factor(stat$direction[relevant_probeset])
	return (slope_probeset)
}

#3rd step - Features selection by random-forest to shrink signatures to most important components
machine_learning <- function(eset,druglist,signatures,disjoint=FALSE,min=20,max=100,
                             baseline='DMSO', # DMSO, Media or Others
                             cor_threshold=0.75) {
  output<-list()
  
  if(disjoint){
    # Remove overlapping probesets from signatures
    Dname<-c()
    for(i in 1:length(signatures)){
      Dname<-append(Dname,names(signatures[[i]]))
    }
    Sname<-table(Dname)
    Duplicate<-names(Sname[Sname>1])
    for(i in 1:length(signatures)){
      signatures[[i]] = signatures[[i]][!(names(signatures[[i]]) %in% Duplicate)]
    }
    
  }
  for(i in 1:length(druglist)){
    ML_result<-loop_ML(eset,druglist[i],druglist,signatures[[druglist[i]]],min,max,baseline,cor_threshold)
    output<-append(output,list(ML_result))
    names(output)[i]<-druglist[i]
    print(paste(i,' - ',druglist[i],sep=''))
  }
  cat('End Machine learning\n')
  return(output)    
}

loop_ML <- function(eset,drug,druglist,signature,Smin,Smax,baseline,cor_threshold){
  require(caret)
  require(mlbench)
  require(Hmisc)
  require(randomForest)
  require(dendextend)
  require(e1071)
  set.seed(100)
  
  subset<-get.array.subset.exprset(eset,"Time",6)
  if(baseline=='DMSO'){
    subset<-get.array.subset.exprset(subset,"Drug", c(drug,"DMSO"))
    subset<-get.array.subset.exprset(subset,"Concentration", c("IC20","NONE"))
  } else if(baseline=='Media'){
    subset<-get.array.subset.exprset(subset,"Drug", c(drug,"Media"))
    subset<-get.array.subset.exprset(subset,"Concentration", c("IC20","NONE"))
  } else if(baseline=='Others'){
    subset1<-get.array.subset.exprset(subset,"Drug", setdiff(druglist,drug))
    subset1<-get.array.subset.exprset(subset1,"Concentration","IC20")
    subset2<-get.array.subset.exprset(subset,"Drug", drug)
    subset2<-get.array.subset.exprset(subset2,"Concentration",c("IC20","IC20_LOW"))
    subset<-Biobase::combine(subset1,subset2)
  }
  info<-as.factor(pData(subset)$Drug==drug)
  
  # Remove correlated features
  correlationMatrix <- cor(t(exprs(eset)[names(signature),]))
  highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=cor_threshold)
  signature<-signature[-highlyCorrelated]
  
  if(Smin<min(summary(signature))/2){
    Smin=round(min(summary(signature))/2,0)
  }
  
  if(Smax>min(summary(signature))){
    Smax=min(summary(signature))
  }
  if(Smax<Smin){Smax=Smax*2}
  
  # Prenormalization with center and scaling
  dat<-data.frame(t(exprs(subset)[names(signature),]))
  normalization<-preProcess(dat)
  x <- predict(normalization, dat)
  x <- as.data.frame(x)
  
  # Rank features by importance
  range <- seq(Smin,Smax,round(Smax/10))
  print(info)

  ctrl <- rfeControl(functions = rfFuncs,
                     method = "cv",
                     number=10,
                     verbose = FALSE)
    
  # Start random forest
  rfProfile <- rfe(x,info,
                   sizes = range,
                   rfeControl = ctrl,
                   allowParallel = TRUE)
   
  optimal_signature<-predictors(rfProfile)
  return (signature[optimal_signature])
}

#4th step - Subpopulation targeting - Generate interaction matrix
population_targeting <- function(eset,druglist,signatures,within.sig=TRUE)
{
	cross.kill<-data.frame(matrix(data=0,ncol=length(druglist), nrow=length(druglist), dimnames=list(druglist,druglist)))
	cross.survival<-cross.kill
  
	for(i in 1:length(druglist)) {
    targeting<-similarity(eset,druglist[i],druglist,signatures[[i]],within.sig,timepoint)
    cross.kill[druglist[i],rownames(targeting)]=targeting$kill
    cross.survival[druglist[i],rownames(targeting)]=targeting$survival
    print(paste(i,' - ',druglist[i],sep=''))
	}
  
  # Combination interaction matrices
  output<-append(list(cross.kill),list(cross.survival))
  names(output)<-c('kill','survival')
  return (output)
}

cosine <- function(x,y) {return(x %*% y / sqrt(x%*%x * y%*%y))}

similarity<-function(eset,drug,druglist,signature,within.sig,timepoint) {

  require(simpleaffy)
  set.seed(1234)
  
	eset_drug<-get.array.subset.exprset(eset, "Concentration","IC20")
	eset_drug<-get.array.subset.exprset(eset_drug, "Time",timepoint)
	
  # Matrix expression #
  data<-exprs(eset_drug)
  
  if(within.sig){data<-data[names(signature),]}
  
  # row normalize std
  data.mean <- apply(data[,],1,mean,na.rm=T)
  data.sd <- apply(data[,],1,sd,na.rm=T)
  norm_data<-as.data.frame((data-data.mean)/data.sd)
  
  # Format signature
  DOWN_genes<-names(signature[signature=='d'])
  UP_genes<-names(signature[signature=='u'])

  yDOWN <- as.numeric(rownames(norm_data) %in% DOWN_genes)
  yUP <- abs(yDOWN-1)
  
  DOWN_SCORE<--apply(norm_data,2,cosine,y=yDOWN)
	UP_SCORE<--apply(norm_data,2,cosine,y=yUP)
  
  output<-data.frame(row.names=pData(eset_drug)$Drug,kill=scale(DOWN_SCORE),survival=scale(UP_SCORE),stringsAsFactors = FALSE)
	return (output)
}

assymetric_killing<-function(interaction,population){
  # Find drug combinations that are asymetric in terms of killing
  symetry<-sign(interaction$kill)
  diag(symetry)<-0
  hierarchy <- matrix(data=0,ncol=ncol(interaction$kill), nrow=nrow(interaction$kill),dimnames=list(rownames(interaction$kill),colnames(interaction$kill)))
  rownames(hierarchy)<-rownames(symetry)
  colnames(hierarchy)<-colnames(symetry)
  for(i in 1:nrow(symetry)){
    for(j in 1:ncol(symetry)){
      if(symetry[i,j]==1 && symetry[j,i]==-1){hierarchy[i,j]=1}
    }
  }
  
  am.graph<-new("graphAM", adjMat=as.matrix(hierarchy), edgemode="directed")
  plot(am.graph,"dot", attrs = list(node = list(fillcolor = "lightblue",fontsize=15,shape='circle',width=1,height=0.5),
                                    edge = list(arrowsize=0.5)))
  
  # Social network Analysis :)
  require(sna)
  # calculate reachibility
  elimination<-reachability(as.matrix(hierarchy))
  drug_spectrum<-rowSums(elimination)
  names(drug_spectrum)<-rownames(interaction$kill)
  hierarchy%*%population
  
  prediction<-function(eset,druglist,interaction,hierarchy,subpopulation=TRUE,spectrum=TRUE,scale=TRUE,top=10){
    druglist<-colnames(targeting_matching[[1]])
    interaction<-targeting_matching[['targeting']]
    if(scale){interaction<-scale(interaction)}
    final_prediction<-data.frame(drug1=character(),drug2=character(),rank=numeric())
    
    if(!subpopulation) {
      for(i in 1:(length(druglist)-1)) {
        for(j in (i+1):length(druglist)) {
          drug_combination<-interaction[druglist[i],druglist[j]]+interaction[druglist[j],druglist[i]]
          final_prediction<-rbind(final_prediction,data.frame(drug1=druglist[i],drug2=druglist[j],rank=drug_combination))
        }	
      }
    } else {
      # Hierarchy based scoring
      targeting_matrix<- hierarchy_subpopulation(targeting_matching,eset,subpopulations,graph_display=FALSE)
      for(i in 1:(length(druglist)-1)) {
        for(j in (i+1):length(druglist)) {
          if(!spectrum) {
            # Hierarchy based on drug and subpopulation at timepoint 3
            drug_combination <- interaction[druglist[i],druglist[j]] * children_killing(druglist[j],targeting_matrix,interaction) + interaction[druglist[j],druglist[i]] * children_killing(druglist[i],targeting_matrix,interaction)
          } else {
            drug_combination <- interaction[druglist[i],druglist[j]] * children_killing(druglist[j],targeting_matrix,interaction) * (1/targeting_matrix$spectrum[[druglist[j]]]) + 
              interaction[druglist[j],druglist[i]] * children_killing(druglist[i],targeting_matrix,interaction) * (1/targeting_matrix$spectrum[[druglist[i]]])
          }
          final_prediction<-rbind(final_prediction,data.frame(drug1=druglist[i],drug2=druglist[j],rank=drug_combination))
        }
      }
    }
    
    rank<-order(final_prediction$rank,decreasing=TRUE)
    final_prediction$rank<-rank
    final_prediction<-final_prediction[with(final_prediction, order(rank)), ]
    rank<-ncbi_scoring(final_prediction)
    cat(paste('Top combination prediction: ',network_scoring(final_prediction,druglist,top),sep=""))
  }
  
}
children_killing<- function(drug,targeting_matrix,interaction) {
  
  population<-targeting_matrix[['drug']]
  children<-colnames(population)[population[drug,]>0]
  round<-0
  efficiency<-c(0)
  past_children<-vector()
  repeat{
    if(any(is.na(children))) { break }
    if(round>3) {break}
    efficiency<-append(efficiency,as.vector(unlist(interaction[children,])))
    past_children<-append(past_children,children)
    children<-colnames(population)[population[children,]>0]
    children<-setdiff(children,past_children)
    round<-round+1
  }
  
  return (max(efficiency))
}

hierarchy_subpopulation <- function(data,eset,subpopulations,graph_display=FALSE) {
  library(Rgraphviz)
  #Use surviving population from knn at t=24hr
  clustering<-subpopulations$clustering
  names(clustering)<-pData(eset)[names(subpopulations$clustering),]$Drug
  temp<-names(table(clustering))[table(clustering)>1]
  
  drug_population<-matrix(data=0,nrow=14,ncol=length(temp))
  colnames(drug_population)<-temp
  rownames(drug_population)<-rownames(data[[1]])
  for(i in 1:length(temp)){
    drug_population[unique(names(clustering)[clustering==temp[i]]),i]<-1
  }
  
  # Rows are parent nodes, Cols are children nodes
  hierarchy<-as.matrix(data[['des']])
  hierarchy[,]<-0	
  # Find drug combinations that are asymetric
  sym_test<-sign(as.matrix(data[['des']]))
  sym_test[sym_test<0]<-0
  for(i in 1:length(sym_test[,1])) {
    for(j in 1:length(sym_test[1,])) {
      if (sym_test[i,j]>sym_test[j,i]) {hierarchy[i,j]=1;}
    }
  }
  
  #Combine
  output<-change_matrix_dimension(hierarchy,drug_population)
  
  am.graph<-new("graphAM", adjMat=output, edgemode="directed")
  plot(am.graph,"dot", attrs = list(node = list(fillcolor = "lightblue",fontsize=15,shape='rectangle',width=1,height=0.5),
                              edge = list(arrowsize=1)))
  
# 	graph <- graph.adjacency(as.matrix(output),mode="directed",weighted=T)
# 	mst <- minimum.spanning.tree(graph)
# 	if(graph_display) {
# 		plot.igraph(mst,layout=layout.fruchterman.reingold,vertex.size=20,
# 			vertex.label=V(mst)$name,	
# 			main='Hierarchy of subpopulations\nbased on non-symetric killing',
# 			vertex.label.color="black",
# 			edge.label.color="black"
# 	)
# 	}
# 	#Return graph as matrix with weight
  x<-list(hierarchy)
  importance<-rowSums(output)
  spectrum<-as.vector(drug_population%*%importance)
  names(spectrum)<-rownames(drug_population)
  
  x<-append(x,list(spectrum))

  names(x)<-c('drug','spectrum')
	return (x)
}

# Extra function
factor_dataframe <- function(dataframe){
    class.data  <- sapply(dataframe, class)
    factor.vars <- class.data[class.data == "factor"]
    for (colname in names(factor.vars)) {
       dataframe[,colname] <- as.character(dataframe[,colname])
    }
    return (dataframe)
}

change_matrix_dimension<-function(M,V){
  #M Original matrix
  #V Rule matrix
  
  output<-matrix(data=NA,nrow=ncol(V),ncol=ncol(V))
  rownames(output)<-colnames(V)
  colnames(output)<-colnames(V)
  
  for(i in 1:ncol(V)){
    for(j in 1:ncol(V)){
      dimension1<-rownames(V)[V[,i]>0]
      dimension2<-rownames(V)[V[,j]>0]
      output[i,j]<-sum(M[dimension1,dimension2])
    }
  }
  return(output)  
}