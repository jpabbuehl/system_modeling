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

knn_subpopulation<-function(eset,metric='manhattan'){
  subpopulation<-list()
  ICeset<-get.array.subset.exprset(eset,"Concentration", "IC20")  
  #Fct for KNN clustering based on silhouette
  OptK <- function(x, ns=3,nk=14,title,...) {
    asw<-numeric(nk-ns)
    for (k in ns:(nk-1)) {
      asw[k-ns+1] <- pam(x, k,metric=metric)$silinfo$avg.width
      k.best <- which.max(asw)+ns-1
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

knn_process<-function(eset,druglist,clustering,timepoint){
  clust<-data.frame(drug=pData(eset)[names(clustering[[timepoint]]),1],population=clustering[[timepoint]])
  clust<-clust[match(clust[,1],druglist),]
  clust[,2]<-as.factor(clust[,2])
  population=matrix(data=0,ncol=length(levels(clust[,2])), nrow=length(druglist), dimnames=list(druglist,1:length(levels(clust[,2]))))
  for(i in 1:length(levels(clust[,2]))){
    population[,i]=as.numeric(clust$population==i)
  }
  return (population)  
}

# 1st step - Find differentially expressed genes
time_filter <- function(
eset, #From initialize_eset function
druglist,
ctrl='DMSO', # Either DMSO, Media or Time (6hr used as control)
threshold=.95, # above which is interesting
timepoint=1 # Either 1,2 or 12
)
{
  output<-list()
  
	for(i in 1:length(druglist)){
		time_result <-loop_time_filter(eset,druglist[i],ctrl,threshold,timepoint)
		output<-append(output,list(time_result))
    names(output)[i]<-druglist[i]
    #print(paste(i,' - ',druglist[i],sep=''))
	}
	cat('End Time Filtering\n')
	return (output)
}

loop_time_filter <- function(eset,drug_test,ctrl,threshold,timepoint) {
  if(ctrl=='Time'){
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
  if(time==12){
    output1<-rep(1,length(rownames(dat1)))
    names(output1)<-rownames(dat1)
    dat2<-dat2[!(rownames(dat2) %in% rownames(dat1))]
    output2<-rep(2,length(rownames(dat2)))
    names(output2)<-rownames(dat2)
    output<-append(output1,output2)
  }
  else if(time==1){
    output<-rep(1,length(rownames(dat1)))
    names(output)<-rownames(dat1)
  }
  else if(time==2){
    output<-rep(2,length(rownames(dat2)))
    names(output)<-rownames(dat2)
  } else{
    stop('Selection parameters wrongs (1,2,12)')
  }
  return(output)
}

#2nd step - Probeset-Drug concentration model fitting
model_fit <- function(
eset,druglist,relevant_probeset,
timepoint=1,
drug_flag=TRUE # Use drug data instead of assuming IC vector c(0,0.02,0.2)
)
{
  output<-list()
	for(i in 1:length(druglist)){
		fit_result<-loop_model_fit(eset,druglist[i],names(relevant_probeset[[druglist[i]]]),
                               drug_flag,timepoint)
		output<-append(output,list(fit_result))
		names(output)[i]<-druglist[i]
		#print(paste(i,' - ',druglist[i],sep=''))
	}
  #cat('End Model Fitting\n')
	return(output)
}

loop_model_fit <- function(eset,drug,probeset,drug_flag,timepoint) {
	
	# Libraries and options
		
	set.seed(1234)
  
  # Preparation Data
	if(timepoint==1){timepoint=6}
	if(timepoint==2){timepoint=12}
	if(timepoint==3){timepoint=24}
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
    #print(paste(i,' - ',druglist[i],sep=''))
  }
  #cat('End Machine learning\n')
  return(output)    
}

loop_ML <- function(eset,drug,druglist,signature,Smin,Smax,baseline,cor_threshold){
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

  ctrl <- rfeControl(functions = rfFuncs,
                     method = "cv",
                     number=10,
                     verbose = FALSE)
    
  # Start random forest
  rfProfile <- rfe(x,info,
                   sizes = range,
                   rfeControl = ctrl,
                   allowParallel = FALSE)
   
  optimal_signature<-caret::predictors(rfProfile)
  return (signature[optimal_signature])
}

#4th step - Subpopulation targeting - Generate interaction matrix
population_targeting <- function(eset,druglist,signatures,within.sig=TRUE,timepoint)
{
	cross.kill<-data.frame(matrix(data=0,ncol=length(druglist), nrow=length(druglist), dimnames=list(druglist,druglist)))
	cross.survival<-cross.kill
  
	for(i in 1:length(druglist)) {
    targeting<-similarity(eset,druglist[i],druglist,signatures[[i]],within.sig,timepoint)
    cross.kill[druglist[i],rownames(targeting)]=targeting$kill
    cross.survival[druglist[i],rownames(targeting)]=targeting$survival
    #print(paste(i,' - ',druglist[i],sep=''))
	}
  
  # Combination interaction matrices
  output<-append(list(cross.kill),list(cross.survival))
  names(output)<-c('kill','survival')
  return (output)
}

cosine <- function(x,y) {return(x %*% y / sqrt(x%*%x * y%*%y))}

similarity<-function(eset,drug,druglist,signature,within.sig,timepoint) {

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

assymetric_killing<-function(interaction,population,plot=FALSE){
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
  
  if(plot){
    am.graph<-new("graphAM", adjMat=as.matrix(hierarchy), edgemode="directed")
    plot(am.graph,"dot", attrs = list(node = list(fillcolor = "lightblue",fontsize=15,shape='circle',width=1,height=0.5),
                                      edge = list(arrowsize=0.5)))
  }
  # calculate reachibility
  elimination<-reachability(as.matrix(hierarchy))
  drug_spectrum<-rowSums(elimination)
  names(drug_spectrum)<-rownames(interaction$kill)
  output<-hierarchy%*%population
  return (output)
}


prediction<-function(druglist,interaction,hierarchy,drug_data,hierarchy_flag=TRUE,drug_flag=FALSE) {
  
  final_prediction<-data.frame(drug1=character(),drug2=character(),rank=numeric(),category=numeric())
  x<-drug_data$rank
  normalized=1-normalized
  names(normalized)<-rownames(drug_data)
  
  if(!hierarchy_flag) {
    for(i in 1:(length(druglist)-1)) {
      for(j in (i+1):length(druglist)) {
        if(drug_flag){
          drug_combination<-(interaction$kill[druglist[i],druglist[j]]*normalized[names(normalized)==druglist[i]]+interaction$kill[druglist[j],druglist[i]]*normalized[names(normalized)==druglist[j]])-
            (interaction$survival[druglist[i],druglist[j]]+interaction$survival[druglist[j],druglist[i]])
        } else {
        drug_combination<-(interaction$kill[druglist[i],druglist[j]]+interaction$kill[druglist[j],druglist[i]])-
          (interaction$survival[druglist[i],druglist[j]]+interaction$survival[druglist[j],druglist[i]])
        }
        final_prediction<-rbind(final_prediction,data.frame(drug1=druglist[i],drug2=druglist[j],rank=drug_combination,category=1))
      }	
    }
  } else {
    # Hierarchy based scoring  
    for(i in 1:(length(druglist)-1)) {
      for(j in (i+1):length(druglist)) {
        killing<-apply(hierarchy[c(druglist[i:j]),],2,max)
        if(drug_flag){
          drug_combination<-(interaction$kill[druglist[i],druglist[j]]*normalized[names(normalized)==druglist[i]]+interaction$kill[druglist[j],druglist[i]]*normalized[names(normalized)==druglist[j]])-
            (interaction$survival[druglist[i],druglist[j]]+interaction$survival[druglist[j],druglist[i]])
        } else {
          drug_combination<-(interaction$kill[druglist[i],druglist[j]]+interaction$kill[druglist[j],druglist[i]])-
            (interaction$survival[druglist[i],druglist[j]]+interaction$survival[druglist[j],druglist[i]])
        }
        final_prediction<-rbind(final_prediction,data.frame(drug1=druglist[i],drug2=druglist[j],rank=drug_combination,category=sum(killing)))
       }
    }
  }
  
  rank<-order(final_prediction$category,final_prediction$rank,decreasing=TRUE)
  final_prediction<-final_prediction[rank,]
  final_prediction$rank<-seq(1,nrow(final_prediction),1)
  return (final_prediction)
}

dream_ranking <- function (prediction) {
  
  seed<-10000
  data<-read.table(file='inputs/drug_synergy_data_IC20.txt',sep='\t',header=TRUE)
  colnames(data) <- c('drug1','drug2','eob','eob_error')
  unique_drugs <- unique(data$drug1)
  drugpair <- paste(data$drug1,'_',data$drug2,sep="")
  idx <-order(data$eob)
  eob <- data$eob[idx]
  eob_error <-data$eob_error[idx]
  drugpair <- drugpair[idx]
  
  p_matrix <- probability_matrix(eob,eob_error)
  colnames(prediction) <-c("drug1","drug2","rank")
  
  drugpair_pred<-vector()
  # Deal with prediction
  for(i in 1:length(drugpair)) {
    current_compounds<- as.vector(unlist(strsplit(drugpair[i],'_')))
    if(length(unique(current_compounds))>1) {
      match1<-grepl(current_compounds[1],prediction$drug1)
      match2<-grepl(current_compounds[2],prediction$drug2)
      
      # No intersection, switch columns
      if(!sum(match1*match2)) {
        match1<-grepl(current_compounds[2],prediction$drug1)
        match2<-grepl(current_compounds[1],prediction$drug2)	
      }
      
      if(sum(match1*match2)) {
        match_pred<-match1*match2
        class(match_pred)<-'logical'
        rank_pred<-as.numeric(as.vector(unlist(prediction$rank[match_pred])))
        drugpair_pred<-append(drugpair_pred,rank_pred)
      } else {stop(i)}
    }
  }
  if(length(drugpair_pred)!=length(drugpair)) {stop('error prediction size')}
  
  cindex_nulldist<-vector()
  # Computed weighted cindex,null dist and pvalue
  weighted_cindex <- concordance(drugpair_pred,1:length(drugpair),p_matrix)
  for (i in 1:seed) {
    cindex_nulldist[i] <- concordance(sample(length(drugpair)),1:length(drugpair),p_matrix)
  }
  pv_cindex <- length(cindex_nulldist[cindex_nulldist>=weighted_cindex])/length(cindex_nulldist)
  #cat('c-index=',weighted_cindex,'\n')
  #cat('p-value=',pv_cindex,'\n')
  
  #Ranking test
  ranking<-read.table(file='inputs/ranking.txt',sep=";",header=TRUE)
  final_rank<-ranking[ranking$cindex<weighted_cindex,][1,2]
  #cat('Ranking=',final_rank,'\n')
  return(final_rank)
}

probability_matrix <- function (x,x_std) {
  x <- as.vector(unlist(x))
  n <- length(x)
  p_matrix <- rep(0,n)
  X <- matrix(x,n,n)
  X <- X-t(X)
  X_std <- matrix(x_std,n,n)
  X_std <- sqrt(X_std^2 + t(X_std)^2)
  p_matrix <- .5*(1 + erf(X/X_std))
  return (p_matrix)
}

concordance <- function (x, y,p_matrix) {
  x <- as.vector(unlist(x))
  y <- as.vector(unlist(y))
  n <- length(x)
  X <- matrix(x,n,n)
  Y <- matrix(y,n,n)
  C <- sign(X-t(X))==sign(Y-t(Y))
  
  C <- C*(1-t(p_matrix)) + (1-C)*t(p_matrix)
  C <- sum(sum(C[lower.tri(C)]))/n/(n-1)*2;
  return (C)
}

top_ranking<-function(analysis,top){
  gs <- read.table(file='inputs/drug_synergy_data_IC20.txt',sep='\t',header=TRUE,stringsAsFactors = FALSE)
  colnames(gs) <- c('drug1','drug2','eob','eob_error')
  gs <- gs[1:top,1:2]
  analysis <- analysis[1:top,1:2]
  score=0
  
  for(i in 1:nrow(analysis)){
    for(j in 1:nrow(gs)){
      if((gs[i,1]==analysis[j,1] && gs[i,2]==analysis[j,2]) || (gs[i,1]==analysis[j,2] && gs[i,2]==analysis[j,1])){
        score=score+1
      }
    }
  }
  return (score)
}