initialize_eset <- function( filename='inputs/input.txt', metadata='inputs/info.txt',drugfile='inputs/drug_input.txt'){
library(Biobase)
library(simpleaffy)

### Prepare ExpressionSet with expression values and pData ###
data <- read.table(file=filename,sep="\t", header=TRUE)
data<-data[-c(1,2,3),]
rownames(data)<-data$AffyID
affy_gene<-data[,c("AffyID","Genename")] 
data$AffyID<-NULL
data$Genename<-NULL
varInfo<-read.table(metadata,header=TRUE,colClasses="character",stringsAsFactors = FALSE)
rownames(varInfo)<-varInfo$Sample_ID
varInfo$Sample_ID<-NULL
data<-as.matrix(data)
class(data)<-"numeric"
phenoData <- new("AnnotatedDataFrame", data=varInfo)
eset<-new("ExpressionSet", exprs=data, phenoData=phenoData)
cat('ExpressionSet ready\n')
return(eset)
}

initialize_drug <- function(eset) {
### Prepare vector with all drug names, except media and DMSO ###
druglist<-unique(as.vector(eset$Drug))
druglist<-druglist[-which(druglist == "Media")]
druglist<-druglist[-which(druglist == "DMSO")]
druglist<-unique(druglist)
return(druglist)
}

ncbi_scoring <- function (prediction) {
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
  cat('c-index=',weighted_cindex,'\n')
  cat('p-value=',pv_cindex,'\n')
  
  #Ranking test
  ranking<-read.table(file='inputs/ranking.txt',sep=";",header=TRUE)
  final_rank<-ranking[ranking$cindex<weighted_cindex,][1,2]
  cat('Ranking=',final_rank,'\n')
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
  library(VGAM)
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

top_combination_graph<-function(eset,rank_prediction,top=10,landscape_all,druglist) {
	data<-factor_dataframe(rank_prediction)
	par(mfrow=c((top%/%3)+1,3))
	for(i in 1:top) {
		drugA<-landscape_all[[data[i,1]]]
		drugB<-landscape_all[[data[i,2]]]
		combined<-landscape_niveling(drugA,drugB)
		landscape_plot(combined,paste('Rank ',i,': ',data[i,1],' & ',data[i,2],sep=""),data[i,],druglist,col_line=c(2,6))
	}
}

combination_graph<-function(eset,rank_prediction,druglist,fit_result) {
	data<-factor_dataframe(rank_prediction)
	last<-length(data[,1])
	first<-landscape_definition(eset,data[1,1],druglist,fit_result)
	second<-landscape_definition(eset,data[1,2],druglist,fit_result)
	third<-landscape_definition(eset,data[last,1],druglist,fit_result)
	fourth<-landscape_definition(eset,data[last,2],druglist,fit_result)
		
	best_combination<-landscape_niveling(first,second)
	worst_combination<-landscape_niveling(third,fourth)
	
	par(mfrow=c(2,3))
	landscape_plot(first,paste('Drug A: ',data[1,1],sep=""),data[1,1],druglist,col_line=c(2,2))
	landscape_plot(second,paste('Drug B: ',data[1,2],sep=""),data[1,2],druglist,col_line=c(6,6))
	landscape_plot(best_combination,paste('Best predicted combination: ',data[1,1],' & ',data[1,2],sep=""),data[1,],druglist,col_line=c(2,6))
	landscape_plot(third,paste('Drug C: ',data[last,1],sep=""),data[last,1],druglist,col_line=c(2,2))
	landscape_plot(fourth,paste('Drug D: ',data[last,2],sep=""),data[last,2],druglist,col_line=c(6,6))
	landscape_plot(worst_combination,paste('Worst predicted combination: ',data[last,1],' & ',data[last,2],sep=""),data[last,],druglist,col_line=c(2,6))
}

network_scoring<-function(rank_prediction,druglist,top=10) {
  library(igraph)
  library(Matrix)
  set.seed(3952)
  predicted_data<-prediction_matrix(rank_prediction,druglist)
  data<-experimental_matrix(druglist)
  
  for(i in 1:length(druglist)) {
    for(j in 1:length(druglist)) {
      if(predicted_data[i,j]>top && i!=j) {predicted_data[i,j]=0;predicted_data[j,i]=0}
      if(data[i,j]>top && i!=j) {data[i,j]=0;data[j,i]=0}	
    }
  }
  
  data[is.na(data)]<-0
  predicted_data[is.na(predicted_data)]<-0
  as.vector(data*predicted_data)
  top10_overlap<-nnzero(data*predicted_data)/2
  return (top10_overlap)
}

prediction_matrix <- function(prediction,druglist) {
	output<-matrix(nrow=14,ncol=14)
	rownames(output)<-druglist
	colnames(output)<-druglist
	prediction<-factor_dataframe(prediction)
	for(i in 1:length(prediction[,1])) {
		output[prediction[i,1],prediction[i,2]]<-prediction[i,3]
		output[prediction[i,2],prediction[i,1]]<-prediction[i,3]
	}
	return (output)
}

experimental_matrix <- function(druglist) {
	output<-matrix(nrow=14,ncol=14)
	rownames(output)<-druglist
	colnames(output)<-druglist
	data<-factor_dataframe(read.table(file='inputs/drug_synergy_data_IC20.txt',sep='\t',header=TRUE))
	colnames(data) <- c('drug1','drug2','rank','eob_error')
	data$rank<-1:length(data[,1])
	data$eob_error<-NULL
	
	for(i in 1:length(data[,1])) {
		output[data[i,1],data[i,2]]<-data[i,3]
		output[data[i,2],data[i,1]]<-data[i,3]
	}
	return (output)
}

inverse_eset_fct <- function(eset) {
	input<-exprs(eset)
	samples_mean<-apply(input,2,mean)
	output<-matrix(nrow=length(input[,1]),ncol=length(input[1,]))
	rownames(output)<-rownames(input)
	colnames(output)<-colnames(input)
	output[is.na(output)]<-0
	#Per columns, per rows
	for(j in 1:length(input[1,])) {
		for(i in 1:length(input[,1])) {
			diff_to_mean<-abs(samples_mean[j]-input[i,j])
			if(input[i,j]>samples_mean[j]) { output[i,j]<-samples_mean[j]-diff_to_mean}
			else {output[i,j]<-samples_mean[j]+diff_to_mean}
		}
	}
	exprs(eset)<-output
	return (eset)
}

rank_data<-function() {
  data<-read.table(file='inputs/drug_synergy_data_IC20.txt',sep='\t',header=TRUE)
  colnames(data) <- c('drug1','drug2','rank','eob_error')
  idx <-order(data$rank)
  data$rank<-idx
  data$eob_error<-NULL
  return(data)
}

synergy_scoring<-function(rank_prediction,parameters=FALSE) {
	data<-rank_data()
	drugpair <- paste(data$drug1,'_',data$drug2,sep="")
	drugpair_pred<-vector()
	# Deal with rank_prediction
	for(i in 1:length(drugpair)) {
		current_compounds<- as.vector(unlist(strsplit(drugpair[i],'_')))
		if(length(unique(current_compounds))>1) {
		  match1<-grepl(current_compounds[1],rank_prediction$drug1)
		  match2<-grepl(current_compounds[2],rank_prediction$drug2)
		  
		  # No intersection, switch columns
		  if(!sum(match1*match2)) {
			match1<-grepl(current_compounds[2],rank_prediction$drug1)
			match2<-grepl(current_compounds[1],rank_prediction$drug2)	
		  }
		  
		  if(sum(match1*match2)) {
			match_pred<-match1*match2
			class(match_pred)<-'logical'
			rank_pred<-as.numeric(as.vector(unlist(rank_prediction$rank[match_pred])))
			drugpair_pred<-append(drugpair_pred,rank_pred)
		  } else {stop(i)}
		}
	}
	match_prediction<-cbind(drugpair_pred[1:24],drugpair_pred[1:24]<25)
	final<-cbind(data[1:24,],match_prediction)
	colnames(final)[c(4,5)]<-c("predicted_rank","predicted_synergy")
	detected<-sum(final$predicted_synergy)
	final$predicted_synergy<-sub(1,'yes',final$predicted_synergy)
	final$predicted_synergy<-sub(0,'no',final$predicted_synergy)
	synergy_detected<-(detected/24*100)
	if(parameters) { return (synergy_detected)}
	else {return (final)}
}
