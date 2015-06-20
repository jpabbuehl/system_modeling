# additional approach - First determine the number of subpopulations by subtypes the samples
#Clustering KNN rank selection
rm_outlier_replicate <- function(eset) {
  # out of 3 replicates, keep only the best 2 correlated
  data<-cor(exprs(eset))
  replicate_info<-pData(eset)
  time<-unique(eset$Time)
  druglist<-initialize_drug(eset)
  concentration<-c("IC20","IC20_LOW")
  array_to_keep<-vector()
  for(i in 1:length(time)){
    for(j in 1:length(druglist)){
      for(z in 1:length(concentration)){
        arrays_to_compare<-rownames(subset(replicate_info , Time==time[i] & Drug==druglist[j] & Concentration==concentration[z]))
        cor_arrays<-data[arrays_to_compare,arrays_to_compare]
        diag(cor_arrays)<-0
        array_to_keep<-append(array_to_keep,rownames(
          cor_arrays[which(cor_arrays == max(cor_arrays), arr.ind = TRUE)[,1],]
        ))
      }
    }
  }

  # Add dmso data
  for(i in 1:length(time)){
    arrays_to_compare<-rownames(subset(replicate_info , Time==time[i] & Drug=='DMSO'))
    cor_arrays<-data[arrays_to_compare,arrays_to_compare]
    diag(cor_arrays)<-0
    array_to_keep<-append(array_to_keep,rownames(
      cor_arrays[which(cor_arrays == max(cor_arrays), arr.ind = TRUE)[,1],]
    ))
  }
  
  filtered_eset<-eset[, sampleNames(eset) %in% array_to_keep]
  return (filtered_eset)
}

knn_subpopulation<-function(eset){
  library(cluster)
  ICeset<-get.array.subset.exprset(eset,"Concentration", "IC20")
  subpopulation<-list()
  
  #Fct for KNN clustering based on silhouette
  OptK <- function(x, nk=14,title,...) {
    asw <- numeric(nk) 
    for (k in 2:nk) {  
      asw[k] <- pam(x, k, ...)$silinfo$avg.width
      k.best <- which.max(asw)
    }
    plot(1:nk, asw, type="s", main=title,
      xlab="K (number of clusters)", ylab = "mean silhouette width")
      
    axis(1, k.best, paste("best",k.best,sep="\n"), col="red", col.axis="red")
    return(pam(x, k.best, ...))
  }
  
  par(mfrow=c(1,3))
  time<-unique(as.vector(unlist(pData(ICeset)$Time)))
  for(t in 1:length(time)){
    subset<-get.array.subset.exprset(ICeset,"Time", time[t])
    clust <- OptK(t(exprs(subset)), 14, title=paste("Clustering Optimization for ",time[t],"hr",sep=""))
    subpopulation<-append(subpopulation,list(clust))
  }
  
  subpopulation<-subpopulation[[3]]
  return (subpopulation)
}

# 1st step - Find differentially expressed genes with Betr package
betr_filter <- function(
eset, #From initialize_eset function
ctrl='DMSO', # Either DMSO, Media or Time (6hr used as control)
threshold=.95,
report=TRUE
)
{
	library(betr)
	druglist<-unique(pData(eset)$Drug)
  druglist<-druglist[!grepl('DMSO',unique(pData(eset)$Drug))]
	if(report) {betr_report<-data.frame(drug=NA,detected=NA);}
	output<-list()
	
	for(i in 1:length(druglist)){
		betr_result <-loop_betr_filter(eset,druglist[i],ctrl,threshold)
		output<-append(output,list(betr_result))
		
		if(report) {
			betr_report<-rbind(betr_report,c(druglist[i],length(betr_result)))
			#betr_graph(eset,druglist[i],betr_result,25) Not working yet
		}
	}
	
	if(report) {
		betr_report<-betr_report[-1,]
		betr_report$detected<-as.numeric(betr_report$detected)
		betr_report$percent<- apply(betr_report,1,function(row) as.numeric(format(as.numeric(row[2])*100/length(featureNames(eset)),digits=2)))
		write.table(betr_report,file='betr_report.tab', sep="\t")	
	}
  names(output)<-druglist
	cat('End Betr Filtering\n')
	return (output)
}

loop_betr_filter <- function(eset,drug_test,ctrl,threshold) {
	if(grepl(ctrl,'Time')) {
		twoCondition_flag=FALSE 
		eset <- get.array.subset.exprset(eset, "Concentration", c("IC20"))
		eset <- get.array.subset.exprset(eset, "Drug", drug_test)
		condition_data <-NULL
		timepoint_data <-pData(eset)$Time
		timepoint_data<-sub(6,0,timepoint_data)
		timepoint_data<-sub(12,1,timepoint_data)
		timepoint_data<-sub(24,2,timepoint_data)
		timepoint_data<-as.numeric(timepoint_data)
		replicate_data<-pData(eset)$Replicate		
	} else {
		twoCondition_flag=TRUE
		eset <- get.array.subset.exprset(eset, "Concentration", c("IC20","NONE"))
		eset <- get.array.subset.exprset(eset, "Drug", c(ctrl,drug_test))
		condition_data <-pData(eset)$Drug
		timepoint_data <-as.numeric(pData(eset)$Time)
		timepoint_data<-sub(6,0,timepoint_data)
		timepoint_data<-sub(12,1,timepoint_data)
		timepoint_data<-sub(24,2,timepoint_data)
		timepoint_data<-as.numeric(timepoint_data)
    replicate_data<-rep(NA,length(timepoint_data))
    for(i in 0:2){
      replicate_data[grepl("DMSO",condition_data) & timepoint_data==i]<-seq(1,2)
	  	replicate_data[!grepl("DMSO",condition_data) & timepoint_data==i]<-seq(1,2)
    }
 	}
	options(warn=-1) # Dont show warnings message, most cases they are simply variance equal null, so ignore them
	prob <- betr(eset, cond=condition_data,timepoint=timepoint_data,replicate=replicate_data, alpha=0.05, twoCondition=twoCondition_flag,verbose=FALSE)
	cat(drug_test," vs ",ctrl," - ",length(which(threshold<prob)), "genes detected\n")
	probe_id<-rownames(as.data.frame(which(threshold<prob)))
	return (probe_id)
}

#2nd step - Probeset-Drug concentration model fitting
model_fit <- function(
eset,
betr_result,
timepoint=6,
drug_flag=TRUE, # Use drug data instead of assuming IC vector c(0,0.02,0.2)
report=TRUE,
choice='E2' # List: E2, Williams, Marcus, M, ModM
)
{
  druglist<-unique(pData(eset)$Drug)
  druglist<-druglist[!grepl('DMSO',unique(pData(eset)$Drug))]
  n<-length(druglist)
  output<-list()
  
	if(report) {model_report<-data.frame(drug=character(n),detected=numeric(n),stringsAsFactors=FALSE)}

	for(i in 1:length(druglist)){
		fit_result<-loop_model_fit(eset,druglist[i],as.vector(unlist(betr_result[druglist[i]])),
                               drug_flag,timepoint,choice)
		output<-append(output,list(fit_result))	  
		if(report) {model_report[i,]<-c(druglist[i],length(fit_result[,1]))}
	}
	names(output) <-druglist
	if(report) {write.table(model_report,file='model_report.tab', sep="\t")}
	return(output)
}

loop_model_fit <- function(eset,drug,probeset,drug_flag,timepoint,choice) {
	
	# Libraries and options
	library(IsoGene)
	library(drc)
  library(tcltk)			
	set.seed(1234)
  
  # Preparation Data
	eset<-get.array.subset.exprset(eset,"Drug", c("DMSO",drug))
	eset<-get.array.subset.exprset(eset,"Time", 6)

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
  
	Expression<-as.data.frame(exprs(eset[probeset,])) 
  
  if(grepl(''))
	IsoSAM.obj<-IsoTestSAM(IC, y=Expression, fudge="pooled", niter=100,FDR=0.05, stat=choice) 
  significant_probes<-as.vector(unlist(IsoSAM.obj[[1]]$Probe.ID))		
	slope<-IsoGenem(IC,Expression[significant_probes,])$direction
  output<-data.frame(weight=rep(1,length(slope)),direction=slope)

	cat(drug,' - ',length(probeset),'probesets - ',nrow(output),' remaining\n')
	return (output)
}

#3rd step - Features selection by random-forest to shrink signatures to most important components
signature_optimization <- function(drug,fit_result,no_overlap=FALSE,optimize=TRUE,min_comp=30,max_comp=150) {
  
  library(caret)
  library(mlbench)
  library(Hmisc)
  library(randomForest)
  library(dendextend)
  
    subset<-get.array.subset.exprset(eset,"Drug", c(drug,"DMSO"))
    subset<-get.array.subset.exprset(subset,"Concentration", c("IC20","NONE"))
    replicate<-exprs(subset)
    replicate<-replicate[rownames(fit_result[[drug]]),]
    n<-length(rownames(fit_result[[drug]]))
    replicate<-data.frame(t(replicate))
    info<-as.factor(pData(subset)$Drug)
    
    # Prenormalization with center and scaling
    normalization<-preProcess(replicate)
    x <- predict(normalization, replicate)
    x <- as.data.frame(x)
    
    subsets <- seq(min_comp,max_comp,20)
    set.seed(10)
    
    ctrl <- rfeControl(functions = rfFuncs,
                       method = "repeatedcv",
                       number=10, repeats=10,
                       verbose = FALSE)
    
    # Start random forest
    rfProfile <- rfe(x,info,
                     sizes = subsets,
                     rfeControl = ctrl,
                     allowParallel = TRUE)
   
    optimal_signature<-predictors(rfProfile)  
    groupCodes<-pData(subset)$Drug
    rownames(replicate)<-make.unique(groupCodes)
    par(mfrow=c(1,2))
    plot(hclust(dist(replicate)))
    plot(hclust(dist(replicate[,optimal_signature])))
    optimal_signature<-sub("X","",optimal_signature)
    if(nrow(fit_result[[drug]][optimal_signature,])==length(optimal_signature) &
         abs(table(fit_result[[drug]][optimal_signature,])[1]-table(fit_result[[drug]][optimal_signature,])[2])<length(optimal_signature)*.8) {
      fit_result[[drug]]<-fit_result[[drug]][optimal_signature,]
    } else {
      cat("Problem with signatures")
    }
    cat(paste(drug,' signature - ',n,' parameters reduced to ',length(optimal_signature),' parameters',sep=""))
return (fit_result)
}

signature_formating <- function(fit_result,drug) {
  temp<-fit_result[[drug]]
  output<-list(rownames(temp[temp$direction=='d',]))
  output<-append(output,list(rownames(temp[temp$direction=='u',])))
  names(output)<-c('DES','DSS')
  return (output)  
}

#4th step - Subpopulation targeting - Generate interaction matrix
population_targeting <- function(eset,signature_dataset,within.sig=FALSE,nresmpl=100,pval=.05)
{
  druglist<-unique(pData(eset)$Drug)
  druglist<-druglist[!grepl('DMSO',unique(pData(eset)$Drug))]
  
	interaction.des<-data.frame(matrix(data=0,ncol=length(druglist), nrow=length(druglist), dimnames=list(druglist,druglist)))
  interaction.dss<-interaction.des
  
	for(i in 1:length(druglist)) {
			signature<-signature_formating(signature_dataset,druglist[i])				
      targeting<-similarity(eset,signature,within.sig,nresmpl)
      
      # Process DES
      effect<-subset(targeting,des_pval<pval)
      samples<-unique(as.vector(unlist(effect$drug)))
      for(j in 1:length(samples)) {
        interaction.des[druglist[i],samples[j]]<-mean(subset(effect,grepl(samples[j],drug))$des)  
      }
      
      # Process DSS
			effect<-subset(targeting,dss_pval<pval)
			samples<-unique(as.vector(unlist(effect$drug)))
			for(j in 1:length(samples)) {
			  interaction.dss[druglist[i],samples[j]]<-mean(subset(effect,grepl(samples[j],drug))$dss)  
			}
	}
  
  # Combination interaction matrices
  des<--interaction.des
  des[des<0]<-0
  dss<-interaction.dss
  dss[dss<0]<-0
  interaction<-des+dss
  
  output<-list(interaction)
  output<-append(output,list(interaction.des))
  output<-append(output,list(interaction.dss))
  names(output)<-c('targeting','des','dss')
  return (output)
}
  

similarity<-function(eset,signature,within.sig,nresmpl) {
	# Advanced setting
	set.seed(1234)

	eset_drug<-get.array.subset.exprset(eset, "Concentration","IC20")
	eset_drug<-get.array.subset.exprset(eset_drug, "Time","6")
	
  # Matrix expression #
  data<-exprs(eset_drug)
  
  if(within.sig){data<-data[as.vector(unlist(signature)),]}
  
  # row normalize std
  data.mean <- apply(data[,],1,mean,na.rm=T)
  data.sd <- apply(data[,],1,sd,na.rm=T)
  norm_data<-as.data.frame((data-data.mean)/data.sd)

	cos.sim <- function(x) {return(sum(x)/sqrt(sum(x^2))*sqrt(length(x)))}
  
  des<-1-apply(norm_data[signature[['DES']],],2,cos.sim)
	dss<-1-apply(norm_data[signature[['DSS']],],2,cos.sim)
  
  # Compute permutation distance for DES
  nominal.p.des<-rep(NA,ncol(norm_data))
	nominal.p.dss<-nominal.p.des
  
  for(i in 1:ncol(norm_data)){
    rnd.feature.matrix.DES<-matrix(0,nrow=length(signature[['DES']]),ncol=nresmpl)
    rnd.feature.matrix.DSS<-matrix(0,nrow=length(signature[['DSS']]),ncol=nresmpl)
    for (p in 1:nresmpl){
      rnd.feature.matrix.DES[,p]<-sample(norm_data[,i],length(signature[['DES']]),replace=F)
      rnd.feature.matrix.DSS[,p]<-sample(norm_data[,i],length(signature[['DSS']]),replace=F)
    }
    perm.des<-abs(1-apply(rnd.feature.matrix.DES,2,cos.sim))
    perm.dss<-abs(1-apply(rnd.feature.matrix.DSS,2,cos.sim))
    
    stat_rank.des<-rank(-c(abs(des[i]),perm.des))
    stat_rank.dss<-rank(-c(abs(dss[i]),perm.dss))
    
    nominal.p.des[i]<-stat_rank.des[1]/length(stat_rank.des)
    nominal.p.dss[i]<-stat_rank.dss[1]/length(stat_rank.dss)
  }

	# MCT correction
	BH.FDR.DES<-nominal.p.des*ncol(norm_data)/rank(-nominal.p.des)
	BH.FDR.DES[BH.FDR.DES>1]<-1

	BH.FDR.DSS<-nominal.p.dss*ncol(norm_data)/rank(-nominal.p.dss)
	BH.FDR.DSS[BH.FDR.DSS>1]<-1
  
	output<-data.frame(name=sampleNames(eset_drug),
                     drug=pData(eset_drug)$Drug,
                     des=des,
                     dss=dss,
                     des_pval=BH.FDR.DES,
	                   dss_pval=BH.FDR.DSS
                     )
	return (output)
}

prediction<-function(eset,subpopulations,targeting_matching,subpopulation=TRUE,spectrum=TRUE,scale=TRUE,top=10){
  
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