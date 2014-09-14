# additional approach - First determine the number of subpopulations by subtypes the samples
nmf_test <- function(
eset,
sd=0.9
)
{
  library(NMF)
  library(fastICA)
  time=unique(eset$Time)
  druglist=initialize_drug(eset)
  ICeset<-eset[,eset@phenoData@data$Concentration=='IC20']
  
  #test ranking 42 subtypes for different SD
  set.seed(1234)
  data<-exprs(ICeset)
  SampleVar <- apply(data[1:nrow(data),],1,var)
  sd=.5
  test<-list()
  for(x in 0:9) { 
    idxHighVar <- SampleVar > quantile(SampleVar,probs = (sd+.05*x)
    data_highSD <- data[idxHighVar,]
    cl_sample <- nmf(data_highSD,rank=42,method="lee",seed='random',.opt="-v+P",nrun=30)
    test<-append(test,cl_sample)  
  }
  
  
  # Timewise clustering
  # Determine the reduction of heterogeneity over time, ideally 9 at t1, 7 at t2 and 3 at t3
  # Min 2, Max 14 
  set.seed(1234)
  clustering_data<-list()
  clustering_control<-list()
  
  for(t in 1:length(time)){
    subset<-ICeset[,ICeset@phenoData@data$Time==time[t]]
    data<-exprs(subset)
    cat("Removing genes with low sample variance...\n")
    SampleVar <- apply(data[1:nrow(data),],1,var)
    idxHighVar <- SampleVar > quantile(SampleVar,probs = sd)
    data_highSD <- data[idxHighVar,]
    
    max_ranking <- length(druglist)
    data_random<-randomize(data_highSD)
  
    for(i in 2:max_ranking){
      cat(paste("Evaluating ",i," metagroups\n",sep=""))
      cl_sample <- nmf(data_highSD,rank=i,method="lee",seed='random',.opt="-v+P",nrun=30)
      cl_random <- nmf(data_random,rank=i,method="lee",seed='random',.opt="-v+P",nrun=30)
      clustering_data<-append(clustering_data,list(cl_sample))
      clustering_control<-append(clustering_control,list(cl_random))
      }
    } # End ranking evaluation
  } # End Timewise clustering
  
  save(file='timewise clustering data.rda',clustering_data)
  save(file='timewise clustering ctrl.rda',clustering_control)
  
  # Drugwise clustering
  # Determine the drug effect on heterogeneity over time
  # Min 1, max 3
  set.seed(1234)
  clustering_data2<-list()
  clustering_control2<-list()

  for(t in 1:length(druglist)){
    subset<-ICeset[,ICeset@phenoData@data$Drug==druglist[t]]
    data<-exprs(subset)
    cat("Removing genes with low sample variance...\n")
    SampleVar <- apply(data[1:nrow(data),],1,var)
    idxHighVar <- SampleVar > quantile(SampleVar,probs = sd)
    data_highSD <- data[idxHighVar,]
    
    max_ranking <- length(time)
    data_random<-randomize(data_highSD)
    
    for(i in 2:max_ranking){
      cat(paste("Evaluating ",i," metagroups\n",sep=""))
      cl_sample <- nmf(data_highSD,rank=i,method="lee",seed='random',.opt="-v+P",nrun=30)
      cl_random <- nmf(data_random,rank=i,method="lee",seed='random',.opt="-v+P",nrun=30)
      clustering_data2<-append(clustering_data2,list(cl_sample))
      clustering_control2<-append(clustering_control2,list(cl_random))
      }
    } # End ranking evaluation
  } # End Drugwise clustering
  
  save(file='drugwise clustering data.rda',clustering_data2)
  save(file='drugwise clustering ctrl.rda',clustering_control2) 
  
  # Process all data, make table with 3 rows (time) and 14 columns (drug)
  # What if replicates are not robust ?
  # Timewise clustering: Colors
  # Drugwise clustering: 
  dRSS_data<-vector()
  dRSS_ctrl<-vector()
  
  for(i in 1:(length(clustering_data)-1)){
    dRSS_data<-append(dRSS_data,residuals(clustering_data[[i]])-residuals(clustering_data[[i+1]]))
    dRSS_ctrl<-append(dRSS_ctrl,residuals(clustering_control[[i]])-residuals(clustering_control[[i+1]]))
  }
  
}

# 1st step - Find differentially expressed genes with Betr package
betr_filter <- function(
eset, #From initialize_eset function
ctrl='DMSO', # Either DMSO, Media or Time (6hr used as control)
threshold=.95,
report=FALSE
)
{
	library(betr)
	druglist<-initialize_drug(eset)
	if(report) {betr_report<-data.frame(drug=NA,detected=NA);}
	output<-list()
	
	for(i in 1:length(druglist)){
		betr_result <-loop_betr_filter(eset,druglist[i],ctrl,threshold)
		output<-c(output,list(betr_result))
		
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
		replicate_data<-pData(eset)$Replicate
 	}
	options(warn=-1) # Dont show warnings message, most cases they are simply variance equal null, so ignore them
	prob <- betr(eset, cond=condition_data,timepoint=timepoint_data,replicate=replicate_data, alpha=0.05, twoCondition=twoCondition_flag,verbose=FALSE)
	cat(drug_test," vs ",ctrl," - ",length(which(threshold<prob)), "genes detected\n")
	probe_id<-rownames(as.data.frame(which(threshold<prob)))
	return (probe_id)
}

#2nd new step - Probeset-drug concentration model fitting from basis NMF factor
model_fit2<- function(
eset,
nmf_result,
model='DRC',
pval=.05,
size=100,
drug_data=FALSE,
report=TRUE
)
{
  druglist<-initialize_drug(eset)
  time<-unique(eset$Time)
  subset<-eset[,eset@phenoData@data$Time==time[1]]
  
  output<-list()
  
  # Generate dataframe of basis vectors of probeset (rows) for all drugs (columns) at timepoint 6hr
  
  for(i in 1:length(druglist)){
    # First, find the subtype of each drug, extract basis vector and sort by contribution
    # Check just at first timepoint
    
    drug_subset<-subset[,subset@phenoData@data$Drug==druglist[i]]
    fit_result<-loop_model_fit(drug_subset,nmf_basis,model,pval,size,drug_data)
    output<-append(output,list(fit_result))
  }
  return (output)
}
  
  
#2nd step - Probeset-Drug concentration model fitting
model_fit <- function(
eset,
betr_result, # List with same order as druglist
model='linear', # Linear, DRC, NLS, IsoGene
pval_anova=.05, # Pvalue for model fitting
drug_data=FALSE, # Use drug data instead of assuming IC vector c(0,0.02,0.2)
report=TRUE
)
{
	druglist<-initialize_drug(eset)
	if(report) {model_report<-data.frame(drug=NA,detected=NA)}
	output<-list()
	total<-length(druglist)
	pb <- winProgressBar(title = "progress bar", min = 0,max = total, width = 300)	
	for(i in 1:length(druglist)){
	  current<-betr_result[i]
		fit_result<-loop_model_fit(eset,druglist[i],current,model,pval_anova,drug_data) # Fit result is dataframe ?
		output<-c(output,list(fit_result))
	  setWinProgressBar(pb, i, title=paste( round(i/total*100, 0),"% done"))
	  
		if(report) {
			model_report<-rbind(model_report,c(druglist[i],length(fit_result[,1])))
			#fit_graph(eset,druglist[i],pattern_result,25) Not working yet
		}
	}
	names(output) <-druglist
	
	if(report) {
		model_report<-model_report[-1,]
		write.table(model_report,file='model_report.tab', sep="\t")	
	}
	close(pb)
	return(output)
}

loop_model_fit <- function(eset,drug,probeset,model,pval_anova,drug_flag) {
	
	# Libraries and options
	library(IsoGene)
	library(drc)
	options(show.error.messages = FALSE)				
	set.seed(1234)
	
	# Preparation Output
	output<-data.frame(A=NA,B=NA,C=NA)
	class(output)<-"numeric"
	
	# Preparation Data
	eset<-get.array.subset.exprset(eset,"Drug", c("DMSO",drug))
	timelist<-unique(as.vector(eset$Time))
	probeset<-as.vector(unlist(probeset))
	threshold<-.75
	
	# Used if NLS model, Isogene or DRC
	stat<-data.frame(probeset=NULL,time=NULL,quality=NULL,slope=NULL)
	
	# Prepare dosage vectors with data instead of IC
	if(drug_flag) {
		temp_data<-read.table(file='inputs/drug_input.txt',sep='\t',header=TRUE)
		drug_data<-as.vector(temp_data[,2])
		names(drug_data)<-temp_data[,1]
		experimental_ic<-as.numeric(drug_data[drug])
	}
	
	if(grepl('IsoGene',model)) {
		for(j in 1:length(timelist)) {
			loop_subset<-get.array.subset.exprset(eset,"Time", timelist[j])
			IC<-loop_subset$Concentration
			for(z in 1:length(IC)) {
				if(drug_flag) {
					if(IC[z]=='NONE') {IC[z]=0}
					if(IC[z]=='IC20') {IC[z]=experimental_ic}	
					if(IC[z]=='IC20_LOW') {IC[z]=experimental_ic*.1}
				} else {
					if(IC[z]=='NONE') {IC[z]=0}
					if(IC[z]=='IC20') {IC[z]=0.2}	
					if(IC[z]=='IC20_LOW') {IC[z]=0.02}
				}
			}
			class(IC)<-"numeric"
			Expression<-as.data.frame(exprs(loop_subset[probeset,]))
			rawpvalue<-IsoRawp(IC,Expression,niter=1000,progressBar=FALSE)
			isodata<-rawpvalue[[2]]
			signif<-IsoTestBH(isodata,FDR=0.05,type="BY",stat="ModifM")
			signif<-subset(isodata,ModM<0.05)
			for(i in 1:length(signif[,1])) {
				current<-as.vector(unlist(Expression[signif[i,1],]))
				m<-vector() # Median vector
				slope<-0
				for(xxx in 1:3) {
					m<-append(m,median(current[grepl(unique(sort(IC))[xxx],IC)]))
				}
			
				if(m[1]>m[2] && m[2]>m[3]) {slope<-1}
				if(m[1]<m[2] && m[2]<m[3]) {slope<--1}
				adding_row<-list(signif[i,1],as.numeric(j),'0',as.numeric(slope))
				stat<-factor_dataframe(rbind(stat,adding_row))
				colnames(stat)<-c('probeset','time','quality','slope')
			}
		} # End time loop		
	# Remove useless data
	stat<-stat[stat$slope!=0,] 
	} # End Isogene condition
	else {
		# All others method proceed probeset by probeset
		for(i in 1:length(probeset)){
			dose_effect<-NULL
			for(j in 1:length(timelist)) {
		
				loop_subset<-get.array.subset.exprset(eset,"Time", timelist[j])
				IC<-loop_subset$Concentration
				Expression<-as.vector(exprs(loop_subset[probeset[i],]))
			
				for(z in 1:length(IC)) {
					if(drug_flag && !grepl('DRC',model)) {
						if(IC[z]=='NONE') {IC[z]=0}
						if(IC[z]=='IC20') {IC[z]=experimental_ic}	
						if(IC[z]=='IC20_LOW') {IC[z]=experimental_ic*.1}
					} else {
						if(IC[z]=='NONE') {IC[z]=0}
						if(IC[z]=='IC20') {IC[z]=0.2}	
						if(IC[z]=='IC20_LOW') {IC[z]=0.02}
					}
				}
				class(IC)<-"numeric"
				class(Expression)<-"numeric" 

				# Option linear
				if(grepl('linear',model)){
					mm_fit<-lm(Expression~IC)
					pvalue<-as.vector(unlist(anova(mm_fit)["Pr(>F)"]))[1]
					if(pvalue<pval_anova) {
							slope<-coef(mm_fit)[2]
							if(slope>0) {dose_effect[j]<- 1}
							else {dose_effect[j]<- -1}
					} else {dose_effect[j]<- 0}
				}
				
				# Option linear
				if(grepl('NLS',model)) {
					mm_fit<-itfit<-try(nls(Expression ~ cbind(1,(IC/(c+IC))),start=list(c=0.1),data=as.data.frame(cbind(IC,Expression)),alg = "plinear"),silent=TRUE)
					# Test convergence
					if(!inherits(itfit,'try-error')) {
						if(coef(mm_fit)[1]>0 && coef(mm_fit)[1]<1) {
							rss<-sum(residuals(mm_fit)^2)
							tss<-0
							for(xxx in 1:3) {
								m<-mean(Expression[grepl(unique(IC)[xxx],IC)])
								tss<-tss+sum((Expression[grepl(unique(IC)[xxx],IC)]-m)^2)
							}
							final<-1-rss/tss
							# Coef 3 is Vmax-Vmin, unlike DRC
							if(coef(mm_fit)[3]>0){slope=1} else {slope=-1}
							
							adding_row<-list(probeset[i],as.numeric(j),as.numeric(final),as.numeric(slope))
							stat<-factor_dataframe(rbind(stat,adding_row))
							colnames(stat)<-c('probeset','time','quality','slope')
						}
					}	
				}
				
				# Option DRC
				if(grepl('DRC', model)) {
					mm_fit<-itfit<-try(drm(Expression~IC,fct=MM.3(),control=drmc(noMessage=TRUE)),silent=TRUE)
					# Test convergence
					if(!inherits(itfit,'try-error')) {
						if(as.numeric(coef(mm_fit)[1])<as.numeric(coef(mm_fit)[2])) {slope=1} else {slope=-1}
						adding_row<-list(probeset[i],as.numeric(j),as.numeric(noEffect(mm_fit)[3]),as.numeric(slope))
						stat<-factor_dataframe(rbind(stat,adding_row))
						colnames(stat)<-c('probeset','time','pval','slope')
					}
				}
				
			} # End loop time
			
			# For linear model (if not, output is overwritten later...)
			final<-c(dose_effect[1],dose_effect[2],dose_effect[3])
			output<-rbind(output,final)
			rownames(output)[nrow(output)]<-probeset[i]
			
		} # End loop probeset
	} # End Algorithm

	# Process stats from different models and save into output
	
	# Linear - Filter out probeset without any fitting at any timepoint
	if(grepl('linear',model)) {
		output<-output[-1,] # Eliminate output entry
		disc<-vector()
		for(i in 1:length(output[,1])) {
			if(sum(!grepl(0,output[i,]))<1) {
				disc<-append(disc,i)
			}
		}
		output<-output[-disc,]
	}
	
	# NLS - Perform quantile threshold to get balanced signature
	if(grepl('NLS',model)) {
		probeset_up<-list()
		probeset_down<-list()
		for( t in 1:3) {
			current<-stat[stat$time==t,]
			
			current_up<-current[current$slope==1,]
			thrs_up<-as.numeric(quantile(current_up$quality,p=threshold))
			probeset_up<-c(probeset_up,list(current_up[current_up$quality>thrs_up,'probeset']))
			
			current_down<-current[current$slope==-1,]
			thrs_down<-as.numeric(quantile(current_down$quality,p=threshold))
			probeset_down<-c(probeset_down,list(current_down[current_down$quality>thrs_down,'probeset']))
		}
		
		all_probes<-unique(append(as.vector(unlist(probeset_up)),as.vector(unlist(probeset_down))))
		output<-as.data.frame(matrix(0,nrow=length(all_probes),ncol=3))
		colnames(output)<-c("A","B","C")
		rownames(output)<-all_probes
		
		for( t in 1:3) {
			output[as.vector(unlist(probeset_up[t])),t]<-1
			output[as.vector(unlist(probeset_down[t])),t]<--1
		}
	}
	
	if(grepl('DRC',model)) {
		probeset_up<-list()
		probeset_down<-list()
		for( t in 1:3) {
			filtered<-stat[stat$time==t,]
      current<-subset(filtered,pval<pval_anova/length(filtered[,1])) # BH Correction
			probeset_up<-c(probeset_up,list(current[current$slope==1,'probeset']))		
			probeset_down<-c(probeset_down,list(current[current$slope==-1,'probeset']))
		}
		
		all_probes<-unique(append(as.vector(unlist(probeset_up)),as.vector(unlist(probeset_down))))
		output<-as.data.frame(matrix(0,nrow=length(all_probes),ncol=3))
		colnames(output)<-c("A","B","C")
		rownames(output)<-all_probes
		
		for( t in 1:3) {
			output[as.vector(unlist(probeset_up[t])),t]<-1
			output[as.vector(unlist(probeset_down[t])),t]<--1
		}
	} # End DRC code
	
	if(grepl('IsoGene',model)) {
		probeset_up<-list()
		probeset_down<-list()
		for( t in 1:3) {
			current<-stat[stat$time==t,]
			probeset_up<-c(probeset_up,list(current[current$slope==1,'probeset']))		
			probeset_down<-c(probeset_down,list(current[current$slope==-1,'probeset']))
		}
		
		all_probes<-unique(append(as.vector(unlist(probeset_up)),as.vector(unlist(probeset_down))))
		output<-as.data.frame(matrix(0,nrow=length(all_probes),ncol=3))
		colnames(output)<-c("A","B","C")
		rownames(output)<-all_probes
		
		for( t in 1:3) {
			output[as.vector(unlist(probeset_up[t])),t]<-1
			output[as.vector(unlist(probeset_down[t])),t]<--1
		}
	
	} # End IsoGene code

	cat(drug,' - ',length(probeset),'probesets - ',length(output[,1]),' remaining\n')
	return (output)
}

#3rd step - Subpopulations definitions
capture_signature <- function(drug,fit_result,timepoint,permanent_effect,distinct_signature) {
	druglist<-initialize_drug(eset)
	DSS<-as.vector(NULL) # Up genes -->  increased expression of genes with increasing drug concentration --> Subpopulation spared
	DES<-as.vector(NULL) # Down genes -->  decreased expression of genes with increasing drug concentration --> Subpopulation eliminated
	all_probes<-vector()
	
  for(i in 1:length(druglist)){
		if(grepl(drug,druglist[i])) {
			convert<-as.data.frame(fit_result[i])
			if(nrow(convert)>0) {
  			current<-convert[,timepoint]
  			names(current)<-rownames(convert)
  			for(j in 1:length(current)) {
  				probeset<-names(current)[j]
  				
  				# Check if not equal 0
  				if(current[j]!=0) {
  					# No time-restriction of probeset
  					if(!permanent_effect || timepoint==3) {
  						if(current[j]>0) {DSS<-append(DSS,probeset)}
  						else {DES<-append(DES,probeset)}
  					} else {	
  					# Permanent effect over time
  					
  						#Detect 1st timepoint with any correlation
  						flag<-FALSE
  						for(t in 1:3) {
  							if(convert[probeset,t]!=0) {
  								flag<-TRUE
  								correlation_sign<-convert[probeset,t]
  								break
  							}
  						}
  						
  						#Check if next timepoints are all correlated the same way
  						if(flag) {
  							flag2<-TRUE
  							for(s in (t+1):3) {
  								if(convert[probeset,s]!=correlation_sign) {flag2<-FALSE;break}}
  						}
  						
  						#Sustained drug-probeset correlation over time eg. 0 1 1
  						if(flag2) {
  							if(current[j]>0) {DSS<-append(DSS,probeset)}
  							else {DES<-append(DES,probeset)}
  						}	
  					} 
				  }
			  }
		  }
		}
    all_probes<-append(all_probes,rownames(fit_result[[i]]))
	}
	output<-list()
  if(distinct_signature) {
    all_probes<-all_probes[duplicated(all_probes)]
    output<-c(output,list(setdiff(DSS,all_probes)),list(setdiff(DES,all_probes))) 
  } else {
	  output<-c(output,list(DSS),list(DES))
  }
	names(output)<-c('DSS','DES')
	return (output)
}

#4th step - Subpopulation targeting - Generate interaction matrix
population_targeting <- function(eset,model,permanent_effect=FALSE,time_matching=FALSE,drug_specific=FALSE,within.sig=FALSE,distinct_signature=FALSE,inverted_score=FALSE,
limit_num_features=5,
limit_FDR=0.05,
report=FALSE,
inverted_eset=NULL
)
{
	druglist<-initialize_drug(eset)
	interaction<-as.data.frame(matrix(ncol=length(druglist), nrow=length(druglist), dimnames=list(druglist,druglist)))
	interaction[is.na(interaction)]<-0 
	if(report) {
		interaction_report_DES_t1<-list()
		iot_DES<-as.data.frame(matrix(ncol=3, nrow=length(druglist), dimnames=list(druglist,c('t1','t2','t3'))))
		iot_DES[is.na(iot_DES)]<-0
		for(i in 1:length(druglist)) {
			interaction_report_DES_t1[[i]]<-iot_DES
		}
		names(interaction_report_DES_t1)<-druglist
		interaction_report_DES_t2<-interaction_report_DES_t1
		interaction_report_DES_t3<-interaction_report_DES_t1
		interaction_report_DSS_t1<-interaction_report_DES_t1
		interaction_report_DSS_t2<-interaction_report_DES_t1
		interaction_report_DSS_t3<-interaction_report_DES_t1
	}
	
	for(i in 1:length(druglist)) {
		for(t in 1:3) {
			signature<-capture_signature(druglist[i],model,t,permanent_effect,distinct_signature)
			
			if(length(signature[['DSS']])>limit_num_features && length(signature[['DES']])>limit_num_features){
							
				targeting<-NTPez(eset,druglist[i],druglist,signature,t,time_matching,drug_specific,within.sig)
				targeting<-factor_dataframe(subset(targeting,FDR<limit_FDR))
        
				complement_drug<-unique(targeting$drug)
				
				if(length(complement_drug)>0) {
					for(z in 1:length(complement_drug)) {
						complement_targeting<-subset(targeting,grepl(complement_drug[z],drug))
						
						#DES=decreased expression with increasing drug concentration 2 --> Subpopulation eliminated by drug1 --> good combination if DES matching in other treated samples
						DES<-median(subset(complement_targeting,grepl('DES',label))[,'distance'])
						if(is.na(DES) || length(signature[['DES']])<limit_num_features) {DES<-0}
						
						#DSS=increased expression with increasing drug concentration 1 --> Subpopulation spared by drug1 --> not useful if DSS matching in other treated samples
						DSS<-median(subset(complement_targeting,grepl('DSS',label))[,'distance'])
						if(is.na(DSS) || length(signature[['DES']])<limit_num_features) {DSS<-0}
						
						if(inverted_score) {
							interaction[druglist[i],complement_drug[z]]<-interaction[druglist[i],complement_drug[z]]+DES
						} else {
							interaction[druglist[i],complement_drug[z]]<-interaction[druglist[i],complement_drug[z]]+DES-DSS
						}
						if(report) {
						
							time_matching2<-unique(complement_targeting$timepoint)
							for(o in 1:length(time_matching2)) {
								complement_targeting2<-subset(complement_targeting,timepoint==time_matching2[o])
								
								DES<-median(subset(complement_targeting2,grepl('DES',label))[,'distance'])
								if(is.na(DES) || length(signature[['DES']])<limit_num_features) {DES<-0}
								DSS<-median(subset(complement_targeting2,grepl('DSS',label))[,'distance'])
								if(is.na(DSS) || length(signature[['DES']])<limit_num_features) {DSS<-0}
								
								if(t==1) {
								interaction_report_DES_t1[[druglist[i]]][complement_drug[z],o]<-DES
								interaction_report_DSS_t1[[druglist[i]]][complement_drug[z],o]<-DSS
								}
								if(t==2) {
								interaction_report_DES_t2[[druglist[i]]][complement_drug[z],o]<-DES
								interaction_report_DSS_t2[[druglist[i]]][complement_drug[z],o]<-DSS
								}
								if(t==3) {
								interaction_report_DES_t3[[druglist[i]]][complement_drug[z],o]<-DES
								interaction_report_DSS_t3[[druglist[i]]][complement_drug[z],o]<-DSS
								}
							}
						}			
					}
				}
				
				if(inverted_score) {
					# Parameter Inverted_score on !
					# DSS substracted only when found overpresent in inverted Eset
					targeting<-NTPez(inverted_eset,druglist[i],druglist,signature,t,time_matching,drug_specific,within.sig)
					targeting<-factor_dataframe(subset(targeting,FDR<limit_FDR))
					complement_drug<-unique(targeting$drug)
				
					if(length(complement_drug)>0) {
						for(z in 1:length(complement_drug)) {
							complement_targeting<-subset(targeting,grepl(complement_drug[z],drug))
							
							#DES=decreased expression with increasing drug concentration 2 --> Subpopulation eliminated by drug1 --> good combination if DES matching in other treated samples
							DES<-median(subset(complement_targeting,grepl('DES',label))[,'distance'])
							if(is.na(DES) || length(signature[['DES']])<limit_num_features) {DES<-0}
							
							#DSS=increased expression with increasing drug concentration 1 --> Subpopulation spared by drug1 --> not useful if DSS matching in other treated samples
							DSS<-median(subset(complement_targeting,grepl('DSS',label))[,'distance'])
							if(is.na(DSS) || length(signature[['DES']])<limit_num_features) {DSS<-0}
							
							interaction[druglist[i],complement_drug[z]]<-interaction[druglist[i],complement_drug[z]]-DSS
						}
					}
				}	   
			}
		}
	}
	
	#Make final 3x3 matrix
	if(report) {
		save(interaction_report_DES_t1,file='DES interaction report at t1.rda')
		save(interaction_report_DES_t2,file='DES interaction report at t2.rda')
		save(interaction_report_DES_t3,file='DES interaction report at t3.rda')
		save(interaction_report_DSS_t1,file='DSS interaction report at t1.rda')
		save(interaction_report_DSS_t2,file='DSS interaction report at t2.rda')
		save(interaction_report_DSS_t3,file='DSS interaction report at t3.rda')
	}	
	return (interaction)
}

NTPez<-function(eset,drug,druglist,signature,timepoint,time_matching,drug_specific,within.sig) {
	# Advanced setting
	set.seed(7392854)
	nresmpl<-1000

	# Data timepoint matching with pData(eset)
	if(timepoint==1) {timepoint=6}
	if(timepoint==2) {timepoint=12}
	if(timepoint==3) {timepoint=24}
	
	class_matrix<-data.frame(feature=NULL,group=NULL)
	### design signatures matrix ###
	for(i in 1:length(signature)) {
		class_name<-names(signature[i])
		feature_class<-as.vector(unlist(signature[i]))
		to_merge<-data.frame(feature=feature_class,group=rep(i,length(feature_class)))
		class_matrix<-rbind(class_matrix,to_merge)
	}
	class_matrix<-factor_dataframe(class_matrix) #Remove factor definition within dataframe

	#Format expression data from eset (restriction to IC20 data, no Media and DMSO data
	#	two parameters:
	#	- time_matching: if TRUE, DES and DSS at a specific timepoint will be evaluated against all samples at the same timepoint
	#	- drug_specific: if TRUE, differential expression of each drugs against all drugs will be calculated prior analysis
	eset_drug<-get.array.subset.exprset(eset, "Concentration","IC20")
	if(time_matching) {
		eset_drug<-get.array.subset.exprset(eset_drug, "Time",timepoint)
	}

	if(drug_specific) {
		eset_drug<-get.array.subset.exprset(eset_drug,"Drug",druglist)
		library(limma)
		group<-factor(pData(eset_drug)[,'Drug'])
		design<-model.matrix(~0+group)
		colnames(design)<-sub("group","",colnames(design))
		fit<-lmFit(eset_drug,design)
		
		#Loop contrast matrix
		druglist2<-setdiff(druglist,drug)
		for(i in 1:length(druglist2)) {
			comparaison<-paste(druglist2[i],'-',drug,sep="")
			comparaison<-as.vector(comparaison)
			contrast<-makeContrasts(contrasts=comparaison,levels=design)
			fit2<-contrasts.fit(fit,contrast)
			fit2<-eBayes(fit2)
			tab<-topTable(fit2,coef=1,adjust="fdr",number=nrow(fit2))[,c("ID","logFC")]
			extract<-tab$logFC
			names(extract)<-tab$ID
			extract<-extract[order(names(extract))]
			if(i==1) {
				exp.dataset<-data.frame(id=names(extract))
			}
			exp.dataset[,druglist2[i]]<-extract
		}
		rownames(exp.dataset)<-exp.dataset$id
		exp.dataset$id<-NULL 
	} else {
		samples_drug<-rownames(subset(pData(eset_drug),grepl(drug,Drug)))
		druglist2<-druglist
		eset_drug<-get.array.subset.exprset(eset_drug,"Drug",druglist2)
		exp.dataset<-exprs(eset_drug)
		samples<-setdiff(colnames(exp.dataset),samples_drug)
	}
			
	num.features<-length(unique(class_matrix$feature))
	num.cls<-length(signature)
	num.samples<-length(colnames(exp.dataset))
	# row normalize std
	normed.exp.dataset<-exp.dataset
	if(drug_specific) {
	exp.mean <- apply(exp.dataset[,],1,mean,na.rm=T)
	exp.sd <- apply(exp.dataset[,],1,sd,na.rm=T)	
	} else {
	exp.mean <- apply(exp.dataset[,samples],1,mean,na.rm=T)
	exp.sd <- apply(exp.dataset[,samples],1,sd,na.rm=T)
	}
	#for(i in 1:length(exp.sd)) { if(exp.sd[i]==0) {exp.sd[i]=0.00577350}} required if sample nb is low
	normed.exp.dataset<-as.data.frame((exp.dataset-exp.mean)/exp.sd)

	# Shrink eset to features of class_matrix only
	exp.dataset.extract<-normed.exp.dataset[class_matrix$feature,]	 

	#### Start algorithm #### 
	# Design Group-specific templates
	for (i in 1:num.cls){
		temp.temp<-as.numeric(as.vector(class_matrix$group))
		temp.temp[temp.temp!=i]<-0
		temp.temp[temp.temp==i]<-1
		eval(parse(text=paste("temp.",i,"<-temp.temp",sep="")))
	}

	# compute distance and p-value #
	predict.label<-vector(length=num.samples,mode="numeric")
	dist.to.template<-vector(length=num.samples,mode="numeric")
	dist.to.cls1<-vector(length=num.samples,mode="numeric")
	rnd.feature.matrix<-matrix(0,nrow=num.features,ncol=nresmpl)
	perm.dist.vector<-vector(length=nresmpl*num.cls,mode="numeric")
	nominal.p<-vector(length=num.samples,mode="numeric")
	BH.FDR<-vector(length=num.samples,mode="numeric")
	Bonferroni.p<-vector(length=num.samples,mode="numeric")

	for (i in 1:num.samples){
      
			current.sample <- as.numeric(as.vector(unlist(exp.dataset.extract[,i])))

			# compute original distance
			orig.dist.to.all.temp <- vector(length=num.cls,mode="numeric")
			
			# Compute distance to all templates
			for (o in 1:num.cls){
				eval(parse(text=paste("current.temp <- temp.",o,sep="")))
				orig.dist.to.all.temp[o]<-sum(current.temp*current.sample)/(sqrt(sum(current.temp^2))*sqrt(sum(current.sample^2)))
			}
				
		# find nearest neighbors (2 classes)
		 if (num.cls==2){           
			  if (orig.dist.to.all.temp[1]>=orig.dist.to.all.temp[2]){
				predict.label[i]<-1
				dist.to.template[i]<-1-orig.dist.to.all.temp[1]
				dist.to.cls1[i]<--(orig.dist.to.all.temp[1]+1)
			  }
			  if (orig.dist.to.all.temp[1]<orig.dist.to.all.temp[2]){
				predict.label[i]<-2
				dist.to.template[i]<-1-orig.dist.to.all.temp[2]
				dist.to.cls1[i]<-orig.dist.to.all.temp[2]+1
			  }
		 } else {
			for (o in 1:num.cls){       # find nearest neighbor (>2 classes)
				if (is.na(orig.dist.to.all.temp[o])!=T){
					if (orig.dist.to.all.temp[o]==max(orig.dist.to.all.temp,na.rm=T)){
						predict.label[i]<-o
						dist.to.template[i]<-1-orig.dist.to.all.temp[o]
						dist.to.cls1[i]<-(1-orig.dist.to.all.temp[o])+o
					}
				}
			}
		}

		# permutation test
		if (!within.sig){     
		# generate resampled features from all probes in normalized eset
			for (p in 1:nresmpl){
				rnd.feature.matrix[,p]<-sample(normed.exp.dataset[,i],num.features,replace=F)
			}
		} else {
			# generate resampled features from only features in normalized eset
			for (p in 1:nresmpl){
				rnd.feature.matrix[,p]<-sample(exp.dataset.extract[,i],num.features,replace=F)
			}
		}

		# compute distance to all templates
		for (res in 1:num.cls){
			eval(parse(text=paste("temp.resmpl<-temp.",res,sep="")))

			prod.sum<-apply(t(t(rnd.feature.matrix)*temp.resmpl),2,sum)
			data.sq.sum<-apply(rnd.feature.matrix^2,2,sum)
			temp.sq.sum<-sum(temp.resmpl^2)

			perm.dist.vector[(1+(nresmpl*(res-1))):(nresmpl*res)]<-(1-(prod.sum/(sqrt(data.sq.sum)*sqrt(temp.sq.sum))))
		}

		# compute nominal p-value
		combined.stats.rank<-rank(c(dist.to.template[i],perm.dist.vector))
		nominal.p[i]<-combined.stats.rank[1]/length(combined.stats.rank)
	} #End loop samples

	# MCT correction
	BH.FDR<-nominal.p*num.samples/rank(nominal.p)
	Bonferroni.p<-nominal.p*num.samples
	BH.FDR[BH.FDR>1]<-1
	Bonferroni.p[Bonferroni.p>1]<-1
		
	# Format output results
	if(drug_specific) {
		sample_names<-colnames(exp.dataset.extract)
		for(i in 1:num.cls) {
			predict.label<-sub(paste(i),names(signature)[i],predict.label)
		}
		output<-data.frame(name=sample_names,drug=sample_names,label=predict.label,distance=dist.to.template,FDR=BH.FDR,pval=Bonferroni.p)
	} else {
		sample_names<-sampleNames(eset_drug)
		sample_drug<-pData(eset_drug)[sample_names,'Drug']
		sample_timepoint<-pData(eset_drug)[sample_names,'Time']
		for(i in 1:num.cls) {
			predict.label<-sub(paste(i),names(signature)[i],predict.label)
		}
		output<-data.frame(name=sample_names,drug=sample_drug,timepoint=sample_timepoint,label=predict.label,distance=dist.to.template,FDR=BH.FDR,pval=Bonferroni.p)
	}
	return (output)
}

prediction<-function(interaction,drug_concentration=FALSE,partial_hierarchy=FALSE,full_hierarchy=FALSE,parameters_testing=FALSE){
	druglist<-colnames(interaction)
	final_prediction<-data.frame(drug1=character(),drug2=character(),rank=numeric())
	if(drug_concentration) {
		temp_data<-read.table(file='inputs/drug_input.txt',sep='\t',header=TRUE)
		drug_data<-as.vector(temp_data[,2])
		names(drug_data)<-temp_data[,1]
	}
	library(memisc)
	
	if(!partial_hierarchy && !full_hierarchy) {
		for(i in 1:(length(druglist)-1)) {
			for(j in (i+1):length(druglist)) {
				if(drug_concentration) {
					drug_combination<-(1/drug_data[druglist[j]]) * interaction[druglist[i],druglist[j]] + (1/drug_data[druglist[i]]) * interaction[druglist[j],druglist[i]]
				} else {
					drug_combination<-interaction[druglist[i],druglist[j]]+interaction[druglist[j],druglist[i]]
				}
				final_prediction<-rbind(final_prediction,data.frame(drug1=druglist[i],drug2=druglist[j],rank=drug_combination))
			}	
		}
	} else {
	  # Hierarchy based scoring
		targeting_matrix<- collateral_targeting_matrix(interaction,druglist)
		for(i in 1:(length(druglist)-1)) {
			for(j in (i+1):length(druglist)) {
				if(drug_concentration && full_hierarchy) {
					#Case Full Hierarchy + Drug data
					drug_combination<-max((1/drug_data[druglist[j]]) * interaction[druglist[i],druglist[j]] * sum(targeting_matrix[druglist[j],]),0)+max((1/drug_data[druglist[i]]) * interaction[druglist[j],druglist[i]] * sum(targeting_matrix[druglist[i],]),0)
				} else if(!drug_concentration && full_hierarchy) {
					#Case Full Hierarchy
					drug_combination<-max(interaction[druglist[i],druglist[j]] * sum(targeting_matrix[druglist[j],]),0)+max(interaction[druglist[j],druglist[i]] * sum(targeting_matrix[druglist[i],]),0)
				} else if(drug_concentration && partial_hierarchy && sign(interaction[i,j])==sign(interaction[j,i])) {
					#Case Partial Hierarchy, drug data, symetric killing
					drug_combination<-max((1/drug_data[druglist[j]]) * interaction[druglist[i],druglist[j]] * sum(targeting_matrix[druglist[j],]),0)+max((1/drug_data[druglist[i]]) * interaction[druglist[j],druglist[i]] * sum(targeting_matrix[druglist[i],]),0)
				} else if(!drug_concentration && partial_hierarchy && sign(interaction[i,j])==sign(interaction[j,i])) {
					#Case Partial Hierarchy, symetric killing
					drug_combination<-max(interaction[druglist[i],druglist[j]] * sum(targeting_matrix[druglist[j],]),0)+max(interaction[druglist[j],druglist[i]] * sum(targeting_matrix[druglist[i],]),0)	
				} else if(drug_concentration && partial_hierarchy && sign(interaction[i,j])!=sign(interaction[j,i])) {
					#Case Partial Hierarchy, drug data, asymetric killing
					drug_combination<-max((1/drug_data[druglist[j]]) * interaction[druglist[i],druglist[j]] * max(targeting_matrix[druglist[j],]),0)+max((1/drug_data[druglist[i]]) * interaction[druglist[j],druglist[i]] * max(targeting_matrix[druglist[i],]),0)
				} else if(!drug_concentration && partial_hierarchy && sign(interaction[i,j])!=sign(interaction[j,i])) {
					#Case Partial Hierarchy, asymetric killing
					drug_combination<-max(interaction[druglist[i],druglist[j]] * max(targeting_matrix[druglist[j],]),0)+max(interaction[druglist[j],druglist[i]] * max(targeting_matrix[druglist[i],]),0)	
				}
				final_prediction<-rbind(final_prediction,data.frame(drug1=druglist[i],drug2=druglist[j],rank=drug_combination))
			}
		}
	}
	
	rank<-order(final_prediction$rank,decreasing=TRUE)
	final_prediction$rank<-rank
	final_prediction<-final_prediction[with(final_prediction, order(rank)), ]
	rank<-ncbi_scoring(final_prediction)
	if(parameters_testing) {
		return (rank)
	}
	else {
		print(rank)
		return (final_prediction)
	}
}

collateral_targeting_matrix<- function (interaction,druglist) {
	population<-hierarchy_subpopulation(interaction)
	interaction_population<-interaction/population
	interaction_population[!is.finite(as.matrix(interaction_population))]<-0
	population_output<-interaction_population
	for(i in 1:length(druglist)) {
		# Get populations that are killed by a drug, then merge rows
		update_row<-interaction_population[i,]
		next_round<-update_row!=0
		correction<-rowSums(interaction_population[next_round,])!=0
		next_round<-rep(FALSE,length(interaction_population[,1]))
		names(next_round)<-rownames(interaction_population)
		next_round[names(correction)]<-correction
		
		repeat{
			if(sum(next_round)==0) { break }
			
			neighbours<-interaction_population[next_round,]
			neighbours<-neighbours[!duplicated(neighbours),]
			neighbours<-neighbours[which(rowSums(neighbours) > 0),] 
			update_row<-rbind(update_row,colSums(neighbours))
			
			next_round<-colSums(neighbours)!=0
			correction<-rowSums(interaction_population[next_round,])!=0
			next_round<-rep(FALSE,length(interaction_population[,1]))
			names(next_round)<-rownames(interaction_population)
			next_round[names(correction)]<-correction
		}
		population_output[i,]<-colSums(update_row)
	}
	#Maximize diagonal values
	for(i in 1:length(druglist)) {
		population_output[i,i]<-max(population_output)
	}
	return (population_output)
}

hierarchy_subpopulation <- function(data,graph_display=FALSE) {
	library(igraph)
	# Rows are parent nodes, Cols are children nodes
	hierarchy<-data
	hierarchy[,]<-0
	weight<-hierarchy

	# Find drug combinations that are asymetric
	sym_test<-sign(data)
	for(i in 1:length(sym_test[,1])) {
		for(j in 1:length(sym_test[1,])) {
			if (sym_test[i,j]>sym_test[j,i]) {hierarchy[i,j]=1;weight[i,j]=data[i,j]}
		}
	}
	weight<-round(weight,2)
	graph <- graph.adjacency(as.matrix(weight),mode="directed",weighted=T)
	set.seed(112) 
	mst <- minimum.spanning.tree(graph)
	if(graph_display) {
		plot.igraph(mst,layout=layout.fruchterman.reingold,vertex.label.dist=-0.2,vertex.size=5,
			vertex.label=V(mst)$name,	
			main='Hierarchy of subpopulations\nbased on non-symetric killing',
			vertex.label.color="black",
			edge.label=E(mst)$weight,
			edge.label.color="black"
	)
	}
	#Return graph as matrix with weight
	return (as.matrix(get.adjacency(mst,attr='weight')))
}

parameters_optimization<-function(eset,model) {
  
  # Inputs
  inverted_eset<-inverse_eset_fct(eset)
  druglist<-initialize_drug(eset)
  
  # Matrix with all parameters to test
  parameters<-binary_matrix(2000,7)
  parameters<-parameters[!duplicated(parameters),]
  colnames(parameters)<-c('permanent_effect','time_matching','drug_specific','within.sig','drug_concentration','distinct_signature','inverted_score')
  
  # Outputs parameters_optimization
  ncbi_estimation_default<-vector()
  ncbi_estimation_partial_hierarchy<-vector()
  ncbi_estimation_full_hierarchy<-vector()
  
  network_estimation_default<-vector()
  network_estimation_partial_hierarchy<-vector()
  network_estimation_full_hierarchy<-vector()
  
  synergy_estimation_default<-vector()
  synergy_estimation_partial_hierarchy<-vector()
  synergy_estimation_full_hierarchy<-vector()
 
  total<-length(parameters[,1])
  pb <- winProgressBar(title = "progress bar", min = 0,max = total, width = 300)	
  for(i in 1:total) {
    CTT<-parameters[i,]
    print(CTT)
    interaction<-population_targeting(eset,model,permanent_effect=CTT[1],time_matching=CTT[2],drug_specific=CTT[3],within.sig=CTT[4],distinct_signature=CTT[6],inverted_score=CTT[7],inverted_eset=inverted_eset)

	#NCBI Scoring based - Closer to 1, the best
    ncbi_estimation_default<-append(ncbi_estimation_default,prediction(interaction,drug_concentration=CTT[5],partial_hierarchy=FALSE,full_hierarchy=FALSE,parameters_testing=TRUE))
	  ncbi_estimation_partial_hierarchy<-append(ncbi_estimation_partial_hierarchy,prediction(interaction,drug_concentration=CTT[5],partial_hierarchy=TRUE,full_hierarchy=FALSE,parameters_testing=TRUE))
  	ncbi_estimation_full_hierarchy<-append(ncbi_estimation_full_hierarchy,prediction(interaction,drug_concentration=CTT[5],partial_hierarchy=FALSE,full_hierarchy=TRUE,parameters_testing=TRUE))
    
	#Network Scoring based - Closer to 10, the best
	network_estimation_default<-append(network_estimation_default,network_scoring(prediction(interaction,drug_concentration=CTT[5],partial_hierarchy=FALSE,full_hierarchy=FALSE,parameters_testing=FALSE),druglist))
	network_estimation_partial_hierarchy<-append(network_estimation_partial_hierarchy,network_scoring(prediction(interaction,drug_concentration=CTT[5],partial_hierarchy=TRUE,full_hierarchy=FALSE,parameters_testing=FALSE),druglist))
	network_estimation_full_hierarchy<-append(network_estimation_full_hierarchy,network_scoring(prediction(interaction,drug_concentration=CTT[5],partial_hierarchy=FALSE,full_hierarchy=TRUE,parameters_testing=FALSE),druglist))
	
    #Synergy Scoring based - Closer to 100, the best
	synergy_estimation_default<-append(synergy_estimation_default,synergy_scoring(prediction(interaction,drug_concentration=CTT[5],partial_hierarchy=FALSE,full_hierarchy=FALSE,parameters_testing=FALSE),parameters=TRUE))
	synergy_estimation_partial_hierarchy<-append(synergy_estimation_partial_hierarchy,synergy_scoring(prediction(interaction,drug_concentration=CTT[5],partial_hierarchy=TRUE,full_hierarchy=FALSE,parameters_testing=FALSE),parameters=TRUE))
	synergy_estimation_full_hierarchy<-append(synergy_estimation_full_hierarchy,synergy_scoring(prediction(interaction,drug_concentration=CTT[5],partial_hierarchy=FALSE,full_hierarchy=TRUE,parameters_testing=FALSE),parameters=TRUE))
	setWinProgressBar(pb, i, title=paste( round(i/total*100, 0),"% done"))
 }
  
  
  output<-cbind(parameters,ncbi_estimation_default,ncbi_estimation_partial_hierarchy,ncbi_estimation_full_hierarchy)
  output<-cbind(output,network_estimation_default,network_estimation_partial_hierarchy,network_estimation_full_hierarchy)
  output<-cbind(output,synergy_estimation_default,synergy_estimation_partial_hierarchy,synergy_estimation_full_hierarchy)
  return(output)
}

binary_matrix <- function(m, n) {
  tmp <- sample.int(2L, size = m * n, replace = TRUE) - 1L
  dim(tmp) <- c(m, n)
  tmp
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

