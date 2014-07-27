# Report Betr Graph - To implement
betr_graph <- function(eset,drug,probe_id,limit_display,vsctrl=FALSE) {
	random<-sample(1:length(probe_id),limit_display,replace=F)
	control_probe<-setdiff(rownames(eset),probe_id)
	
	if(vsctrl) {
	#DMSO data
		i<-random[1]
		ctrl_probe<-control_probe[1]
		
		baseline<-mean(exprs(eset[probe_id[i],eset$Drug=="Media"]))
		dmso<-mean(exprs(eset[probe_id[i],eset$Drug=="DMSO" & eset$Time==6])) 
		dmso<-append(dmso,mean(exprs(eset[probe_id[i],eset$Drug=="DMSO" & eset$Time==12])))
		dmso<-append(dmso,mean(exprs(eset[probe_id[i],eset$Drug=="DMSO" & eset$Time==24])))
		probe_value<-mean(exprs(eset[probe_id[i],eset$Drug==drug & eset$Time==6])) 
		probe_value<-append(probe_value,mean(exprs(eset[probe_id[i],eset$Drug==drug & eset$Time==12]))) 
		probe_value<-append(probe_value,mean(exprs(eset[probe_id[i],eset$Drug==drug & eset$Time==24])))	
		dmso<-append(baseline,dmso,1)	
		probe_value<-append(baseline,probe_value,1)
		adjusted_probe_value<-(log2(probe_value)-log2(dmso))/log2(dmso)
		
		# now ctrl
		baseline<-mean(exprs(eset[ctrl_probe,eset$Drug=="Media"]))
		dmso<-mean(exprs(eset[ctrl_probe,eset$Drug=="DMSO" & eset$Time==6])) 
		dmso<-append(dmso,mean(exprs(eset[ctrl_probe,eset$Drug=="DMSO" & eset$Time==12])))
		dmso<-append(dmso,mean(exprs(eset[ctrl_probe,eset$Drug=="DMSO" & eset$Time==24])))
		probe_value<-mean(exprs(eset[ctrl_probe,eset$Drug==drug & eset$Time==6])) 
		probe_value<-append(probe_value,mean(exprs(eset[ctrl_probe,eset$Drug==drug & eset$Time==12]))) 
		probe_value<-append(probe_value,mean(exprs(eset[ctrl_probe,eset$Drug==drug & eset$Time==24])))	
		dmso<-append(baseline,dmso,1)	
		probe_value<-append(baseline,probe_value,1)
		adjusted_probe_value<-rbind(adjusted_probe_value,(log2(probe_value)-log2(dmso))/log2(dmso))
		

		library(ggplot2)
	library(reshape2)
	betr_data<-adjusted_probe_value
	rownames(betr_data)<-c("Detected","Rejected")
	colnames(betr_data)<-c("0","6","12","24")

	
	betr_graph<-melt(betr_data)
	windowsFonts(myriad=windowsFont("Myriad pro"))

	
	graph_out<-ggplot(betr_graph,aes(x=Var2,y=value,colour=Var1,group=Var1)) + geom_line(size=2) + geom_point(size=5) +coord_cartesian(ylim = c(-.06,.11)) 
	graph_out<-graph_out+xlab("Time") + ylab("Log-ratio of expression") + ggtitle('Time-course expression profiles\nof two illustrative genes')
	graph_out <-graph_out+theme_set( theme_bw( base_family= "myriad"))
	graph_out <-graph_out+theme(axis.text.x=element_text(size=18))+theme(axis.text.y=element_text(size=19))+theme(axis.title.y=element_text(size=20))+theme(axis.title.x=element_text(size=20))
	graph_out<-graph_out+theme(plot.title = element_text(size =30))
	graph_out<-graph_out+ theme(legend.position = c(.15, .9))+scale_colour_brewer(name = "Probeset",palette="Set1")+theme(legend.text = element_text(size = 20),legend.title=element_text(size=24))
	print(graph_out)
	
	} else {
	for(z in 1:limit_display) {
		i<-random[z]
		baseline<-mean(exprs(eset[probe_id[i],eset$Drug=="Media"]))
		baseline_sd<-sd(as.vector(exprs(eset[probe_id[i],eset$Drug=="Media"])))
		
		dmso<-mean(exprs(eset[probe_id[i],eset$Drug=="DMSO" & eset$Time==6])) 
		dmso<-append(dmso,mean(exprs(eset[probe_id[i],eset$Drug=="DMSO" & eset$Time==12])))
		dmso<-append(dmso,mean(exprs(eset[probe_id[i],eset$Drug=="DMSO" & eset$Time==24])))
		dmso_sd<-sd(as.vector(exprs(eset[probe_id[i],eset$Drug=="DMSO" & eset$Time==6])))
		dmso_sd<-append(dmso_sd,sd(as.vector(exprs(eset[probe_id[i],eset$Drug=="DMSO" & eset$Time==12]))))
		dmso_sd<-append(dmso_sd,sd(as.vector(exprs(eset[probe_id[i],eset$Drug=="DMSO" & eset$Time==24]))))
		
		probe_value<-mean(exprs(eset[probe_id[i],eset$Drug==drug & eset$Time==6])) 
		probe_value<-append(probe_value,mean(exprs(eset[probe_id[i],eset$Drug==drug & eset$Time==12]))) 
		probe_value<-append(probe_value,mean(exprs(eset[probe_id[i],eset$Drug==drug & eset$Time==24])))
		probe_value_sd<-sd(as.vector(exprs(eset[probe_id[i],eset$Drug==drug & eset$Time==6]))) 
		probe_value_sd<-append(probe_value_sd,sd(as.vector(exprs(eset[probe_id[i],eset$Drug==drug & eset$Time==12])))) 
		probe_value_sd<-append(probe_value_sd,sd(as.vector(exprs(eset[probe_id[i],eset$Drug==drug & eset$Time==24]))))
		
		dmso<-append(baseline,dmso,1)
		dmso_sd<-append(baseline_sd,dmso_sd,1)
		
		probe_value<-append(baseline,probe_value,1)
		probe_value_sd<-append(baseline_sd,probe_value_sd,1)
		
		adjusted_probe_value<-(log2(probe_value)-log2(dmso))/log2(dmso)
		#Wrong, dont know how to calculate this SD
		adjusted_probe_value_sd<-(log2(probe_value_sd)-log2(dmso_sd))/log2(dmso_sd)
		
		if(z>1) {
			betr_data<-rbind(betr_data,adjusted_probe_value)
			betr_data_sd<-rbind(betr_data_sd,adjusted_probe_value_sd)
		} else {
			betr_data<-rbind(adjusted_probe_value)
			betr_data_sd<-rbind(adjusted_probe_value_sd)
		}
		rownames(betr_data)[z]<-probe_id[i]
		rownames(betr_data_sd)[z]<-probe_id[i]
	}

	library(ggplot2)
	library(reshape2)
	
	colnames(betr_data)<-c("0","6","12","24")
	colnames(betr_data_sd)<-c("0","6","12","24")
	### Write PDF with graph in the folder, no output ###
	
	betr_graph<-melt(betr_data)
	betr_graph_sd<-melt(betr_data_sd)
	colnames(betr_graph_sd)<-colnames(betr_graph)
	
	graph_out<-ggplot(betr_graph,aes(x=Var2,y=value,colour=Var1,group=Var1)) + geom_line(size=1) + geom_point(size=3)
	graph_out<-graph_out+xlab("Hours") + ylab("Gene Expression") + ggtitle(paste(limit_display,' probes differently expressed (',length(probe_id),' total)\n',drug, ' Drug - BETR Filtering',sep="")) + theme(plot.title=element_text(size=14))
	print(graph_out)
	}
}

# Graph fig 1B
fit_result_boxplot <- function(drug_list,fit_result,ratio=FALSE) {

	DSS_A<-vector()
	DSS_B<-vector()
	DSS_C<-vector()
	DES_A<-vector()
	DES_B<-vector()
	DES_C<-vector()
	
	for(i in 1:length(fit_result)) {
		current<-as.data.frame(fit_result[[i]])
		colnames(current)<-c('A','B','C')
		DSS_A<-append(DSS_A,length(subset(current,A==1)[,1]))
		DSS_B<-append(DSS_B,length(subset(current,B==1)[,2]))
		DSS_C<-append(DSS_C,length(subset(current,C==1)[,3]))
		DES_A<-append(DES_A,length(subset(current,A==-1)[,1]))
		DES_B<-append(DES_B,length(subset(current,B==-1)[,2]))
		DES_C<-append(DES_C,length(subset(current,C==-1)[,3]))
	}
	
	if(ratio) {
		data<-cbind(DSS_A/DES_A,DSS_B/DES_B,DSS_C/DES_C)
		name<-c('6hr','12hr','24hr')
		boxplot(data,names=name,ylab = "Number of probesets in DSS/DES")
		title("Distribution of numbers of probesets in normalized DSS/DES")
	} else {
	data<-cbind(DES_A,DSS_A,DES_B,DSS_B,DES_C,DSS_C)
	name<-c('DES/6hr','DSS/6hr','DES/12hr','DSS/12hr','DES/24hr','DSS/24hr')
	boxplot(data,names=name,ylab = "Number of probesets")
	title("Distribution of numbers of probesets per DES and DSS")
	}
}

# Interaction matrix display, fig2A
interaction_graph <- function(data,limit_scale=c(-2.1,2.1),sign=FALSE) {
	if(sign) {
		library(corrplot)
		limit_scale=c(-1,1)
		data<-sign(data)
		output<-data
		output[!is.na(output)]<-0
		for(i in 1:length(data[,1])) {
			for(j in 1:length(data[1,])) {
				if (data[i,j]==data[j,i]) {output[i,j]=1;output[j,i]=1}
				else {output[i,j]=-1;output[j,i]=-1}
			}
		}
		col1 <- colorRampPalette(c("blue","white","red"))
		corrplot(output, method = "color",type="upper",col=col1(3))
	}
	
	else {
	
	library(ggplot2)
	library(reshape2)
	output_order<-rownames(data)[order(rownames(data),decreasing=TRUE)]
	data<-as.matrix(data)
	data<-melt(data)
	data<-within(data,Var1<-factor(Var1,levels=output_order))
	graph_out <- ggplot(data, aes(Var2, Var1)) + geom_tile(aes(fill = value),colour = "white") + scale_fill_gradientn(limits=limit_scale,colours=c("blue","white","red"))
	base_size <- 18
	graph_out <- graph_out + theme_grey(base_size = base_size) + labs(x = "Drug-Specific Signatures",y = "Drug-treated samples") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
	graph_out <- graph_out + theme(axis.ticks = element_blank(), axis.text.x = element_text(size = base_size, angle = 330, hjust = 0, colour = "black"), axis.text.y = element_text(size =base_size, colour = "black"))
	graph_out <- graph_out + ggtitle('Drug-Drug interaction matrix') + theme(axis.text.x=element_text(size=base_size*0.9),axis.text.y=element_text(size=base_size*0.9),title=element_text(size=base_size*2))
	print(graph_out)
	}
}



# to work, must run population_targeting with report=ON
timepoint_graph <- function() {

load(interaction_report_DES_t1,file='DES interaction report at t1.rda')
load(interaction_report_DES_t2,file='DES interaction report at t2.rda')
load(interaction_report_DES_t3,file='DES interaction report at t3.rda')
load(interaction_report_DSS_t1,file='DSS interaction report at t1.rda')
load(interaction_report_DSS_t2,file='DSS interaction report at t2.rda')
load(interaction_report_DSS_t3,file='DSS interaction report at t3.rda')

#Normalize by rows (eg. timepoint) for any drug
for( z in 1:length(druglist)) {
	interaction_report_DES_t1[[druglist[z]]]<-prop.table(as.matrix(interaction_report_DES_t1[[druglist[z]]]),1)
	interaction_report_DES_t1[[druglist[z]]][is.na(interaction_report_DES_t1[[druglist[z]]])]<-0
	interaction_report_DSS_t1[[druglist[z]]]<-prop.table(as.matrix(interaction_report_DSS_t1[[druglist[z]]]),1)
	interaction_report_DSS_t1[[druglist[z]]][is.na(interaction_report_DSS_t1[[druglist[z]]])]<-0
			
	interaction_report_DES_t2[[druglist[z]]]<-prop.table(as.matrix(interaction_report_DES_t2[[druglist[z]]]),1)
	interaction_report_DES_t2[[druglist[z]]][is.na(interaction_report_DES_t2[[druglist[z]]])]<-0
	interaction_report_DSS_t2[[druglist[z]]]<-prop.table(as.matrix(interaction_report_DSS_t2[[druglist[z]]]),1)
	interaction_report_DSS_t2[[druglist[z]]][is.na(interaction_report_DSS_t2[[druglist[z]]])]<-0
			
	interaction_report_DES_t3[[druglist[z]]]<-prop.table(as.matrix(interaction_report_DES_t3[[druglist[z]]]),1)
	interaction_report_DES_t3[[druglist[z]]][is.na(interaction_report_DES_t3[[druglist[z]]])]<-0
	interaction_report_DSS_t3[[druglist[z]]]<-prop.table(as.matrix(interaction_report_DSS_t3[[druglist[z]]]),1)
	interaction_report_DSS_t3[[druglist[z]]][is.na(interaction_report_DSS_t3[[druglist[z]]])]<-0	
}			
	#Rows are timepoints used to generate DES or DSS
	#Columns are timepoints of drug-treated samples with matches
		
	DES_time_report<-matrix(ncol=3,nrow=3)
	rownames(DES_time_report)<- c('6hr signature','12hr signature','24hr signature')
	colnames(DES_time_report)<- c('6hr treated samples','12hr treated samples','24hr treated_samples')
	DSS_time_report<-DES_time_report
		
	for(i in 1:3) {
			any_drug_DSS_t1<-vector()
			any_drug_DES_t1<-vector()
			any_drug_DSS_t2<-vector()
			any_drug_DES_t2<-vector()
			any_drug_DSS_t3<-vector()
			any_drug_DES_t3<-vector()
			
			any_drug_DES_t1<-append(any_drug_DES_t1,interaction_report_DES_t1[[druglist[z]]][,i])
			any_drug_DSS_t1<-append(any_drug_DSS_t1,interaction_report_DSS_t1[[druglist[z]]][,i])
			any_drug_DES_t2<-append(any_drug_DES_t2,interaction_report_DES_t2[[druglist[z]]][,i])
			any_drug_DSS_t2<-append(any_drug_DSS_t2,interaction_report_DSS_t2[[druglist[z]]][,i])
			any_drug_DES_t3<-append(any_drug_DES_t3,interaction_report_DES_t3[[druglist[z]]][,i])
			any_drug_DSS_t3<-append(any_drug_DSS_t3,interaction_report_DSS_t3[[druglist[z]]][,i])
			
			DES_time_report[1,i] <-median(any_drug_DES_t1)
			DES_time_report[2,i] <-median(any_drug_DES_t2)
			DES_time_report[3,i] <-median(any_drug_DES_t3)
			DSS_time_report[1,i] <-median(any_drug_DSS_t1)
			DSS_time_report[2,i] <-median(any_drug_DSS_t2)
			DSS_time_report[3,i] <-median(any_drug_DSS_t3)
	}
		
	# Plot 3x3 matrix for DES and DSS for any drugs
	library(ggplot2)
	library(reshape2)
	base_size <- 20
	limit_scale<-c(0,1)
	
	DES_output_order<-rev(rownames(DES_time_report))
	DES_time_report<-melt(DES_time_report)
	DES_time_report<-within(DES_time_report,Var1<-factor(Var1,levels=DES_output_order))
	DES_time_report$Var2<-relevel(DES_time_report$Var2, '6hr treated samples')
	
	DESgraph <- ggplot(DES_time_report, aes(Var2, Var1)) + geom_tile(aes(fill = value),colour = "black") + scale_fill_gradientn(limits=limit_scale,colours=c("white","red"),breaks=limit_scale)
	DESgraph <- DESgraph + theme_grey(base_size = base_size) + labs(x = "Samples",y = "Signatures") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
	DESgraph <- DESgraph + theme(axis.ticks = element_blank(), axis.text.x = element_text(size = base_size, angle = 330, hjust = 0, colour = "black"), axis.text.y = element_text(size =base_size, colour = "black"))
	DESgraph<- DESgraph + ggtitle('DES - Definition vs Matching') + theme(axis.text.x=element_text(size=base_size*0.9),axis.text.y=element_text(size=base_size*0.9),title=element_text(size=base_size*2))

	DSS_output_order<-rev(rownames(DSS_time_report))
	DSS_time_report<-melt(DSS_time_report)
	DSS_time_report<-within(DSS_time_report,Var1<-factor(Var1,levels=DSS_output_order))
	DSS_time_report$Var2<-relevel(DSS_time_report$Var2, '6hr treated samples')
	
	DSSgraph <- ggplot(DSS_time_report, aes(Var2, Var1)) + geom_tile(aes(fill = value),colour = "black") + scale_fill_gradientn(limits=limit_scale,colours=c("white","red"),breaks=limit_scale)
	DSSgraph <- DSSgraph + theme_grey(base_size = base_size) + labs(x = "Samples",y = "Signatures") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
	DSSgraph <- DSSgraph + theme(axis.ticks = element_blank(), axis.text.x = element_text(size = base_size, angle = 330, hjust = 0, colour = "black"), axis.text.y = element_text(size =base_size, colour = "black"))
	DSSgraph<- DSSgraph + ggtitle('DSS - Definition vs Matching') + theme(axis.text.x=element_text(size=base_size*0.9),axis.text.y=element_text(size=base_size*0.9),title=element_text(size=base_size*2))

	multiplot(DESgraph,DSSgraph,cols=2)		
}

#GGplot side by side
multiplot <- function(..., plotlist=NULL, cols) {
    require(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # Make the panel
    plotCols = cols                          # Number of columns of plots
    plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols

    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x, y)
        viewport(layout.pos.row = x, layout.pos.col = y)

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
        curRow = ceiling(i/plotCols)
        curCol = (i-1) %% plotCols + 1
        print(plots[[i]], vp = vplayout(curRow, curCol ))
    }

}

graph_model_fit<-function(eset,drug,probes,timepoint) {
	#Probes, 1st up, 2nd down
	library(Biobase)
	library(simpleaffy)
	library(drc)
	
	subset<-get.array.subset.exprset(eset, "Time", timepoint)
	subset<-get.array.subset.exprset(subset, "Drug" , c(drug,'DMSO'))
	
	IC<-subset$Concentration
	Expression_up<-as.vector(exprs(subset[probes[1],]))
	Expression_down<-as.vector(exprs(subset[probes[2],]))
	
	for(z in 1:length(IC)) {
			if(IC[z]=='NONE') {IC[z]=0.00}
			if(IC[z]=='IC20') {IC[z]=0.2}	
			if(IC[z]=='IC20_LOW') {IC[z]=0.02}

	}
	class(IC)<-"numeric"
	class(Expression_up)<-"numeric" 
	class(Expression_down)<-"numeric" 

	mm_fit_up<-itfit<-try(drm(Expression_up~IC,fct=MM.3(),control=drmc(noMessage=TRUE)),silent=TRUE)
	mm_fit_down<-itfit<-try(drm(Expression_down~IC,fct=MM.3(),control=drmc(noMessage=TRUE)),silent=TRUE)
	coef_up<-round(as.vector(coef(mm_fit_up)),3)	
	coef_down<-round(as.vector(coef(mm_fit_down)),3)	
	
	text="Shifted Michaelis-Menten Model"
	text_up=paste(text,"\nParameters: C=",coef_up[1]," D=",coef_up[2]," E=",coef_up[3],sep="")	
	text_down=paste(text,"\nParameters: C=",coef_down[1]," D=",coef_down[2]," E=",coef_down[3],sep="")
	
	plot(mm_fit_down,lwd=2,type="all",main="Probeset classified in DES",xlab="Dosage",ylab="Probeset Intensity",legend=TRUE,legendText=text_down)
	plot(mm_fit_up,lwd=2,type="all",main="Probeset classified in DSS",xlab="Dosage",ylab="Probeset Intensity",legend=TRUE,legendText=text_up)
}

landscape_definition<-function(eset,drug,druglist,fit_result) {
  #output is a matrix 16x3 from 0 (alive) to 1 (death)
  output<-matrix(rep(-2,14*3),nrow=14,ncol=3)
  drug<-as.vector(unlist(drug))

  for(t in 1:3) {
    signature<-capture_signature(drug,fit_result,t,permanent_effect=FALSE,distinct_signature=FALSE) 
	targeting<-factor_dataframe(NTPez(eset,drug,druglist,signature,t,time_matching=FALSE,drug_specific=FALSE,within.sig=FALSE))
	for(i in 1:length(druglist)){
		if(!grepl(drug,druglist[i])) {
			#Keep relevant matches
			case<-subset(targeting,grepl(druglist[i],drug))
			case<-subset(case,FDR<0.05)
			
			if(length(case[,1])>0) {
				# Get most repeated Label first
				label<-names(which.max(table(case$label)))
				distance<-median(case$distance)
			} else {
				#Neither DES or DSS are significant
				distance<-0
				label<-'DES'
			}
			
			if(grepl('DSS',label)) {distance<--distance}# Anticorrelated
		
		} else {
			case<-subset(targeting,grepl(druglist[i],drug))
			case<-subset(case,FDR<0.05)
			distance<-median(case$distance)
		}    
      output[i,t]<-distance 
    }
  }
  #Scale between 0 and 1, per column (=timepoint)
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  output<-apply(output,2,range01)
  #Expand matrix for nice plot
  val_min<-min(output)
  val_max<-max(output)
  nice_output<-matrix(rep(val_min,30*8),nrow=30,ncol=8)
  nice_output[c(1,30),]<-val_max
	
	for(i in seq(2, 28, 2)) {
		nice_output[i:(i+1),(3:4)]<-output[i/2,1]
		nice_output[i:(i+1),(5:6)]<-output[i/2,2]
		nice_output[i:(i+1),(7:8)]<-output[i/2,3]
	}
	return (nice_output)
}

landscape_plot <- function(z,title,drug,druglist,col_line=c(2,6),colors=colorRampPalette( c("white","green") )){
	x <- seq(0,29,length.out=nrow(z))
	y <- seq(0,7,length.out=ncol(z))
	## getting the value of the midpoint
    zz <- (z[-1,-1] + z[-1,-ncol(z)] + z[-nrow(z),-1] + z[-nrow(z),-ncol(z)])/4
    ## calculating the breaks
    breaks <- hist(zz, plot=FALSE)$breaks
    ## cutting up zz
    cols <- colors(length(breaks)-1)
    zzz <- cut(zz, breaks=breaks, labels=cols)
    ## plotting
    persp(x,y,z,r=4,main=title,xlab='Subpopulations',ticktype='simple',ylab='Timepoint',zlab='Death',col=as.character(zzz),theta=140,phi=30,expand=0.2)->res
    
	## Draw lines for baseline drug
	index<-match(drug,druglist)*2
	
	for(i in 1:length(index)) {
		lines (trans3d(x=index[i], y, z = max(z[index[i],]), pmat = res), col = col_line[i],lw=4,lty=1)
		lines (trans3d(x=index[i]-1, y, z = max(z[index[i],]), pmat = res), col = col_line[i],lw=4,lty=1)
	}

	
    ## return breaks and colors for the legend
    list(breaks=breaks, colors=cols) 
}

landscape_niveling<- function(landscape,drug1,drug2,druglist,hierarchy_targeting=FALSE,hierarchy=NULL) {
	# Find coordinates drug1/drug2
	for(i in 1:length(druglist)) {
		if(grepl(drug1,druglist[i])) {x<-i}
		if(grepl(drug2,druglist[i])) {y<-i}
	}
	if(hierarchy_targeting) {
		new_landscape1<-hierarchy_landscape(drug1,landscape[[x]],hierarchy)
		new_landscape2<-hierarchy_landscape(drug2,landscape[[y]],hierarchy)
		return (pmax(new_landscape1,new_landscape2))
	} else {	
	return (pmax(landscape[[x]],landscape[[y]]))
	}
}

hierarchy_landscape<-function(drug,landscape_drug,hierarchy) {
# Build hierarchy data for drug
hierarchy_drug<-hierarchy[drug,]
next_round<-hierarchy[drug,]!=0
repeat {
	if(sum(next_round)==0) { break }
	if(sum(next_round)==1) {
		hierarchy_drug<-hierarchy_drug+hierarchy[next_round,]
		next_round<-hierarchy[next_round,]!=0
	}
	if(sum(next_round)>1) {
		hierarchy_drug<-hierarchy_drug+colSums(hierarchy[next_round,])
		next_round<-hierarchy[next_round,]!=0
	}
}
	
# Set anything to 1 or zero
hierarchy_drug[hierarchy_drug==0]<-0
hierarchy_drug[hierarchy_drug!=0]<-1
hierarchy_drug[drug]<-2

#What was done for landscape definition
temp<-c(0)
for(i in 1:length(hierarchy_drug)) {
	temp<-append(temp,hierarchy_drug[i])
	temp<-append(temp,hierarchy_drug[i])
}
temp<-append(temp,0)
names(temp)<-NULL
hierarchy_drug<-temp
landscape_drug[hierarchy_drug==1,(3:8)]<-max(landscape_drug[hierarchy_drug==2,(3:4)])
return (landscape_drug)
}

landscape_definition_alldrug<-function(eset,druglist,fit_result) {
	output<-list()
	for(i in 1:length(druglist)){
		land<-landscape_definition(eset,druglist[i],druglist,fit_result)
		output<-cbind(output,list(land))
	}
	return(output)	
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

network_graph<-function(rank_prediction,druglist,top=10,title='Network') {
  library(igraph)
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
  
  graph1 <- graph.adjacency(as.matrix(predicted_data),weighted=T,mode="undirected",diag=FALSE)
  graph2 <- graph.adjacency(as.matrix(data),weighted=T,mode="undirected",diag=FALSE)
  layout_final <- layout.fruchterman.reingold(graph2, niter=10000, area=30*vcount(graph2)^2)
  par(mfrow=c(1,2)) 
    plot.igraph(graph1,
                rescale=TRUE,
                vertex.size=45,
                vertex.label.family="serif",
                edge.label.family="Palatino",
                vertex.label=V(graph1)$name,
                vertex.color='white',
                vertex.frame.color='black',
                vertex.label.color='black',	
                vertex.label.font=0.5,
                vertex.shape="rectangle",
                layout=layout_final,
                main='Predicted combination',
                vertex.label.color="black",
                edge.color="black",
                edge.width=3,
                edge.arrow.size=.2,
                edge.curved=TRUE
    )
    
    plot.igraph(graph2,
                rescale=TRUE,
                vertex.size=45,
                vertex.label.family="serif",
                edge.label.family="Palatino",
                vertex.label=V(graph2)$name,
                vertex.color='white',
                vertex.frame.color='black',
                vertex.label.color='black',	
                vertex.label.font=0.5,
                vertex.shape="rectangle",
                layout=layout_final,
                main='Experimental results',
                vertex.label.color="black",
                edge.color="black",
                edge.width=3,
                edge.arrow.size=.2,
                edge.curved=TRUE
    )
}

hclust_eset<-function(eset,druglist,topsd=FALSE,threshold=500) {
 library(simpleaffy)
 library(ape)
 eset<-get.array.subset.exprset(eset,"Concentration","IC20")
 eset<-get.array.subset.exprset(eset,"Drug", druglist)
 timelist<-c(6,12,24)
 feat<-length(featureNames(eset))
 
 par(mfrow=c(1,3))
 for(t in 1:length(timelist)) {
 output<-matrix(nrow=feat,ncol=length(druglist))
 rownames(output)<-featureNames(eset)
 colnames(output)<-druglist
	current<-get.array.subset.exprset(eset,"Time", timelist[t])
	for(i in 1:length(druglist)) {
		temp<-exprs(get.array.subset.exprset(current,"Drug", druglist[i]))
		output[,i]<-apply(temp, 1, median)
	}
	
	if(topsd){
		filtersd<-apply(output,1,sd)
		ranksd<-rank(-filtersd)
		for(x in 1:length(ranksd)) {
			if(ranksd[x]<=threshold) {ranksd[x]=TRUE} else {ranksd[x]=FALSE}
		}
		output<-output[indexsd,]
	}
	d<-dist(t(output), method = "euclidean")
	fit<-hclust(d,method="ward")
	plot(as.phylo(fit),label.offset=1,type='cladogram',cex=1.2)
	title(paste("Drug-treated samples at ",timelist[t],"hr",sep=""),cex.main=2)
 }
}