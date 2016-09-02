#Author: mw201608
calc_msbb_array_deg=function(eset,trait=c('CDR','Braak','CERAD','PLQ_Mn','NPrSum','NTrSum'),n.permute){
#eset, an expressionSet object
#n.permute, number of permutations
	library(limma)
	library(Biobase)
	trait=match.arg(trait)
	#
	stages=c('Normal','Low','High')
	AD=stagingAD(eset,trait)
	eset=eset[,!is.na(AD)]
	AD=AD[!is.na(AD)]
	AD=stages[AD+1]
	stages=stages[stages %in% names(table(AD)[table(AD) >= 3])]
	eset=eset[,AD %in% stages]
	AD=AD[AD %in% stages]
	f <- ordered(AD, levels=stages)
	design <- model.matrix(~0+f)
	colnames(design) <- stages
	fit <- lmFit(eset, design)
	if(all(c('Normal','Low','High') %in% stages)){
		iStages=1
		cont <- makeContrasts(LowvsNormal="Low-Normal",HighvsNormal="High-Normal",HighvsLow="High-Low",levels=design)
	}else if(all(c('Normal','Low') %in% stages)){
		iStages=2
		cont <- makeContrasts(LowvsNormal="Low-Normal",levels=design)
	}else if(all(c('Normal','High') %in% stages)){
		iStages=3
		cont <- makeContrasts(HighvsNormal="High-Normal",levels=design)
	}else if(all(c('Low','High') %in% stages)){
		iStages=4
		cont <- makeContrasts(HighvsLow="High-Low",levels=design)
	}
	fit2 <- contrasts.fit(fit, cont)
	fit2 <- eBayes(fit2)
	sig1=sig2=sig3=NULL
	if(iStages==1 || iStages==2) sig1=cbind(topTable(fit2, coef='LowvsNormal', adjust="BH",n= Inf),contrast="Low-Normal")
	if(iStages==1 || iStages==3) sig2=cbind(topTable(fit2, coef='HighvsNormal', adjust="BH",n= Inf),contrast="High-Normal")
	if(iStages==1 || iStages==4) sig3=cbind(topTable(fit2, coef='HighvsLow', adjust="BH",n= Inf),contrast="High-Low")
	if(n.permute>0){
		set.seed(12345)
		tperm=lapply(1:n.permute,function(i){
			f0 <- sample(f,length(f),FALSE)
			design <- model.matrix(~0+f0)
			colnames(design) <- stages
			fit0 <- lmFit(eset, design)
			if(iStages==1){
				cont <- makeContrasts(LowvsNormal="Low-Normal",HighvsNormal="High-Normal",HighvsLow="High-Low",levels=design)
			}else if(iStages==2){
				cont <- makeContrasts(LowvsNormal="Low-Normal",levels=design)
			}else if(iStages==3){
				cont <- makeContrasts(HighvsNormal="High-Normal",levels=design)
			}else if(iStages==4){
				cont <- makeContrasts(HighvsLow="High-Low",levels=design)
			}
			fit0 <- contrasts.fit(fit0, cont)
			fit0 <- eBayes(fit0)
			fit0$t
		})
		tperm=do.call(rbind,tperm)
		if(iStages==1 || iStages==2) sig1$FDR=getFDR(abs(sig1$t),abs(tperm[,'LowvsNormal']))
		if(iStages==1 || iStages==3) sig2$FDR=getFDR(abs(sig2$t),abs(tperm[,'HighvsNormal']))
		if(iStages==1 || iStages==4) sig3$FDR=getFDR(abs(sig3$t),abs(tperm[,'HighvsLow']))
	}
	sig=rbind(sig1,sig2,sig3)
	return(sig)
}
#group samples into low, medium and severe categories
stagingAD=function(eset,trait){
	if(trait=='CDR'){
		#0: 0; 1: 0.5,1,2; 2: 3~5
		AD=ifelse(eset$CDR>0,1,0)
		AD=AD+ifelse(eset$CDR>2,1,0)
	}else if (trait=='Braak'){
		#0: 0~2; 1: 3~4; 2: 5-6 #0 is merged with 1-2 because there is none or no more than two 0 in each region
		AD=ifelse(eset$Braak>2,1,0)
		AD=AD+ifelse(eset$Braak>5,1,0)
	}else if(trait=='CERAD'){
		#0: 1; 1: 3-4; 2: 2 #1=normal, 2=definite AD, 3=Probable AD, 4=possible AD
		AD=ifelse(eset$CERAD>1,1,0)
		AD=AD+ifelse(eset$CERAD==2,1,0)
	}else if(trait=='PLQ_Mn'){
		#0: 0; 1: 1~9; 2: > 9
		AD=ifelse(eset$PLQ_Mn>0,1,0)
		AD=AD+ifelse(eset$PLQ_Mn>9,1,0)
	}else if(trait=='NPrSum'){
		#0: 0; 1: 1~17; 2: > 17
		AD=ifelse(eset$NPrSum>0,1,0)
		AD=AD+ifelse(eset$NPrSum>17,1,0)
	}else if(trait=='NTrSum'){
		#0: 0~2; 1: 2~10; 2: > 10
		AD=ifelse(eset$NTrSum>0,1,0)
		AD=AD+ifelse(eset$NTrSum>10,1,0)
	}
	AD
}
getFDR=function(x,x0){
#x, a vector of non-negative test statistics from original data.
#x0, a vector of non-negative test statistics from permutation
	if(any(x<0) || any(x0<0)) stop('Input must be non-negative\n')
	n1=length(x)
	n2=length(x0)
	x0=sort(x0)
	order.x=order(x)
	rank.x=rep(0L,n1)
	rank.x[order.x]=1:n1
	p=rep(0,n1)
	p[order.x]=n2-findInterval(x[order.x],x0)
	p=(p/n2)/((n1+1-rank.x)/(n1))
	p[p>1]=1
	if(n1==1) return(p)
	for(i in 1:(n1-1)) if(p[order.x[i+1]] > p[order.x[i]]){p[order.x[i+1]]=p[order.x[i]]} #make sure the larger statistics will have a FDR no bigger than that of the smaller statistics
	p
}
