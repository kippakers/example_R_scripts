data_out<-NULL
#this loops through every gene/tissue combination to examine correlations
for (f in 1:nrow(gene_trait)) {
	#for each gene
	gene<-as.character(gene_trait[f,"ensembl2"])
	symbol<-as.character(gene_trait[f,"gene"])
	trait<-as.character(gene_trait[f,"trait"])
	#for each tissue
	for (tissue in list(AOR, Blood, FC, LIV, MAM, MP, SF, SKLM, VAF)) {
		#if the tissue is expressed
		if(gene %in% rownames(tissue)){ 
			tissuename<-gsub("_.*", "", colnames(tissue)[1])
			###top % vs bottom %
			bottom10=names(sort(tissue[gene,])[1:(ncol(tissue)*0.1)])
			bottom10=gsub(".*_", "", bottom10)
			top10=names(sort(tissue[gene,], decreasing=T)[1:(ncol(tissue)*0.1)])
			top10=gsub(".*_", "", top10)
			results<-t.test(clin[which(clin$starnet.ID %in% bottom10),trait], clin[which(clin$starnet.ID %in% top10),trait])
			# effect of going from bottom10% to top10%
			fc<-log2(results$estimate[2]/results$estimate[1])
			### Linear Modeling:
			#merge clinical and gene expression data
			tempExpr<-as.data.frame(tissue[gene,])
			rownames(tempExpr)<-gsub(".*_", "", rownames(tempExpr))
			colnames(tempExpr)<-gene
			if (trait =="BMI") {
				tempdf<-merge(clin[,c("starnet.ID", trait, "Age", "Sex")], tempExpr, by.x="starnet.ID", by.y=0)
				tempdf2<-merge(tempdf, covariates[which(covariates$tissue==tissuename),], by.x="starnet.ID", by.y="patient")
				linearmodsimple<-summary(lm(tempdf2[,2]~tempdf2[,5]))
				if (length(unique(tempdf2$lab)) > 1) {
					linearmodcomplex<-summary(lm(tempdf2[,2]~tempdf2[,5]+tempdf2[,3]+tempdf2[,4]+tempdf2[,8]))
				}
				else {
					linearmodcomplex<-summary(lm(tempdf2[,2]~tempdf2[,5]+tempdf2[,3]+tempdf2[,4]))
				}
			}
			else {
				tempdf<-merge(clin[,c("starnet.ID", trait, "Age", "Sex", "BMI")], tempExpr, by.x="starnet.ID", by.y=0)
				tempdf2<-merge(tempdf, covariates[which(covariates$tissue==tissuename),], by.x="starnet.ID", by.y="patient")
				linearmodsimple<-summary(lm(tempdf2[,2]~tempdf2[,6]))
				if (length(unique(tempdf2$lab)) > 1) {
					linearmodcomplex<-summary(lm(tempdf2[,2]~tempdf2[,6]+tempdf2[,3]+tempdf2[,4]+tempdf2[,5]+tempdf2[,9]))
				}
				else {
					linearmodcomplex<-summary(lm(tempdf2[,2]~tempdf2[,6]+tempdf2[,3]+tempdf2[,4]+tempdf2[,5]))
				}
			}
			linearmodsum<-cbind(lmSimpEst=linearmodsimple$coefficients[2,1],lmSimpP=linearmodsimple$coefficients[2,4], lmFullEst=linearmodcomplex$coefficients[2,1], lmFullP=linearmodcomplex$coefficients[2,4], intercept=linearmodsimple$coefficients[1,1], adjR2=linearmodsimple$adj.r.squared)		
			#add it up
			data_out<-rbind(data_out,data.frame(cbind(gene, symbol, trait, Tissue=tissuename, subsetfc=fc, subsetpval=results$p.value, linearmodsum), row.names=f))
		}
	}
}

#for every gene/tissue : expression value, + low, high, stdv. 
gene_data<-NULL
for (gene in unique(data_out$gene)) {
	for (tissue in list(AOR, Blood, FC, LIV, MAM, MP, SF, SKLM, VAF)) {
		tissuename<-gsub("_.*", "", colnames(tissue)[1])
		if (gene %in% rownames(tissue)) { 
			mysum<-summary(tissue[gene,])
			mysd<-sd(tissue[gene,])
			gene_data<-rbind(gene_data, c(gene=gene, Tissue=tissuename, mysum, stdv=mysd))
		}
		else { 
			mysum<-rep("NA", 7)
			gene_data<-rbind(gene_data, c(gene=gene, Tissue=tissuename, mysum))
		}
	}
}
gene_data<-as.data.frame(gene_data)

##generate images (linear):
for (f in 1:nrow(data_out)) {
	mygene=as.character(data_out[f,"gene"])
	mytissue=data_out[f,"Tissue"]
	mytrait=as.character(data_out[f,"trait"])
	pdfname=paste("plots/starnet", gene_converter[mygene,], mytissue, mytrait, "simple", "pdf", sep=".")
	pdf(pdfname)
	print(plotGenelm(mygene, get(as.character(mytissue)), mytrait, data_out[f,"intercept"], data_out[f,"lmSimpEst"], data_out[f, "adjR2"]))
	dev.off()
}

#generate images (boxplot):
for (f in 1:nrow(data_out)) {
	mygene=as.character(data_out[f,"gene"])
	mytissue=data_out[f,"Tissue"]
	mytrait=as.character(data_out[f,"trait"])
	pdfname=paste("plots/starnet", gene_converter[mygene,], mytissue, mytrait, "10pctboxplot", "pdf", sep=".")
	pdf(pdfname)
	print(plotGeneBox(mygene, get(as.character(mytissue)), mytrait))
	dev.off()
}


library(RColorBrewer)
plotGene<-function(gene, tissue, trait) {
	localenv <- environment()
	tempExpr<-t(as.data.frame(tissue[gene,]))
	tissuename<-gsub("_.*", "", colnames(tissue)[1])
	rownames(tempExpr)<-gsub(".*_", "", rownames(tempExpr))
	colnames(tempExpr)<-gene
	tempdf<-unique(merge(clin[,c("starnet.ID", trait)], tempExpr, by.x="starnet.ID", by.y=0))
	breaks<-unique(quantile(tempdf[,gene],probs=seq(0,1,by=0.1)))
	tempdf$deciles = cut(tempdf[,gene],breaks=breaks,include.lowest=TRUE)
	colours=brewer.pal(name="RdYlGn", n=nlevels(tempdf$deciles))
	trait<-traitConverter[trait,2]
	qplot(x=tempdf[, 3], y = tempdf[, 2], colour = tempdf[,4], environment = localenv) + geom_point() + labs(y=trait, x=gene, title=tissuename)+ scale_color_manual(values=colours, guide=F)
}
plotGenelm<-function(gene, tissue, trait, myintercept, myestimate, rsquared) {
	tempExpr<-t(as.data.frame(tissue[gene,]))
	tissuename<-gsub("_.*", "", colnames(tissue)[1])
	rownames(tempExpr)<-gsub(".*_", "", rownames(tempExpr))
	colnames(tempExpr)<-gene
	tempdf<-merge(clin[,c("starnet.ID", trait)], tempExpr, by.x="starnet.ID", by.y=0)
	mytext<-paste("Adj.R^2 =", rsquared, sep=" ")
	qplot(x=tempdf[, 3], y = tempdf[, 2]) + geom_point() +geom_abline(intercept=as.numeric(as.character(myintercept)), slope=as.numeric(as.character(myestimate)))+labs(y=trait, x=gene, title=tissuename) + annotate("text", x=max(tempdf[,3], na.rm=T), y=min(tempdf[,2], na.rm=T), label=mytext, size=3 )
}
plotGeneBox<-function(gene, tissue, trait) {
	localenv <- environment()
	tempExpr<-t(as.data.frame(tissue[gene,]))
	tissuename<-gsub("_.*", "", colnames(tissue)[1])
	rownames(tempExpr)<-gsub(".*_", "", rownames(tempExpr))
	colnames(tempExpr)<-gene
	tempdf<-unique(merge(clin[,c("starnet.ID", trait)], tempExpr, by.x="starnet.ID", by.y=0))
	breaks<-unique(quantile(tempdf[,gene],probs=seq(0,1,by=0.1)))
	tempdf$deciles = cut(tempdf[,gene],breaks=breaks,include.lowest=TRUE)
	trait<-traitConverter[trait,2]
	qplot(x=tempdf[which(tempdf$deciles==levels(tempdf$deciles)[1] | tempdf$deciles==levels(tempdf$deciles)[10]),"deciles"], y=tempdf[which(tempdf$deciles==levels(tempdf$deciles)[1] | tempdf$deciles==levels(tempdf$deciles)[10]),2], environment = localenv)+geom_boxplot()+labs(y=trait, x=gene, title=tissuename)

} 
plotGeneResids<-function(gene, tissue, trait) {
	localenv <- environment()
	tempExpr<-t(as.data.frame(tissue[gene,]))
	tissuename<-gsub("_.*", "", colnames(tissue)[1])
	rownames(tempExpr)<-gsub(".*_", "", rownames(tempExpr))
	colnames(tempExpr)<-gene
	tempdf<-unique(merge(clin[,c("starnet.ID", trait, "Age", "Sex", "BMI")], tempExpr, by.x="starnet.ID", by.y=0))	
	tempdf2<-merge(tempdf, covariates[which(covariates$tissue==tissuename),], by.x="starnet.ID", by.y="patient")
	tempdf2<-na.omit(tempdf2)
	if (length(unique(tempdf2$lab)) > 1) {
		linearmodcomplex<-(lm(tempdf2[,2]~tempdf2[,3]+tempdf2[,4]+tempdf2[,5]+tempdf2[,9]))
	}
	else {
		linearmodcomplex<-(lm(tempdf2[,2]~tempdf2[,3]+tempdf2[,4]+tempdf2[,5]))
	}
	tempdf2$resids<-residuals(linearmodcomplex)
	breaks<-unique(quantile(tempdf2[,gene],probs=seq(0,1,by=0.1)))
	tempdf2$deciles = cut(tempdf2[,gene],breaks=breaks,include.lowest=TRUE)
	colours=brewer.pal(name="RdYlGn", n=nlevels(tempdf2$deciles))
	trait<-traitConverter[trait,2]
	qplot(x=tempdf2[, gene], y = tempdf2[, "resids"], colour = tempdf2[,"deciles"], environment = localenv) + geom_point() + labs(y=paste("Age, Sex, BMI, and Site Ajusted", trait, "Residuals", sep=" "), x=gene, title=tissuename)+ scale_color_manual(values=colours, guide=F)
}


makepics<-function(gene, tissue, trait) {
	tissuename<-gsub("_.*", "", colnames(tissue)[1])
	#pdf(paste(gene, trait, tissuename, "basic", "pdf", sep="."))
	#pdf(paste(gene, trait, tissuename, "pdf", sep="."))
	p1<-plotGene(gene, tissue, trait)
	ggsave(p1, filename=paste(gene, trait, tissuename, "basic", "pdf", sep="."))
	#pdf(paste(gene, trait, tissuename, "boxplots", "pdf", sep="."))
	p2<-plotGeneBox(gene, tissue, trait)
	ggsave(p2, filename=paste(gene, trait, tissuename, "boxplots", "pdf", sep="."))
	#dev.off()
	#pdf(paste(gene, trait, tissuename, "adjusted", "pdf", sep="."))
	p3<-plotGeneResids(gene, tissue, trait)
	ggsave(p3, filename=paste(gene, trait, tissuename, "adjusted", "pdf", sep="."))
	#dev.off()
}

