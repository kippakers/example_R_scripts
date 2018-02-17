load("blca.round2.forBojan.RData")
wideScreen()
library(plyr)
library(survival)
library(GGally)
library(ggplot2)
library(rms)
library(MASS)

## Divide Samples--Pick one of the strategies
#   strategy 1: supplement batch 1 samples
testSamples<-sampleSheet$ID[which(sampleSheet$Set=="test")] 
trainSamples<-sampleSheet$ID[which(sampleSheet$Set=="train")]
#	strategy 2: mix/randomize batch1/2
set.seed(1234) ; testSamples<-sample(sampleSheet$ID, 194) 
trainSamples<-sampleSheet$ID[! sampleSheet$ID %in% testSamples]

#shuffling samples
sampleSheet$Set2 <- sampleSheet$ID %in% testSamples

#reduce measures:
PatientsWithEpitope<-apply(EpMeasuresClass1all, 2, function(x) sum(x>0.5))
Percent<-round(nrow(EpMeasuresClass1all)*0.5)
EpMeasuresClass1<-EpMeasuresClass1all[,which(PatientsWithEpitope>Percent)]

PatientsWithEpitope<-apply(EpMeasuresClass2all, 2, function(x) sum(x>0.5))
Percent<-round(nrow(EpMeasuresClass2all)*0.5)
EpMeasuresClass2<-EpMeasuresClass2all[,which(PatientsWithEpitope>Percent)]

#separate training data
trainData<-merge(EpMeasuresClass1[which(rownames(EpMeasuresClass1) %in% trainSamples),] , blcaPatient, by.x=0, by.y="ID")
colnames(trainData)[1]<-"ID"
trainDatac2<-merge(EpMeasuresClass2[which(rownames(EpMeasuresClass2) %in% trainSamples),] , blcaPatient, by.x=0, by.y="ID")
colnames(trainDatac2)[1]<-"ID"

# Class 1: perform first pass ep measure association testing:
coxResultslog<-data.frame()
for (f in (colnames(EpMeasuresClass1))) {
	myformula<-as.formula(paste('Surv(Overall.Survival..Months., OS.Status)~log(', f, ")+ SimpleStageNumeric + Diagnosis.Age + Pass2quant", collapse=""))
 	mycox<-cph(myformula, data=trainData, x=T, y=T)
 	coxResultslog<-rbind(coxResultslog, c(mycox$loglik[2], mycox$stats)) 
}
colnames(coxResultslog)[1]<-c("loglik")
colnames(coxResultslog)[2:12]<-names(mycox$stats)
coxResultslog$model<-(colnames(EpMeasuresClass1))
coxResultslog$EpitopeCount<-colSums(trainData[,coxResultslog$model])-(0.5*nrow(trainData))
coxResultslog$NumPatients<-apply(trainData[,coxResultslog$model], 2, function(x) sum(x>0.5))
head(coxResultslog)
# Class 1: select models, perform bootstrapping
mymodels<-head(coxResultslog[order(coxResultslog$Dxy, decreasing=T),"model"], n=100)
coxResultsReduced<-data.frame()
for (model in mymodels) {
        myformula<-as.formula(paste('Surv(Overall.Survival..Months., OS.Status)~log(', model, ')+ SimpleStage + Diagnosis.Age + Pass2quant', collapse=""))
        mycox<-cph(myformula, data=trainData, x=T, y=T)
        myvalid<-validate(mycox, B=200, method="boot")
        coxResultsReduced<-rbind(coxResultsReduced, c(mycox$loglik[2], mycox$coefficients[1],  myvalid[1,1:5] ))
}
colnames(coxResultsReduced)<-c("loglik", "coef", "index.orig", "training", "test", "optimism", "index.corrected")
coxResultsReduced$model<-mymodels
coxResultsReduced$EpitopeCount<-colSums(trainData[,coxResultsReduced$model])-(0.5*nrow(trainData))
coxResultsReduced$NumPatients<-apply(trainData[,mymodels], 2, function(x) sum(x>0.5))
coxResultsReduced$RankIndex.c<-rank(-coxResultsReduced$index.corrected)
coxResultsReduced$RankLogLik<-rank(-coxResultsReduced$loglik)
head(coxResultsReduced[order(coxResultsReduced$RankIndex.c),], n=20)
# Class 2: perform first pass ep measure association testing:
coxResultslogc2<-data.frame()
for (f in (colnames(EpMeasuresClass2))) {
myformula<-as.formula(paste('Surv(Overall.Survival..Months., OS.Status)~log(', f, ")+ SimpleStageNumeric + Diagnosis.Age + Pass2quant", collapse=""))
 mycox<-cph(myformula, data=trainDatac2, x=T, y=T)
 coxResultslogc2<-rbind(coxResultslogc2, c(mycox$loglik[2], mycox$stats)) }
colnames(coxResultslogc2)[1]<-c("loglik")
colnames(coxResultslogc2)[2:12]<-names(mycox$stats)
head(coxResultslogc2)
coxResultslogc2$model<-(colnames(EpMeasuresClass2))
coxResultslogc2$EpitopeCount<-colSums(trainDatac2[,coxResultslogc2$model])-(0.5*nrow(trainDatac2))
coxResultslogc2$NumPatients<-apply(trainDatac2[,coxResultslogc2$model], 2, function(x) sum(x>0.5))
head(coxResultslogc2)
# Class 2: select models, perform bootstrapping
mymodelsc2<-head(coxResultslogc2[order(coxResultslogc2$Dxy, decreasing=T),"model"], n=100)

coxResultsReducedc2<-data.frame()
for (model in mymodelsc2) {
        myformula<-as.formula(paste('Surv(Overall.Survival..Months., OS.Status)~log(', model, ")+ SimpleStage + Diagnosis.Age + Pass2quant", collapse=""))
        mycox<-cph(myformula, data=trainDatac2, x=T, y=T)
        myvalid<-validate(mycox, B=200, method="boot")
        coxResultsReducedc2<-rbind(coxResultsReducedc2, c(mycox$loglik[2], mycox$coefficients[1],  myvalid[1,1:5] ))
}
colnames(coxResultsReducedc2)<-c("loglik", "coef", "index.orig", "training", "test", "optimism", "index.corrected")
coxResultsReducedc2$model<-mymodelsc2
coxResultsReducedc2$EpitopeCount<-colSums(trainDatac2[,coxResultsReducedc2$model])-(0.5*nrow(trainDatac2))
coxResultsReducedc2$NumPatients<-apply(trainDatac2[,mymodelsc2], 2, function(x) sum(x>0.5))
coxResultsReducedc2$RankIndex.c<-rank(-coxResultsReducedc2$index.corrected)
coxResultsReducedc2$RankLogLik<-rank(-coxResultsReducedc2$loglik)
head(coxResultsReducedc2[order(coxResultsReducedc2$RankIndex.c),], n=20)
head(coxResultsReduced[order(coxResultsReduced$RankIndex.c),], n=20)
mutationModel<-cph(Surv(Overall.Survival..Months., OS.Status)~log(Mutations+0.5) + SimpleStage + Diagnosis.Age + Pass2quant, data=trainData, x=T, y=T)
mutationModel


#select top epitope measures:
class1EC<-head(coxResultsReduced[order(coxResultsReduced$RankIndex.c),"model"], n=1)
class2EC<-head(coxResultsReducedc2[order(coxResultsReducedc2$RankIndex.c),"model"], n=1)
temp<-data.frame(Class2EC=trainDatac2[,colnames(trainDatac2) %in% class2EC], ID=trainDatac2$ID)
temp2<-merge(trainData, temp, all.x=T)
temp2$Class2EC[is.na(temp2$Class2EC)]<- 0.5
colnames(temp2)[which( colnames(temp2) %in% class1EC)]<-"Class1EC"
trainData<-temp2


testModel<-cph(Surv(Overall.Survival..Months., OS.Status)~log(Class1EC) + log(Class2EC) + SimpleStageNumeric + Diagnosis.Age + EpsNormSelfc1 , data=trainData, x=T, y=T)
testModel2<-cph(Surv(Overall.Survival..Months., OS.Status)~log(Mutations+.5) + Diagnosis.Age + HLA.B , data=trainData, x=T, y=T)


#data frames for model building
mycols4<-c("Class1EC", "Class2EC", "Overall.Survival..Months.", "OS.Status", "Mutations","SimpleStageNumeric","Diagnosis.Age","Person.Gender","CTLA4","CD4","CD19","CD274","CD8A","FOXP3","PDCD1","ImmuneReads","ImmuneClones","EpsNormSelfc1","EpsNormSelfc2","HLA.A","HLA.B","HLA.C","DQA1","DQB1","DRB1", "clones","Pass2quant")
myfullform<-Surv(Overall.Survival..Months., OS.Status) ~ log(Class1EC) + log(Class2EC) + SimpleStageNumeric + Diagnosis.Age + Person.Gender + CTLA4 + CD4 + CD19  + CD274 + CD8A + FOXP3 + PDCD1+log((ImmuneClones+0.5)/(ImmuneReads+0.5)) + EpsNormSelfc1 + EpsNormSelfc2 + HLA.A +HLA.B + HLA.C + DQA1 + DQB1 + DRB1 + clones + Pass2quant 
mutfullform<-Surv(Overall.Survival..Months., OS.Status) ~ log(Mutations + 0.5) + SimpleStageNumeric + Diagnosis.Age + Person.Gender + CTLA4 + CD4 + CD19 + CD274 + CD8A + FOXP3 + PDCD1 +log((ImmuneClones+0.5)/(ImmuneReads+0.5)) +  HLA.A +HLA.B + HLA.C + DQA1 + DQB1 + DRB1 + clones+ Pass2quant 
mutfullformInt<-Surv(Overall.Survival..Months., OS.Status) ~ (log(Mutations + 0.5) + SimpleStageNumeric + Diagnosis.Age + Person.Gender + CTLA4 + CD4 + CD19 + CD8A + FOXP3 + PDCD1 + log((ImmuneClones+0.5)/(ImmuneReads+0.5))+  HLA.A +HLA.B + HLA.C + DQA1 + DQB1 + DRB1 + clones + Pass2quant)^2
myfullformInt<-Surv(Overall.Survival..Months., OS.Status) ~ (log(Class1EC) + log(Class2EC)  + SimpleStageNumeric + Diagnosis.Age + Person.Gender + CTLA4 + CD4 + CD19 + CD274 + CD8A + FOXP3 + PDCD1+log((ImmuneClones+0.5)/(ImmuneReads+0.5)) + EpsNormSelfc1 + EpsNormSelfc2 + HLA.A +HLA.B + HLA.C + DQA1 + DQB1 + DRB1 + clones + Pass2quant)^2
myreducedform <- Surv(Overall.Survival..Months., OS.Status) ~ log(Class1EC) + log(Class2EC) + SimpleStageNumeric + Diagnosis.Age + Person.Gender + clones + EpsNormSelfc1 + EpsNormSelfc2 + log(ImmuneReads + 0.5)
mutreducedform <-  Surv(Overall.Survival..Months., OS.Status) ~ log(Mutations +0.5) + SimpleStageNumeric + Diagnosis.Age + Person.Gender + clones  + log(ImmuneReads + 0.5)
#model building
EpsCox <- stepAIC(cph(Surv(Overall.Survival..Months., OS.Status) ~ log(Class1EC)+ log(Class2EC) + SimpleStageNumeric + Diagnosis.Age, data=na.omit(trainData[,mycols4]), surv=T), scope=myfullform, direction="both", trace=0)
EpsCoxInt <- stepAIC(cph(Surv(Overall.Survival..Months., OS.Status) ~ log(Class1EC) + log(Class2EC) + SimpleStageNumeric + Diagnosis.Age, data=na.omit(trainData[,mycols4]), surv=T), scope=myfullformInt,direction="both", trace=0)
MutCox <- stepAIC(cph(Surv(Overall.Survival..Months., OS.Status) ~ log(Mutations + 0.5) + SimpleStageNumeric + Diagnosis.Age, data=na.omit(trainData[,mycols4]), surv=T), scope=mutfullform, direction="both", trace=0)
MutCoxInt <- stepAIC(cph(Surv(Overall.Survival..Months., OS.Status) ~ log(Mutations + 0.5) + SimpleStageNumeric + Diagnosis.Age, data=na.omit(trainData[,mycols4]), surv=T), scope=mutfullformInt, direction="both", trace=0)


pec_test <- pec( object = list(epitopes = EpsCox, mutations = MutCox), 
				formula = Surv(Overall.Survival..Months., OS.Status)~1, data =  na.omit(trainData[,mycols4]),
				B = 100, verbose = T,keep.index = T, keep.matrix = T, splitMethod = 'Boot632plus')


form <- Surv(Overall.Survival..Months., OS.Status)~ Set + SimpleStageNumeric + Diagnosis.Age
cph(form, data = merge(blcaPatient, sampleSheet), x = T, y = T , surv = T)

###testing 
#set up data frame
testData<-merge(EpMeasuresClass1[which(rownames(EpMeasuresClass1) %in% testSamples),] , blcaPatient, by.x=0, by.y="ID")
colnames(testData)[1]<-"ID"
testDatac2<-merge(EpMeasuresClass2[which(rownames(EpMeasuresClass2) %in% testSamples),] , blcaPatient, by.x=0, by.y="ID")
colnames(testDatac2)[1]<-"ID"
temp<-data.frame(Class2EC=testDatac2[,colnames(testDatac2) %in% class2EC], ID=testDatac2$ID)
temp2<-merge(testData, temp, all.x=T)
temp2$Class2EC[is.na(temp2$Class2EC)]<- 0.5
colnames(temp2)[which( colnames(temp2) %in% class1EC)]<-"Class1EC"
testData<-temp2

# compare models
cph(formula=formula(EpsCox), data=testData)
cph(formula=formula(MutCox), data=testData)

###Random Forests
rfs_epitope <- rfsrc(myreducedform, data = na.omit(trainData[,mycols4]),forest = TRUE, ntree = 2000, importance=T, tree.err=T, nsplit = 1)
rfs_mutation <- rfsrc(mutreducedform, data = na.omit(trainData[,mycols4]),forest = TRUE, ntree = 2000, importance=T, tree.err=T, nsplit = 1)

#fix samples with negative followup time
testData$Overall.Survival..Months.[which(testData$Overall.Survival..Months. <0)] <- 0 

#predict with rf modles
rfs_epitope_predictions <- predict(rfs_epitope,  na.omit(testData[,mycols4]), importance = T, tree.err = T)
rfs_mutation_predictions <- predict(rfs_mutation,  na.omit(testData[,mycols4]), importance = T, tree.err = T)

#plot RF features
gg_md_epitope <- gg_minimal_depth(rfs_epitope_predictions) 
gg_md_plot[[i]] <- plot( gg_md )
gg_minmal_vimp_plot[[i]] <- plot( gg_minimal_vimp(gg_md) )

#Compare models with PEC
pec_rf <- pec( object = list(epitopes = rfs_epitope, mutations = rfs_mutation), 
				formula = Surv(Overall.Survival..Months., OS.Status)~1, data =  na.omit(testData[,mycols4]),
				verbose = T,keep.index = T, keep.matrix = T, splitMethod = 'none')

#variable importance
sf_variable_importance <-plot.variable(fit_master2[['D']]$classification$rfsrc$pred_on_B_J_class, xvar.names = c('PC_gene1',  
			'PC_gene2'), partial = T, smooth.lines = T, main = 'WC Temporal Cluster D on DC') 
			#'Sum.reads', 'UMI.per.clone','PC_gene2', 'PC_exon1', 'PC_exon2'), partial = T, smooth.lines = T)
#
plot(fit_master2[['B']]$classification$rfsrc$pred_on_B_J_class)
#
rsf_variable_importance <-plot.variable(fit_master2[['B']]$classification$rfsrc$pred_on_B_J_class, xvar.names = c('PC_gene2',  
			'PC_gene3'), partial = T, smooth.lines = T, main = 'WC Temporal Cluster B on Janssen Blood') 
#
plot(fit_master2[['17']]$classification$rfsrc$pred_on_PC_class)
#
rsf_variable_importance <-plot.variable(fit_master2[['17']]$classification$rfsrc$pred_on_PC_class, xvar.names = c('PC_gene1',  
			'PC_gene2', 'PC_exon2', 'median_score_gene'), partial = T, smooth.lines = T, main = 'WC 17 days on PC') 


rsfF <- rfsrc(formfull, data = na.omit(clindat[idx,idxc]), 
                    forest = TRUE, ntree = 2000, importance=T, tree.err=T, nsplit = 0)
rsf_variable_importance <-plot.variable(rsfF, xvar.names = c('clones',  'Mutation.Count', 'ImmuneReads'), partial = T, smooth.lines = T, surv.type = 'surv')

