load("blca.round2.forBojan.RData")
wideScreen()
library(plyr)
library(survival)
library(GGally)
library(ggplot2)
library(rms)
library(MASS)
load("July2017/fixed.round1.round2.epitopes.RData")

#fix the numeric stage:
blcaPatient$SimpleStage<-(factor(blcaPatient$SimpleStage, levels=c("StageI+II", "Stage III", "Stage IV")))
blcaPatient$SimpleStageNumeric<-as.numeric(blcaPatient$SimpleStage)
#fix time less than 0
blcaPatient[which(blcaPatient$Overall.Survival..Months. < 0),"Overall.Survival..Months."]<-0


#HLA Expression
HLA_thresh<-ddply(individualEpitopes[which(individualEpitopes$HLALocus %in% c("A", "B", "C")),], .(ID), summarise, noHLAcut=sum(ic50 <=500)+0.5, HLAcut25=sum(ic50<=500 & HLApercentileCancer >0.10 )+0.5)
test_survival<-merge(HLA_thresh, blcaPatient)
survModel1<-cph(Surv(Overall.Survival..Months., OS.Status)~log(noHLAcut)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
survModel2<-cph(Surv(Overall.Survival..Months., OS.Status)~log(HLAcut25)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
mutationModel<-cph(Surv(Overall.Survival..Months., OS.Status)~log(Mutations + 0.5)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
#results: mutations >>> allHLA > only expressed hla

#Location
Loc_filter<-ddply(individualEpitopes[which(individualEpitopes$HLALocus %in% c("A", "B", "C")),], .(ID), summarise, noLocCut=sum(ic50 <=500)+0.5, LocPMS=sum(ic50<=500 & Location %in% c("plasma membrane", "secreted"))+0.5)
test_survival<-merge(Loc_filter, blcaPatient)
survModel1<-cph(Surv(Overall.Survival..Months., OS.Status)~log(noLocCut)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
survModel2<-cph(Surv(Overall.Survival..Months., OS.Status)~log(LocPMS)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
mutationModel<-cph(Surv(Overall.Survival..Months., OS.Status)~log(Mutations + 0.5)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
#results: mutations >>> plasma membrane +secreted ~ nocut

#VAF
VAF_filter<-ddply(individualEpitopes[which(individualEpitopes$HLALocus %in% c("A", "B", "C")),], .(ID), summarise, noVAFcut=sum(ic50 <=500)+0.5, VAFlow=sum(ic50 <=500 & VAFfraction < 0.25, na.rm=T)+0.5, VAFhigh=sum(ic50 <=500 & VAFfraction > 0.75, na.rm=T)+0.5)
test_survival<-merge(VAF_filter, blcaPatient)
survModel1<-cph(Surv(Overall.Survival..Months., OS.Status)~log(noVAFcut)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
survModel2<-cph(Surv(Overall.Survival..Months., OS.Status)~log(VAFlow)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
survModel3<-cph(Surv(Overall.Survival..Months., OS.Status)~log(VAFhigh)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
mutationModel<-cph(Surv(Overall.Survival..Months., OS.Status)~log(Mutations + 0.5)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
#results mutations > nocut ~ vaflow > vafhigh

#Gene Expression
GE_filter<-ddply(individualEpitopes[which(individualEpitopes$HLALocus %in% c("A", "B", "C")),], .(ID), plyr::summarise, noGEcut=sum(ic50 <=500)+0.5, GEexpr=sum(ic50 <=500 & logCPM > -2.5, na.rm=T)+0.5, GEmoderate=sum(ic50 <=500 & logCPM > 2.5, na.rm=T)+0.5, GEhigh=sum(ic50 <=500 & logCPM > 5, na.rm=T)+0.5)
cor(GE_filter[,-1])
test_survival<-merge(GE_filter, blcaPatient)
survModel1<-cph(Surv(Overall.Survival..Months., OS.Status)~log(noGEcut)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
survModel2<-cph(Surv(Overall.Survival..Months., OS.Status)~log(GEexpr)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
survModel3<-cph(Surv(Overall.Survival..Months., OS.Status)~log(GEmoderate)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
survModel4<-cph(Surv(Overall.Survival..Months., OS.Status)~log(GEhigh)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
mutationModel<-cph(Surv(Overall.Survival..Months., OS.Status)~log(Mutations + 0.5)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
#results mutations > nocut ~ (2.5) ~ (-2.5) > (5)

#combine GE, VAF
GEVAF_filter<-ddply(individualEpitopes[which(individualEpitopes$HLALocus %in% c("A", "B", "C")),], .(ID), plyr::summarise, nocut=sum(ic50 <=500)+0.5, GEVAF=sum(ic50 <=500 & logCPM > 2.5 & VAFfraction < 0.25, na.rm=T)+0.5)
cor(GEVAF_filter[,-1])
test_survival<-merge(GEVAF_filter, blcaPatient)
survModel1<-cph(Surv(Overall.Survival..Months., OS.Status)~log(nocut)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
survModel2<-cph(Surv(Overall.Survival..Months., OS.Status)~log(GEVAF)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
mutationModel<-cph(Surv(Overall.Survival..Months., OS.Status)~log(Mutations + 0.5)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
#results: the combo works! still not as good as mutations alone.  

#examine Mutations
mutationsDF<-unique(individualEpitopes[, c("ID", "chr", "pos", "logCPM", "VAF", "VAFfraction", "Location")])
#VAF filter
VAF_filter<-ddply(mutationsDF, .(ID), summarise, VAFlow=sum(VAFfraction < 0.25, na.rm=T)+0.5, VAFhigh=sum(VAFfraction > 0.75, na.rm=T)+0.5, nocut=sum(VAFfraction < 1.25, na.rm=T)+0.5)
test_survival<-merge(VAF_filter, blcaPatient)
mutationModel<-cph(Surv(Overall.Survival..Months., OS.Status)~log(nocut)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
survModel1<-cph(Surv(Overall.Survival..Months., OS.Status)~log(VAFlow)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
survModel2<-cph(Surv(Overall.Survival..Months., OS.Status)~log(VAFhigh)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival)
#results: filtering VAF innefective

mutation_filter<-ddply(mutationsDF, .(ID), plyr::summarise, nocut=sum(pos > -1)+0.5, GE=sum(logCPM > 2.5, na.rm=T)+0.5, VAF=sum(VAFfraction < 0.25, na.rm=T)+0.5, GEVAF=sum(logCPM > 2.5 & VAFfraction < 0.25, na.rm=T)+0.5, GEVAFLOC=sum(logCPM > 2.5 & VAFfraction < 0.25 & Location %in% c("plasma membrane", "secreted"), na.rm=T)+0.5, LOC=sum(Location %in% c("plasma membrane", "secreted"), na.rm=T)+0.5, GELOC=sum(logCPM > 2.5 & Location %in% c("plasma membrane", "secreted"), na.rm=T)+0.5, LOCVAF=sum(VAFfraction < 0.25 & Location %in% c("plasma membrane", "secreted"), na.rm=T)+0.5) 
cor(mutation_filter[,-1])
test_survival_muts<-merge(mutation_filter, blcaPatient)
mutationModel_muts<-cph(Surv(Overall.Survival..Months., OS.Status)~log(nocut)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival_muts)
survModel1_muts<-cph(Surv(Overall.Survival..Months., OS.Status)~log(GE)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival_muts)
survModel2_muts<-cph(Surv(Overall.Survival..Months., OS.Status)~log(VAF)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival_muts)
survModel3_muts<-cph(Surv(Overall.Survival..Months., OS.Status)~log(GEVAF)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival_muts)
survModel4_muts<-cph(Surv(Overall.Survival..Months., OS.Status)~log(LOC)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival_muts)
survModel5_muts<-cph(Surv(Overall.Survival..Months., OS.Status)~log(GEVAFLOC)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival_muts)
survModel6_muts<-cph(Surv(Overall.Survival..Months., OS.Status)~log(GELOC)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival_muts)

# results: GEVAFLOC > GEVAF ~ nocut ~ GE~LOC > VAF
mutationResults<-as.data.frame(rbind(mutationModel_muts$stats, survModel1_muts$stats, survModel2_muts$stats, survModel3_muts$stats,survModel4_muts$stats, survModel5_muts$stats, survModel6_muts$stats))
mutationResults$test<-c("mutations", "GEhigh", "VAFlow", "GEVAF", "LocPMS", "GEVAFLOC", "GELOC")
mutationResults
colnames(mutationResults)[3]<-"ModelLR"
ggplot(mutationResults, aes(x=test, y=Dxy)) + geom_point(size=4)

#epitopes comparison
epitopes_filter<-ddply(individualEpitopes[which(individualEpitopes$HLALocus %in% c("A", "B", "C")),], .(ID), plyr::summarise, nocut=sum(ic50 <=500)+0.5, GE=sum(ic50 <=500 & logCPM > 2.5, na.rm=T)+0.5, VAF=sum(ic50<=500 & VAFfraction < 0.25, na.rm=T)+0.5, GEVAF=sum(ic50 <=500 & logCPM > 2.5 & VAFfraction < 0.25, na.rm=T)+0.5, GEVAFLOC=sum(ic50 <=500 & logCPM > 2.5 & VAFfraction < 0.25 &  Location %in% c("plasma membrane", "secreted"), na.rm=T)+0.5, LOC=sum(ic50 <=500 &  Location %in% c("plasma membrane", "secreted"), na.rm=T)+0.5, GELOC=sum(ic50 <=500 & logCPM > 2.5 & Location %in% c("plasma membrane", "secreted"), na.rm=T)+0.5, LOCVAF=sum(ic50 <=500 & VAFfraction < 0.25 &  Location %in% c("plasma membrane", "secreted"), na.rm=T)+0.5)
cor(epitopes_filter[,-1])
test_survival<-merge(epitopes_filter, blcaPatient)
mutationModel<-cph(Surv(Overall.Survival..Months., OS.Status)~log(nocut)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival,surv=T)
survModel1<-cph(Surv(Overall.Survival..Months., OS.Status)~log(GE)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival,surv=T)
survModel2<-cph(Surv(Overall.Survival..Months., OS.Status)~log(VAF)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival,surv=T)
survModel3<-cph(Surv(Overall.Survival..Months., OS.Status)~log(GEVAF)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival,surv=T)
survModel4<-cph(Surv(Overall.Survival..Months., OS.Status)~log(LOC)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival,surv=T)
survModel5<-cph(Surv(Overall.Survival..Months., OS.Status)~log(GEVAFLOC)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival,surv=T)
survModel6<-cph(Surv(Overall.Survival..Months., OS.Status)~log(GELOC)+ SimpleStageNumeric + Diagnosis.Age, data=test_survival,surv=T)
epitopeResults<-as.data.frame(rbind(mutationModel$stats, survModel1$stats, survModel2$stats, survModel3$stats,survModel4$stats, survModel5$stats, survModel6$stats))
epitopeResults$test<-c("mutations", "GEhigh", "VAFlow", "GEVAF", "LocPMS", "GEVAFLOC", "GELOC")
epitopeResults
colnames(epitopeResults)[3]<-"ModelLR"
ggplot(epitopeResults, aes(x=test, y=Dxy)) + geom_point(size=4)

mutationResults$measure<-"mutations"
epitopesResults$measure<-"epitopes"
myresults<-rbind(epitopeResults, mutationResults)
ggplot(myresults, aes(x=test, y=Dxy, color=measure)) + geom_point(size=4)
#results: by Dxy, epitopes wins out.  By LR, usually mutations wins, but GEVAFLOC eps wins out.  Overall GEVAFLOC always the best.
#export:
myresults_blca<-myresults
save(myresults_blca, file="July2017/eps_muts_blca.RData")
test_survival_muts_blca<-test_survival_muts
test_survival_eps_blca<-test_survival
save(test_survival_muts_blca,test_survival_eps_blca, file="July2017/eps_muts_blca_surv.RData")

#pec testing:
library(pec)
library(doMC)
registerDoMC()
testframe<-na.omit(test_survival[,c("nocut","GE","VAF", "GEVAF","LOC", "GEVAFLOC","GELOC", "Overall.Survival..Months.","OS.Status", "SimpleStageNumeric","Diagnosis.Age" )])
pec_train_eps <- pec( object = list( GE = survModel1, GELOC=survModel6, GELVAFLOC=survModel5, unfiltered=mutationModel ), 
				formula = Surv(Overall.Survival..Months., OS.Status)~1, data = testframe ,
				B = 100, verbose = T,keep.index = T, keep.matrix = T, splitMethod = 'Boot632plus')

testframe<-na.omit(test_survival_muts[,c("nocut","GE","VAF", "GEVAF","LOC", "GEVAFLOC","GELOC", "Overall.Survival..Months.","OS.Status", "SimpleStageNumeric","Diagnosis.Age" )])
pec_train_muts <- pec( object = list( GE = survModel1, GELOC=survModel6, GELVAFLOC=survModel5, unfiltered=mutationModel ), 
				formula = Surv(Overall.Survival..Months., OS.Status)~1, data = testframe ,
				B = 100, verbose = T,keep.index = T, keep.matrix = T, splitMethod = 'Boot632plus')
				
#


#patient pairing
#list of ids:
myPatients<-blcaPatient$ID[(blcaPatient$ID %in% individualEpitopes$ID)]
common_alleles<-unique(individualEpitopes[,c("ID", "allele", "HLALocus")]
sort(table(common_alleles[which(common_alleles$HLALocus=="A"),"allele"]))
#start with A*0101 and A*0202
length(unique(individualEpitopes[which(individualEpitopes$allele %in% c("HLA-A*0201", "HLA-A*0101")),"ID"]))
myPatients2<-unique(individualEpitopes[which(individualEpitopes$allele %in% c("HLA-A*0201", "HLA-A*0101")),"ID"])
pair_frame<-data.frame(index=1:(210*209))
rownumber<-0
for (f in 1:210) { for (g in 1:210) { if (f==g) next ; rownumber<-rownumber+1 ; pair_frame$sample1[rownumber]<-myPatients2[f] ; pair_frame$sample2[rownumber]<-myPatients2[g] }}
for (z in 1:nrow(pair_frame)) {
	pair_frame[z, "ks.test"]<-ks.test(individualEpitopes[which(individualEpitopes$allele %in% c("HLA-A*0201", "HLA-A*0101") & individualEpitopes$ID==pair_frame[z,"sample1"]),"ic50"], individualEpitopes[which(individualEpitopes$allele %in% c("HLA-A*0201", "HLA-A*0101") & individualEpitopes$ID==pair_frame[z,"sample2"]),"ic50"], alternative='less')$p.value	
	print(z)
}
myPatients2df<-data.frame(sample1=myPatients2)
for (x in 1:nrow(myPatients2df)) {
	myPatients2df[x,"sample1alleles"]<-length(unique(individualEpitopes[which(individualEpitopes$allele %in% c("HLA-A*0201", "HLA-A*0101") & individualEpitopes$ID==myPatients2df[x,"sample1"]),"allele"]))
}
pair_frame2<-merge(pair_frame, myPatients2df)
colnames(myPatients2df)<-c("sample2", "sample2alleles")
pair_frame3<-merge(pair_frame2, myPatients2df)
for (x in 1:nrow(pair_frame3)) {
	pair_frame3[x,"dVDJ"]<-blcaPatient[which(blcaPatient$ID==pair_frame3[x,"sample1"]),"ImmuneReads"]-blcaPatient[which(blcaPatient$ID==pair_frame3[x,"sample2"]),"ImmuneReads"]
}
ggplot(pair_frame3, aes(x=dVDJ, y=-log10(ks.test))) + geom_point()
ggplot(pair_frame3[which(pair_frame3$dVDJ >= 0),], aes(x=dVDJ, y=-log10(ks.test))) + geom_point(alpha=0.4) + scale_x_log10()
ggplot(pair_frame3[which(pair_frame3$dVDJ >0),], aes(x=log10(dVDJ), y=-log10(ks.test))) + geom_point(alpha=0.1) + stat_density2d(aes(fill=..level..), geom='polygon')


save.image("blca.RData")