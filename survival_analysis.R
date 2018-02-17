library(survival)
library(GGally)
library(scales)
library(ggplot2)

clindat<-PatientInfo[which(PatientInfo$Study.ID=="skcm_tcga"),]
#cox proportional hazard models:
fullmodel<-coxph(Surv(Overall.Survival..Months., Overall.Survival.Status=="DECEASED")~Mutation.Count + SimpleStage + Diagnosis.Age + CTLA4 + Year.Cancer.Initial.Diagnosis, data=clindat)
reducedmodel<-coxph(Surv(Overall.Survival..Months., Overall.Survival.Status=="DECEASED")~Mutation.Count + SimpleStage + Diagnosis.Age + CTLA4, data=clindat)
#likelihood ratio test: NOTE all models must have same number of data points.  
anova(fullmodel, reducedmodel)

#kaplan-meier curve:
pdf("km_plot.pdf")
#set your cut point for using a continuous variable.  
medianVal<-median(clindat$Mutation.Count, na.rm=T)
#create the plots
ggsurv(survfit(Surv(Overall.Survival..Months., Overall.Survival.Status=="DECEASED")~Mutation.Count < medianVal, data=clindat))
ggsurv(survfit(Surv(Overall.Survival..Months., Overall.Survival.Status=="DECEASED")~Mutation.Count < medianVal, data=clindat), CI=T)
# a nicer looking 95% CI.  apparently imperfect: see: http://stackoverflow.com/questions/33874909/how-do-i-add-shading-and-color-to-the-confidence-intervals-in-ggplot-2-generated
ggsurv(survfit(Surv(Overall.Survival..Months., Overall.Survival.Status=="DECEASED")~Mutation.Count < medianVal, data=clindat), CI=F) + geom_ribbon(aes(ymin=low,ymax=up,fill=group),alpha=0.3)
dev.off() 