# MDC.RData should have your expression data. (ie dat1 and dat2 in the function)
# I use this script by creating a submit lsf script, then submitting it ~200 times).  

library(coexpp)
library(plyr)
load("MDC.RData")

#change the seed so different threads are different. 
t<-as.numeric(Sys.time()) 
set.seed((t - floor(t)) * 1e8 );

modDiffCon<-NULL
modDiffCon2<-NULL

permuteFDR <- 
function(dat1=aordatExpr, dat2=mamdatExpr, permutations=1000, modules=moduleColors) {
        modDiffCon<-data.frame(row.names=names(table(modules)))
        modDiffCon2<-data.frame(row.names=names(table(modules)))
        bigdat<-rbind(dat1,dat2)
        for (i in 1:permutations) {
                print(i)
                tempaorAdjacency<-adjacency(dat1[,sample(ncol(dat1))], type="signed", power=softPower)
                tempmamAdjacency<-adjacency(dat2[,sample(ncol(dat2))], type="signed", power=softPower)
                tempIC_1<-intramodularConnectivity(tempaorAdjacency, modules)
                tempIC_2<-intramodularConnectivity(tempmamAdjacency, modules)
                tempIC_1$module<-modules
                tempIC_2$module<-modules
                tempdat<-ddply(tempIC_2, .(module), summarize, sum(kWithin))[,2]/ddply(tempIC_1, .(module), summarize, sum(kWithin))[,2]
                modDiffCon[,i]<-tempdat

                tempaorAdjacency<-adjacency(bigdat[sample(nrow(bigdat), size=nrow(dat1)),], type="signed", power=softPower)
                tempmamAdjacency<-adjacency(bigdat[sample(nrow(bigdat), size=nrow(dat2)),], type="signed", power=softPower)
                tempIC_1<-intramodularConnectivity(tempaorAdjacency, modules)
                tempIC_2<-intramodularConnectivity(tempmamAdjacency, modules)
                tempIC_1$module<-modules
                tempIC_2$module<-modules
                tempdat<-ddply(tempIC_2, .(module), summarize, sum(kWithin))[,2]/ddply(tempIC_1, .(module), summarize, sum(kWithin))[,2]
                modDiffCon2[,i]<-tempdat
        }
	rownames(modDiffCon)<-names(table(modules))
        colnames(modDiffCon)<-paste("MDC1", colnames(modDiffCon), sep=".")
        colnames(modDiffCon2)<-paste("MDC2", colnames(modDiffCon2), sep=".")
        rownames(modDiffCon2)<-names(table(modules))
        MDC<-cbind(modDiffCon, modDiffCon2)
        return(MDC)
}

myIndex<-round(runif(1, 1, 1e6))
myOutfileName<-paste("output", myIndex, "txt", sep=".")
print(myOutfileName)
mdc<-permuteFDR(dat1=multiExpr[[2]]$data, dat2=multiExpr[[1]]$data, permutations=5, modules=moduleColors)

write.table(mdc, file=paste("MDC", myOutfileName, sep="."), quote=F, sep="\t")

