library(WGCNA)
library(gdata)
library(gplots)
options(stringsAsFactors = FALSE);
load("Consensus-NetworkConstruction-man.RData")
load("tomplotData.RData")

selectColors<-c("green")
mamselectdiss<-mamDiss[which(moduleColors %in% selectColors), which(moduleColors %in% selectColors)]
consTreeSelectMAM<-hclust(as.dist(mamselectdiss), method="average")
aorselectdiss<-aorDiss[which(moduleColors %in% selectColors), which(moduleColors %in% selectColors)]
consTreeSelectAOR<-hclust(as.dist(aorselectdiss), method="average")

#split matrix: MAM on the upper triangle. 
m1<-mamselectdiss[consTreeSelectMAM$order, consTreeSelectMAM$order]^4
m2<-aorselectdiss[consTreeSelectMAM$order, consTreeSelectMAM$order]^4
splitselectDiss<-m1
lowerTriangle(splitselectDiss)<-lowerTriangle(m2)

#cb palette
my_palette<-colorRampPalette(c("#D55E00", "#9ad0f3"))
png("Plots/TOMs/split_greenCon_mamDendro_CB.png", height=1000, width=1000, res=200)
heatmap.2(splitselectDiss, dendrogram="none", Rowv=F, Colv=F, labRow="", labCol="", trace="none", col=my_palette(250))
dev.off()
png("Plots/TOMs/split_greenCon_mamDendro.png", height=1000, width=1000, res=200)
heatmap.2(splitselectDiss, dendrogram="none", Rowv=F, Colv=F, labRow="", labCol="", trace="none")
dev.off()


#other modules:
for (f in c("greenyellow", "lightyellow", "pink","antiquewhite4")) {
	selectColors<-f
	mamselectdiss<-mamDiss[which(moduleColors %in% selectColors), which(moduleColors %in% selectColors)]
	consTreeSelectMAM<-hclust(as.dist(mamselectdiss), method="average")
	aorselectdiss<-aorDiss[which(moduleColors %in% selectColors), which(moduleColors %in% selectColors)]
	consTreeSelectAOR<-hclust(as.dist(aorselectdiss), method="average")
	#split matrix: MAM on the upper triangle.  this puts mam on the  
	m1<-mamselectdiss[consTreeSelectMAM$order, consTreeSelectMAM$order]^4
	m2<-aorselectdiss[consTreeSelectAOR$order, consTreeSelectAOR$order]^4
	splitselectDiss<-m1
	lowerTriangle(splitselectDiss)<-lowerTriangle(m2)
	#cb palette
	my_palette<-colorRampPalette(c("#D55E00", "#9ad0f3"))
	png(paste("Plots/split2_", f, "Con_mamDendro_CB.png", sep=""), height=1000, width=1000, res=200)
	heatmap.2(splitselectDiss, dendrogram="none", Rowv=F, Colv=F, labRow="", labCol="", trace="none", col=my_palette(250))
	dev.off()
	png(paste("Plots/split2_", f, "Con_mamDendro.png", sep=""), height=1000, width=1000, res=200)
	heatmap.2(splitselectDiss, dendrogram="none", Rowv=F, Colv=F, labRow="", labCol="", trace="none")
	dev.off()
}
###Full Plot
m2<-aorDiss[consTreeMAM$order, consTreeMAM$order]^4
m1<-mamDiss[consTreeMAM$order, consTreeMAM$order]^4
splitDiss<-m1
lowerTriangle(splitDiss)<-lowerTriangle(m2)
png("Plots/TOMs/split_all_mamDendro_CB.png", height=2000, width=2000, res=150)
heatmap.2(splitDiss, dendrogram="none", Rowv=F, Colv=F, labRow="", labCol="", trace="none")
dev.off()


