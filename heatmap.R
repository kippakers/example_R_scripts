library(method)
library(methods)
library(edgeR)
library(gplots)
#red/green color pallete
my_palette <- colorRampPalette(c("red", "black", "green"))(n = 20)

#load your data:

dat<- # place here your samples x Genes matrix
phenos<- #place here a samples x Phenos matrix.  should match dat


#assign colors to your samples.  Just a guess of how this should look.  This is saying if the sample PC1 is less than -10 color it blue, if not, orange.  
mycolors<-ifelse(phenos$PC1 < -10, "blue", "orange")
#draw a heatmap:
heatmap.2(dat, Rowv = T, scale = 'none', Colv = T, trace = 'none', srtCol = 45, cexRow = 0.5, cexCol=0.5, margins= c(7, 7), col = my_palette, labCol=colnames(dat), colCol=mycolors)

