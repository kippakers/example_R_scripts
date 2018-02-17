library(limma)
library(edgeR)

counts <- read.table( 'hb_ucsc.counts.matrix2', row.names = 1 , header=T)
phenotype <- read.table( 'hb.phenotype_mod', row.names = 1, header=T)

#counts <- counts[ , order(colnames(counts))]

dim(phenotype)
dim(counts)
summary(phenotype)
head(phenotype)
summary(phenotype$Phenotype)

#exclude xenograft PDX, R-PDX for now

phenotype_NT_T_R <- (phenotype[phenotype$Phenotype != "PDX" & phenotype$Phenotype != "R-PDX",])
counts_NT_T_R <- counts[, rownames(phenotype_NT_T_R)]
dim(counts_NT_T_R)

summary(phenotype_NT_T_R$Phenotype)
phenotype_NT_T_R$Phenotype <- factor(phenotype_NT_T_R$Phenotype, c("NT", "R", "T"))
summary(phenotype_NT_T_R$Phenotype)

dge <- DGEList(counts = counts_NT_T_R)
dge <- dge[ rowSums(cpm(dge$counts) > 1) >20,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

v_NT_T_R <- voom( dge, model.matrix( ~as.factor(Phenotype), data = phenotype_NT_T_R) )

dim(v_NT_T_R)

##fit_NT_T_R <- eBayes(lmFit( v_NT_T_R, model.matrix( ~as.factor(Phenotype) + as.factor(Machine) + as.factor( FlowcellID), data = phenotype_NT_T_R) ) )
#under-determined FlowcellID

#empty MACHINE covariate
fit_NT_T_R <- eBayes(lmFit( v_NT_T_R, model.matrix( ~as.factor(Phenotype) + as.factor(Machine), data = phenotype_NT_T_R) ) )
summary(decideTests(fit_NT_T_R))

#minimal model suffices: test for contrasts above logFC=0
fit_NT_T_R <- eBayes(lmFit( v_NT_T_R, model.matrix( ~as.factor(Phenotype) , data = phenotype_NT_T_R) ) )
summary(decideTests(fit_NT_T_R))

##playing around

#topTable( fit_NT_T_R, coef = 2, num = 20, sort.by = "p", p.value = .05)
#topTable( fit_NT_T_R, coef = 3, num = 20, sort.by = "p", p.value = .05)
#topTable( fit_NT_T_R, coef = 3, num = Inf, sort.by = "p", p.value = .05)["IGF2",]
#topTable( fit_NT_T_R, coef = 3, num = Inf, sort.by = "p", p.value = .05)["H19",]
#topTable( fit_NT_T_R, coef = 2, num = Inf, sort.by = "p", p.value = .05)["H19",]
#topTable( fit_NT_T_R, coef = 2, num = Inf, sort.by = "p", p.value = .05)["IGF2",]
#topTable( fit_NT_T_R, coef = 2, num = Inf, sort.by = "p", p.value = .1)["IGF2",]
#topTable( fit_NT_T_R, coef = 2, num = Inf, sort.by = "p", p.value = 1)["IGF2",]
#topTable( fit_NT_T_R, coef = 2, num = Inf, sort.by = "p", p.value = 1)["IGF1",]

dim(topTable( fit_NT_T_R, coef = 2, num = Inf, sort.by = "p", p.value = .05))
dim(topTable( fit_NT_T_R, coef = 3, num = Inf, sort.by = "p", p.value = .05))

#now test empirical p-values against a nonzero abs(logFC) threshold of 1 for the R comparison, since the group is small (4)
#and there is a strong signal otherwise in coef=3 (epigenetic misregulation confers such a strong DEG signal?)

fit_NT_T_R2 <- treat(lmFit( v_NT_T_R , model.matrix( ~as.factor(Phenotype), data = phenotype_NT_T_R) ) , lfc = 1)
summary(decideTests(fit_NT_T_R2))

#

dim(topTreat( fit_NT_T_R2, coef = 2, num = Inf, p.value = .05))
dim(topTreat( fit_NT_T_R2, coef = 3, num = Inf, p.value = .05))


#write and plot results


write.table( topTreat( fit_NT_T_R2, coef = 2, num = Inf, p.value = .05), 'DEG_topTreat_lfc_1__R_q_05', row.names = T, col.names = T, sep = "\t", quote = F)
write.table( topTreat( fit_NT_T_R2, coef = 3, num = Inf, p.value = .05), 'DEG_topTreat_lfc_1__NT_q_05', row.names = T, col.names = T, sep = "\t", quote = F)

write.table( topTable( fit_NT_T_R, coef = 3, num = Inf, sort.by = "p", p.value = .05), 'DEG_topTable_lfc_0__R_q_05', row.names = T, col.names = T, sep = "\t", quote = F)
write.table( topTable( fit_NT_T_R, coef = 2, num = Inf, sort.by = "p", p.value = .05), 'DEG_topTable_lfc_0__NT_q_05', row.names = T, col.names = T, sep = "\t", quote = F)



##slightly modify volcanoplot from limma to better display outliers#####


volcanoplot2 <- function (fit, coef = 1, highlight = 0, names = fit$genes$ID,
    xlab = "Log Fold Change", ylab = "Log Odds", pch = 16, cex = 0.35,
    ...)
{
    if (!is(fit, "MArrayLM"))
        stop("fit must be an MArrayLM")
    if (is.null(fit$lods))
        stop("No B-statistics found, perhaps eBayes() not yet run")
    x <- as.matrix(fit$coef)[, coef]
    y <- as.matrix(fit$lods)[, coef]
    plot(x, y, xlab = xlab, ylab = ylab, pch = pch, cex = cex,
        ...)
    if (highlight > 0) {
        if (is.null(names))
            names <- 1:length(x)
        names <- as.character(names)
        o <- order(x, decreasing = TRUE)
        on <- order(-x, decreasing=TRUE)
        i <- o[1:highlight]
        text(x[i], y[i], labels = substring(names[i], 1, 8),
            cex = 0.5, col = "blue", pos=3)
        i2 <- on[1:round(highlight)]
        text(x[i2], y[i2], labels = substring(names[i2], 1, 8),
            cex = 0.5, col = "red", pos = 3)
    }
    invisible()
}   

##################################


pdf('HB_basic_model.pdf')
#voomfit
v_NT_T_R <- voom( dge, model.matrix( ~as.factor(Phenotype), data = phenotype_NT_T_R) , plot = T)
#qq plots
qqt(fit_NT_T_R$t[,2],df=fit_NT_T_R$df.residual + fit_NT_T_R$df.prior)
abline(0,1)
#MDS plots
plotMDS(v_NT_T_R, col = as.numeric(phenotype_NT_T_R$Phenotype), labels = colnames( v_NT_T_R), cex = 0.6, main = "voom HB: model.matrix( ~as.factor(Phenotype), data = phenotype_NT_T_R) ) ",cex.main=0.7,cex.lab=0.8 )
legend( "bottomright", pch = 20, col = c("black","red", "green"), legend = levels(phenotype_NT_T_R$Phenotype) , bty = 'n', cex = .75)
#volcanoplots using modified volcanoplots function
volcanoplot2( fit_NT_T_R, highlight = 55, names = rownames(fit_NT_T_R), cex = 0.1, cex.lab = 0.8, pch = 16, main = "fit_NT_T_R: NT vs T contrast", coef = 3)
volcanoplot2( fit_NT_T_R, highlight = 55, names = rownames(fit_NT_T_R), cex = 0.1, cex.lab = 0.8, pch = 16, main = "fit_NT_T_R: R vs (T & NT) contrast", coef = 2)
#volcanoplot of R condition on treated fit
volcanoplot2( eBayes(fit_NT_T_R2), highlight = 55, names = rownames(fit_NT_T_R2), cex = 0.1, cex.lab = 0.8, pch = 16, main = "fit_NT_T_R2: R vs (T & NT) contrast, treat min_lfc=1", coef = 2)
#MA plots
plot(fit_NT_T_R,coef=2,main = "MA plot for R vs (NT & T) contrast")
plot(fit_NT_T_R,coef=3, main = "MA plot for T vs (NT & R) contrast")
dev.off()

#savehistory('HB_gene_level.RHistory')

#corfit <- duplicateCorrelation(v_full, model.matrix(~as.factor( Phenotype) , data = phenotype), block = phenotype$PatientID)
#fit_block_full <- eBayes(lmFit( v_full, model.matrix(~as.factor( Phenotype) , data = phenotype), block = phenotype$PatientID, correlation = corfit$consensus.correlation))
#fit_block_full <- eBayes(lmFit( voom(dge, model.matrix(~as.factor(Phenotype), data = phenotype), correlation = corfit$consensus.correlation, block = phenotype$PatientID), 
#model.matrix(~as.factor( Phenotype) , data = phenotype), block = phenotype$PatientID, correlation = corfit$consensus.correlation))
#cor(topTable(fit_full, num = nrow(counts), sort.by = "none")$adj.P.Val, topTable(fit_block_full, num = nrow(counts), sort.by = "none")$adj.P.Val  )
#fit_block_full_double <- eBayes(lmFit( voom(dge, model.matrix(~as.factor(Phenotype), data = phenotype), correlation = corfit$consensus.correlation, block = phenotype$PatientID), 
#model.matrix(~as.factor( Phenotype) , data = phenotype), block = phenotype$PatientID, correlation = corfit$consensus.correlation))
#fit_block_full_single <- eBayes(lmFit( v_full, model.matrix(~as.factor( Phenotype) , data = phenotype), block = phenotype$PatientID, correlation = corfit$consensus.correlation))
#cor(topTable(fit_block_full_double, num = nrow(counts), sort.by = "none")$adj.P.Val, topTable(fit_block_full_single, num = nrow(counts), sort.by = "none")$adj.P.Val  )




####Full analysis with patientID blocking

counts <- read.table( 'hb_ucsc.counts.matrix2', row.names = 1 , header=T)
counts <- counts[ , order(colnames(counts))]
phenotype <- read.table( 'hb.phenotype_mod', row.names = 1, header=T)
summary(phenotype)

dge <- DGEList(counts = counts)
dge <- dge[rowSums(cpm(dge$counts) > 1) > 20,,keep.lib.size = F]
dge <- calcNormFactors(dge)

design_Phenotype <- model.matrix( ~as.factor(Phenotype), data = phenotype)

v_full <- voom( dge, design_Phenotype, plot = T)
fit_full <- eBayes(lmFit( v_full, design_Phenotype   ))

corfit <- duplicateCorrelation(v_full,design_Phenotype, block = phenotype$PatientID)

fit_full_block_1 <- eBayes(lmFit( v_full, design_Phenotype, block = phenotype$PatientID, correlation = corfit$consensus.correlation))

summary(decideTests(fit_full_block_1))

#check effect of voom itself using blocking
v_full_block <- voom( dge, design_Phenotype, block = phenotype$PatientID, correlation = corfit$consensus.correlation)
fit_full_block <- eBayes(lmFit( v_full_block, design_Phenotype, block = phenotype$PatientID, correlation = corfit$consensus.correlation   ))

cor(topTable(fit_full_block, num = nrow(counts), sort.by = "none", coef = 4)$adj.P.Val, topTable(fit_full, num = nrow(counts), sort.by = "none", coef = 4)$adj.P.Val  )

summary(decideTests(fit_full))


