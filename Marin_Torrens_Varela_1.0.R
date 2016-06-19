library(SummarizedExperiment)
library(edgeR)
library(geneplotter)
library(limma)
sc<-readRDS("seLIHC.rds")
dge<-DGEList(counts=assays(sc)$counts,genes=mcols(sc))
assays(sc)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5) #transform data into "Counts per Million reads"

ord <- order(dge$sample$lib.size/1e6)
barplot(dge$sample$lib.size[ord]/1e6, las=1, ylab="Millions of reads",
        xlab="Samples", col=c("blue", "red")[(sc$type[ord] == "tumor") + 1])
legend("topleft", c("tumor", "normal"), fill=c("red", "blue"), inset=0.01)


#select only samples with more than 40 million reads
libmask<-dge$samples$lib.size>40000000 
sc<-sc[,libmask]
dge<-dge[,libmask]

ord <- order(dge$sample$lib.size/1e6)
barplot(dge$sample$lib.size[ord]/1e6, las=1, ylab="Millions of reads",
        xlab="Samples", col=c("blue", "red")[(sc$type[ord] == "tumor") + 1])
legend("topleft", c("tumor", "normal"), fill=c("red", "blue"), inset=0.01)

#select only samples that are paired
masknames<-rep(FALSE,dim(sc)[2])
count=1
for  (i in sc$bcr_patient_barcode){
  count2=0
  for (j in sc$bcr_patient_barcode) {
    if (identical(i,j) && !is.na(i)){
      count2=count2+1
      if (count2==2){ #each paired sample appears twice, we should only count it once
        masknames[count]=TRUE
      }
    }
  }
  count=count+1
}


#filter by paired data
sc<-sc[,masknames]
dge<-dge[,masknames]


par(mfrow=c(1, 2))
multidensity(as.list(as.data.frame(assays(sc[, sc$type == "tumor"])$logCPM)),
             xlab="log 2 CPM", legend=NULL, main="Tumor samples", las=1)
multidensity(as.list(as.data.frame(assays(sc[, sc$type == "normal"])$logCPM)),
             xlab="log 2 CPM", legend=NULL, main="Normal samples", las=1)

#plot average expresion for each gene
avgexp <- rowMeans(assays(sc)$logCPM)
hist(avgexp, xlab="log2 CPM", main="", las=1)
abline(v=1, col="red", lwd=2)


#eliminate genes with expresion lower than 1.
mask <- avgexp > 1
sc <- sc[mask, ]
dge <- dge[mask, ]


dge <- calcNormFactors(dge)

#MA-plots

assays(sc)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5)
par(mfrow=c(6, 5), mar=c(1, 1, 1, 1))
sctmp <- sc[, sc$type == "tumor"]
dgetmp <- dge[, sc$type == "tumor"]
for (i in 1:ncol(sctmp)) {
  A <- rowMeans(assays(sctmp)$logCPM)
  M <- assays(sctmp)$logCPM[, i] - A
  samplename <- as.character(sctmp$bcr_patient_barcode[i])
  print (samplename)
  smoothScatter(A, M, main=samplename, las=1)
  abline(h=0, col="blue", lwd=2)
  lo <- lowess(M ~ A)
  lines(lo$x, lo$y, col="red", lwd=2)
}
#Sample with name "TCGA-ES-A2HT" seems aberrant according to the MA PLOT, so we elimate it.
namemask2<-sc$bcr_patient_barcode!="TCGA-ES-A2HT"
sc<-sc[,namemask2]
dge<-dge[,namemask2]

#Repeat the MA-plot to check if the sample is really out.
par(mfrow=c(6, 5), mar=c(1, 1, 1, 1))
sctmp <- sc[, sc$type == "normal"]
dgetmp <- dge[, sc$type == "normal"]
for (i in 1:ncol(sctmp)) {
  A <- rowMeans(assays(sctmp)$logCPM)
  M <- assays(sctmp)$logCPM[, i] - A
  samplename <- as.character(sctmp$bcr_patient_barcode[i])
  print (samplename)
  smoothScatter(A, M, main=samplename, las=1)
  abline(h=0, col="blue", lwd=2)
  lo <- lowess(M ~ A)
  lines(lo$x, lo$y, col="red", lwd=2)
}

#Studying posible batch effects
tss <- substr(colnames(sc), 6, 7)
table(tss)
center <- substr(colnames(sc), 27, 28)
table(center)
plate <- substr(colnames(sc), 22, 25)
table(plate)
portionanalyte <- substr(colnames(sc), 18, 20)
table(portionanalyte)
samplevial <- substr(colnames(sc), 14, 16)
table(samplevial)

#is our outcome of interest (tumor or normal) well distributed according to the tss?
table(data.frame(TYPE=sc$type, TSS=tss))

#Note: As the data is paired, a desequilibrium is not possible


#Build a hierarchical clustering of the samples to identify possible batches.
par(mfrow=c(1,1))
logCPM <- cpm(dge, log=TRUE, prior.count=3)
d <- as.dist(1-cor(logCPM, method="spearman"))
sampleClustering <- hclust(d)
batch <- as.integer(factor(tss))
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
names(batch) <- colnames(sc)
outcome <- paste(substr(colnames(sc), 9, 12), as.character(sc$type), sep="-")
names(outcome) <- colnames(sc)

sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col=as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)
plot(sampleDendrogram, main="Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch)), levels(factor(tss))), fill=sort(unique(batch)))


#MDS will help to identify batches
plotMDS(dge, labels=outcome, col=batch)
legend("bottomleft", paste("Batch", sort(unique(batch)), levels(factor(tss))),
       fill=sort(unique(batch)), inset=0.05)

#Sample with name "TCGA-BC-A10X" seems aberrant according to the Hierarchical clustering, so we elimate it.
namemask2<-sc$bcr_patient_barcode!="TCGA-BC-A10X"
sc<-sc[,namemask2]
dge<-dge[,namemask2]

#first analysis without taking any batch into account
library(sva)
mod <- model.matrix(~ sc$type, colData(sc))
mod0 <- model.matrix(~1, colData(sc))
pv <- f.pvalue(assays(sc)$logCPM, mod, mod0)

par(mfrow=c(1,2),mar=c(5,4,4,5))

hist(pv, las=1, main="p-value distribution without sva", xlab="p-value",xlim=c(0,1), las=1,ylim = c(0,7000)) #distribution of pvalues
sv <- sva(assays(sc)$logCPM, mod, mod0) #analysis of surrogate variables.
#do again the analysis taking the founded surrogates into acount
modsv <- cbind(mod, sv$sv)
mod0sv <- cbind(mod0, sv$sv)
pvsv <- f.pvalue(assays(sc)$logCPM, modsv, mod0sv)
hist(pvsv, las=1,main="p-value distribution with sva",xlab="p-value",xlim=c(0,1), las=1,ylim = c(0,7000))
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
## new part 2

# Fold changes
FDRcutoff <- 0.05
tumorExp <- rowMeans(logCPM[, sc$type == "tumor"])
normalExp <- rowMeans(logCPM[, sc$type == "normal"])
par(mfrow=c(1,2), mar=c(4, 4, 4, 4))
plot(normalExp, tumorExp, xlab = "Normal", ylab = "Tumor", pch = ".")
lines(c(1:15),c(1:15), col="red")
plot((normalExp + tumorExp)/2, tumorExp - normalExp, xlab = "Gene mean expression", ylab = "Fold change", pch = ".")
lines(c(1:15),rep(0,15), col="red")

# top 10 DE
log2fc <- tumorExp - normalExp #as tumorExp and normalExp are in log scale, their substraction is equal to compute the foldchange.
ranking <- sort(abs(log2fc), decreasing = TRUE, index.return = TRUE)
#ranking$x has the list of differences ordered by values, ranking$ix only has the indexes ordered according to the difference values.
#head(data.frame(Log2FC = log2fc[ranking$x], FC = 2^log2fc[ranking$x]), n = 10)
sc_fd <- sc[head(ranking$ix,n=10),] 
dge_fd <- dge[head(ranking$ix,n=10),] 

# Adjust for unknown covariates
sc$bcr_patient_barcode<-droplevels(sc$bcr_patient_barcode)
par(mfrow=c(1,1))
mod_sv2<-model.matrix(~sc$type+sc$bcr_patient_barcode,data=colData(sc)) # paired data
mod0_sv2 <- model.matrix(~sc$bcr_patient_barcode, colData(sc))
v<-voom(dge,mod_sv2,plot=TRUE)
sv2 <- sva(v$E, mod = mod_sv2, mod0 = mod0_sv2)
designsv <- cbind(mod_sv2, sv2$sv)
colnames(designsv) <- c(colnames(designsv)[1:29], paste0("SV", 1:sv2$n))
fit4<-lmFit(v,designsv)
fit4<-eBayes(fit4)
res4 <- decideTests(fit4, p.value = FDRcutoff)
#summary(res4)
tt4 <- topTable(fit4, coef = 2, n = Inf)
DEgenes <- rownames(tt4)[tt4$adj.P.Val < FDRcutoff]
length(DEgenes) #number of DE genes.

par(mfrow = c(1, 2), mar = c(4, 5, 2, 2))
hist(tt4$P.Value, xlab = "P-values", main = "", las = 1)
qqt(fit4$t[, 2], df = fit4$df.prior + fit4$df.residual, main = "", pch = ".", cex = 3)
abline(0,1,lwd=2,col="red")


#boxplot of best p-values:
par(mfrow=c(2,5))
for (i in 1:10){
  boxplot(cpm(dge$counts[rownames(tt4)[i],], log=TRUE, prior.count=0.5) ~ sc$type,main=tt4$symbol[i],ylim=c(6,18),las=2, col=c('lightskyblue3','lightsalmon'))
}

#volcano plot
par(mfrow = c(1, 1))

plot(tt4$logFC, -log10(tt4$P.Value), xlab = "Log2 fold-change", ylab = "-log10 P-value", 
     pch = ".", cex = 5, col = grey(0.75), cex.axis = 1.2, cex.lab = 1.5, las = 1)
points(tt4[tt4$adj.P.Val < 0.05, "logFC"], -log10(tt4[tt4$adj.P.Val < 0.05, "P.Value"]), pch = ".", 
       cex = 5, col = "red")
abline(h = -log10(max(tt4[tt4$adj.P.Val < 0.05, "P.Value"])), col = grey(0.5), lty = 2)


# Diagnostic plots (MA-plot)
par(mfrow=c(1, 1), mar=c(2, 2, 2, 2))
top7 <- order(fit4$lods[, 2], decreasing = TRUE)[1:7]
limma::plotMA(fit4, coef = 2, status = rownames(fit4$lods) %in% DEgenes, legend = FALSE,
              main = "Model 2", hl.pch = 46, hl.cex = 4, bg.pch = 46, bg.cex = 3, las = 1) 
text(fit4$Amean[top7], fit4$coef[top7, 2], fit4$genes$symbol[top7], cex = 0.5, pos = 4)

# Factorial design
newfac <- factor(paste(sc$gender, sc$type, sep = "."))
head(newfac)
design <- model.matrix(~0 + newfac, colData(sc))
head(design, n = 3)
fit6 <- lmFit(v, design)
cont.matrix <- makeContrasts(men = newfacMALE.normal - newfacMALE.tumor, women =newfacFEMALE.normal - newfacFEMALE.tumor , levels = design)
fit6 <- contrasts.fit(fit6, cont.matrix)
fit6 <- eBayes(fit6)
ttconmales <- topTable(fit6, coef = "men", n = Inf)
ttconfemales <- topTable(fit6, coef = "women", n = Inf)
res <- decideTests(fit6, p.value = 0.05)
vennDiagram(res,circle.col=c("red","blue"))


##### Functional annotations
#Create a data.frame object with gene metadata (chromosome and symbol)
library(GOstats)
geneUniverse=rownames(sc)
params <- new("GOHyperGParams", geneIds=DEgenes, universeGeneIds=geneUniverse,
              annotation="org.Hs.eg.db", ontology="BP",
              pvalueCutoff=0.05, testDirection="over")
hgOver <- hyperGTest(params)
htmlReport(hgOver, file = "gotests.html")
#browseURL("gotests.html")

# Conditional test
conditional(params) <- TRUE
hgOverCond <- hyperGTest(params)
htmlReport(hgOver, file = "gotests_cond.html")
goresults <- summary(hgOverCond)


#filter the previous results by a minimum value on the Count and Size 
goresults <- goresults[goresults$Size >= 10 & goresults$Count >= 5, ]
goresults <- goresults[order(goresults$OddsRatio, decreasing = TRUE), ]
head(goresults)


# extract the genes that enrich each GO term and paste it to the result
geneIDs <- geneIdsByCategory(hgOverCond)[goresults$GOBPID]
geneSYMs <- sapply(geneIDs, function(id) select(org.Hs.eg.db, columns = "SYMBOL", key = id, 
                                                keytype = "ENTREZID")$SYMBOL)
geneSYMs <- sapply(geneSYMs, paste, collapse = ", ")
goresults <- cbind(goresults, Genes = geneSYMs)
rownames(goresults) <- 1:nrow(goresults)

library(xtable)
xtab <- xtable(goresults, align = "l|c|r|r|r|r|r|p{3cm}|p{3cm}|")
print(xtab, file = "goresults.html", type = "html")

###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

#########       GENE SETS ANALYSIS       ##########################################################################

library(GSEABase)
library(GSVAdata)
data(c2BroadSets)
#select only three gene collections: KEGG,REACTOME and BIOCARTA
c2BroadSets <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
                             grep("^REACTOME", names(c2BroadSets)), grep("^BIOCARTA", names(c2BroadSets)))]
gsc <- GeneSetCollection(c(c2BroadSets))
gsc <- mapIdentifiers(gsc, AnnoOrEntrezIdentifier(metadata(sc)$annotation)) # This is not necessary becase our genes are already in Entrex IDs

Im <- incidence(gsc) #Build a matrix indicating which genes are expressed in each gene set.
Im <- Im[, colnames(Im) %in% rownames(sc)] #discard genes that are not in our study group

#discard genes from our original study set that are not included in the c2BroadSets library
#reducing the number of genes diminishes multiple test correction problems.

sctemp <- sc[colnames(Im), ] 
dgetemp <- dge[colnames(Im), ]

# We start the simple GSEA analysis by doing classical paired DE without calling any gene DE:
#We do this test again because now we have less genes -> improves multiple testing correction

design<-model.matrix(~sctemp$type+sctemp$bcr_patient_barcode,data=colData(sctemp)) # These step are the same we did before
v <- voom(dge, design, plot = FALSE)
mod0 <- model.matrix(~sctemp$bcr_patient_barcode, data=colData(sctemp))
sv <- sva(v$E, mod = design, mod0 = mod0)   
design <- cbind(design, sv$sv)
fit<-lmFit(v,design)
fit<-eBayes(fit)
tt <- topTable(fit, coef = 2, n = Inf)

# qq-plot
qq <- qqnorm(tt$t)
abline(0, 1,col="red")
chroutliers <- tt$chr[abs(tt$t) > 10]
text(qq$x[abs(qq$y) > 10], qq$y[abs(qq$y) > 10], chroutliers, pos = 4)

# remove gene sets with less than 5 genes and more than 200
Im <- Im[rowSums(Im) >= 5 && rowSums(Im) <= 200 , ]

#calculate the Z statistic
#Get the t-statistic value of the genes in the incidence matrix
tGSgenes <- tt[match(colnames(Im), rownames(tt)), "t"] #tt[corresponding positions of (the genes in the incidence matrix) in tt table,select the "t" column]

zS <- sqrt(rowSums(Im)) * (as.vector(Im %*% tGSgenes)/rowSums(Im))
#zS=media de t-statistico para o gen set de interes*numero de genes no gene set
#rowSums(Im) number of genes present in one gene set
#as.vector(Im %*% tGSgenes) gets only the t-staticstic of the genes present in the gene set of interes.(Im has 0 or 1, so if the gene is not present, its t-statisc gets * by 0)
#/rowSums(Im) we divide by the number of genes to obtain the mean of t-statitcs for that gene set.

# qq plot of gene set Z-scores (gives a quick overview of how many promising gene sets could be DE)
qqnorm(zS)
abline(-1.5, 0,col="red")
chroutliers <- tt$chr[abs(tt$t) > 10]
text(qq$x[abs(qq$y) > 10], qq$y[abs(qq$y) > 10], chroutliers, pos = 4)

# first few gene sets with largest Z-score:
rnkGS <- sort(abs(zS), decreasing = TRUE)

# scatter plots for the sets with larger Z scores
plotGS <- function(se, gs, pheno, ...) {
  l <- levels(colData(se)[, pheno])
  idxSamples1 <- colData(se)[, pheno] == l[1]
  idxSamples2 <- colData(se)[, pheno] == l[2]
  exps1 <- rowMeans(assays(se)$logCPM[gs, idxSamples1])
  exps2 <- rowMeans(assays(se)$logCPM[gs, idxSamples2])
  rng <- range(c(exps1, exps2))
  plot(exps1, exps2, pch = 21, col = "black", bg = "black", xlim = rng, ylim = rng, 
       xlab = "logCPM(normal)", ylab = "logCPM(tumor)", ...)
  abline(a = 0, b = 1, lwd = 2, col = "red")
}

#pick the top 2 best ranked sets, and get gene names present in them (that is: names(rnkGS)[1] and names(rnkGS)[2]):
#from colnames(Im) (which is a list of gene identifiers) pick only the ones which are present in the number one gene set (names(rnkGS)[1]) 
genesGS1 <- colnames(Im)[which(Im[names(rnkGS)[1], ] == 1)] 
genesGS2 <- colnames(Im)[which(Im[names(rnkGS)[2], ] == 1)]

#plot the gene expression of the genes in the best ranked gene set spliting by type
par(mfrow = c(1, 2), mar = c(4, 5, 3, 4))
plotGS(sctemp, genesGS1, "type", main = names(ful)[1], cex.lab = 2, las = 1)
plotGS(sctemp, genesGS2, "type", main = names(rnkGS)[2], cex.lab = 2, las = 1)

# z-test
pv <- pmin(pnorm(zS), 1 - pnorm(zS)) 
#pnorm(x) returns the probability that a normally distributed RANDOM number is smaller than x!
#notice that once you have passed the mean (by default 0) the probability of finding a randonly generated number bigger than the 
#mean keeps increasing despite that you are getting far apart from the mean.
#that is why we need to pick the minimun between pnorm(zS) and 1 - pnorm(zS).

sum(pv < 0.05)
pvadj <- p.adjust(pv, method = "fdr")
pvadj<-sort(pvadj)
DEgs <- names(pvadj)[which(pvadj < 0.05)] #get the gene set names that have a pvalue < 0.05


# Chi square test, main difference is that now we do not lose those gene set that have a equal number of
#of upregulated and downregulated genes (before: we calculated the mean of t-statistic which can result in zero. Now: sum((x - mean(x))^2) )
library(Category)
Im2 <- Im
Im2 <- Im2[rowSums(Im2) >= 20, ]
tGSgenes2 <- tt[match(colnames(Im2), rownames(tt)), "t"] 
xS <- applyByCategory(tGSgenes2, Im2, function(x) (sum((x - mean(x))^2) - (length(x) - 1))/(2 * (length(x) - 1)))
rnkGS <- sort(abs(xS), decreasing = TRUE) #results from Chi-Square test
pv <- pmin(pnorm(xS), 1 - pnorm(xS))
pvadj2 <- p.adjust(pv)
pvadj2<-sort(pvadj2)
DEgsByScale <- names(pvadj2)[which(pvadj2 < 0.05)]

#get genes from best ranked gene sets according to Chi-Square punctuation.
topchi1genes <- colnames(Im2)[which(Im2[names(pvadj2)[1], ] == 1)]
topchi2genes <- colnames(Im2)[which(Im2[names(pvadj2)[2], ] == 1)]
par(mfrow = c(1, 2))
plotGS(sc, topchi1genes, "type", main = names(pvadj2[1]), cex.lab = 2, las = 1)
plotGS(sc, topchi2genes, "type", main = names(pvadj2[2]), cex.lab = 2, las = 1)

# union z-score and Chi square with smallest p-values
full_pv <- pvadj2
for (i in seq(length(pvadj))) {
  if ((names(pvadj[i]) %in% names(pvadj2))==FALSE) {
    full_pv <- c(full_pv,pvadj[i])
  } else if (((names(pvadj[i]) %in% names(pvadj2))==TRUE) && (pvadj2[[names(pvadj[i])]]<pvadj[[i]])) {
    full_pv[(names(pvadj[i]))] <- pvadj[i]
  }
}

full_pv<-sort(full_pv)

#get genes from best ranked gene sets according to both tests union.
topgs1genes <- colnames(Im)[which(Im[names(full_pv)[1], ] == 1)]
topgs2genes <- colnames(Im)[which(Im[names(full_pv)[2], ] == 1)]

plotGS(sc, topgs1genes, "type", main = names(full_pv[1]), cex.lab = 2, las = 1)
plotGS(sc, topgs2genes, "type", main = names(full_pv)[2], cex.lab = 2, las = 1)


# Overlaps
library(GSVA)
gsov <- computeGeneSetsOverlap(gsc[fullDEgs], rownames(sctemp)) # put to 0 genes not present in dataset from sets
trimask <- upper.tri(gsov)
rnkOv <- data.frame(gs1 = row(gsov)[trimask], gs2 = col(gsov)[trimask], ov = gsov[trimask])
rnkOv <- rnkOv[order(rnkOv$ov, decreasing = TRUE), ] # order values (decreasing)
rnkOv$gs1 <- rownames(gsov)[rnkOv$gs1]
rnkOv$gs2 <- rownames(gsov)[rnkOv$gs2]
sum(rnkOv$ov == 1)  ## how many pairs of gene sets are identical?
sum(rnkOv$ov < 0.2)  ## how many pairs of gene sets share less than 20% of the genes?


#GSVA analysis #########################################################3

library(GSVA)
#GSexpr <- gsva(assays(sc)$counts, gsc, rnaseq = TRUE, min.sz = 10, max.sz = 300, verbose = FALSE,method="ssgsea")$es.obs
GSexpr <- gsva(assays(sc)$counts, gsc, rnaseq = TRUE, min.sz = 25, max.sz = 300, verbose = FALSE,kernel=FALSE)
dim(GSexpr) # gene sets x enrichment score (ES)


mod<-model.matrix(~sc$type+sc$bcr_patient_barcode,data=colData(sc))
mod0 <- model.matrix(~sc$bcr_patient_barcode, colData(sc))
svaobj <- sva(GSexpr,mod= mod,mod0= mod0) # Not working
modSVs <- cbind(mod, svaobj$sv)
fit <- lmFit(GSexpr, modSVs)
fit <- eBayes(fit)
ttgs <- topTable(fit, coef = 2, n = Inf)
DEgs2 <- rownames(ttgs[ttgs$adj.P.Val < 0.05, , drop = FALSE])

par(mfrow=c(2,5))
for (i in 1:10){
  boxplot(GSexpr[DEgs2[i], ] ~ sc$type,main=DEgs2[i],las=2, cex.main=0.6, col=c('lightskyblue3','lightsalmon'))
}