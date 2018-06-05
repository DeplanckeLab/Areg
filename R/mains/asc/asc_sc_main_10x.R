###### Initialise variables and load functionality #####
source("plot-util.R")
source("util-asc.R")

#### Initialise path #####
mypath=inpath="/10x/"
myname="mASC10x-allData-NoXistKrt-review-FF"

#### remove all cells with total nr. of genes under 1000
genescut=1000
#### load libraries
library(SC3); library(cellrangerRkit); library(beanplot); library(scater)
library(M3Drop); library(scran); library(org.Mm.eg.db)

#### set of genes of interest; for plotting
goi=unique(c(c(allmarkers),mychoiceMarkerAll))


#### Output plots into a .pdf
pdf(paste(mypath,myname,".",genescut,"genecut.scran.final.Pub.pdf",sep=""),width=7,height=7)

#### load the CellRanger results
cellranger_pipestance_path="/10x/mASC1-lsf/"
gbm = load_cellranger_matrix(cellranger_pipestance_path)
countmatrix=as.matrix(exprs(gbm))
rownames(countmatrix)=fData(gbm)[,1];colnames(countmatrix)=pData(gbm)[,1]

#### Select Xist positive and Krt18/19 negative cells  ####
xisti=which(countmatrix[fData(gbm)[which(fData(gbm)[,2]%in%"Xist"),1],]>1)
krt=intersect(union(which(countmatrix[fData(gbm)[which(fData(gbm)[,2]%in%"Krt18"),1],]<1),
	which(countmatrix[fData(gbm)[which(fData(gbm)[,2]%in%"Krt19"),1],]<1)),which(countmatrix[fData(gbm)[which(fData(gbm)[,2]%in%"Epcam"),1],]<1))
selected=intersect(xisti, krt) #1896 cells no Xist; 1869 cells no Krt/Epcam + Xist 
countmatrix=countmatrix[,selected]

#### Get some info about the data: nr. of genes/cell, nr. UMI/cell & filter
genesPerCell=apply(countmatrix,2,function(x) length(which(x>0)))
okmatrix=countmatrix[apply(countmatrix,1,function(x) sum(x))>=10, which(genesPerCell>=2)]
umiPerCell=apply(countmatrix,2,function(x) sum(x))
genesPerCell=apply(countmatrix,2,function(x) length(which(x>0)))
# hard filtering if necessary?
# okmatrix=countmatrix[,which(genesPerCell>=genescut)]; okmatrix=countmatrix[,which(genesPerCell<=10000)]
okmatrix=okmatrix[apply(okmatrix,1,function(x) sum(x))>=10,]; okmatrix=okmatrix[apply(okmatrix,1,function(x) length(which(x>0)))>=2,]
genesPerCell=apply(okmatrix,2,function(x) length(which(x>0)))
dim(countmatrix); dim(okmatrix) #countmatrix 27998  1869 # okmatrix 13597  1869

#### Process cell/gene expression matrix & filter using scran
sce <- newSCESet(countData=okmatrix); sce <- calculateQCMetrics(sce)
hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="", breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="", breaks=20, col="grey80", ylab="Number of cells")
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE) #6
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE) #38
sce <- sce[,!(libsize.drop | feature.drop )]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), Remaining=ncol(sce))
ave.counts <- rowMeans(counts(sce)); # keep <- ave.counts >= 1 #additional filtering?
hist(log10(ave.counts), breaks=100, main="", col="grey80", xlab=expression(Log[10]~"average count"))
abline(v=log10(1), col="blue", lwd=2, lty=2)
plotQC(sce, type = "highest-expression", n=50) 

#### estimate normalisation factors 
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusters)
plot(sizeFactors(sce), sce$total_counts/1e6, log="xy", ylab="Library size (millions)", xlab="Size factor")
#1896 Xist cells, 1858 after filtering,13597 features 1831 cells; G1 cells #1831   23    4

#### predict cell cycle using cyclone and only keep G1 predicted cells 
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assigned <- cyclone(sce, pairs=mm.pairs)
table(assigned$phases) #1804 (G1)   23 (G2M)   4 (S) final
cellcycle=assigned$phases; names(cellcycle)=colnames(sce)
sce <- sce[,assigned$phases=="G1"]
temp=convertTo(sce, type="monocle"); temp.saved=temp
write.table(exprs(temp), paste(inpath,myname,".bam-DGE.4000.fil.",genescut,"genes10DGE2Cells.scranNorm.G1only.arrayExpress.txt",sep=""),
		row.names=T, col.names=T,quote=F,sep="\t")

##### plot some info for the final retained cells and gene values
par(mfrow=c(2,2))
hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="", breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="", breaks=20, col="grey80", ylab="Number of cells")
libsize=cbind(sce$total_counts/1e6,sce$total_features);rownames(libsize)=sampleNames(temp)
colnames(libsize)=c("Libsize","Genenr")

#### output raw for Figure
write.table(libsize,paste(path,"FigS2ab-raw.10xlibsize.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

#### obtain DE genes with M3Drop
scrancounts=convertTo(sce, type="edgeR")
fits <- M3DropDropoutModels(exprs(temp))
data.frame(MM=fits$MMFit$SAr, Logistic=fits$LogiFit$SAr,DoubleExpo=fits$ExpoFit$SAr) #Log best fit but MM close fit
DE_genes_FDR <- M3DropDifferentialExpression(exprs(temp),mt_method="fdr", mt_threshold=0.05)
DE_genes <- M3DropDifferentialExpression(exprs(temp),mt_method="fdr", mt_threshold=0.0001) 
large_DE_genes <- M3DropDifferentialExpression(exprs(temp),mt_method="fdr", mt_threshold=1)
topDE_genes <- M3DropDifferentialExpression(exprs(temp),mt_method="fdr", mt_threshold=0.0000001)
write.table(DE_genes, paste(inpath,myname,".bam-DGE.4000.fil.",genescut,"genes10DGE2Cells.scranNorm.DEgenes.txt",sep=""),
	row.names=T,col.names=T,quote=F,sep="\t") #9
write.table(topDE_genes, paste(inpath,myname,".bam-DGE.4000.fil.",genescut,"genes10DGE2Cells.scranNorm.topDEgenes.txt",sep=""),
	row.names=T,col.names=T,quote=F,sep="\t")

#### obtain HVG genes 
var.fit <- trendVar(sce, trend="loess", span=0.4,use.spikes=F)
var.out <- decomposeVar(sce, var.fit)
top.hvgs <- order(var.out$bio, decreasing=TRUE)
decomp=var.out; fit=var.fit
plot(decomp$mean, decomp$total, xlab="Mean log-expression", ylab="Variance")
o <- order(decomp$mean); lines(decomp$mean[o], decomp$tech[o], col="red", lwd=2); points(fit$mean, fit$var, col="red", pch=16)
null.dist <- correlateNull(ncol(sce))
hvg.out <- var.out[which(var.out$FDR <= 0.1 & var.out$bio >= 0.3),] 
hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] #### 223 genes #### 1942 genes
write.table(hvg.out, paste(inpath,myname,".bam-DGE.4000.fil.",genescut,"genes10DGE2Cells.scranNorm.hvg.txt",sep=""),
	row.names=T,col.names=T,quote=F,sep="\t") #223 #209
	
#### high variability and top DE from M3Drop for tSNE ####
mymatrix.log=log1p(exprs(temp)[unique(c(rownames(hvg.out),rownames(DE_genes))),]) ### 631 genes 550 M3Drop genes 223 HVG genes #1804 cells
write.table(mymatrix.log, paste(inpath,myname,".bam-DGE.4000.fil.",genescut,"genes10DGE2Cells.scranNorm.mymatrixlog.txt",sep=""),
	row.names=T,col.names=T,quote=F,sep="\t") #642 1831 #631 1804

#### Use SC3 for clustering & marker gene analysis	
mysceset = newSCESet(countData = exprs(temp)[which(rownames(exprs(temp))%in%rownames(mymatrix.log)),]); 
mysceset.est <- sc3_estimate_k(mysceset)
myk=mysceset.est@sc3$k_estimation #13 clusters 
mysceset <- calculateQCMetrics(mysceset); mysceset <- sc3(mysceset, ks = c(2:13),n_cores=1)
scclusters.2=pData(mysceset)[,"sc3_2_clusters"]; scclusters.3=pData(mysceset)[,"sc3_3_clusters"]; scclusters.4=pData(mysceset)[,"sc3_4_clusters"]
scclusters.5=pData(mysceset)[,"sc3_5_clusters"]; scclusters.6=pData(mysceset)[,"sc3_6_clusters"]; scclusters.10=pData(mysceset)[,"sc3_10_clusters"]; 
scclusters.7=pData(mysceset)[,"sc3_7_clusters"]; scclusters.8=pData(mysceset)[,"sc3_8_clusters"]; scclusters.9=pData(mysceset)[,"sc3_9_clusters"]
scclusters.11=pData(mysceset)[,"sc3_11_clusters"]; scclusters.12=pData(mysceset)[,"sc3_12_clusters"]; scclusters.13=pData(mysceset)[,"sc3_13_clusters"]
names(scclusters.3)=names(scclusters.4)=names(scclusters.5)=names(scclusters.6)=rownames(pData(mysceset)) 
names(scclusters.7)=names(scclusters.8)=names(scclusters.9)=names(scclusters.10)=rownames(pData(mysceset))
names(scclusters.11)=names(scclusters.12)=names(scclusters.13)=rownames(pData(mysceset))
#myk=mysceset.est@sc3$k_estimation; if(myk>10) {myk=10}; if(myk<2) {myk=2} #### 

#set nr. clusters to 4 & plot cluster stability and silhouette for distinct nr. of k
myk=4
sc3_plot_cluster_stability(mysceset, k = 2); sc3_plot_cluster_stability(mysceset, k = 3); sc3_plot_cluster_stability(mysceset, k = 4); 
sc3_plot_cluster_stability(mysceset, k = 5); sc3_plot_cluster_stability(mysceset, k = 6); sc3_plot_cluster_stability(mysceset, k = 7)
sc3_plot_cluster_stability(mysceset, k = 8);sc3_plot_cluster_stability(mysceset, k = 9); sc3_plot_cluster_stability(mysceset, k = 10)
sc3_plot_cluster_stability(mysceset, k = 11);sc3_plot_cluster_stability(mysceset, k = 12); sc3_plot_cluster_stability(mysceset, k = myk)
sc3_plot_silhouette(mysceset, k = 2); sc3_plot_silhouette(mysceset, k = 3); sc3_plot_silhouette(mysceset, k = 4); sc3_plot_silhouette(mysceset, k = 5); 
sc3_plot_silhouette(mysceset, k = 6); sc3_plot_silhouette(mysceset, k = 7); 
sc3_plot_silhouette(mysceset, k = 8);sc3_plot_silhouette(mysceset, k = 9); sc3_plot_silhouette(mysceset, k = 10)
sc3_plot_silhouette(mysceset, k = 11);sc3_plot_silhouette(mysceset, k = 12); sc3_plot_silhouette(mysceset, k = 13)
scclusters.myk=pData(mysceset)[,"sc3_2_clusters"]
if (myk==3){ scclusters.myk=pData(mysceset)[,"sc3_3_clusters"]}  else { 
	if (myk==4){ scclusters.myk=pData(mysceset)[,"sc3_4_clusters"]} else{
		if (myk==5){ scclusters.myk=pData(mysceset)[,"sc3_5_clusters"]} else{
			if (myk==6){ scclusters.myk=pData(mysceset)[,"sc3_6_clusters"]} else{
				if (myk==7){ scclusters.myk=pData(mysceset)[,"sc3_7_clusters"]}  else{
						if (myk==8){ scclusters.myk=pData(mysceset)[,"sc3_8_clusters"]} else{
						if (myk==9){ scclusters.myk=pData(mysceset)[,"sc3_9_clusters"]} else{
							if (myk==10){ scclusters.myk=pData(mysceset)[,"sc3_10_clusters"]}}}}}}}}


# Plots for distinct nr. of clusters
mydata=scale(mymatrix.log); wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares, scaled")
#### output raw for Figure
names(wss)=c(1:15)
write.table(wss,paste(path,"FigS8k-raw.ClusterNr.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
#### Save data for future reference 												
save(scclusters.myk, myk,mysceset.est,mysceset, file=paste(mypath,myname,".",genescut,"genecut.SC3.scran.2.RData",sep=""))
##### save(scclusters.myk, myk,mysceset.est,mysceset, file=paste(mypath,myname,".",genescut,"genecut.SC3.scran.fullData.RData",sep="")) #### this data contains the full 2,919 cells

### color-code according to different cluster nrs. 
filcols=c(brewer.pal(9,"Pastel1"),brewer.pal(8,"Pastel2"))
marker_genes_sc_myk <- M3DropGetMarkers(exprs(temp), scclusters.myk)
sc3.myk=filcols[scclusters.myk]; names(sc3.myk)=rownames(pData(mysceset))
sc3.3=filcols[scclusters.3]; names(sc3.3)=rownames(pData(mysceset)); sc3.4=filcols[scclusters.4]; names(sc3.4)=rownames(pData(mysceset))
sc3.5=filcols[scclusters.5]; names(sc3.5)=rownames(pData(mysceset)); sc3.6=filcols[scclusters.6]; names(sc3.6)=rownames(pData(mysceset)); 
sc3.7=filcols[scclusters.7]; names(sc3.7)=rownames(pData(mysceset)); sc3.8=filcols[scclusters.8]; names(sc3.8)=rownames(pData(mysceset))
sc3.9=filcols[scclusters.9]; names(sc3.10)=rownames(pData(mysceset)); sc3.10=filcols[scclusters.10]; names(sc3.10)=rownames(pData(mysceset))
sc3.11=filcols[scclusters.11]; names(sc3.11)=rownames(pData(mysceset)); sc3.12=filcols[scclusters.12]; names(sc3.12)=rownames(pData(mysceset))
sc3.13=filcols[scclusters.13]; names(sc3.13)=rownames(pData(mysceset))

#### get marker genes for k=4
P.top.sc3.20=list(NA)
P.top.sc3.100=list(NA)
P.top.sc3.200=list(NA)
P.top.sc3.500=list(NA)
for(k in 1:myk) {
	P.top.sc3.20[[k]]=head(marker_genes_sc_myk[marker_genes_sc_myk$Group==k,],20)
	P.top.sc3.100[[k]]=head(marker_genes_sc_myk[marker_genes_sc_myk$Group==k,],100)
	P.top.sc3.200[[k]]=head(marker_genes_sc_myk[marker_genes_sc_myk$Group==k,],200)
	P.top.sc3.500[[k]]=head(marker_genes_sc_myk[marker_genes_sc_myk$Group==k,],500)
}

### t-SNE of the cells, 2D representation 
rtsne_out.5 <- Rtsne(as.matrix(t(mymatrix.log)),theta=0.0,perplexity=20)
##rtsne_out.5 <- Rtsne(as.matrix(t(mymatrix.log)),theta=0.0000001,perplexity=10)

rtsne_out=rtsne_out.5
goi=unique(c(c(allmarkers),mychoiceMarkerAll,elementMetadata(spikedgenes84.GR[rownames(hvg.out)])$geneName,
elementMetadata(genes84.GR[rownames(topDE_genes)])$geneName))
#### Save the t-SNE for future reference
#save(sce,rtsne_out.5,topDE_genes,DE_genes,hvg.out,file=paste(mypath,myname,".",genescut,"genecut.hvgAndTopM3Drop.touse.2.Rtsne",sep=""))
load(paste(mypath,myname,".",genescut,"genecut.hvgAndTopM3Drop.touse.Rtsne",sep=""))

#### plot the t-SNE representation
mymatrix.withlab.big=log1p(exprs(temp))
par(mfrow=c(3,3))
toplot=mymatrix.withlab.big; rownames(toplot)[dim(toplot)[1]]="chrchrtdRFP"
rtsne_out=rtsne_out.5
#### plot all the distinct clusterings, for different nrs. of k
plot(rtsne_out$Y,pch=19,col=sc3.myk[colnames(toplot)],xlab="dim 1",ylab="dim 2"); legend("topright",legend=c(1:myk),col=filcols[1:myk],lwd=2)
plot(rtsne_out$Y,pch=19,col=sc3.3[colnames(toplot)],xlab="dim 1",ylab="dim 2"); legend("topright",legend=c(1:3),col=filcols[1:3],lwd=2)
plot(rtsne_out$Y,pch=19,col=sc3.4[colnames(toplot)],xlab="dim 1",ylab="dim 2");legend("topright",legend=c(1:4),col=filcols[1:4],lwd=2)
plot(rtsne_out$Y,pch=19,col=sc3.5[colnames(toplot)],xlab="dim 1",ylab="dim 2"); legend("topright",legend=c(1:5),col=filcols[1:5],lwd=2)
plot(rtsne_out$Y,pch=19,col=sc3.6[colnames(toplot)],xlab="dim 1",ylab="dim 2"); legend("topright",legend=c(1:6),col=filcols[1:6],lwd=2)
plot(rtsne_out$Y,pch=19,col=sc3.7[colnames(toplot)],xlab="dim 1",ylab="dim 2"); legend("topright",legend=c(1:7),col=filcols[1:7],lwd=2)
plot(rtsne_out$Y,pch=19,col=sc3.8[colnames(toplot)],xlab="dim 1",ylab="dim 2"); legend("topright",legend=c(1:8),col=filcols[1:8],lwd=2)
plot(rtsne_out$Y,pch=19,col=sc3.9[colnames(toplot)],xlab="dim 1",ylab="dim 2"); legend("topright",legend=c(1:9),col=filcols[1:9],lwd=2)
plot(rtsne_out$Y,pch=19,col=sc3.10[colnames(toplot)],xlab="dim 1",ylab="dim 2"); legend("topright",legend=c(1:10),col=filcols[1:10],lwd=2)
plot(rtsne_out$Y,pch=19,col=sc3.13[colnames(toplot)],xlab="dim 1",ylab="dim 2"); legend("topright",legend=c(1:13),col=filcols[1:13],lwd=2)

#### color according to gene expression, for marker genes and gois
toexclude=NULL
rownames(toplot)=elementMetadata(spikedgenes84.GR[rownames(toplot),])$geneName
mymax=max(mymatrix.withlab.big)*2/3
goi.plot=unique(c(unlist(goi),
	elementMetadata(spikedgenes84.GR[unlist(rownames(topDE_genes))])$geneName))
plotMarker=sort(goi.plot[which(goi.plot%in%rownames(toplot))]); plotMarker=plotMarker[which(!plotMarker%in%toexclude)]
choosecollist=list(c("ivory2","darkgrey", "black"),c("ivory2","orange", "red"),c("ivory2","lightblue", "blue"),c("ivory2","lightgreen", "darkgreen"))
for(j in 1:length(choosecollist)){
	for (i in 1:length(plotMarker)) {
	        choosecols=choosecollist[[j]]; plotcol=valuesToColorsAbs(as.numeric(toplot[plotMarker[i],]),choosecols,colnames(toplot),0,mymax)
		plot(rtsne_out$Y,pch=19, main=plotMarker[i],col=plotcol,xlab="dim 1",ylab="dim 2")}
	}
### color according to Xist expression 
choosecols=choosecollist[[1]]
plotcol=valuesToColorsAbs(log1p(countmatrix[getEnsgID(genes84.GR,"Xist"),colnames(toplot)]),choosecols,colnames(toplot),0,max(log1p(countmatrix))*2/3)
plot(rtsne_out$Y,pch=19, main="Xist",col=plotcol,xlab="dim 1",ylab="dim 2")

### (pre)-adipogenic genes tested in the manuscript; Xist and Mki67 in the case of the full dataset analysis (Fig. S8j)
geneIds=c("Fabp4","Cd34","Zfp423","Prdm16","Ppargc1a","Ebf2","Cd24a","Pdgfra","Pdgfrb","Adipoq","Wt1","Dlk1","Pparg","Xist","Mki67")
geneIds=geneIds[which(geneIds%in%rownames(toplot))]
fig1j=cbind(rtsne_out$Y,sc3.myk[colnames(toplot)], t(toplot[geneIds,]))
#### add batch and expression values for all these genes
colnames(fig1j)[1:3]=c("x","y","SC3.4")
#### output raw for Figure
write.table(fig1j,paste(path,"Fig1egS2h-raw.tSNE10x.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
#write.table(fig1j,paste(path,"FigS8j-raw.tSNE10xAllData.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

#### compare the markers with the C1-data derived population markers
load("mASC-M3Drop-ASC.htseq-bow-mmseqEnsg84-runCut-20161121-161121.208.Markers.RData")
m.10x=lapply(1:myk,function(x) rownames(head(marker_genes_sc_myk[marker_genes_sc_myk$Group==x,],100)))
mymatrix=matrix(data=0,ncol=6,nrow=myk)
for(i in 1:myk){
	mymatrix[i,]=sapply(list(p1.m,p2.m,p3.m,p1.m.hvg,p2.m.hvg,p3.m.hvg),function(x) length(which(m.10x[[i]]%in%x)))
}
par(mfrow=c(2,2))
barplot(table(scclusters.myk),main="Cells in Clusters - 10x")
rownames(mymatrix)=paste("C",1:myk,sep="")
colnames(mymatrix)=c("P1","P2","P3","P1.hvg","P2.hvg","P3.hvg")
barplot(t(mymatrix[,1:3]),col=c("darkgreen","red","blue"),main="M3Drop") #### Markers of the clusters derived from the M3Drop genes
barplot(t(mymatrix[,4:6]),col=c("darkgreen","red","blue"),main="HVG")  #### Markers of the clusters derived from the HVG genes
#### output raw for Figure
write.table(table(scclusters.myk),paste(path,"FigS2d-raw.CellsInCluster10x.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(mymatrix,paste(path,"Fig1fS2e-raw.ClusterOverlaps.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

#### and plot the expression of FACS marker genes (use distinct versions of the data with regards to normalisation)
countmatrix.10x=log1p(as.matrix(exprs(gbm)))
countmatrix.10x.interest=lapply(list("Actb","Ptprc","Pecam1","Cd34","Ly6a","Itgb1","Xist","Mki67","Adipoq","Retn","Cidec"),function(x) countmatrix.10x[getEnsgID(genes84.GR,x),])
countmatrix.filtered.interest=lapply(list("Actb","Ptprc","Pecam1","Cd34","Ly6a","Itgb1","Xist","Mki67","Adipoq","Retn","Cidec"),function(x) log1p(countmatrix[getEnsgID(genes84.GR,x),]))
final.interest=lapply(list("Actb","Actb","Actb","Cd34","Ly6a","Itgb1","Xist","Mki67","Adipoq","Retn","Retn"),function(x) log1p(exprs(temp)[getEnsgID(genes84.GR,x),]))
final.interest[[2]]=final.interest[[3]]=final.interest[[11]]=rep(0,length(colnames(exprs(temp))))
boxplot(countmatrix.10x.interest,names=c("Actb","Ptprc","Pecam1","Cd34","Ly6a","Itgb1","Xist","Mki67","Adipoq","Retn","Cidec"),col=c("grey","grey","grey","salmon3","salmon3","salmon3","yellow","yellow","red","red","red"))
boxplot(countmatrix.filtered.interest,names=c("Actb","Ptprc","Pecam1","Cd34","Ly6a","Itgb1","Xist","Mki67","Adipoq","Retn","Cidec"),col=c("grey","grey","grey","salmon3","salmon3","salmon3","yellow","yellow","red","red","red"))
boxplot(final.interest,names=c("Actb","Ptprc","Pecam1","Cd34","Ly6a","Itgb1","Xist","Mki67","Adipoq","Retn","Cidec"),col=c("grey","grey","grey","salmon3","salmon3","salmon3","yellow","yellow","red","red","red"))
temp1=as.data.frame(countmatrix.10x.interest);colnames(temp1)=c("Actb","Ptprc","Pecam1","Cd34","Ly6a","Itgb1","Xist","Mki67","Adipoq","Retn","Cidec")
#### output raw for Figure
write.table(temp1,paste(path,"FigS8i-raw.MarkerExpression10x.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F) #### this is generated on the total dataset, with n=2919 cells (adjust at input)
genevals=as.data.frame(final.interest)
colnames(genevals)=c("Actb","Ptprc","Pecam1","Cd34","Ly6a","Itgb1","Xist","Mki67","Adipoq","Retn","Cidec")
write.table(genevals,paste(path,"FigS2c-raw.Markergenes10x.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F) #### this is generated on the 1,804 cells

##### calculate the lineage scores #####
myclusters=scclusters.myk;names(myclusters)=colnames(exprs(temp))
for(j in c("adipo","stem","osteo","chondro","endo","immuno","hemato")) {
	temp1=log1p(countmatrix.10x[getEnsgID(spikedgenes84.GR,gsub(" ","",t(mymarkers[j,]))),])
		adiposcore=apply(temp1,2,function(x) sum(x));
	boxplot(lapply(1:myk,function(x) adiposcore[names(myclusters)[which(myclusters%in%x)]]),main=j,ylab="score")
	beanplot(lapply(1:myk,function(x) adiposcore[names(myclusters)[which(myclusters%in%x)]]),main=j,ylab="score")
	plot(density(adiposcore),main=j)
}
#### Output raw for Figure - Plot the scores for adipo & stemness
j="adipo"; temp1=log1p(countmatrix.10x[getEnsgID(spikedgenes84.GR,gsub(" ","",t(mymarkers[j,]))),])
adiposcore=apply(temp1,2,function(x) sum(x));
j="stem"; temp1=log1p(countmatrix.10x[getEnsgID(spikedgenes84.GR,gsub(" ","",t(mymarkers[j,]))),])
stemscore=apply(temp1,2,function(x) sum(x));
test=lapply(1:myk,function(x) adiposcore[names(myclusters)[which(myclusters%in%x)]])
temp=unlist(test); nrel=sapply(test,function(x) length(x))
tempcat=c(rep("P1",nrel[1]),rep("P2",nrel[2]),rep("P3",nrel[3]),rep("P4",nrel[4]))
temp=cbind(temp,unlist(lapply(1:myk,function(x) stemscore[names(myclusters)[which(myclusters%in%x)]])),
	tempcat);colnames(temp)=c("AdipoScore","StemScore","Clusters")
write.table(temp,paste(path,"FigS2f-raw.AdipoScore10x.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

##### merge the clusters G1 & G4 and get the markers; repeat analysis as above
myclusters.merged=myclusters
myclusters.merged[which(myclusters%in%4)]=1
marker_genes_sc_myk_merged <- M3DropGetMarkers(exprs(temp), myclusters.merged)
sc3.myk.merged=filcols[myclusters.merged]; names(sc3.myk)=rownames(pData(mysceset))

P.top.sc3.20.m=list(NA)
P.top.sc3.100.m=list(NA)
P.top.sc3.200.m=list(NA)
P.top.sc3.500.m=list(NA)
for(k in 1:myk) {
	P.top.sc3.20.m[[k]]=head(marker_genes_sc_myk_merged[marker_genes_sc_myk_merged$Group==k,],20)
	P.top.sc3.100.m[[k]]=head(marker_genes_sc_myk_merged[marker_genes_sc_myk_merged$Group==k,],100)
	P.top.sc3.200.m[[k]]=head(marker_genes_sc_myk_merged[marker_genes_sc_myk_merged$Group==k,],200)
	P.top.sc3.500.m[[k]]=head(marker_genes_sc_myk_merged[marker_genes_sc_myk_merged$Group==k,],500)
}
p1.10x.m=rownames(head(marker_genes_sc_myk_merged[marker_genes_sc_myk_merged$Group==1,],100))
p2.10x.m=rownames(head(marker_genes_sc_myk_merged[marker_genes_sc_myk_merged$Group==2,],100))
p3.10x.m=rownames(head(marker_genes_sc_myk_merged[marker_genes_sc_myk_merged$Group==3,],100))
p4.10x.m=rownames(head(marker_genes_sc_myk_merged[marker_genes_sc_myk_merged$Group==4,],100))
m.10x.m=lapply(1:4,function(x) rownames(head(marker_genes_sc_myk_merged[marker_genes_sc_myk_merged$Group==x,],100)))

mymatrix=matrix(data=0,ncol=6,nrow=4)
for(i in 1:4){
	mymatrix[i,]=sapply(list(p1.m,p2.m,p3.m,p1.m.hvg,p2.m.hvg,p3.m.hvg),function(x) length(which(m.10x.m[[i]]%in%x)))
}
par(mfrow=c(2,2))
barplot(table(myclusters.merged),main="Cells in Clusters merged - p10") #1018  664  122
rownames(mymatrix)=paste("C",1:4,sep="")
colnames(mymatrix)=c("P1","P2","P3","P1.hvg","P2.hvg","P3.hvg")
barplot(t(mymatrix[,1:3]),col=c("darkgreen","red","blue"),main="M3Drop")
barplot(t(mymatrix[,4:6]),col=c("darkgreen","red","blue"),main="HVG")
#### output raw for Figure
write.table(mymatrix,paste(path,"FigS2e-raw.ClusterOverlaps10xmerged.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
dev.off()

#### additional data plots 
pdf(paste(mypath,myname,".",genescut,"genecut.scran.final2.additions.Pub.pdf",sep=""),width=7,height=7)
par(mfrow=c(2,2))

#### Plot Silhouette values
sc3_plot_silhouette(mysceset, k = 13)
mys=c(1, 0.9, 0.77, 0.8, 0.7, 0.57, 0.62, 0.55, 0.52, 0.49, 0.53, 0.53, 0.52)
plot(2:14, mys, type="b", xlab="Number of Clusters", ylab="Silhouette",xlim=c(0,15))
#### Output raw for Figures
mysil=cbind(1:13, mys); colnames(mysil)=c("k","Silhouette")
write.table(mysil,paste(path,"Figs8l-raw.Silhouette10x.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

#### Check which lineage genes are among the marker genes
plotMarker=unique(c("cdtdRFP","Dlk1","Ppara","Pdgfra","Pdgfrb","Cd24a","Zfp423","Adipoq","Wt1","Ebf2",allmarkers))
genes84.GR[rownames(DE_genes[which(rownames(DE_genes)%in%getEnsgID(spikedgenes84.GR,plotMarker)),])] ### Fabp4 as only one
genes84.GR[rownames(hvg.out[which(rownames(hvg.out)%in%getEnsgID(spikedgenes84.GR,plotMarker)),])] ### Fabp4 as only one
hvg.marker=elementMetadata(genes84.GR[rownames(hvg.out[which(rownames(hvg.out)%in%getEnsgID(spikedgenes84.GR,plotMarker)),])])$geneName
DE.marker=elementMetadata(genes84.GR[rownames(DE_genes[which(rownames(DE_genes)%in%getEnsgID(spikedgenes84.GR,plotMarker)),])])$geneName
#hvg "Fabp4" "Sdc1"  "Ogn"   "Ebf2"  "Cd36"  "Thbd" 
#DE "Cd44"  "Thbd"  "Fabp4" "Vcam1" "Cd36"  "Sox5"  "Cd24a" "Aoc3"  "Sdc1"  "Runx2"
elementMetadata(genes84.GR[rownames(P.top.sc3.100[[1]])[which(rownames(P.top.sc3.100[[1]])%in%getEnsgID(spikedgenes84.GR,plotMarker))]])$geneName #"Ebf2"
elementMetadata(genes84.GR[rownames(P.top.sc3.100[[2]])[which(rownames(P.top.sc3.100[[2]])%in%getEnsgID(spikedgenes84.GR,plotMarker))]])$geneName
elementMetadata(genes84.GR[rownames(P.top.sc3.100[[3]])[which(rownames(P.top.sc3.100[[3]])%in%getEnsgID(spikedgenes84.GR,plotMarker))]])$geneName 
elementMetadata(genes84.GR[rownames(P.top.sc3.100[[4]])[which(rownames(P.top.sc3.100[[4]])%in%getEnsgID(spikedgenes84.GR,plotMarker))]])$geneName #"Ly6a" "Ogn"  "Cd34"
elementMetadata(genes84.GR[p1.m[which(p1.m%in%getEnsgID(spikedgenes84.GR,plotMarker))]])$geneName #"Ly6a" "Thbd"
elementMetadata(genes84.GR[p2.m[which(p2.m%in%getEnsgID(spikedgenes84.GR,plotMarker))]])$geneName #"Fabp4" "Sdc1"  "Aoc3"  "Pparg"
elementMetadata(genes84.GR[p3.m[which(p3.m%in%getEnsgID(spikedgenes84.GR,plotMarker))]])$geneName 


#### get correlation with Fabp4, Sca1 and Cd34
countmatrix.k.cpm=exprs(temp)
#### Fabp4
mycors.p=c(sapply(c("Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ebf2","Adipoq"),
	function(x) cor.test(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Fabp4")),]),
		log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c(x))[1],]),method="spearman")$p.val))
mycors=c(sapply(c("Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ebf2","Adipoq"),
		function(x) cor(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Fabp4")),]),
		log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c(x))[1],]),method="spearman")))
cor.table=cbind(mycors.p,mycors)
write.table(cor.table,
	paste(mypath,"QualPlots-M3Drop-",myname,".corTable.txt",sep=""),row.names=T,col.names=T,sep="\t",quote=F)
cor.table=cbind(p.adjust(mycors.p),mycors)
#### Output raw for Figures
write.table(cor.table,paste(path,"Figs1h-Fabp4-raw.CorTable10x.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

barplot(sort(-log10(mycors.p),decreasing=T),main="-logp.Cor.Fabp4",las=2)
barplot(mycors[names(sort(-log10(mycors.p),decreasing=T))],main="Cor.Fabp4",las=2)

#### Cd34
mycors.p.cd34=c(sapply(c("Fabp4","Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ebf2","Adipoq"),
	function(x) cor.test(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Cd34")),]),
		log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c(x))[1],]),method="spearman")$p.val))
mycors.cd34=c(sapply(c("Fabp4","Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ebf2","Adipoq"),
		function(x) cor(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Cd34")),]),
		log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c(x))[1],]),method="spearman")))
cor.table.cd34=cbind(mycors.p.cd34,mycors.cd34)
cor.table.cd34=cbind(p.adjust(mycors.p.cd34),mycors.cd34)
#### Output raw for Figures
write.table(cor.table.cd34,paste(path,"Figs6h-Cd34-raw.CorTable10x.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

##### Sca1 and Cd29 (not included)
mycors.p.sca1=c(sapply(c("Fabp4","Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ebf2","Adipoq"),
	function(x) cor.test(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Ly6a")),]),
		log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c(x))[1],]),method="spearman")$p.val))
mycors.sca1=c(sapply(c("Fabp4","Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ebf2","Adipoq"),
		function(x) cor(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Ly6a")),]),
		log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c(x))[1],]),method="spearman")))
cor.table.sca1=cbind(mycors.p.sca1,mycors.sca1)

mycors.p.itgb1=c(sapply(c("Fabp4","Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ebf2","Adipoq"),
	function(x) cor.test(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Itgb1")),]),
		log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c(x))[1],]),method="spearman")$p.val))
mycors.itgb1=c(sapply(c("Fabp4","Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ebf2","Adipoq"),
		function(x) cor(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Itgb1")),]),
		log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c(x))[1],]),method="spearman")))
cor.table.itgb1=cbind(mycors.p.itgb1,mycors.itgb1)		
par(mfrow=c(2,2))
write.table(cor.table.cd34, paste(mypath,"QualPlots-M3Drop-",myname,".corTable.Cd34.txt",sep=""),row.names=T,col.names=T,sep="\t",quote=F)
write.table(cor.table.sca1, paste(mypath,"QualPlots-M3Drop-",myname,".corTable.Sca1.txt",sep=""),row.names=T,col.names=T,sep="\t",quote=F)
write.table(cor.table.itgb1, paste(mypath,"QualPlots-M3Drop-",myname,".corTable.Cd29.txt",sep=""),row.names=T,col.names=T,sep="\t",quote=F)

barplot(sort(-log10(mycors.p.cd34),decreasing=T),main="-logp.Cor.Cd34",las=2)
barplot(mycors.cd34[names(sort(-log10(mycors.p.cd34),decreasing=T))],main="Cor.Cd34",las=2)
barplot(sort(-log10(mycors.p.sca1),decreasing=T),main="-logp.Cor.Sca1",las=2)
barplot(mycors.sca1[names(sort(-log10(mycors.p.sca1),decreasing=T))],main="Cor.Sca1",las=2)
barplot(sort(-log10(mycors.p.itgb1),decreasing=T),main="-logp.Cor.Itgb1",las=2)
barplot(mycors.itgb1[names(sort(-log10(mycors.p.itgb1),decreasing=T))],main="Cor.Itgb1",las=2)

#### Plot the expression of GOIs across the different populations ####
countmatrix.10x.interest=lapply(list("Adam12","Aoc3","Fabp4","Peg3","Pparg","Cd55","Il13ra1","Creb5","F3","Abcg1","Meox2","Cd34","Ly6a","Itgb1","Pdgfra","Pdgfrb"),function(x) 
	lapply(1:myk ,function(y) countmatrix.10x[getEnsgID(genes84.GR,x),names(myclusters[which(myclusters%in%y)])]))
countmatrix.filtered.interest=lapply(list("Adam12","Aoc3","Fabp4","Peg3","Pparg","Cd55","Il13ra1","Creb5","F3","Abcg1","Meox2","Cd34","Ly6a","Itgb1","Pdgfra","Pdgfrb"),function(x) 
	lapply(1:myk ,function(y)  log1p(countmatrix[getEnsgID(genes84.GR,x),names(myclusters[which(myclusters%in%y)])])))
final.interest=lapply(list("Adam12","Aoc3","Fabp4","Peg3","Pparg","Cd55","Il13ra1","Creb5","F3","Abcg1","Meox2","Cd34","Ly6a","Itgb1","Pdgfra","Pdgfrb"),function(x) 
	lapply(1:myk ,function(y)   log1p(exprs(temp)[getEnsgID(genes84.GR,x),names(myclusters[which(myclusters%in%y)])])))

beancols=list(c("green",rep("black",6)),c("gray",rep("black",6)),c("blue",rep("black",6)),c("darkgreen",rep("black",6)),c("orange",rep("black",6)),c("red",rep("black",6)))
beancols=list(c("darkgreen",rep("black",4)),c("green",rep("black",4)),c("red",rep("black",4)),c("blue",rep("black",4)))

ccol=c("green","gray","blue","darkgreen","orange","red")
ccol=c("green","red","blue","darkgreen")
mygoi=list("Adam12","Aoc3","Fabp4","Peg3","Pparg","Cd55","Il13ra1","Creb5","F3","Abcg1","Meox2","Cd34","Ly6a","Itgb1","Pdgfra","Pdgfrb")
for(i in 1:length(countmatrix.10x.interest)) {
	boxplot(countmatrix.10x.interest[[i]],main=mygoi[i],col=ccol)
	beanplot(countmatrix.10x.interest[[i]],main=mygoi[i],col=beancols)
	boxplot(countmatrix.filtered.interest[[i]],main=mygoi[i],col=c("green","gray","blue","darkgreen","orange","red"))
	beanplot(countmatrix.filtered.interest[[i]],main=mygoi[i],col=beancols)
	boxplot(final.interest[[i]],main=mygoi[i],col=c("green","gray","blue","darkgreen","orange","red"))
	beanplot(final.interest[[i]],main=mygoi[i],col=beancols)
}
#### Get the p-values for the expression ####
mypvals=matrix(data=1,nrow=16,ncol=6)
for(i in 1:length(countmatrix.10x.interest)) {
	temp=c(wilcox.test(final.interest[[i]][[1]],final.interest[[i]][[2]])$p.val,wilcox.test(final.interest[[i]][[1]],final.interest[[i]][[3]])$p.val,
		wilcox.test(final.interest[[i]][[1]],final.interest[[i]][[4]])$p.val,wilcox.test(final.interest[[i]][[2]],final.interest[[i]][[3]])$p.val,
		wilcox.test(final.interest[[i]][[2]],final.interest[[i]][[4]])$p.val,wilcox.test(final.interest[[i]][[3]],final.interest[[i]][[4]])$p.val)
	mypvals[i,]=p.adjust(temp,n=96,method="BH")
}
rownames(mypvals)=mygoi
colnames(mypvals)=c("1-2","1-3","1-4","2-3","2-4","3-4")
write.table(mypvals,paste(mypath,myname,".",genescut,"genecut.scran.final.additions2.pvals.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
write.table(mypvals,paste(path,"FigS2g-raw.ExpressionPvals.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

#### and do densities ####
for(j in c("adipo","stem","osteo","endo","immuno","hemato","chondro")) {
	temp1=log1p(countmatrix.10x[getEnsgID(spikedgenes84.GR,gsub(" ","",t(mymarkers[j,]))),])
		adiposcore=apply(temp1,2,function(x) sum(x));
	boxplot(lapply(1:myk,function(x) adiposcore[names(myclusters)[which(myclusters%in%x)]]),main=j,ylab="score")
	beanplot(lapply(1:myk,function(x) adiposcore[names(myclusters)[which(myclusters%in%x)]]),main=j,ylab="score")
	plot(density(adiposcore),main=j,xlim=c(0,10))
}
temp1=log1p(countmatrix.10x[getEnsgID(spikedgenes84.GR,gsub(" ","",t(mymarkers[j,]))),])
adiposcore=apply(temp1,2,function(x) sum(x)); adiposcore=lapply(1:myk,function(x) adiposcore[names(myclusters)[which(myclusters%in%x)]])
adipo.pval.10x=c(wilcox.test(adiposcore[[1]],adiposcore[[2]])$p.val,wilcox.test(adiposcore[[1]],adiposcore[[3]])$p.val,
	wilcox.test(adiposcore[[1]],adiposcore[[4]])$p.val,wilcox.test(adiposcore[[2]],adiposcore[[3]])$p.val,
	wilcox.test(adiposcore[[2]],adiposcore[[4]])$p.val,wilcox.test(adiposcore[[3]],adiposcore[[4]])$p.val)
names(adipo.pval.10x)=c("G1-G2","G1-G3","G1-G4","G2-G3","G2-G4","G3-G4")
#### output raw for Figure
write.table(adipo.pval.10x,paste(path,"FigS2f-raw.AdipoScore10xPvals.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


#### Output adipogenic/stem score as well as Marker score for raw Figures (duplicated)
goiMatrix=sapply(final.interest,function(x) unlist(x))
colnames(goiMatrix)=mygoi
cellnr=sapply(final.interest[[1]],function(x) length(x)) #699 664 122 319
myclasses=c(rep("C1",cellnr[1]),rep("C2",cellnr[2]),rep("C3",cellnr[3]),rep("C4",cellnr[4]))
temp1=log1p(countmatrix.10x[getEnsgID(spikedgenes84.GR,gsub(" ","",t(mymarkers["adipo",]))),])
adiposcore=apply(temp1,2,function(x) sum(x)); adp=unlist(lapply(1:myk,function(x) adiposcore[names(myclusters)[which(myclusters%in%x)]]))
temp1=log1p(countmatrix.10x[getEnsgID(spikedgenes84.GR,gsub(" ","",t(mymarkers["stem",]))),])
adiposcore=apply(temp1,2,function(x) sum(x)); stmn=unlist(lapply(1:myk,function(x) stemscore[names(myclusters)[which(myclusters%in%x)]]))
goiMatrix=cbind(adp,stmn,goiMatrix,myclasses)
write.table(goiMatrix,paste(path,"FigS2f-raw.AdipoMarkerScore10x.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
sapply(m.10x,function(x) elementMetadata(genes84.GR[x])$geneName[which(elementMetadata(genes84.GR[x])$geneName%in%c("Cd55",
	"F3","Il13ra1","Aoc3","Abcg1","Adam12","Fabp4","Pparg","Meox2","Creb5","Mgp","Tgfbi"))])

dev.off()
