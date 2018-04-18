###### Initialise variables and load functionality #####
source("plot-util.R")
source("util-asc.R")

mmseqpath="src/mmseq-latest/"
myset="bow-mmseqEnsg84-runCut-20161121"
mmset=myset
mydate="161121"
choosecols=c("ivory3","orange", "red")
choosecolsgenes=c("steelblue3","yellow", "red")
#### Load interesting sets of genes ######
load("genomes/annotation/GeneTypes84.R")
myname="ASC"
exonLength84=width(exons84.GR)
names(exonLength84)=names(exons84.GR)


##### READ MMSEQ FILES AS REQUIRED ######
####### Read the data from the allmmseq directory #######
path="adipo/ETH/mapped/asc/mmseq/"
outpath="adipo/ETH/mapped/asc/mmseq/plots/"
mynames=list.files(path,pattern="gene.mmseq")

##### Read information about samples (decide which cells not to include) #####
shortnames=gsub("[.][.]","",paste(gsub("BD","D",gsub("_",".",substr(mynames,1,9))), substr(mynames,38,39),sep=""))
load("adipo/ETH/mapped/asc/info/CellInfo_MergedASC.RData")
cellinfo=bigtable.ok
cellinfo[which(cellinfo[,"Cell"]%in%"ND"),"Cell"]="NE"; cellinfo[which(cellinfo[,"Cell"]%in%c("2RD","RDE")),"Cell"]="RD"
cellinfo[which(cellinfo[,"Cell"]%in%"R+N"),"Cell"]="RN"; cellinfo[which(cellinfo[,"Cell"]%in%"2R+N"),"Cell"]="RN"
cellinfo[which(cellinfo[,"Cell"]%in%"3N+R"),"Cell"]="RN"; cellinfo[which(cellinfo[,"Cell"]%in%"POS"),"Cell"]="R"
cellinfo[which(cellinfo[,"Cell"]%in%"NEG"),"Cell"]="N"
#### exclude doublets, debris, low alignment cells
exclude=unique(c(which(as.numeric(cellinfo[,"Aligned"])/as.numeric(cellinfo[,"Total"])<=0.4),
	which(as.numeric(cellinfo[,"Aligned"])<=400000),
	which(!cellinfo[,"Qual"]%in%c("","ok","super","ok (weird profile; possibly degraded)")),
	which(!cellinfo[,"Cell"]%in%c("N","R","NE","R?","RD")),grep("N[.]2[.]",rownames(cellinfo)),grep("N[.]1d[.]",rownames(cellinfo)),
	grep("P[.]d[.]",rownames(cellinfo)),grep("P[.]0[.]",rownames(cellinfo))))
cellinfo.ok=cellinfo[-exclude,] #### unique names substr(cellinfo.ok[,"Name"],1,8)
mynames=substr(mysplit(filemap[rownames(cellinfo.ok)],"09-",2,2),1,39)
names(mynames)=gsub("[.][.]","",paste(gsub("BD","D",gsub("_",".",substr(mynames,1,9))), substr(mynames,38,39),sep=""))

##### Read STAR cells ####
# path to the gsnap/htseq results
starpath="adipo/ETH/mapped/asc/STAR/"
#readSTARHTSeqASC(starpath,"cut",mynames,names(mynames),mydate) 
load(paste(outpath,"SSAll-STARHTSeqTags.-",mydate,".RData",sep="")) #htseqgenes and htseqexons #SSAll-STARHTSeqTags.-140717.RData

#countmatrix=mmseqcounts; myname=paste(myname,"mmseq",sep=".")
countmatrix=htseqgenes; myname=paste(myname,"htseq",sep=".")

#### scater inspection of cells - if needed
library(scater); library(destiny) ;library("ROCR")
#example_sceset <- newSCESet(countData = countmatrix); keep_feature <- rowSums(exprs(example_sceset) > 5) > 0
#example_sceset <- example_sceset[keep_feature,]
#example_sceset <- calculateQCMetrics(example_sceset, feature_controls = grep("ERCC-",rownames(example_sceset)))
#pdf(paste(outpath,"QualPlots-Scatter-",myname,"-",myset,"-",mydate,".pdf",sep=""),height=10,width=10)
#plot(example_sceset); plotQC(example_sceset)
#plotPCA(example_sceset); plotTSNE(example_sceset); plotDiffusionMap(example_sceset) 
#### overview of some fun genes ####
#plotExpression(example_sceset,getEnsgID(genes84.GR,c("Actb","Gapdh","Cd34","Ly6a","Itgb1","Ptprc","Pecam1",
#	"Pdgfra","Pdgfrb","Zfp423","Dlk1","Cd24a"))[1:12])
##CD45, Actb, Dlk1, Gapdh, Pecam1, Cd29, Pdgfra, Cd34, Pdgfrb, Cd24, Zfp423, Ly6a
#dev.off()

#### load the population data too 
load("adipo/ETH/mapped/pop/mmseq/plots/SSAll-STARHTSeqTags.-161121.Counts.RData")
genesPerCell.pop=apply(htseqgenes.pop,2,function(x) length(which(x>0)))
genereadPerCell.pop=apply(htseqgenes.pop,2,function(x) sum(x))/1000000
countmatrix.pop.cpm=sapply(1:length(genereadPerCell.pop),function(x) htseqgenes.pop[,x]/genereadPerCell.pop[x])
colnames(countmatrix.pop.cpm)=colnames(htseqgenes.pop)

library(M3Drop);  library(M3DExampleData); library(SC3); .pardefault <- par()
#par(.pardefault)

#### Only the Postitive ####
countmatrix.k=countmatrix[,c(grep("P1",colnames(countmatrix)),grep("P2",colnames(countmatrix)),grep("P3",colnames(countmatrix)))]
cellinfo.k=cellinfo.ok[colnames(countmatrix.k),]#
cellinfo.k[which(cellinfo.k[,"Cell"]%in%"R?"),"Cell"]="R"; cellinfo.k[which(cellinfo.k[,"Cell"]%in%"NE"),"Cell"]="N"
PvsNFr=sapply(c("R","N"), function(x) length(which(cellinfo.k[,"Cell"]%in%x)))
genesPerCell=apply(countmatrix.k,2,function(x) length(which(x>0)))
genereadPerCell=apply(countmatrix.k,2,function(x) sum(x))/1000000
readsPerCell=as.numeric(apply(cellinfo.k,1,function(x) x[7]))/1000000
names(readsPerCell)=names(genesPerCell)
example_sceset = newSCESet(countData = countmatrix.k); example_sceset =calculateQCMetrics(example_sceset, feature_controls = grep("ERCC-",rownames(example_sceset)))
countmatrix.k.cpm=sapply(1:length(genereadPerCell),function(x) countmatrix.k[,x]/genereadPerCell[x])
colnames(countmatrix.k.cpm)=colnames(countmatrix.k)
	
##### make different colors 
filcols=brewer.pal(9,"Pastel1")
# batch
batchcols=colnames(countmatrix.k); names(batchcols)=colnames(countmatrix.k)
batchcols[grep("P1",batchcols)]=filcols[1]; batchcols[grep("P2",batchcols)]=filcols[2]
batchcols[grep("P3",batchcols)]=filcols[3]
# R vs. P
rfpcols=as.vector(cellinfo.k[,"Cell"]); names(rfpcols)=colnames(countmatrix.k);
rfpcols[grep("R",rfpcols)]="red"; rfpcols[grep("N",rfpcols)]="snow3"
# Readnrs.
readnrcols=log2(as.numeric(as.character(cellinfo.k[,7]))); names(readnrcols)=colnames(countmatrix.k);
readnrcols=valuesToColorsAbs(readnrcols,choosecolsgenes,names(readnrcols),min(readnrcols),max(readnrcols))
# Genenrs.
genenrcols=log2(genesPerCell); names(genenrcols)=colnames(countmatrix.k);
genenrcols=valuesToColorsAbs(genenrcols,choosecolsgenes,names(genenrcols),min(genenrcols),max(genenrcols))
# RFP counts
genemax=max(log1p(countmatrix.k.cpm))*2/3
rrfpcols=log1p(countmatrix.k.cpm["chrtdRFP",]); names(rrfpcols)=colnames(countmatrix.k.cpm);
rrfpcols=valuesToColorsAbs(rrfpcols,choosecols,names(rrfpcols),0,genemax)
# Dlk1 counts
prefcols=log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,"Dlk1")[1],]); names(prefcols)=colnames(countmatrix.k.cpm);
prefcols=valuesToColorsAbs(prefcols,choosecols,names(prefcols),0,genemax)
# Fabp4 counts
fabp4cols=log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,"Fabp4")[1],]); names(fabp4cols)=colnames(countmatrix.k.cpm);
fabp4cols=valuesToColorsAbs(fabp4cols,choosecols,names(fabp4cols),0,genemax)
# Cd34 counts
cd34cols=log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,"Cd34")[1],]); names(fabp4cols)=colnames(countmatrix.k.cpm);
cd34cols=valuesToColorsAbs(cd34cols,choosecols,names(cd34cols),0,genemax)
# adipogenic genes
#mypca=PCA(t(scale(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,
#	unique(c(c(GENETYPES84[[9]][[4]],GENETYPES84[[9]][[5]],GENETYPES84[[9]][[6]]),gsub(" ","",t(mymarkers["adipo",]))))),]))),graph=FALSE)
#adipocols=mypca$ind$coord[,1];names(adipocols)=colnames(countmatrix.k);
temp=log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,gsub(" ","",t(mymarkers["adipo",]))),])
adipocols=adiposcore=apply(temp,2,function(x) sum(x));names(adipocols)=colnames(countmatrix.k);
adipocols=valuesToColorsAbs(adipocols,choosecols,names(adipocols),0,30)
# stem genes
#mypca=PCA(t(scale(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,gsub(" ","",t(mymarkers["stem",]))),]))),graph=FALSE)
#stemcols=mypca$ind$coord[,1];names(stemcols)=colnames(countmatrix.k);
temp=log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,gsub(" ","",t(mymarkers["stem",]))),])
stemcols=stemscore=apply(temp,2,function(x) sum(x));names(stemcols)=colnames(countmatrix.k);
stemcols=valuesToColorsAbs(stemcols,choosecols,names(stemcols),0,30)
# osteo genes
#mypca=PCA(t(scale(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,gsub(" ","",t(mymarkers["stem",]))),]))),graph=FALSE)
#stemcols=mypca$ind$coord[,1];names(stemcols)=colnames(countmatrix.k);
temp=log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,gsub(" ","",t(mymarkers["osteo",]))),])
osteocols=osteoscore=apply(temp,2,function(x) sum(x));names(osteocols)=colnames(countmatrix.k);
osteocols=valuesToColorsAbs(osteocols,choosecols,names(osteocols),0,30)
### endo
temp=log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,gsub(" ","",t(mymarkers["endo",]))),])
endocols=endoscore=apply(temp,2,function(x) sum(x));names(endocols)=colnames(countmatrix.k);
endocols=valuesToColorsAbs(endocols,choosecols,names(endocols),0,30)
### Immuno
temp=log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,gsub(" ","",t(mymarkers["immuno",]))),])
immunocols=immunoscore=apply(temp,2,function(x) sum(x));names(immunocols)=colnames(countmatrix.k);
immunocols=valuesToColorsAbs(immunocols,choosecols,names(immunocols),0,30)

colvector=list(batchcols,readnrcols,genenrcols,rfpcols,rrfpcols,prefcols,cd34cols,fabp4cols,adipocols,stemcols,osteocols,endocols,immunocols)
names(colvector)=c("Batch","ReadNrs","GeneNrs","Red","RFP","Pref","Cd34","Fabp4","Adipo","Stem","Osteo","Endo","Immuno")

pdf(paste(outpath,"QualPlots-M3Drop-",myname,"-",myset,"-",mydate,".pdf",sep=""),height=10,width=10)
par(mfrow=c(2,2))

#### General info on the three replicates (in this case P1, P2, P3, each performed on one day)
barplot(PvsNFr/sum(PvsNFr),col=c("red","snow3"),names=c("P","N"),main="Captured Cells")
beanplot(lapply(c("P1","P2","P3"),function(x) genesPerCell[grep(x,names(genesPerCell))]),names=c("P1","P2","P3"),main="Genes/Cell")
beanplot(lapply(c("P1","P2","P3"),function(x) genereadPerCell[grep(x,names(genereadPerCell))]),names=c("P1","P2","P3"),main="Genereads/Cell")
beanplot(lapply(c("P1","P2","P3"),function(x) readsPerCell[grep(x,names(readsPerCell))]),names=c("P1","P2","P3"),main="Reads/Cell")
#### output raw for Figure
temp=unlist(lapply(c("P1","P2","P3"),function(x) genereadPerCell[grep(x,names(genereadPerCell))]))
nrel=sapply(c("P1","P2","P3"), function(x) length(grep(x,names(genereadPerCell))))
tempcat=c(rep("P1",nrel[1]),rep("P2",nrel[2]),rep("P3",nrel[3]))
temp=cbind(temp,tempcat);colnames(temp)=c("Reads/Cell","Batch")
write.table(temp,paste(path,"FigS1a-raw.alignedReads.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F) #74 71 63

#### Comparison with population data 
temp=lapply(c("P1","P2","P3"),function(x) genesPerCell[grep(x,names(genesPerCell))])
temp2.counts=sapply(c("P1","P2","P3"), 
	function(x) log1p(apply(countmatrix.k.cpm[,grep(x,colnames(countmatrix.k.cpm))],1,function(y) mean(y))))
temp2=apply(temp2.counts,2,function(x) length(which(x>0)))
temp3=genesPerCell.pop
bigGenesPerCell=list(temp[[1]],temp[[2]],temp[[3]],temp2,genesPerCell.pop[19:22])
names(bigGenesPerCell)=c("P1","P2","P3","Merged_SC","Pop_SMART")
beanplot(bigGenesPerCell,xlab="Nr. Genes/Sample")
#### output raw for Figure
temp=unlist(bigGenesPerCell); nrel=sapply(bigGenesPerCell,function(x) length(x))
tempcat=c(rep("P1",nrel[1]),rep("P2",nrel[2]),rep("P3",nrel[3]),rep("SC",nrel[4]),rep("Bulk",nrel[5]))
temp=cbind(temp,tempcat);colnames(temp)=c("Genes/Cell","Sample")
write.table(temp,paste(path,"FigS1d-raw.genesSample.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


### subset of genes that show expression 
selgenes.ASC=which(apply(temp2.counts,1,function(x) sum(x))>0);selgenes.Pop=which(apply(countmatrix.pop.cpm,1,function(x) sum(x))>0)
selgenes.Pop=names(selgenes.ASC)[which(names(selgenes.ASC)%in%names(selgenes.Pop))]
selgenes.ASC=which(apply(countmatrix.k.cpm,1,function(x) sum(x))>0)
allcors=lapply(1:3,function(x) sapply(19:22, function(y)
	scatterSmoothPlot(log1p(temp2.counts[,x]),log1p(countmatrix.pop.cpm[,y]),x,y,"spearman")))
names(allcors)=c("P1","P2","P3")
beanplot(allcors,main="SC vs. Pop",ylab="Spearman cor")
#### output raw for Figure
write.table(as.data.frame(allcors),paste(path,"FigS1c-raw.CorSCvsPop.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

##### Single-cell vs. single-cell correlation - all genes or ERCC genes (choose only a subset of cells)
allcors.SC=lapply(seq(1,200,by=30),function(x) sapply(seq(2,200,by=30), function(y)
	scatterSmoothPlot(log1p(countmatrix.k.cpm[,x]),log1p(countmatrix.k.cpm[,y]),x,y,"spearman")))
allcors.SC.ERCC=lapply(seq(1,200,by=30),function(x) sapply(seq(2,200,by=30), function(y)
	scatterSmoothPlot(log1p(countmatrix.k.cpm[grep("ERCC-",rownames(countmatrix.k.cpm)),x]),log1p(countmatrix.k.cpm[grep("ERCC-",rownames(countmatrix.k.cpm)),y]),x,y,"spearman")))
beanplot(list(unlist(allcors.SC),unlist(allcors.SC.ERCC)),main="ASC-cors",names=c("All","ERCC"),ylab="Spearman cor")
#### output raw for Figure
temp=cbind(unlist(allcors.SC),unlist(allcors.SC.ERCC))
colnames(temp)=c("All","ERCC")
write.table(temp,paste(path,"FigS1b-raw.CorASCERCC.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

#### plot the density of the scores based on the lineage markers 
plot(adiposcore,stemscore,main="Adipo vs. Stem",pch=19,col=colvector[[1]]); text(adiposcore,stemscore,labels=substr(names(adiposcore),3,8),col=colvector[[4]],cex.lab=0.6)
plot(adiposcore,osteoscore,main="Adipo vs. Osteo",pch=19,col=colvector[[1]]);text(adiposcore,osteoscore,labels=substr(names(adiposcore),3,8),col=colvector[[4]],cex.lab=0.6)
plot(adiposcore,endoscore,main="Adipo vs. Endo",pch=19,col=colvector[[1]]); text(adiposcore,endoscore,labels=substr(names(adiposcore),3,8),col=colvector[[4]],cex.lab=0.6)
plot(stemscore,endoscore,main="Stem vs. Endo",pch=19,col=colvector[[1]]);text(stemscore,endoscore,labels=substr(names(adiposcore),3,8),col=colvector[[4]],cex.lab=0.6)
plot(density(adiposcore),main="Adipogenic Score",xlim=c(0,50)); plot(density(stemscore),main="Stem Score",xlim=c(0,50))
plot(density(osteoscore),main="Osteo Score",xlim=c(0,50)); plot(density(endoscore),main="Endo Score",xlim=c(0,50))
plot(density(immunoscore),main="Immuno Score",xlim=c(0,50));

##### clean the display
par(.pardefault)
plotQC(example_sceset)
par(.pardefault)

###### Apply M3Drop to get informative (DE) genes across all cells 
Normalized_data <- M3DropCleanData(countmatrix.k, is.counts=T, min_detected_genes=3000)
#dim(Normalized_data$data) # [1] 17287 208
fits <- M3DropDropoutModels(Normalized_data$data)
data.frame(MM=fits$MMFit$SAr, Logistic=fits$LogiFit$SAr,DoubleExpo=fits$ExpoFit$SAr) #MM best fit
data.frame(MM=fits$MMFit$SSr, Logistic=fits$LogiFit$SSr, DoubleExpo=fits$ExpoFit$SSr) #MM best fit
DE_genes <- M3DropDifferentialExpression(Normalized_data$data,mt_method="fdr", mt_threshold=0.05) # 285 stringently DE genes ... 527 at FDR 0.05 741
plotMarker=unique(c("cdtdRFP","Dlk1","Ppara","Pdgfra","Pdgfrb","Cd24a","Zfp423","Adipoq","Wt1","Ebf2",allmarkers))
genes84.GR[rownames(DE_genes[which(rownames(DE_genes)%in%getEnsgID(spikedgenes84.GR,plotMarker)),])] ### Fabp4 as only one among lineage geens
	
heat_out = M3DropExpressionHeatmap(DE_genes$Gene, Normalized_data$data, as.vector(cellinfo.k[colnames(Normalized_data$data),][,"Cell"]))
par(.pardefault)
heat_out = M3DropExpressionHeatmap(DE_genes$Gene, Normalized_data$data, colvector[[1]])
par(.pardefault)

cell_populations_4 = M3DropGetHeatmapCellClusters(heat_out, k=4); cell_populations_3 = M3DropGetHeatmapCellClusters(heat_out, k=3)
marker_genes_4 = M3DropGetMarkers(Normalized_data$data, cell_populations_4); marker_genes_3 = M3DropGetMarkers(Normalized_data$data, cell_populations_3)

mymatrix.log=log1p(Normalized_data$data[rownames(DE_genes),]); #### run SC3 on informative genes only (used in manuscript)
mysceset = newSCESet(countData = Normalized_data$data[rownames(DE_genes),]); 
mysceset <- calculateQCMetrics(mysceset); mysceset <- sc3(mysceset, ks = 3:4,n_cores=1)
mysceset.all = newSCESet(countData = log1p(countmatrix.k)); #### run SC3 on entire matrix
mysceset.all <- calculateQCMetrics(mysceset.all); mysceset.all <- sc3(mysceset.all, ks = 2:6,n_cores=1)
#scclusters.3=sc3_summarise_results(mysceset, k = 3); scclusters.4=sc3_summarise_results(mysceset, k = 4)
scclusters.3=pData(mysceset)[,"sc3_3_clusters"]; scclusters.4=pData(mysceset)[,"sc3_4_clusters"]
names(scclusters.3)=names(scclusters.4)=rownames(pData(mysceset))
scclusters.3.all=pData(mysceset.all)[,"sc3_3_clusters"]; scclusters.4.all=pData(mysceset.all)[,"sc3_4_clusters"]
names(scclusters.3.all)=names(scclusters.4.all)=rownames(pData(mysceset.all))

#### explore clustering #####
sc3_plot_cluster_stability(mysceset, k = 3); sc3_plot_cluster_stability(mysceset, k = 4)
sc3_plot_silhouette(mysceset, k = 3); sc3_plot_silhouette(mysceset, k = 4)
mysceset.est <- sc3_estimate_k(mysceset); mysceset.est2 <- sc3_estimate_k(mysceset.all)
mysceset.est@sc3$k_estimation #4; mysceset.est2@sc3$k_estimation #6
sc3_plot_consensus( mysceset, k = 3, show_pdata = c("cell_id","log10_total_features","sc3_3_clusters", "sc3_3_log2_outlier_score"))
sc3_plot_expression(mysceset, k = 3, show_pdata = c("cell_id", "log10_total_features","sc3_3_clusters", "sc3_3_log2_outlier_score"))
sc3_plot_cluster_stability(mysceset, k = 3); sc3_plot_cluster_stability(mysceset, k = 4)
mysilhouette=c(0.77,0.76,0.54, 0.58,0.49); names(mysilhouette)=2:6

##### get the clusters - either from hierarchical clustering or from SC3 + get the markers for each cluster ####
heat.3=filcols[cell_populations_3]; names(heat.3)=names(cell_populations_3)
heat.4=filcols[cell_populations_4]; names(heat.4)=names(cell_populations_4)
sc3.3=filcols[scclusters.3]; names(sc3.3)=names(cell_populations_3); cID=scclusters.3;names(cID)=names(cell_populations_3)
sc3.4=filcols[scclusters.4]; names(sc3.4)=names(cell_populations_4)
marker_genes_sc_3 <- M3DropGetMarkers(Normalized_data$data, scclusters.3)
marker_genes_sc_4 <- M3DropGetMarkers(Normalized_data$data, scclusters.4)
##### cluster colors 
ccols=list(heat.3, heat.4, sc3.3,sc3.4); names(ccols)=c("heat.3", "heat.4","sc3.3","sc3.4")

### Heatmap & SC3 markers, top 30 and top 50
P3.top=head(marker_genes_3[marker_genes_3$Group==1,],30); P2.top=head(marker_genes_3[marker_genes_3$Group==2,],30)
P1.top=head(marker_genes_3[marker_genes_3$Group==3,],30)
### SC3 
P3.top.sc3=head(marker_genes_sc_3[marker_genes_sc_3$Group==1,],30); P2.top.sc3=head(marker_genes_sc_3[marker_genes_sc_3$Group==2,],30)
P1.top.sc3=head(marker_genes_sc_3[marker_genes_sc_3$Group==3,],30)
### top 50
P3.top.sc3.100=head(marker_genes_sc_3[marker_genes_sc_3$Group==1,],50); P2.top.sc3.100=head(marker_genes_sc_3[marker_genes_sc_3$Group==2,],50)
P1.top.sc3.100=head(marker_genes_sc_3[marker_genes_sc_3$Group==3,],50)
ofinterest= unlist(c(GENETYPES84[[1]],GENETYPES84[[19]],GENETYPES84[[10]],GENETYPES84[[11]],GENETYPES84[[16]],GENETYPES84[[18]],GENETYPES84[[20]]))
P3.top.sc3.100=P3.top.sc3.100[which(rownames(P3.top.sc3.100)%in%ofinterest),]
P2.top.sc3.100=P2.top.sc3.100[which(rownames(P2.top.sc3.100)%in%ofinterest),]
P1.top.sc3.100=P1.top.sc3.100[which(rownames(P1.top.sc3.100)%in%ofinterest),]

#### Generate different matrices for visualisation
mymatrix.withlab.big=log1p(countmatrix.k.cpm)
mymatrix.withlab=log1p(Normalized_data$data[c(rownames(P3.top),rownames(P2.top),rownames(P1.top)),])
mymatrix.withlab2=log1p(Normalized_data$data[c(rownames(P3.top.sc3),rownames(P2.top.sc3),rownames(P1.top.sc3)),])
mymatrix.withlab3=log1p(Normalized_data$data[c(rownames(P3.top.sc3.100),rownames(P2.top.sc3.100),rownames(P1.top.sc3.100)),])
rownames(mymatrix.withlab)[grep("ENS",rownames(mymatrix.withlab))]=elementMetadata(genes84.GR[rownames(mymatrix.withlab)[grep("ENS",rownames(mymatrix.withlab))]])$geneName
rownames(mymatrix.withlab2)[grep("ENS",rownames(mymatrix.withlab2))]=elementMetadata(genes84.GR[rownames(mymatrix.withlab2)[grep("ENS",rownames(mymatrix.withlab2))]])$geneName
rownames(mymatrix.withlab3)[grep("ENS",rownames(mymatrix.withlab3))]=elementMetadata(genes84.GR[rownames(mymatrix.withlab3)[grep("ENS",rownames(mymatrix.withlab3))]])$geneName
rownames(mymatrix.withlab.big)[grep("ENS",rownames(mymatrix.withlab.big))]=elementMetadata(genes84.GR[rownames(mymatrix.withlab.big)[grep("ENS",rownames(mymatrix.withlab.big))]])$geneName
##### Save the data for later processing 
#save(Normalized_data,mysceset,P3.top.sc3.100, P2.top.sc3.100,P1.top.sc3.100,mymatrix.withlab.big,
#	marker_genes_3, marker_genes_4,countmatrix.k, countmatrix,
#	file=paste(outpath,"ASC-saved-M3Drop-",myname,"-",myset,"-",mydate,".18.RData",sep=""))
##### load previously saved data
load(paste(outpath,"ASC-saved-M3Drop-",myname,"-",myset,"-",mydate,".18.RData",sep=""))

##### Plot - diffusion map and t-SNE (used in manuscript); color according to clusters & markers
par(mfrow=c(3,3))
plot(2:6,mysilhouette,main="Silhouette values")
#plotDiffusionMap(mysceset) 
####dif <- DiffusionMap(t(mymatrix.log))
#dif.DE527=dif
#plot(dif)
dif=dif.DE527
toplot=mymatrix.withlab.big
for (i in 1:length(colvector)) { plot(dif, pch = 19, col = colvector[[i]],main =paste("Cats - ",names(colvector)[i]))}
for (i in 1:length(ccols)) {plot(dif, pch = 19, col = ccols[[i]],main =paste("Clusters - ",names(ccols)[i])) }
plot(eigenvalues(dif), ylim = c(0,0.3), pch = 20, xlab = "Diffusion component (DC)", ylab = "Eigenvalue")
mymax=max(mymatrix.withlab.big)*2/3; plotMarker=unique(c(myMarkerAll,mychoiceMarkerAll,allmarkers,rownames(mymatrix.withlab3)))
plotMarker=sort(plotMarker[which(plotMarker%in%rownames(toplot))])
choosecollist=list(c("ivory2","darkgrey", "black"),c("ivory2","orange", "red"),c("ivory2","lightblue", "blue"),c("ivory2","lightgreen", "darkgreen"))
for (j in 1:length(choosecollist)) {
for (i in 1:length(plotMarker)) {
        choosecols=choosecollist[[j]]; plotcol=valuesToColorsAbs(toplot[plotMarker[i],],choosecols,colnames(mymatrix.withlab),0,mymax)
	plot(dif, pch = 19, col = plotcol,main =plotMarker[i])
}}

##### generate t-SNE or load the previous one
#rtsne_out <- Rtsne(as.matrix(t(mymatrix.log)),theta=0.0000001,perplexity=as.integer(dim(mymatrix.log)[2]/6))
#rtsne_out.DE527=rtsne_out
load(paste(outpath,"QualPlots-M3Drop-",myname,"-",myset,"-",mydate,".RtSNE.527.208.RData",sep=""))
rtsne_out=rtsne_out.DE527
toplot=mymatrix.withlab.big
for (i in 1:length(colvector)) {
	plot(rtsne_out$Y,pch=19, main=paste("Cats - ",names(colvector)[i]),col=colvector[[i]][colnames(mymatrix.withlab)],xlab="dim 1",ylab="dim 2")
#	text(rtsne_out$Y, labels=colnames(mymatrix.withlab),col=colvector[[i]][colnames(mymatrix.withlab)],cex=0.2) # cell labels off
}
for (i in 1:length(ccols)) {
	plot(rtsne_out$Y,pch=19, main=paste("Clusters - ",names(ccols)[i]),col=ccols[[i]][colnames(mymatrix.withlab)],xlab="dim 1",ylab="dim 2")
#	text(rtsne_out$Y, labels=colnames(mymatrix.withlab),col=ccols[[i]][colnames(mymatrix.withlab)],cex=0.2)  # cell labels off
}
mymax=max(mymatrix.withlab.big)*2/3
plotMarker=unique(c(myMarkerAll,mychoiceMarkerAll,allmarkers,rownames(mymatrix.withlab3)))
plotMarker=sort(plotMarker[which(plotMarker%in%rownames(toplot))])
choosecollist=list(c("ivory2","darkgrey", "black"),c("ivory2","orange", "red"),c("ivory2","lightblue", "blue"),c("ivory2","lightgreen", "darkgreen"))
for (j in 1:length(choosecollist)) {
for (i in 1:length(plotMarker)) {
        choosecols=choosecollist[[j]]; plotcol=valuesToColorsAbs(toplot[plotMarker[i],],choosecols,colnames(mymatrix.withlab),0,mymax)
        plot(rtsne_out$Y,pch=19, main=plotMarker[i],col=plotcol,xlab="dim 1",ylab="dim 2")
#        text(rtsne_out$Y, labels=colnames(mymatrix.withlab),col=plotcol,cex=0.2)	
}}
#### output raw for Figure
#### gois + add batch and expression values for all these genes; also and RFP+ and Neg.
geneIds=c("Fabp4","Cd34","Zfp423","Prdm16","Ppargc1a","Ebf2","Cd24a","Pdgfra","Pdgfrb","Adipoq","Wt1","Dlk1","Pparg","chrtdRFP")
fig1b=cbind(rtsne_out$Y,colvector[[1]][colnames(mymatrix.withlab)],
	ccols[[3]][colnames(mymatrix.withlab)],ccols[[4]][colnames(mymatrix.withlab)],
	colvector[[4]][colnames(mymatrix.withlab)],
	t(toplot[geneIds,]))
colnames(fig1b)[1:6]=c("x","y","Batch","SC3.3","SC3.4","RFP")
write.table(fig1b,paste(path,"Fig1bS1S8-raw.tSNEC1.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

##### Additionall plotting: heatmap
par(mfrow=c(1,1))
# heatmap
makeprefheatmaprlog(mymatrix.withlab,spikedgenes84.GR,colvals=colvector[[9]],mainname="row.scaled", rowvals="grey")
makeprefheatmaprlog(mymatrix.withlab,spikedgenes84.GR,colvals=ccols[[1]],mainname="row.scaled", rowvals="grey")
makeprefheatmaprlog(mymatrix.withlab,spikedgenes84.GR,colvals=ccols[[3]],mainname="row.scaled", rowvals="grey")
#SC3
makeprefheatmaprlog(mymatrix.withlab2,spikedgenes84.GR,colvals=colvector[[9]],mainname="row.scaled", rowvals="grey")
makeprefheatmaprlog(mymatrix.withlab2,spikedgenes84.GR,colvals=ccols[[1]],mainname="row.scaled", rowvals="grey")
makeprefheatmaprlog(mymatrix.withlab2,spikedgenes84.GR,colvals=ccols[[3]],mainname="row.scaled", rowvals="grey")
#chosen
makeprefheatmaprlog(mymatrix.withlab3,spikedgenes84.GR,colvals=colvector[[9]],mainname="row.scaled", rowvals="grey")
makeprefheatmaprlog(mymatrix.withlab3,spikedgenes84.GR,colvals=ccols[[1]],mainname="row.scaled", rowvals="grey")
makeprefheatmaprlog(mymatrix.withlab3,spikedgenes84.GR,colvals=ccols[[3]],mainname="row.scaled", rowvals="grey")

####### second round of plots on only top30 marker genes for the 3 populations 
par(mfrow=c(3,3))
#dif2 <- DiffusionMap(t(mymatrix.withlab2))
#plot(dif)
dif2=dif.DE527.top30
for (i in 1:length(colvector)) { plot(dif2, 1:2, pch = 19, col = colvector[[i]],main =paste("Cats - ",names(colvector)[i]))}
for (i in 1:length(ccols)) { plot(dif2, 1:2, pch = 19, col = ccols[[i]],main =paste("Clusters - ",names(ccols)[i]))}
mymax=max(mymatrix.withlab.big)*2/3
plotMarker=unique(c(myMarkerAll,mychoiceMarkerAll,allmarkers,rownames(mymatrix.withlab3)))
plotMarker=sort(plotMarker[which(plotMarker%in%rownames(toplot))])
choosecollist=list(c("ivory2","darkgrey", "black"),c("ivory2","orange", "red"),c("ivory2","lightblue", "blue"),c("ivory2","lightgreen", "darkgreen"))
for (j in 1:length(choosecollist)) {
for (i in 1:length(plotMarker)) {
        choosecols=choosecollist[[j]]
        plotcol=valuesToColorsAbs(toplot[plotMarker[i],],choosecols,colnames(mymatrix.withlab),0,mymax)
	plot(dif2, 1:2,pch = 19, col = plotcol,main =plotMarker[i])
}}
plot(eigenvalues(dif2), ylim = 0:1, pch = 19, xlab = "Diffusion component (DC)", ylab = "Eigenvalue")
#dif.DE527.top30=dif2

##### and tSNE
#rtsne_out2 <- Rtsne(as.matrix(t(mymatrix.withlab2)),theta=0.0000001,perplexity=as.integer(dim(mymatrix.log)[2]/6))
#rtsne_out.DE527.top30=rtsne_out2 
rtsne_out2=rtsne_out.DE527.top30
toplot=mymatrix.withlab.big
for (i in 1:length(colvector)) {
	plot(rtsne_out2$Y,pch=19, main=paste("Cats - ",names(colvector)[i]),col=colvector[[i]][colnames(mymatrix.withlab)],xlab="dim 1",ylab="dim 2")
#	text(rtsne_out2$Y, labels=colnames(mymatrix.withlab),col=colvector[[i]][colnames(mymatrix.withlab)],cex=0.2)
}
for (i in 1:length(ccols)) {
	plot(rtsne_out2$Y,pch=19, main=paste("Clusters - ",names(ccols)[i]),col=ccols[[i]][colnames(mymatrix.withlab)],xlab="dim 1",ylab="dim 2")
#	text(rtsne_out2$Y, labels=colnames(mymatrix.withlab),col=ccols[[i]][colnames(mymatrix.withlab)],cex=0.2)
}
#par(mfrow=c(3,3)); 
mymax=max(mymatrix.withlab.big)*2/3
plotMarker=unique(c(myMarkerAll,mychoiceMarkerAll,allmarkers,rownames(mymatrix.withlab3)))
plotMarker=sort(plotMarker[which(plotMarker%in%rownames(toplot))])
choosecollist=list(c("ivory2","darkgrey", "black"),c("ivory2","orange", "red"),c("ivory2","lightblue", "blue"),c("ivory2","lightgreen", "darkgreen"))
for (j in 1:length(choosecollist)) {
for (i in 1:length(plotMarker)) {
        choosecols=choosecollist[[j]]; plotcol=valuesToColorsAbs(toplot[plotMarker[i],],choosecols,colnames(mymatrix.withlab),0,mymax)
        plot(rtsne_out2$Y,pch=19, main=plotMarker[i],col=plotcol,xlab="dim 1",ylab="dim 2")
  #      text(rtsne_out2$Y, labels=colnames(mymatrix.withlab),col=plotcol,cex=0.2)	
}}
##### Save the maps if analysis is run for the first time
#save(rtsne_out.DE527,dif.DE527,dif.DE527.top30,rtsne_out.DE527.top30,
#	file=paste(outpath,"QualPlots-M3Drop-",myname,"-",myset,"-",mydate,".RtSNE.527.208.RData",sep=""))



################ Reviews: second way of selecting genes: HVG, using scran 
#### load library & set-up; assign cell cycle phases; get variance; get HVGs
library(scran)
sce = newSCESet(countData = Normalized_data$data); 
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assigned <- cyclone(sce, pairs=mm.pairs) 
cell.cycle=assigned$phases; names(cell.cycle)=colnames(Normalized_data$data)
sce <- computeSumFactors(sce, clusters=NULL)
var.fit <- trendVar(sce, trend="loess", span=0.4,use.spikes=F)
var.out <- decomposeVar(sce, var.fit); top.hvgs <- order(var.out$bio, decreasing=TRUE)
head(var.out[top.hvgs,])
decomp=var.out; fit=var.fit
plot(decomp$mean, decomp$total, xlab="Mean log-expression", ylab="Variance")
o <- order(decomp$mean)
lines(decomp$mean[o], decomp$tech[o], col="red", lwd=2)
points(fit$mean, fit$var, col="red", pch=16)
hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] #### 1830 genes

#### overlaps between HVG and M3Drop genes 
which(rownames(DE_genes)%in%rownames(hvg.out)) ### 520/527 genes amoung high variability genes

#### analysis on the set of 1830 genes
temp=convertTo(sce, type="monocle")
mymatrix.hvg.log=log1p(exprs(temp)[unique(c(rownames(hvg.out))),])  
mysceset.hvg = newSCESet(countData = exprs(temp)[unique(c(rownames(hvg.out))),]); 
mysceset.hvg <- calculateQCMetrics(mysceset.hvg); mysceset.hvg <- sc3(mysceset.hvg, ks = 3:4,n_cores=1)
scclusters.3.hvg=pData(mysceset.hvg)[,"sc3_3_clusters"]; scclusters.4.hvg=pData(mysceset.hvg)[,"sc3_4_clusters"]
names(scclusters.3.hvg)=names(scclusters.4.hvg)=rownames(pData(mysceset.hvg))

### Hvg P1=P1; P2=P3; P3=P2 	|| 4 clusters: Hvg P1=P1; P2=P2; P3=P4; P4=P3
### Simple P1=P1; P2=P3; P3=P2	|| 4 clusters: Simple P1=P1; P2=P4; P3=P3; P4=P2

#### adjust cluters nrs. & do comparisons 
temp=scclusters.4.hvg; temp[which(scclusters.4.hvg%in%"2")]=4
temp[which(scclusters.4.hvg%in%"3")]=2; temp[which(scclusters.4.hvg%in%"4")]=3
#### compare clustering using all, M3Drop, and complete
variousclusters=cbind(scclusters.3,scclusters.3.all,scclusters.3.hvg,scclusters.4,scclusters.4.all,temp)
write.table(variousclusters,paste(outpath,"ASC-M3Drop-",myname,"-",myset,"-",mydate,".SC3.variousclusters.txt",sep=""),
	quote=F,row.names=T,col.names=T,sep="\t")
heatmap.2(variousclusters,scale="none",density.info="none",notecex=0.8,notecol="black",
	keysize = 1,trace="none",cexRow=0.6,cexCol=0.4,Rowv=F,Colv=F)
overlaps.hvg.3=table(apply(variousclusters,1,function(x) paste(x[c(1,3)],collapse="")))
overlaps.hvg.4=table(apply(variousclusters,1,function(x) paste(x[c(4,6)],collapse="")))
barplot(cbind(overlaps.hvg.3[1:3],overlaps.hvg.3[4:6],overlaps.hvg.3[7:9]),
	beside=F,col=c("darkgreen","blue","red"),names=c("1","3","2"),main="Cluster Overlaps")
barplot(cbind(c(overlaps.hvg.4[1],0,0,overlaps.hvg.4[2]),overlaps.hvg.4[8:11],c(overlaps.hvg.4[3],0,overlaps.hvg.4[4],0),
	c(overlaps.hvg.4[5],0,overlaps.hvg.4[6:7])),
	beside=F,col=c("darkgreen","red","orange","blue"),names=c("1","2.1","2.2","3"),main="Cluster Overlaps")
#### output raw for Figure
write.table(overlaps.hvg.3,paste(path,"FigS8f-raw.ClusterOverlaps.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

##### also check overlaps with marker genes
overlaps.hvg=table(apply(variousclusters,1,function(x) paste(x[c(1,3)],collapse="")))
mysceset.est.hvg <- sc3_estimate_k(mysceset.hvg)
sc3_plot_consensus( mysceset.hvg, k = 3, show_pdata = c("cell_id","log10_total_features","sc3_3_clusters", "sc3_3_log2_outlier_score"))
sc3_plot_expression(mysceset.hvg, k = 3, show_pdata = c("cell_id", "log10_total_features","sc3_3_clusters", "sc3_3_log2_outlier_score"))
sc3_plot_cluster_stability(mysceset.hvg, k = 3); sc3_plot_cluster_stability(mysceset.hvg, k = 4)
##### get the clusters from SC3
sc3.3.hvg=filcols[scclusters.3.hvg]; names(sc3.3.hvg)=names(cell_populations_3)
sc3.4.hvg=filcols[scclusters.4.hvg]; names(sc3.4.hvg)=names(cell_populations_4)
temp=convertTo(sce, type="monocle")
marker_genes_sc_3.hvg <- M3DropGetMarkers(exprs(temp), scclusters.3.hvg)
marker_genes_sc_4.hvg <- M3DropGetMarkers(exprs(temp), scclusters.4.hvg)

#### Marker table: overlaps betwee top 100 markers
nmarkers=100
#### 3 clusters
p1.m=rownames(head(marker_genes_sc_3[marker_genes_sc_3$Group==1,],nmarkers))
p1.m.hvg=rownames(head(marker_genes_sc_3.hvg[marker_genes_sc_3.hvg$Group==1,],nmarkers))
p3.m=rownames(head(marker_genes_sc_3[marker_genes_sc_3$Group==2,],nmarkers))
p3.m.hvg=rownames(head(marker_genes_sc_3.hvg[marker_genes_sc_3.hvg$Group==2,],nmarkers))
p2.m=rownames(head(marker_genes_sc_3[marker_genes_sc_3$Group==3,],nmarkers))
p2.m.hvg=rownames(head(marker_genes_sc_3.hvg[marker_genes_sc_3.hvg$Group==3,],nmarkers))
#### 4 clusters
p1.m.4=rownames(head(marker_genes_sc_4[marker_genes_sc_4$Group==1,],nmarkers))
p1.m.hvg.4=rownames(head(marker_genes_sc_4.hvg[marker_genes_sc_4.hvg$Group==1,],nmarkers))
p3.m.4=rownames(head(marker_genes_sc_4[marker_genes_sc_4$Group==3,],nmarkers))
p3.m.hvg.4=rownames(head(marker_genes_sc_4.hvg[marker_genes_sc_4.hvg$Group==4,],nmarkers))
p2.m.4=rownames(head(marker_genes_sc_4[marker_genes_sc_4$Group==4,],nmarkers))
p2.m.hvg.4=rownames(head(marker_genes_sc_4.hvg[marker_genes_sc_4.hvg$Group==2,],nmarkers))
p4.m.4=rownames(head(marker_genes_sc_4[marker_genes_sc_4$Group==2,],nmarkers))
p4.m.hvg.4=rownames(head(marker_genes_sc_4.hvg[marker_genes_sc_4.hvg$Group==3,],nmarkers))
###### save the markers for future reference
#save(p1.m,p1.m.hvg,p2.m,p2.m.hvg,p3.m,p3.m.hvg,p1.m.4,p1.m.hvg.4,p3.m.4,p3.m.hvg.4,p2.m.4,p2.m.hvg.4,p4.m.4,p4.m.hvg.4,
#	file=paste(outpath,"mASC-M3Drop-",myname,"-",myset,"-",mydate,".208.Markers.RData",sep=""))
##### load markers from previous analysis
load(paste(outpath,"mASC-M3Drop-",myname,"-",myset,"-",mydate,".208.Markers.RData",sep=""))
markeroverlap=c(length(p1.m[which(p1.m%in%p1.m.hvg)]),length(p2.m[which(p2.m%in%p2.m.hvg)]),length(p3.m[which(p3.m%in%p3.m.hvg)]))
barplot(markeroverlap,names=c("P1","P2","P3"),col=c("green","red","blue"),main="Top 100 Marker overlap")
#### raw output for figure
names(markeroverlap)=c("P1","P2","P3")
write.table(markeroverlap,file=paste(path,"Figs8g-raw.MarkeroverlapC1.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

markeroverlap=c(length(p1.m[which(p1.m%in%p1.m.hvg)]),length(p2.m[which(p2.m%in%p2.m.hvg)]),
	length(p3.m[which(p3.m%in%p3.m.hvg)]),length(p4.m[which(p4.m%in%p4.m.hvg)]),length(p4.m[which(p4.m%in%p1.m.hvg)]))
barplot(markeroverlap,names=c("P1","P2","P3","P2.2","P2.2w1"),col=c("green","red","blue","orange","orange"),
	main="Top 100 Marker overlap")

###### Make t-SNE plots with HVG genes instead of M3Drop genes
mymatrix.log.hvg=log1p(exprs(temp)[unique(c(rownames(hvg.out))),])
#rtsne_out.hvg <- Rtsne(as.matrix(t(mymatrix.log.hvg)),theta=0.0000001,perplexity=as.integer(dim(mymatrix.log.hvg)[2]/6))
#rtsne_out.hvg1830=rtsne_out.hvg
#save(rtsne_out.hvg,rtsne_out.hvg1830,marker_genes_sc_4,marker_genes_sc_3,
#	file=paste(outpath,"QualPlots-M3Drop-",myname,"-",myset,"-",mydate,".RtSNE.1830.208.RData",sep=""))
load(paste(outpath,"QualPlots-M3Drop-",myname,"-",myset,"-",mydate,".RtSNE.1830.208.RData",sep=""))
mymatrix.withlab.hvg=mymatrix.log.hvg
rownames(mymatrix.withlab.hvg)[grep("ENS",rownames(mymatrix.withlab.hvg))]=elementMetadata(genes84.GR[rownames(mymatrix.withlab.hvg)[grep("ENS",rownames(mymatrix.withlab.hvg))]])$geneName
toplot=mymatrix.withlab.big
for (i in 1:length(colvector)) {
	plot(rtsne_out.hvg$Y,pch=19, main=paste("Cats - ",names(colvector)[i]),col=colvector[[i]][colnames(mymatrix.withlab.hvg)],xlab="dim 1",ylab="dim 2")
#	text(rtsne_out.hvg$Y, labels=colnames(mymatrix.withlab.hvg),col=colvector[[i]][colnames(mymatrix.withlab.hvg)],cex=0.2)
}
ccols.hvg=ccols;ccols.hvg[[5]]=sc3.3.hvg;ccols.hvg[[6]]=sc3.4.hvg
for (i in 1:length(ccols.hvg)) {
	plot(rtsne_out.hvg$Y,pch=19, main=paste("Clusters - ",names(ccols.hvg)[i]),col=ccols.hvg[[i]][colnames(mymatrix.withlab.hvg)],xlab="dim 1",ylab="dim 2")
#	text(rtsne_out.hvg$Y, labels=colnames(mymatrix.withlab.hvg),col=ccols.hvg[[i]][colnames(mymatrix.withlab.hvg)],cex=0.2)
}
#par(mfrow=c(3,3)); 
mymax=max(mymatrix.withlab.big)*2/3
plotMarker=unique(c(myMarkerAll,mychoiceMarkerAll,allmarkers,rownames(mymatrix.withlab3)))
plotMarker=sort(plotMarker[which(plotMarker%in%rownames(toplot))])
choosecollist=list(c("ivory2","darkgrey", "black"),c("ivory2","orange", "red"),c("ivory2","lightblue", "blue"),c("ivory2","lightgreen", "darkgreen"))
for (j in 1:length(choosecollist)) {
for (i in 1:length(plotMarker)) {
        choosecols=choosecollist[[j]]; plotcol=valuesToColorsAbs(toplot[plotMarker[i],],choosecols,colnames(mymatrix.withlab),0,mymax)
        plot(rtsne_out.hvg$Y,pch=19, main=plotMarker[i],col=plotcol,xlab="dim 1",ylab="dim 2")
#        text(rtsne_out.hvg$Y, labels=colnames(mymatrix.withlab),col=plotcol,cex=0.2)	
}}
par(mfrow=c(1,1))
# heatmap
makeprefheatmaprlog(mymatrix.withlab.hvg,spikedgenes84.GR,colvals=colvector[[9]],mainname="row.scaled", rowvals="grey")
makeprefheatmaprlog(mymatrix.withlab.hvg,spikedgenes84.GR,colvals=ccols[[1]],mainname="row.scaled", rowvals="grey")
makeprefheatmaprlog(mymatrix.withlab.hvg,spikedgenes84.GR,colvals=ccols[[3]],mainname="row.scaled", rowvals="grey")
makeprefheatmaprlog(mymatrix.withlab.hvg,spikedgenes84.GR,colvals=ccols.hvg[[5]],mainname="row.scaled", rowvals="grey")

#### output raw for Figure 
#### add batch and expression values for all these genes; RFP+ and Neg.
geneIds=c("Fabp4","Cd34","Zfp423","Prdm16","Ppargc1a","Ebf2","Cd24a","Pdgfra","Pdgfrb","Adipoq","Wt1","Dlk1","Pparg","chrtdRFP")
figs8h=cbind(rtsne_out.hvg$Y,colvector[[1]][colnames(mymatrix.withlab)],
	ccols.hvg[[3]][colnames(mymatrix.withlab)],ccols.hvg[[5]][colnames(mymatrix.withlab)],
	t(toplot[geneIds,]))
colnames(figs8h)[1:5]=c("x","y","Batch","SC3.3","SC3.HVG.3")
write.table(figs8h,paste(path,"FigS8h-raw.tSNEC1.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
dev.off()


############################################################################################################################
#### Run GO on the DE genes #####
P1.top.sc3.10=head(marker_genes_sc_3[marker_genes_sc_3$Group==1,],10); P3.top.sc3.10=head(marker_genes_sc_3[marker_genes_sc_3$Group==2,],10)
P2.top.sc3.10=head(marker_genes_sc_3[marker_genes_sc_3$Group==3,],10)
P1.top.sc3.TFs=head(marker_genes_sc_3[marker_genes_sc_3$Group==1,],100); P3.top.sc3.TFs=head(marker_genes_sc_3[marker_genes_sc_3$Group==2,],100)
P2.top.sc3.TFs=head(marker_genes_sc_3[marker_genes_sc_3$Group==3,],100)

ofinterest= unlist(c(GENETYPES84[[1]],GENETYPES84[[20]]))
P3.top.sc3.TFs=P3.top.sc3.TFs[which(rownames(P3.top.sc3.TFs)%in%ofinterest),]
P2.top.sc3.TFs=P2.top.sc3.TFs[which(rownames(P2.top.sc3.TFs)%in%ofinterest),]
P1.top.sc3.TFs=P1.top.sc3.TFs[which(rownames(P1.top.sc3.TFs)%in%ofinterest),]
mymatrix.withlabTFs=log1p(Normalized_data$data[unique(c(rownames(P3.top.sc3.10),rownames(P2.top.sc3.10),rownames(P1.top.sc3.10),
	rownames(P3.top.sc3.TFs),rownames(P2.top.sc3.TFs),rownames(P1.top.sc3.TFs)),
	getEnsgID(spikedgenes84.GR,c("Cd55","Abcg1","F3","Adam12","Aoc3","Il13ra1"))),])
rownames(mymatrix.withlabTFs)[grep("ENS",rownames(mymatrix.withlabTFs))]=elementMetadata(spikedgenes84.GR[rownames(mymatrix.withlabTFs)[grep("ENS",rownames(mymatrix.withlabTFs))]])$geneName


pdf(paste(outpath,"QualPlots-M3Drop-",myname,"-",myset,"-",mydate,".exprOverviews.pdf",sep=""),height=10,width=10)
par(mfrow=c(3,3))
#### distribution of marker genes 
temp=log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Actb","Ptprc","Pecam1","Cd34","Ly6a","Itgb1")),])
rownames(temp)[grep("ENS",rownames(temp))]=elementMetadata(genes84.GR[rownames(temp)[grep("ENS",rownames(temp))]])$geneName
temp=temp[c("Actb","Ptprc","Pecam1","Cd34","Ly6a","Itgb1"),]
boxplot(lapply(1:6,function(x) temp[x,]),names=rownames(temp), col=c("snow3","black","black","salmon2","salmon2","salmon2"),las=2,ylab=c("Log. CPM"))
#### output raw for Figure
write.table(t(temp1),paste(path,"FigS1e-raw.C1Marker.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

#### Correlations with Fabp4
mycors.p=c(sapply(c("Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ppargc1a","Ebf2","Adipoq","Wt1"),
	function(x) cor.test(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Fabp4")),]),
		log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c(x))[1],]),method="spearman")$p.val),
		cor.test(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Fabp4")),]),log1p(countmatrix.k.cpm["chrtdRFP",]))$p.val)
mycors=c(sapply(c("Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ppargc1a","Ebf2","Adipoq","Wt1"),
		function(x) cor(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Fabp4")),]),
		log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c(x))[1],]),method="spearman")),
		cor(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Fabp4")),]),log1p(countmatrix.k.cpm["chrtdRFP",])))
names(mycors.p)[12]=names(mycors)[12]="RFP"

cor.table=cbind(mycors.p,mycors)
write.table(cor.table,
	paste(outpath,"QualPlots-M3Drop-",myname,"-",myset,"-",mydate,".corTable.txt",sep=""),row.names=T,col.names=T,sep="\t",quote=F)
#### output raw for Figure
cor.table=cbind(p.adjust(mycors.p),mycors)
write.table(cor.table,paste(path,"FigS1k-raw.C1Correlations.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
barplot(sort(-log10(mycors.p),decreasing=T),main="-logp.Cor.Fabp4",las=2)
barplot(mycors[names(sort(-log10(mycors.p),decreasing=T))],main="Cor.Fabp4",las=2)

#### and correlation with Cd34 #"Fabp4","Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ebf2","Adipoq"
mycors.p.cd34=c(sapply(c("Fabp4","Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ppargc1a","Ebf2","Adipoq","Wt1"),
	function(x) cor.test(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Cd34")),]),
		log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c(x))[1],]),method="spearman")$p.val),
		cor.test(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Cd34")),]),log1p(countmatrix.k.cpm["chrtdRFP",]))$p.val)
mycors.cd34=c(sapply(c("Fabp4","Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ppargc1a","Ebf2","Adipoq","Wt1"),
		function(x) cor(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Cd34")),]),
		log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c(x))[1],]),method="spearman")),
		cor(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Cd34")),]),log1p(countmatrix.k.cpm["chrtdRFP",])))
names(mycors.p.cd34)[12]=names(mycors)[12]="RFP"
cor.table.cd34=cbind(mycors.p.cd34,mycors.cd34)
write.table(cor.table.cd34,
	paste(outpath,"QualPlots-M3Drop-",myname,"-",myset,"-",mydate,".corTable.Cd34.18.txt",sep=""),row.names=T,col.names=T,sep="\t",quote=F)



##### Adipogenic, stemness etc. score
adiposcore.cats=lapply(c(1,3,2),function(x) adiposcore[names(cID[which(cID%in%x)])]) #1 to 2 0.1391, 1 to 3  3.921e-12, 2 to 3 1.686e-05
stemscore.cats=lapply(c(1,3,2),function(x) stemscore[names(cID[which(cID%in%x)])]) #1 to 2 0.001402, 1 to 3   1.889e-06, 2 to 3 0.3144
osteoscore.cats=lapply(c(1,3,2),function(x) osteoscore[names(cID[which(cID%in%x)])]) #1 to 2 0.001402, 1 to 3   1.889e-06, 2 to 3 0.3144
endoscore.cats=lapply(c(1:3),function(x) endoscore[names(cID[which(cID%in%x)])]) #1 to 2 0.0536, 1 to 3   0.05227, 2 to 3 0.005642
immunoscore.cats=lapply(c(1:3),function(x) immunoscore[names(cID[which(cID%in%x)])]) #1 to 2 0.02373, 1 to 3  0.08961, 3 to 2  0.1672

adipo.pval.c1=c(wilcox.test(adiposcore.cats[[1]],adiposcore.cats[[2]])$p.val,wilcox.test(adiposcore.cats[[1]],adiposcore.cats[[3]])$p.val,
	wilcox.test(adiposcore.cats[[2]],adiposcore.cats[[3]])$p.val)
names(adipo.pval.c1)=c("P1-P2","P1-P3","P2-P3")
#### output raw for Figure
write.table(adipo.pval.c1,paste(path,"FigS1i-raw.AdipoScorePvals.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

cd55.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Cd55")),names(cID[which(cID%in%x)])])) #1 to 2 0.1391, 1 to 3  3.921e-12, 2 to 3 1.686e-05
il13ra1.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Il13ra1")),names(cID[which(cID%in%x)])])) #1 to 2 0.1391, 1 to 3  3.921e-12, 2 to 3 1.686e-05
f3.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("F3")),names(cID[which(cID%in%x)])])) #1 to 2 0.1391, 1 to 3  3.921e-12, 2 to 3 1.686e-05
abcg1.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Abcg1")),names(cID[which(cID%in%x)])])) #1 to 2 0.1391, 1 to 3  3.921e-12, 2 to 3 1.686e-05
adam12.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Adam12")),names(cID[which(cID%in%x)])])) #1 to 2 0.1391, 1 to 3  3.921e-12, 2 to 3 1.686e-05
aoc3.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Aoc3")),names(cID[which(cID%in%x)])])) #1 to 2 0.1391, 1 to 3  3.921e-12, 2 to 3 1.686e-05
adipoq.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Adipoq")),names(cID[which(cID%in%x)])]))
retn.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Retn")),names(cID[which(cID%in%x)])])) 
cd34.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Cd34")),names(cID[which(cID%in%x)])])) 
sca1.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Ly6a")),names(cID[which(cID%in%x)])])) 
cd29.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Itgb1")),names(cID[which(cID%in%x)])])) 


#### output raw for Figure #83 96 29
temp1=unlist(adiposcore.cats); nrel=sapply(adiposcore.cats,function(x) length(x))
tempcat=c(rep("P1",nrel[1]),rep("P2",nrel[2]),rep("P3",nrel[3]))
temp1=cbind(temp1,unlist(stemscore.cats),unlist(cd55.cats),
unlist(il13ra1.cats),unlist(f3.cats),unlist(abcg1.cats),unlist(adam12.cats),
unlist(aoc3.cats),unlist(adipoq.cats),unlist(retn.cats),
unlist(cd34.cats),unlist(sca1.cats),unlist(cd29.cats),tempcat);
colnames(temp1)=c("AdipoScore","StemScore","CD55","Il13ra1","F3","Abcg1","Adam12","Aoc3","Adipoq","Retn","Cd34","Sca1","Cd29","Cluster")
write.table(temp1,paste(path,"FigS12S2-raw.AdipoScore.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

#### and p-values 
markerlist=list(cd55.cats,il13ra1.cats,f3.cats,abcg1.cats,adam12.cats,aoc3.cats,adipoq.cats,retn.cats,
cd34.cats,sca1.cats,cd29.cats)
markerpval=matrix(data=NA, ncol=3,nrow=length(markerlist))
for (i in 1:length(markerlist)) {
	markerpval[i,]=c(p.adjust(wilcox.test(markerlist[[i]][[1]],markerlist[[i]][[2]])$p.val,method="BH"),
		p.adjust(wilcox.test(markerlist[[i]][[1]],markerlist[[i]][[3]])$p.val,method="BH"),
		p.adjust(wilcox.test(markerlist[[i]][[2]],markerlist[[i]][[3]])$p.val,method="BH"))
}
colnames(markerpval)=c("P1-P2","P1-P3","P2-P3"); rownames(markerpval)=c("Cd55","Il13ra1","F3","Abcg1", "Adam12","Aoc3","Adipoq","Retn","Cd34","Sca1","Cd29")
#### output raw for Figure
write.table(markerpval,paste(path,"Fig2aS2b-raw.MarkerPvals.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

#### positive vs. negative (data not shown in manuscript); Dlk1/RFP analysis
par(mfrow=c(3,6))
PvsN.cats=sapply(c(1,3,2),function(x) c(length(names(cID[which(cID%in%x)])[which(names(cID[which(cID%in%x)])%in%colnames(countmatrix.k.cpm[,which(cellinfo.k[,"Cell"]%in%"R")]))]),
	length(names(cID[which(cID%in%x)])[which(names(cID[which(cID%in%x)])%in%colnames(countmatrix.k.cpm[,which(cellinfo.k[,"Cell"]%in%"N")]))])))
#### make beanplot with the expression of different genes per RED vs. SNOW
lapply(c("Fabp4","Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ppargc1a","Ebf2","Adipoq","Wt1"),
	function(y) boxplot(lapply(c("R","N"),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,y)[1],
	which(cellinfo.k[,"Cell"]%in%x)])),col=c("red","snow3"),names=c("R","N"),main=y))
### and RFP
boxplot(lapply(c("R","N"),function(x) log1p(countmatrix.k.cpm["chrtdRFP",
	which(cellinfo.k[,"Cell"]%in%x)])),col=c("red","snow3"),names=c("R","N"),main="RFP") ### wilcox p: 3.901e-06
sepvals=lapply(c("Fabp4","Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ppargc1a","Ebf2","Adipoq","Wt1"),
	function(y) lapply(c("R","N"),function(x) as.numeric(log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,y)[1],
	which(cellinfo.k[,"Cell"]%in%x)]))))
pvals=sapply(sepvals,function(x) wilcox.test(x[[1]],x[[2]], alternative="greater")$p.val)
snames=sapply(c("Fabp4","Pdgfra","Pdgfrb","Pparg","Zfp423","Cd24a","Dlk1","Prdm16","Ppargc1a","Ebf2","Adipoq","Wt1"),
	function(y) getEnsgID(spikedgenes84.GR,y)[1])
names(pvals)=elementMetadata(spikedgenes84.GR[snames])$geneName
pvals=sapply(pvals,function(x) min(x,1-x))
barplot(apply(PvsN.cats,2,function(x)x/sum(x)),beside=F,names=c("P1","P2","P3"),col=c("red","snow3"),las=2)

##### plot the scores
par(mfrow=c(3,3))
barplot(-log10(pvals),main="Wilcox.p",las=2)
beancols=list(c("darkgreen",rep("black",3)),c("salmon3",rep("black",3)),c("blue",rep("black",3)))
beanplot(adiposcore.cats,col=beancols,main="Adipo");beanplot(stemscore.cats,col=beancols,main="Stem")
beanplot(osteoscore.cats,col=beancols,main="Osteo");beanplot(endoscore.cats,col=beancols,main="Endo")
beanplot(immunoscore.cats,col=beancols,main="Immuno")
barplot(PvsN.cats,beside=F,names=c("P1","P2","P3"),col=c("red","snow3"))
beanplot(cd55.cats,col=beancols,main="CD55");beanplot(il13ra1.cats,col=beancols,main="IL13RA1")
beanplot(f3.cats,col=beancols,main="F3");beanplot(abcg1.cats,col=beancols,main="ABCG1")
beanplot(adam12.cats,col=beancols,main="ADAM12");beanplot(aoc3.cats,col=beancols,main="AOC3")


### Output the heatmap with marker gene expression
distance = as.dist(1-cor(t(mymatrix.withlabTFs),method="spearman")); cluster = hclust(distance, method = "ward.D2")
distance.Col = dist(t(mymatrix.withlabTFs)); cluster.Col = hclust(distance, method = "ward.D2")
rowvals=valuesToColors(apply(mymatrix.withlabTFs,1,function(x) median(x)),c("ivory","darkred"),names(apply(mymatrix.withlabTFs,1,function(x) median(x))))

heatmap.2(mymatrix.withlabTFs[,names(sort(cID))],scale="none",col=rev(colorRampPalette(c("red","orange","ivory","lightblue","darkblue"))(216)),
		Colv = FALSE, Rowv = as.dendrogram(cluster), density.info="none",notecex=0.8,notecol="black",keysize = 1,trace="none",cexRow=0.8,cexCol=0.4,
                ColSideColors=ccols[[3]][names(sort(cID))],RowSideColors=rowvals[rownames(mymatrix.withlabTFs)])
heatmap.2(mymatrix.withlabTFs[,names(sort(cID))],scale="none",col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
		Colv = FALSE, Rowv = as.dendrogram(cluster), density.info="none",notecex=0.8,notecol="black",keysize = 1,trace="none",cexRow=0.8,cexCol=0.4,
                ColSideColors=ccols[[3]][names(sort(cID))],RowSideColors=rowvals[rownames(mymatrix.withlabTFs)])
heatmap.2(mymatrix.withlabTFs[,names(sort(cID))],scale="row",col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
		Colv = FALSE, Rowv = as.dendrogram(cluster), density.info="none",notecex=0.8,notecol="black",keysize = 1,trace="none",cexRow=0.8,cexCol=0.4,
                ColSideColors=ccols[[3]][names(sort(cID))],RowSideColors=rowvals[rownames(mymatrix.withlabTFs)])

#### output raw for Figure #83 96 29
write.table(mymatrix.withlabTFs[,names(sort(cID))],
	paste(path,"Fig1c-raw.C1Heatmap.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

allgenes.background=rownames(Normalized_data$data)
bigcut=0.000006
P1.top.sc3.GO=head(marker_genes_sc_3[marker_genes_sc_3$Group==1,],200);#P1.top.sc3.GO= P1.top.sc3.GO[which( P1.top.sc3.GO[,3]<bigcut),]
P3.top.sc3.GO=head(marker_genes_sc_3[marker_genes_sc_3$Group==2,],200);#P3.top.sc3.GO= P3.top.sc3.GO[which( P3.top.sc3.GO[,3]<bigcut),]
P2.top.sc3.GO=head(marker_genes_sc_3[marker_genes_sc_3$Group==3,],200);#P2.top.sc3.GO= P2.top.sc3.GO[which( P2.top.sc3.GO[,3]<bigcut),]

###### run the Gene Ontology analysis on genes specific for one population (not included in .man)
golist.sign=list(P1.top.sc3.GO,P2.top.sc3.GO,P3.top.sc3.GO)
names(golist.sign)=c("P1","P2","P3")
SC3.golist=list(NA)
myp=0.001;showgraph=10; 
mainname="ASC"
for (i in 1:length(golist.sign)){
	myIDs=rownames(golist.sign[[i]])#names(spikedgenes84.GR[which(elementMetadata(spikedgenes84.GR)$geneName%in%rownames(golist.sign[[i]][which(golist.sign[[i]][,2]>1),]))])
	if (length(myIDs)>10) {
		SC3.golist[[i]]=runGO.84(myIDs,allgenes.background,outpath,ensembl84,"BP","elimCount",myp,paste(mainname,"GO-pcut-",i,sep="_"),ontoGOs.84,showgraph)
	}
	}
dev.off()

###### Test additional genes per population 
pdf(paste(outpath,"ExprPerPop-M3Drop-",myname,"-",myset,"-",mydate,".exprOverviews.xTra.pdf",sep=""),height=10,width=10)
par(mfrow=c(3,3))
cd44.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Cd44")),names(cID[which(cID%in%x)])])) #1 to 2 0.1391, 1 to 3  3.921e-12, 2 to 3 1.686e-05
il13ra1.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Il13ra1")),names(cID[which(cID%in%x)])])) #1 to 2 0.1391, 1 to 3  3.921e-12, 2 to 3 1.686e-05
il33.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Il33")),names(cID[which(cID%in%x)])])) #1 to 2 0.1391, 1 to 3  3.921e-12, 2 to 3 1.686e-05
il13ra2.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Il13ra2")),names(cID[which(cID%in%x)])])) #1 to 2 0.1391, 1 to 3  3.921e-12, 2 to 3 1.686e-05
####
postn.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Postn")),names(cID[which(cID%in%x)])])) #1 to 2 0.0008842, 1 to 3  5.102e-05, 2 to 3 9.407e-08
fap.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Fap")),names(cID[which(cID%in%x)])])) #1 to 2 0.1398, 1 to 3  0.03143, 2 to 3 0.00175
col1a1.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("Col1a1")),names(cID[which(cID%in%x)])])) #1 to 2 1.819e-06, 1 to 3   0.0001874, 2 to 3  0.3419
fsp1.cats=lapply(c(1,3,2),function(x) log1p(countmatrix.k.cpm[getEnsgID(spikedgenes84.GR,c("S100a4")),names(cID[which(cID%in%x)])])) #1 to 2 0.3204, 1 to 3 0.3963, 2 to 3 0.9407

beanplot(cd44.cats,col=beancols,main="CD44"); boxplot(cd44.cats,col=boxcols,main="CD44")
beanplot(il13ra1.cats,col=beancols,main="Il13RA1"); boxplot(il13ra1.cats,col=boxcols,main="Il13RA1")
beanplot(il13ra2.cats,col=beancols,main="IL13RA2"); boxplot(il13ra2.cats,col=boxcols,main="IL13RA2")
beanplot(il33.cats,col=beancols,main="IL33"); boxplot(il33.cats,col=boxcols,main="IL33")

beanplot(postn.cats,col=beancols,main="POSTN"); boxplot(postn.cats,col=boxcols,main="POSTN")
beanplot(fap.cats,col=beancols,main="FAP"); boxplot(fap.cats,col=boxcols,main="FAP")
beanplot(col1a1.cats,col=beancols,main="COL1A1"); boxplot(col1a1.cats,col=boxcols,main="COL1A1")
beanplot(fsp1.cats,col=beancols,main="FSP1"); boxplot(fsp1.cats,col=boxcols,main="FSP1")
dev.off()

##### write out the cluster markers used for GO analysis #######
for (i in 1:length(SC3.golist)){
	write.table(SC3.golist[[i]],paste(outpath,myname,"-",myset,"-",mydate,".P",i,".annotation.txt",sep=""),
row.names=T,col.names=T,sep="\t",quote=F) }

##### write out the matrix #######
write.table(countmatrix[,colnames(Normalized_data$data)],paste(outpath,myname,"-",myset,"-",mydate,".counts.STAR.txt",sep=""),
	row.names=T,col.names=T,sep="\t",quote=F)
mybatch=substr(colnames(countmatrix[,colnames(Normalized_data$data)]),3,4)
names(mybatch)=colnames(countmatrix[,colnames(Normalized_data$data)])
write.table(mybatch,paste(outpath,myname,"-",myset,"-",mydate,".batches.txt",sep=""),
	row.names=T,col.names=T,sep="\t",quote=F)
	
##### write out the marker genes 
write.table(rownames(P1.top.sc3.GO),paste(path,"P1MarkerGenes-M3Drop-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(rownames(P2.top.sc3.GO),paste(path,"P2MarkerGenes-M3Drop-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(rownames(P3.top.sc3.GO),paste(path,"P3MarkerGenes-M3Drop-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(allgenes.background,paste(path,"BackgroundGenes-M3Drop-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
### IDS
write.table(elementMetadata(spikedgenes84.GR[rownames(P1.top.sc3.GO)])$geneName,paste(path,"P1MarkerGenes-M3Drop-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(spikedgenes84.GR[rownames(P2.top.sc3.GO)])$geneName,paste(path,"P2MarkerGenes-M3Drop-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(spikedgenes84.GR[rownames(P3.top.sc3.GO)])$geneName,paste(path,"P3MarkerGenes-M3Drop-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(spikedgenes84.GR[allgenes.background])$geneName,paste(path,"BackgroundGenes-M3Drop-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),
	row.names=F,col.names=F,quote=F)

#### check batch vs. cluster
scclusters.3.hvg
batchm=cbind(scclusters.4,mybatch)
batchm.tab=table( apply(batchm,1,function(x) paste(x,collapse="")))
pdf(paste(outpath,"Batchplot-M3Drop-",myname,"-",myset,"-",mydate,".exprOverviews.xTra.pdf",sep=""),height=10,width=10)
par(mfrow=c(3,3))
barplot(cbind(batchm.tab[1:3],batchm.tab[4:6],batchm.tab[7:9],batchm.tab[10:12]),beside=F,
	main="Cluster.Batch",col=filcols[1:3],names=c("P1","P4","P3","P2"))
dev.off()

##### write out the DE genes ################
write.table(log1p(Normalized_data$data[rownames(DE_genes),]),paste(outpath,myname,"-",myset,"-",mydate,".norm.DEgenes.STAR.txt",sep=""),
	row.names=T,col.names=T,sep="\t",quote=F)

##### write out the clustering result #######
write.table(cID,paste(outpath,myname,"-",myset,"-",mydate,".clusters.txt",sep=""),
	row.names=T,col.names=T,sep="\t",quote=F)

	
#### Brbseq & Pop analysis
load("adipo/ETH/mapped/brbseq/STAR/plots/Limma-Nov16.BrbseqCounts.RData")
ASC.P3=Normalized_data$data[,names(cID[which(cID%in%"2")])]
ASC.P3.counts=countmatrix.k.cpm[,names(cID[which(cID%in%"2")])]
Brb.P3=normcounts.brbseq[,colnames(normcounts.brbseq)[grep("F3p",colnames(normcounts.brbseq))][grep("[.]0",colnames(normcounts.brbseq)[grep("F3p",colnames(normcounts.brbseq))])]]
Brb.P3.counts=fullCountTable.brbseq[,colnames(fullCountTable.brbseq)[grep("F3p",colnames(fullCountTable.brbseq))][grep("[.]0",colnames(fullCountTable.brbseq)[grep("F3p",colnames(fullCountTable.brbseq))])]]
Brb.ASC=normcounts.brbseq[,colnames(normcounts.brbseq)[grep("DP",colnames(normcounts.brbseq))][grep("[.]0",colnames(normcounts.brbseq)[grep("DP",colnames(normcounts.brbseq))])]]

okrows=rownames(Brb.P3)[which(rownames(Brb.P3)%in%rownames(ASC.P3))]
ASC.P3=ASC.P3[okrows,];Brb.P3=ASC.P3[okrows,]
ASC.P3.mean=log1p(apply(ASC.P3,1,function(x) mean(x)))
ASC.P3.counts.mean=log1p(apply(ASC.P3.counts,1,function(x) mean(x)))

pdf(paste(outpath,"QualPlots-M3Drop-",myname,"-",myset,"-",mydate,".PopPlots.Pub.pdf",sep=""),height=10,width=10)
par(mfrow=c(3,3))
mycors.P3=sapply(1:4,function(x) cor(log1p(Brb.P3[,x]),ASC.P3.mean,method="pearson"))
beanplot(mycors.P3,main="Brbseq vs. SC",ylab="Pearson Cor")
mycors.P3=sapply(1:4,function(x) cor(log1p(Brb.P3[,x]),ASC.P3.mean,method="spearman"))
beanplot(mycors.P3,main="Brbseq vs. SC",ylab="Spearman Cor")
#function(a,b,namea,nameb,mymethod) 
sapply(1:4,function(x) scatterSmoothPlot(log1p(Brb.P3[,x]),ASC.P3.mean,"Brbseq","SingleCell","pearson"))
sapply(1:4,function(x) scatterSmoothPlot(log1p(Brb.P3[,x]),ASC.P3.mean,"Brbseq","SingleCell","spearman"))

### raw for Figure
write.table(mycors.P3,paste(path,"FigS4c-raw.BrbseqSpearman.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F) #74 71 63
mytab=cbind(log1p(Brb.P3[,x]),ASC.P3.mean); colnames(mytab)=c("Brbseq","ASC.P3.Mean")
write.table(mytab,paste(path,"FigS4c-raw.BrbseqSpearmanExample.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F) #74 71 63


#scatterRegPlot
beanplot(apply(Brb.P3.counts,2,function(x) length(which(x>0))),main="Brbseq.P3",xlab="Nr.genes")
barplot(length(which(ASC.P3.counts.mean>0)),main="ASC.P3",xlab="Nr. genes")

#### cors with population data
popbrbgenes=rownames(Brb.ASC)[which(rownames(Brb.ASC)%in%rownames(countmatrix.pop.cpm))]
mycors.ASC=apply(log1p(countmatrix.pop.cpm[popbrbgenes,c(19:22)]),2,function(y)
	sapply(1:6,function(x) cor(log1p(Brb.ASC[popbrbgenes,x]),y,method="pearson")))
beanplot(mycors.P3,main="Brbseq vs. SC",ylab="Pearson Cor")
mycors.P3=sapply(1:4,function(x) cor(log1p(Brb.P3[,x]),ASC.P3.mean,method="spearman"))
beanplot(mycors.P3,main="Brbseq vs. SC",ylab="Spearman Cor")
#function(a,b,namea,nameb,mymethod) 
sapply(1:4,function(x) scatterSmoothPlot(log1p(Brb.P3[,x]),ASC.P3.mean,"Brbseq","SingleCell","pearson"))
sapply(1:4,function(x) scatterSmoothPlot(log1p(Brb.P3[,x]),ASC.P3.mean,"Brbseq","SingleCell","spearman"))


### Single-cell correlations (as above)
allcors=lapply(1:3,function(x) sapply(19:22, function(y)
	scatterSmoothPlot(log1p(temp2.counts[,x]),log1p(countmatrix.pop.cpm[,y]),x,y,"spearman")))
names(allcors)=c("P1","P2","P3")
beanplot(allcors,main="SC vs. Pop",ylab="Spearman cor")
allcors.SC=lapply(seq(1,200,by=30),function(x) sapply(seq(2,200,by=30), function(y)
	scatterSmoothPlot(log1p(countmatrix.k.cpm[,x]),log1p(countmatrix.k.cpm[,y]),x,y,"spearman")))
allcors.SC.ERCC=lapply(seq(1,200,by=30),function(x) sapply(seq(2,200,by=30), function(y)
	scatterSmoothPlot(log1p(countmatrix.k.cpm[grep("ERCC-",rownames(countmatrix.k.cpm)),x]),log1p(countmatrix.k.cpm[grep("ERCC-",rownames(countmatrix.k.cpm)),y]),x,y,"spearman")))
beanplot(list(unlist(allcors.SC),unlist(allcors.SC.ERCC)),main="ASC-cors",names=c("All","ERCC"),ylab="Spearman cor")

dev.off()







