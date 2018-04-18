#bsub -M 4000000 -R  "rusage[mem=4000]" -Is bash
source("plot-util.R")
source("util-asc.R")
library(beanplot)

#### Initialise path #####
#whichlib="50"
whichlib="6"

path=paste("mapped/hsabrbseq/Demultiplexed_Lib",whichlib,"/",sep="")
#path="mapped/hsabrbseq/Lib50ng/"
#mynames=paste("s",1:37,sep="");names(mynames)=mynames
myname="Jan17";  
filcols=brewer.pal(9,"Pastel1")
myset="bow-mmseqEnsg84-runCut-20170110"
mmset=myset
mydate="20170110"


#### Initialise sample details #####
mynr=12;
mynames=c("i04_T0_F3-","i04_T0_F3+","i04_T0_Total","i05_T0_F3-","i05_T0_F3+",
	"i05_T0_Total","i06_T0_F3-" ,"i06_T0_F3+" ,"i06_T0_Total" ,"i08_T0_F3-" , 
	"i08_T0_F3+" ,"i08_T0_Total")
names(mynames)=mynames
mynames.touse=mynames

##### prepare colors
f3Pblues=brewer.pal(9,"Blues");f3Pgreens=brewer.pal(9,"Greens"); f3Mreds=brewer.pal(9,"Reds"); f3Moranges=brewer.pal(9,"Oranges");CtrlMgreys=brewer.pal(9,"Greys")
mycolors=mysplit(mynames,"T0",2,2);names(mycolors)=mynames
mycolors[which(mycolors%in%"_F3+")]=f3Pblues[3]
mycolors[which(mycolors%in%"_F3-")]=f3Mreds[3]
mycolors[which(mycolors%in%"_Total")]=CtrlMgreys[3]
mycolors.unique=unique(mycolors)
names(mycolors.unique)=c("_F3+","_F3-","_Total")
mycols=mycolors



####################################################################
##### Get general info first: read counts, trimming & mapping ######
logpath=paste("mapped/hsabrbseq/Demultiplexed_Lib",whichlib,"/logs/",sep="")
txtfiles=list.files(logpath,pattern=".txt")

####### First get the trimming & mapping stats & plot these ########
#### triming first
trimmed=read.table(paste(logpath,txtfiles[grep("Trim",txtfiles)],sep=""),header=FALSE,fill=T)[1:mynr,]
filenames=mysplit(trimmed[,1],"[.]",4,1)
mytab1=cbind(filenames,as.vector(trimmed[,5]),as.vector(trimmed[,8]),as.vector(trimmed[,31]), gsub("%","",mysplit(as.vector(trimmed[, 15]), "[(]", 2, 2)))
mytab1=cbind(mytab1,as.numeric(mytab1[,3])/as.numeric(mytab1[,2])*100,as.numeric(mytab1[,4])/as.numeric(mytab1[,2])*100)
colnames(mytab1)=c("Sample","Total","Trimmed","LenRemoved","%QCFiltered","%Trimmed","%LenRemoved")
for (i in seq(2,mynr,by=2)) {
	mytab1[i-1,7]=mytab1[i,7]; mytab1[i-1,4]=mytab1[i,4]
}
mytab1=mytab1[which(!is.na(mytab1[,1])),]
rownames(mytab1)=mytab1[,1]


##### Aligned Reads filcols=brewer.pal(9,"Pastel1")
pdf(paste(path,"TrimmedInfo-",myname,".pdf",sep=""),width=10,height=8)
##### Based on trimming and QC Info
par(mfrow=c(1,1))
mytot=as.numeric(mytab1[,2]);mytotKept=as.numeric(mytab1[,2])-as.numeric(mytab1[,3])
names(mytot)=rownames(mytab1)
barplot(sort(mytot)/1000000,col=filcols[4],names=mynames[mytab1[order(mytot),1]],main="Overview",ylab="Read Nrs. (Millions)",xlab="",las=2,cex.names=0.7)
barplot(mytotKept[order(mytot)]/1000000,col=filcols[3],names="",ylab="",xlab="",add=T,las=2)
legend("topleft",c("Kept","Filtered"),col=c(filcols[3],filcols[4]),lwd=2)

barplot(as.numeric(mytab1[,5])[order(mytot)],col=filcols[3],names=mytab1[order(mytot),1],main="QC Filter",ylab="% Bp QC Filter",xlab="",las=2,cex.names=0.7)
barplot(as.numeric(mytab1[,6])[order(mytot)],col=filcols[3],names=mytab1[order(mytot),1],main="Trimmed",ylab="% Reads Trimmed",xlab="",las=2,cex.names=0.7)
barplot(as.numeric(mytab1[,7])[order(mytot)],col=filcols[3],names=mytab1[order(mytot),1],main="Removed Reads",ylab="% Len Rem",xlab="",las=2,cex.names=0.7)

par(mfrow=c(2,2))
for (i in 2: dim(mytab1)[2]) {
	toplot=as.numeric(mytab1[,i])
	if (mean(toplot)>=1000000) {
		toplot=toplot/1000000
	}
	truehist(toplot,main=colnames(mytab1)[i],col=filcols[3],xlab=colnames(mytab1)[i],ylab="Density")
}
#### check read nrs. vs. concentration 
#myconc=as.numeric(as.vector(allinfo[names(mytot),6]))
#scatterNormalPlot(log1p(mytot),myconc,"Log.Reads","Conc","spearman")
dev.off()

##### REad the star alignment info
starpath=paste("mapped/hsabrbseq/Demultiplexed_Lib",whichlib,"/plots/",sep="")
#starpath="mapped/hsabrbseq/Lib50ng/plots/"
starnrs=read.table(paste(logpath,"STAR-align.txt",sep=""),sep="|",fill=T)
#starnrs=starnrs[-c(grep("E10",starnrs[,1]),grep("E11",starnrs[,1]),grep("E12",starnrs[,1])),]
starinfo=matrix(data=NA,nrow=mynr,ncol=9)
k=1
for(i in 1:dim(starinfo)[1]){
	starinfo[i,]=gsub("%","",gsub(" ","",as.vector(starnrs[,2][k:c(k+8)])))
	k=k+10
}
colnames(starinfo)=as.vector(starnrs[,1][1:9])
rownames(starinfo)=mynames[as.vector(starnrs[seq(10,dim(starnrs)[1],by=10),1][1:mynr])]
#rownames(starinfo)=as.vector(starnrs[seq(10,dim(starnrs)[1],by=10),1][1:mynr])
#starinfo=starinfo[as.character(samples[sindex,"typestagerep"]),]
write.table(starinfo,paste(starpath,"STARnrs-",myname,".txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
starinfo[,7]=as.numeric(starinfo[,1])-as.numeric(starinfo[,2])-as.numeric(starinfo[,5])-as.numeric(starinfo[,6])

mycol=filcols
pdf(paste(path,"STARnrs-",myname,".pdf",sep=""),width=9,height=7)
par(mfrow=c(1,1))
barplot(t(apply(starinfo[,c(2,5:7)],2,function(x) as.numeric(x)/1000000))[,order(apply(starinfo[,c(2,5:7)],1,function(x) as.numeric(x)/1000000)[1,])],
	col=mycol[c(2:5)],las=2,ylab=c("Read Numbers (Millions)"),cex.names=0.6,
	names=rownames(starinfo)[order(apply(starinfo[,c(2,5:7)],1,function(x) as.numeric(x)/1000000)[1,])])
legend("topleft",c("Unique","Multiple","Too multiple","Not Aligned"),col=mycol[c(2:5)],lwd=2)
barplot(apply(starinfo[,c(2,5:7)],1,function(x) as.numeric(x)/sum(as.numeric(x)))[,order(apply(starinfo[,c(2,5:7)],1,
	function(x) as.numeric(x)/sum(as.numeric(x)))[1,])],col=mycol[c(2:5)],las=2, 
	ylab=c("Read Numbers (Millions)"),cex.names=0.6)
legend("topleft",c("Unique","Multiple","Too multiple","Not Aligned"),col=mycol[c(2:5)],lwd=2)
par(mfrow=c(2,2))
truehist(t(apply(starinfo[,c(2,5:7)],2,function(x) as.numeric(x)/1000000))[1,],main="Aligned Reads",col=mycol[2],
	xlab="Nr reads",nbins=20,ylab="Cell density")
truehist(apply(starinfo[,c(2,5:7)],1,function(x) as.numeric(x)/sum(as.numeric(x)))[1,],main="Fraction Aligned Reads",col=mycol[2],
	xlab="Fraction Aligned",nbins=20,ylab="Cell density")
dev.off()


####### Read Htseq quantifications now ########
path=paste("mapped/hsabrbseq/Demultiplexed_Lib",whichlib,"/",sep="")
outpath=paste("mapped/hsabrbseq/Demultiplexed_Lib",whichlib,"/plots/",sep="")
#path="mapped/hsabrbseq/Lib50ng/"
#outpath="mapped/hsabrbseq/Lib50ng/plots/"

myfiles=list.files(path)
myfiles=myfiles[which(!myfiles%in%myfiles[grep("plots",myfiles)])]
myfiles=myfiles[which(!myfiles%in%myfiles[grep("Transcriptome",myfiles)])]
myfiles=myfiles[which(!myfiles%in%myfiles[grep("STAR",myfiles)])]
myfiles=myfiles[which(!myfiles%in%myfiles[grep("pdf",myfiles)])]
myfiles=myfiles[which(!myfiles%in%myfiles[grep("logs",myfiles)])]
myfiles=myfiles[which(!myfiles%in%myfiles[c(grep("E10",myfiles),grep("Limma",myfiles),grep("Background",myfiles))])]
###myfiles=myfiles[which(myfiles%in%names(mynames))]
readmatrix=matrix(data=NA,ncol=length(myfiles),nrow=60775)
for(i in 1:length(myfiles)) {
	truefile=list.files(paste(path,myfiles[i],"/",sep=""),pattern="gene")#
	readmatrix[,i]=read.table(paste(path,myfiles[i],"/",truefile,sep=""),sep="\t",header=FALSE)[,2]
}
rownames(readmatrix)=read.table(paste(path,myfiles[i],"/",truefile,sep=""),sep="\t",header=FALSE)[,1]
colnames(readmatrix)=myfiles; htseqcounts=readmatrix
save(htseqcounts,file=paste(outpath,"STAR-HTSeqTags.-",myname,".RData",sep=""))

load(paste(outpath,"STAR-HTSeqTags.-",myname,".RData",sep=""))
#htseqcounts.6=htseqcounts
#save(htseqcounts.6,file="mapped/hsabrbseq/Demultiplexed_Lib6/plots/STAR-HTSeqTags.-Jan17.5.RData")
#load("mapped/hsabrbseq/Demultiplexed_Lib5/plots/STAR-HTSeqTags.-Jan17.5.RData")
write.table(htseqcounts,paste(outpath,"STAR-HTSeqTags.-",myname,".counts.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

htseqinfo=htseqcounts[(c(dim(htseqcounts)[1]-4):dim(htseqcounts)[1]),]
htseqcounts=htseqcounts[-(c(dim(htseqcounts)[1]-4):dim(htseqcounts)[1]),]
htseqinfo=rbind(apply(htseqcounts,2,function(x) sum(x)),htseqinfo)
rownames(htseqinfo)=c("Gene","No_gene","Ambig","lowQ","no_align","Not_unique")

#### add the higher quality data ####
mycol=filcols
pdf(paste(path,"HTSeqnrs-",myname,".pdf",sep=""),width=9,height=7)
par(mfrow=c(1,1))
barplot(t(apply(htseqinfo[c(1:3,6),],1,function(x) as.numeric(x)/1000000)[order(apply(htseqinfo[c(1:3,6),],1,function(x) as.numeric(x)/1000000)[,1]),]),
	col=mycol[c(2:5)],las=2,ylab=c("Read Numbers (Millions)"),cex.names=0.6,names=colnames(htseqinfo)[order(apply(htseqinfo[c(1:3,6),],1,function(x) as.numeric(x)/1000000)[,1])])
legend("topleft",c("Gene","No_Gene","Ambig","Not_unique"),col=mycol[c(2:5)],lwd=2)
barplot(t(apply(htseqinfo[c(1:3,6),],1,function(x) as.numeric(x)/sum(as.numeric(x)))[order(apply(htseqinfo[c(1:3,6),],1,function(x) as.numeric(x)/sum(as.numeric(x)))[,1]),]),
	col=mycol[c(2:5)],las=2, ylab=c("Read Numbers (Millions)"),cex.names=0.6)
legend("topleft",c("Gene","No_Gene","Ambig","Not_unique"),col=mycol[c(2:5)],lwd=2)
par(mfrow=c(2,2))
truehist(t(apply(htseqinfo[c(1:3,6),],1,function(x) as.numeric(x)/1000000))[1,],main="Genes",col=mycol[2],
	xlab="Nr reads",nbins=20,ylab="Cell density")
truehist(apply(htseqinfo[c(1:3,6),],1,function(x) as.numeric(x)/sum(as.numeric(x)))[,1],main="Fraction Genes",col=mycol[2],
	xlab="Fraction Aligned",nbins=20,ylab="Cell density")
dev.off()


mycols=mycolors
#### make some overview plots: PCA, heatmaps, clustering, ... ####
pdf(paste(path,"STAR-",myname,".overviewplots.pdf",sep=""),width=7,height=7)
par(mfrow=c(2,2))
genesPerSample.1=apply(htseqcounts,2,function(x) length(which(x>0)))
genesPerSample.5=apply(htseqcounts,2,function(x) length(which(x>4)))
countPerSample=apply(htseqcounts,2,function(x) log1p(sum(x)))
truehist(genesPerSample.1,main="genesPerSample.1",ylab="Nr. genes");truehist(countPerSample,main="countPerSample",ylab="log SumCounts")
par(mfrow=c(1,1))
barplot(sort(genesPerSample.1),cex.names=0.4,las=2,ylab="Nr. genes");barplot(sort(countPerSample),cex.names=0.4,las=2,ylab="log SumCounts")
par(mfrow=c(2,2))
scatterNormalPlot(genesPerSample.1,countPerSample,"genePerSample","countPerSample","spearman")
scatterNormalPlot(genesPerSample.1,mytot[names(genesPerSample.1)]/1000000,"genePerSample","totalReads","spearman")
scatterNormalPlot(countPerSample,mytot[names(genesPerSample.1)]/1000000,"countPerSample","totalReads","spearman")

mymatrix=htseqcounts
tonorm=apply(htseqcounts,2,sum)/1000000
mymatrix=apply(mymatrix,2,function(x) x/(sum(x)/1000000))
mymatrix=mymatrix[which(apply(mymatrix,1,function(x) length(which(x>0))>2)),]
truehist(apply(mymatrix,1,function(x) mean(x[which(x>0)])),main="Mean Expression")
mymatrix=mymatrix[which(apply(mymatrix,1,function(x) sum(x))>10),]
mymatrix.log=log1p(mymatrix)
#makeCorPlots(path,allfiles.mmseqnorm.mod[[1]],myset) 
par(mfrow=c(1,1))
mymatrixOverviewPlots(mymatrix.log,mycols,paste(myname,"-all",sep=""),bootnr=20)
##### Redo removing outlier samples #####
#outliers=c("DN.td.4","F3m.td.4","F3p.td.4","DN.tb.4","F3m.tb.4","F3p.tb.4","F3m.ps.1",
#	"F3m.td.1","F3m.tb.1","F3mABCG1m.tb.2","F3pABCG1p.td.2","F3pABCG1p.tb.2","DN.nt.1","F3m.nt.1","DN.td.1",
#	"F3m.td.1","F3p.td.1","DN.tb.1","F3m.tb.1","F3p.tb.1")
#mymatrix.log=mymatrix.log[,which(!colnames(mymatrix.log)%in%outliers)]
#par(mfrow=c(1,1))
#mymatrixOverviewPlots(mymatrix.log,mycols,paste(myname,"-all-nooutliers",sep=""),bootnr=20)
dev.off()


###### Process using limma #####
#fullCountTable=htseqcounts[,which(!colnames(htseqcounts)%in%outliers)]
# F3pABCG1p.ps.3
#torem=-which(colnames(htseqcounts)%in%outliers)
torem=c(1:12)
fullCountTable=htseqcounts[,torem]
fullCountTable = round(fullCountTable[apply(fullCountTable,1,function(x) any(round(x)>0)),] ) # only keep features with counts
dge <- DGEList(counts=fullCountTable, genes=rownames(fullCountTable))
dge <- calcNormFactors(dge); isexpr <- rowSums(cpm(dge) > 2) >= 3
y <- voom(dge[isexpr,],plot=F)
lev <- c(unique(gsub("[+]","p",gsub("-","m",gsub("_","",mysplit(colnames(y$E),"T0",2,2)[torem]))))) 
lev2 <- c(unique(colnames(y$E)[torem]))
f <- factor(gsub("[+]","p",gsub("-","m",gsub("_","",mysplit(colnames(y$E),"T0",2,2)[torem]))), levels=lev)
celltype <- factor(colnames(y$E), levels=lev2)
batch=mysplit(colnames(y$E),"T0",2,1)[torem]; names(batch)=colnames(fullCountTable)
design <- model.matrix(~0+f)
colnames(design) <- levels(f)
##colnames(design) <- lev; corfit <- duplicateCorrelation(y$E,design,block=Mice)
##corfit$consensus fit <- lmFit(y$E,design,block=Mice,correlation=corfit$consensus)

mydesign=cbind(rep(1:length(celltype)),as.character(f),batch)
colnames(mydesign)[1]="sample"; mydesign=data.frame(mydesign)

library(sva)
mycounts=y$E
#mod = model.matrix(~f, data=mydesign); mod0 = model.matrix(~1,data=mydesign)
modcombat = model.matrix(~1, data=mydesign)
mycounts.combat = ComBat(dat=mycounts, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
fit.combat = lmFit(mycounts.combat,design)
fit <- lmFit(y$E,design)

cm <- makeContrasts(F3mvsF3p = F3m-F3p,
	DPvsF3m = F3m-Total,
	DPvsF3p = F3p-Total,
	levels=design) #,
	
fit2 <- contrasts.fit(fit, cm); fit2 <- eBayes(fit2)
fit2.c <- contrasts.fit(fit.combat, cm); fit2.c <- eBayes(fit2.c)

#### check cor. between F3 and genes #####
cor.F3=apply(mycounts,1,function(x) 
	cor(mycounts[getEnsgID(hsagenes84.GR,"F3"),],x,method="spearman"))
toprint=cbind(names(sort(cor.F3,decreasing=T)),
	elementMetadata(hsagenes84.GR[names(sort(cor.F3,decreasing=T))])$geneName,sort(cor.F3,decreasing=T))
write.table(toprint,paste(outpath,"Limma-corF3-",myname,".txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
save(fullCountTable,dge,fit, fit2,mycounts.combat,fit2.c,cor.F3,file=paste(outpath,"Limma-",myname,".RData",sep=""))

load(paste(outpath,"Limma-",myname,".RData",sep=""))

mycounts=mycounts.combat

#mycounts[getEnsgID(hsagenes84.GR,"ADIPOQ"),grep("F3p",colnames(mycounts))]
mycounts[getEnsgID(hsagenes84.GR,"ADIPOQ"),grep("F3-",colnames(mycounts))]
mycounts[getEnsgID(hsagenes84.GR,"CIDEC"),grep("F3-",colnames(mycounts))]
mycounts[getEnsgID(hsagenes84.GR,"FASN"),grep("F3-",colnames(mycounts))]
mycounts[getEnsgID(hsagenes84.GR,"LPL"),grep("F3-",colnames(mycounts))]

mycounts[getEnsgID(hsagenes84.GR,"F3"),grep("ps",colnames(mycounts))]
mycounts[getEnsgID(hsagenes84.GR,"ABCG1"),grep("ps",colnames(mycounts))]

#### F3 comparisons
pcut=0.1;fccut=2;usefit=fit2.c;selnr=20000
#DPvsF3m.0.de=topTable(usefit, adjust="BH",coef="DPvsF3m.0",p.value=pcut,number=selnr,lfc=log2(fccut))
#DPvsF3p.0.de=topTable(usefit, adjust="BH",coef="DPvsF3p.0",p.value=pcut,number=selnr,lfc=log2(fccut))
F3mvsF3p.de=topTable(usefit, adjust="BH",coef="F3mvsF3p",p.value=pcut,number=selnr,lfc=log2(fccut))
DPvsF3m.de=topTable(usefit, adjust="BH",coef="DPvsF3m",p.value=pcut,number=selnr,lfc=log2(fccut))
DPvsF3p.de=topTable(usefit, adjust="BH",coef="DPvsF3p",p.value=pcut,number=selnr,lfc=log2(fccut))

F3mvsF3p.de[getEnsgID(hsagenes84.GR,"F3"),]
F3mvsF3p.de[getEnsgID(hsagenes84.GR,"ABCG1"),]

F3mvsF3p.de[getEnsgID(hsagenes84.GR,"ADIPOQ"),]
F3mvsF3p.de[getEnsgID(hsagenes84.GR,"FASN"),]
F3mvsF3p.de[getEnsgID(hsagenes84.GR,"LPL"),]
F3mvsF3p.de[getEnsgID(hsagenes84.GR,"CIDEC"),]

allDEs=list(F3mvsF3p.de,DPvsF3m.de,DPvsF3p.de)

names(allDEs)=c("F3mvsF3p.de","DPvsF3m.de","DPvsF3p.de")

results <- decideTests(fit2,adjust.method="BH",p.value=pcut,lfc=log2(fccut))
results2 <- decideTests(fit2.c,adjust.method="BH",p.value=pcut,lfc=log2(fccut))

mycols.saved=mycols;mycounts=mycounts.combat

myboot=100
##### vennDiagramms with the results ######### make some overview plots: PCA, heatmaps, clustering, ... ####
pdf(paste(path,"STAR-limma-",myname,".overviewplots.combat.pdf",sep=""),width=7,height=7)
par(mfrow=c(1,1))

###### plots based on y$E values
mymatrixOverviewPlots(mycounts,mycols,paste(myname,"-all",sep=""),bootnr=myboot)

######## plots based on differential genes
mymatrix.log=mycounts[unique(unlist(lapply(allDEs,function(x) rownames(x)))),]; par(mfrow=c(1,1))
mymatrixOverviewPlots(mymatrix.log,mycols,paste(myname,"-all-de",sep=""),bootnr=myboot)

#pdf(paste(path,"STAR-limma-",myname,".overviewplots.combat.v2.de.pdf",sep=""),width=7,height=7)

mymatrix.log=mycounts
vennDiagram(results)
vennDiagram(results2)
barplot(sapply(allDEs, function(x) dim(x)[1]),main="Nr. of DE regions - FDR 5% FC2",ylab="Nr. genes",las=2,cex.names=0.7)
dev.off()

#### output raw for Figure
write.table(results2,paste(path,"FigS6d-raw.hsaVenns.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

colorder=c("i04_T0_F3-","i05_T0_F3-","i06_T0_F3-","i08_T0_F3-",
	"i04_T0_F3+","i05_T0_F3+","i06_T0_F3+","i08_T0_F3+",
	"i04_T0_Total","i05_T0_Total","i06_T0_Total","i08_T0_Total")
colorder=colorder[which(colorder%in%colnames(mymatrix.log))]

allgenes.background=rownames(mycounts);myp=0.001;showgraph=10; 
write.table(allgenes.background,paste(path,"BackgroundGenes-Limma-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(hsagenes84.GR[which(names(hsagenes84.GR)%in%allgenes.background)])$geneName,paste(path,"BackgroundGenes-Limma-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),
	row.names=F,col.names=F,quote=F)


pdf(paste(path,"STAR-limma-",myname,".heatmaps.Pub.pdf",sep=""),width=7,height=7)
#### 1. Validation of SC results & higher adipogenesis in F3- vs. DP? 
# Plot expression of chosen adipogenic genes
plotgenes=toupper(unique(c("F3","Mgp","Meox2","Abcg1","Dpp4","Akr1c18","Cd55","Il13ra1","Sparcl1","Adam12","Aoc3","Cd34","Tnfa","Il6","Cd34","Cd31",
	"Pecam1","Mcam","Vcam1","Pparg","Fabp4","Adipoq","Cidec","Cd36","Retn","Fasn","Dkk3",
	"Fabp12","Lpl","Alpl","Fmo2","Il6",gsub(" ","",unique(c(apply(mymarkers[c(1:6,9:10),],2,function(x) as.character(x))))))))
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(hsagenes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Interest.Genes",
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])

# Plot expression of endothelial genes
plotgenes=toupper(unique(c(c(gsub(" ","",t(mymarkers[c("endo","adipo","osteo","chondro"),]))))))
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(hsagenes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Interest.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])

# Plot expression of interest genes:
daniel=elementMetadata(hsagenes84.GR[c("ENSG00000173546", "ENSG00000143248", "ENSG00000175084", "ENSG00000113721")])$geneName
tgfbgenes="GO:0007179"; fgfgenes="GO:0008543"; vegfgenes="GO:0038084"
adipogenes=c("GO:0045444", "GO:0060612")
hhgenes=c("GO:0007224")

tgfbgenes=getBM( c("external_gene_name","ensembl_gene_id"),filters=c("go_id"),values=tgfbgenes,mart = ensembl84hg)[,1]
fgfgenes=getBM( c("external_gene_name","ensembl_gene_id"),filters=c("go_id"),values=fgfgenes,mart = ensembl84hg)[,1]
vegfgenes=getBM( c("external_gene_name","ensembl_gene_id"),filters=c("go_id"),values=vegfgenes,mart = ensembl84hg)[,1]
adipogenes=getBM( c("external_gene_name","ensembl_gene_id"),filters=c("go_id"),values=adipogenes,mart = ensembl84hg)[,1]
hhgenes=getBM( c("external_gene_name","ensembl_gene_id"),filters=c("go_id"),values=hhgenes,mart = ensembl84hg)[,1]

# Plot expression of Daniels genes
plotgenes=toupper(daniel)
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(hsagenes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Daniel.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])

# Plot expression of tgfb genes
plotgenes=toupper(tgfbgenes)
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(hsagenes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Tgfb.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])

# Plot expression of fgf genes
plotgenes=toupper(fgfgenes)
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(hsagenes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Fgf.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])

# Plot expression of vegf genes
plotgenes=toupper(vegfgenes)
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(hsagenes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Vegf.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])

# Plot expression of hh genes
plotgenes=toupper(hhgenes)
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(hsagenes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Hh.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])

# Plot expression of adipo genes
plotgenes=toupper(adipogenes)
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(hsagenes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Adipo.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])

dev.off()


pdf(paste(path,"STAR-limma-",myname,".heatmaps.Pub.final.pdf",sep=""),width=7,height=7)
# Plot expression of most F3+ genes
plotmatrix=mycounts[rownames(F3mvsF3p.de[which(F3mvsF3p.de[,1]<0),]),colorder]
mainname=myname
myIDsF3pstr.GO=runGO.84(rownames(plotmatrix),allgenes.background,outpath,ensembl84hg,
	"BP","elimCount",myp,paste(mainname,"GO-pcut-F3p-",i,sep="_"),ontoGOs.hg84,showgraph)
write.table(rownames(plotmatrix),paste(path,"F3p-Limma-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(hsagenes84.GR[which(names(hsagenes84.GR)%in%rownames(plotmatrix))])$geneName,
	paste(path,"F3p-Limma-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),row.names=F,col.names=F,quote=F)
temp=plotmatrix[,colorder]; rownames(temp)=elementMetadata(hsagenes84.GR[which(names(hsagenes84.GR)%in%rownames(plotmatrix))])$geneName
write.table(temp,paste(path,"F3p-Limma-",myname,"-",myset,"-",mydate,".details.txt",sep=""),
	row.names=T,col.names=T,quote=F)
rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F, main="F3+ high",
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
#### output raw for Figure
write.table(plotmatrix,paste(path,"FigS6g-raw.hsaAregPHeatmap.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


#### 1. Validation of SC results & higher adipogenesis in F3- vs. DP? 
# Plot expression of most F3- genes
plotmatrix=mycounts[rownames(F3mvsF3p.de[which(F3mvsF3p.de[,1]>0),]),colorder]

myIDsF3mstr.GO=runGO.84(rownames(plotmatrix),allgenes.background,outpath,ensembl84hg,
	"BP","elimCount",myp,paste(mainname,"GO-pcut-F3m-",i,sep="_"),ontoGOs.hg84,showgraph)
write.table(rownames(plotmatrix),paste(path,"F3m-Limma-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(hsagenes84.GR[which(names(hsagenes84.GR)%in%rownames(plotmatrix))])$geneName,
	paste(path,"F3m-Limma-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),row.names=F,col.names=F,quote=F)
temp=plotmatrix[,colorder]; rownames(temp)=elementMetadata(hsagenes84.GR[which(names(hsagenes84.GR)%in%rownames(plotmatrix))])$geneName
write.table(temp,paste(path,"F3m-Limma-",myname,"-",myset,"-",mydate,".details.txt",sep=""),
	row.names=T,col.names=T,quote=F)

rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F, main="F3- high",
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])

#### select only the top 60
if (dim(plotmatrix)[1]>60) {
	plotmatrix=mycounts[rownames(F3mvsF3p.de[which(F3mvsF3p.de[,1]>0),]),colorder]
	plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%rownames(F3mvsF3p.de[order(F3mvsF3p.de[,5]),])[1:100]),colorder]
	plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%rownames(F3mvsF3p.de[order(F3mvsF3p.de[,1],decreasing=T),])[1:100]),colorder]
	rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
	pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F, main="F3- high",
		clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])

	plotmatrix=mycounts[rownames(F3mvsF3p.de[which(F3mvsF3p.de[,1]>0),]),colorder]
	plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%rownames(F3mvsF3p.de[order(F3mvsF3p.de[,5]),])[1:35]),colorder]
#	plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%rownames(F3mvsF3p.de[order(F3mvsF3p.de[,1],decreasing=T),])[1:30]),colorder]
	rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
	pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F, main="F3- high",
		clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])

	#### in the heatmap
	plotmatrix=mycounts[rownames(F3mvsF3p.de[which(F3mvsF3p.de[,1]>0),]),colorder]
	plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%rownames(F3mvsF3p.de[order(F3mvsF3p.de[,5]),])[1:100]),colorder]
	plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%rownames(F3mvsF3p.de[order(F3mvsF3p.de[,1],decreasing=T),])[1:50]),colorder]
	rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
	pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F, main="F3- high",
		clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
	#### output raw for Figure
	write.table(plotmatrix,paste(path,"Fig3e-raw.hsaAregMHeatmap.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
	
	plotmatrix=mycounts[rownames(F3mvsF3p.de[which(F3mvsF3p.de[,1]>0),]),colorder]
	plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%rownames(F3mvsF3p.de[order(F3mvsF3p.de[,5]),])[1:100]),colorder]
#	plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%rownames(F3mvsF3p.de[order(F3mvsF3p.de[,1],decreasing=T),])[1:25]),colorder]
	rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
	pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F, main="F3- high",
		clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
	
}

dev.off()

my5up=as.character(read.table("/home/pschwali/myarchive/adipo/ETH/mapped/hsabrbseq/Demultiplexed_Lib5/F3p-Limma-Jan17-bow-mmseqEnsg84-runCut-20170110-20170110.ID.txt",head=F)[,1])
my6up=as.character(read.table("/home/pschwali/myarchive/adipo/ETH/mapped/hsabrbseq/Demultiplexed_Lib6/F3p-Limma-Jan17-bow-mmseqEnsg84-runCut-20170110-20170110.ID.txt",head=F)[,1])

my5down=as.character(read.table("/home/pschwali/myarchive/adipo/ETH/mapped/hsabrbseq/Demultiplexed_Lib5/F3m-Limma-Jan17-bow-mmseqEnsg84-runCut-20170110-20170110.ID.txt",head=F)[,1])
my6down=as.character(read.table("/home/pschwali/myarchive/adipo/ETH/mapped/hsabrbseq/Demultiplexed_Lib6/F3m-Limma-Jan17-bow-mmseqEnsg84-runCut-20170110-20170110.ID.txt",head=F)[,1])

#### check Fishers for pathways
my5up[which(my5up%in%tgfbgenes)];# SERPINE1 
my5up[which(my5up%in%fgfgenes)]; #"PDE1C" "FGF1"
my5up[which(my5up%in%vegfgenes)]; my5up[which(my5up%in%daniel)]
my6up[which(my6up%in%tgfbgenes)]; my6up[which(my6up%in%fgfgenes)]; my6up[which(my6up%in%vegfgenes)]; my6up[which(my6up%in%daniel)]
my5down[which(my5down%in%tgfbgenes)]; my5down[which(my5down%in%fgfgenes)]; my5down[which(my5down%in%vegfgenes)]; my5down[which(my5down%in%daniel)]
my6down[which(my6down%in%tgfbgenes)]; my6down[which(my6down%in%fgfgenes)]; #"PRKAR2B" "KLB"
my6down[which(my6down%in%vegfgenes)]; my6down[which(my6down%in%daniel)] #"CSPG4"


pdf(paste(path,"STAR-limma-",myname,".hsaAndmmu.pdf",sep=""),width=7,height=7)
allgenes.background=rownames(mycounts);myp=0.001;showgraph=10; 

#### Check overlap with mouse #####

mousepath="/home/pschwali/myarchive/adipo/ETH/mapped/brbseq/STAR/InterestGenes/"
mouse.F3p.D0s=toupper(as.vector(read.table(paste(mousepath,"F3p-D0s-Limma-Nov16-bow-mmseqEnsg84-runCut-20161121-161121.ID.txt",sep=""),head=F)[,1]))
mouse.F3p.5h=toupper(as.vector(read.table(paste(mousepath,"F3p-5hs-Limma-Nov16-bow-mmseqEnsg84-runCut-20161121-161121.ID.txt",sep=""),head=F)[,1]))
mouse.F3p.24h=toupper(as.vector(read.table(paste(mousepath,"F3p-24hs-Limma-Nov16-bow-mmseqEnsg84-runCut-20161121-161121.ID.txt",sep=""),head=F)[,1]))

#### Common genes mouse & human
plotgenes=my5up[which(my5up%in%c(mouse.F3p.D0s,mouse.F3p.5h,mouse.F3p.24h))]
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(hsagenes84.GR,plotgenes))],colorder]
myIDsF3MouseHuman.GO=runGO.84(rownames(plotmatrix),allgenes.background,outpath,ensembl84hg,
	"BP","elimCount",myp,paste(mainname,"GO-pcut-F3pMouseAndHuman-",i,sep="_"),ontoGOs.hg84,showgraph)
write.table(rownames(plotmatrix),paste(path,"F3pMouseAndHuman-Limma-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(hsagenes84.GR[which(names(hsagenes84.GR)%in%rownames(plotmatrix))])$geneName,
	paste(path,"F3pMouseAndHuman-Limma-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),row.names=F,col.names=F,quote=F)
temp=plotmatrix[,colorder]; rownames(temp)=elementMetadata(hsagenes84.GR[which(names(hsagenes84.GR)%in%rownames(plotmatrix))])$geneName
write.table(temp,paste(path,"F3pMouseAndHuman-Limma-",myname,"-",myset,"-",mydate,".details.txt",sep=""),
	row.names=T,col.names=T,quote=F)

rownames(plotmatrix)=elementMetadata(hsagenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="F3pMouseAndHuman.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
dev.off()







	
