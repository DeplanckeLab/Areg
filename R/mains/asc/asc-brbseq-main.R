#bsub -M 4000000 -R  "rusage[mem=4000]" -Is bash
source("plot-util.R")
source("util-asc.R")

library(beanplot)

#### Initialise path #####
path="mapped/brbseq/"

#### Read the annotations ####
infopath="data/brbseq/"
samples=read.csv(paste(infopath,"ASCbrbseq_nov16.csv",sep=""))
rownames(samples)=samples[,1]
samples=cbind(samples,mysplit(as.character(samples[,5]),"_",2,1),mysplit(as.character(samples[,5]),"_",2,2),
	 gsub("DP1","DP",substr(mysplit(as.character(samples[,3]),"_",2,1),1,3)),
	 mysplit(mysplit(as.character(samples[,3]),"_",2,2),"[.]",2,2))
samples=cbind(samples,paste(samples[,9],samples[,10],sep="."))
colnames(samples)=c("ID","Nr","Cells","Nr.cells","Type","Buffer","Replicate","Time","CellType","Time2","CellTypeTime")
samples=as.matrix(samples)
samples[26,3]="DP11_3.0"; samples[40,3]="DP11_4.0"

##### prepare colors
f3Pblues=brewer.pal(9,"Blues"); f3Mreds=brewer.pal(9,"Reds"); CtrlMgreys=brewer.pal(9,"Greys")
mycolors=unique(as.character(samples[,11])); names(mycolors)=unique(samples[,11])
mycolors["DP.0"]=CtrlMgreys[3]; mycolors["F3-.0"]=f3Mreds[3];mycolors["F3+.0"]=f3Pblues[3]
mycolors["F3-.5h"]=f3Mreds[5];mycolors["F3+.5h"]=f3Pblues[5]
mycolors["F3-.24h"]=f3Mreds[7];mycolors["F3+.24h"]=f3Pblues[7]
mycolors["DP.d8"]=CtrlMgreys[9]; mycolors["F3-.d8"]=f3Mreds[9];mycolors["F3+.d8"]=f3Pblues[9]
mycolors.unique=mycolors
mycolors=mycolors[ as.character(samples[,"CellTypeTime"])]; names(mycolors)=samples[,3]
mycols=mycolors

#### Initialise sample details #####
mynames=samples[,3]; names(mynames)=samples[,1] #### Set the names
myname="Nov16"; mynr=52; mynames.touse=mynames 
filcols=brewer.pal(9,"Pastel1")
myset="bow-mmseqEnsg84-runCut-20161121"
mmset=myset
mydate="161121"

####################################################################
##### Get general info first: read counts, trimming & mapping ######
logpath="mapped/brbseq/logs/"
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
rownames(mytab1)=mynames[rownames(mytab1)]
mytab1=mytab1[as.character(samples[,"Cells"]),]


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
dev.off()

##### Read the star alignment info
starpath="mapped/brbseq/STAR/plots/"
starnrs=read.table(paste(logpath,"STAR-align.txt",sep=""),sep="|",fill=T)
starnrs=starnrs[-c(grep("E10",starnrs[,1]),grep("E11",starnrs[,1]),grep("E12",starnrs[,1])),]
starinfo=matrix(data=NA,nrow=mynr,ncol=9)
k=1
for(i in 1:dim(starinfo)[1]){
	starinfo[i,]=gsub("%","",gsub(" ","",as.vector(starnrs[,2][k:c(k+8)])))
	k=k+10
}
colnames(starinfo)=as.vector(starnrs[,1][1:9])
rownames(starinfo)=mynames[as.vector(starnrs[seq(10,dim(starnrs)[1],by=10),1][1:mynr])]
starinfo=starinfo[as.character(samples[,"Cells"]),]
write.table(starinfo,paste(starpath,"STARnrs-",myname,".txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
starinfo[,7]=as.numeric(starinfo[,1])-as.numeric(starinfo[,2])-as.numeric(starinfo[,5])-as.numeric(starinfo[,6])

mycol=filcols
pdf(paste(path,"STARnrs-",myname,".pdf",sep=""),width=9,height=7)
par(mfrow=c(1,1))
barplot(t(apply(starinfo[,c(2,5:7)],2,function(x) as.numeric(x)/1000000))[,order(apply(starinfo[,c(2,5:7)],1,function(x) as.numeric(x)/1000000)[1,])],
	col=mycol[c(2:5)],las=2,ylab=c("Read Numbers (Millions)"),cex.names=0.6,names=rownames(starinfo)[order(apply(starinfo[,c(2,5:7)],1,function(x) as.numeric(x)/1000000)[1,])])
legend("topleft",c("Unique","Multiple","Too multiple","Not Aligned"),col=mycol[c(2:5)],lwd=2)
barplot(apply(starinfo[,c(2,5:7)],1,function(x) as.numeric(x)/sum(as.numeric(x)))[,order(apply(starinfo[,c(2,5:7)],1,
	function(x) as.numeric(x)/sum(as.numeric(x)))[1,])],col=mycol[c(2:5)],las=2, ylab=c("Read Numbers (Millions)"),cex.names=0.6)
legend("topleft",c("Unique","Multiple","Too multiple","Not Aligned"),col=mycol[c(2:5)],lwd=2)
par(mfrow=c(2,2))
truehist(t(apply(starinfo[,c(2,5:7)],2,function(x) as.numeric(x)/1000000))[1,],main="Aligned Reads",col=mycol[2],
	xlab="Nr reads",nbins=20,ylab="Cell density")
truehist(apply(starinfo[,c(2,5:7)],1,function(x) as.numeric(x)/sum(as.numeric(x)))[1,],main="Fraction Aligned Reads",col=mycol[2],
	xlab="Fraction Aligned",nbins=20,ylab="Cell density")
dev.off()

#### output raw for Figure
temp=apply(starinfo[,c(2,5:7)],2,function(x) as.numeric(x)/1000000)
rownames(temp)=rownames(starinfo)
write.table(temp,paste(path,"FigS4c-raw.BrbseqMapping.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


####### Read Htseq quantifications now ########
path="mapped/brbseq/STAR/"
outpath="mapped/brbseq/STAR/plots/"
myfiles=list.files(path)
myfiles=myfiles[which(!myfiles%in%myfiles[grep("plots",myfiles)])]
myfiles=myfiles[which(!myfiles%in%myfiles[grep("Transcriptome",myfiles)])]
myfiles=myfiles[which(!myfiles%in%myfiles[grep("STAR",myfiles)])]
myfiles=myfiles[which(!myfiles%in%myfiles[grep("pdf",myfiles)])]
myfiles=myfiles[which(!myfiles%in%myfiles[c(grep("E10",myfiles),grep("E11",myfiles),grep("E12",myfiles))])]
myfiles=myfiles[which(!myfiles%in%myfiles[grep("genes",myfiles)])]
readmatrix=matrix(data=NA,ncol=length(myfiles),nrow=47829)
for(i in 1:length(myfiles)) {
	truefile=list.files(paste(path,myfiles[i],"/",sep=""),pattern="gene")#
	readmatrix[,i]=read.table(paste(path,myfiles[i],"/",truefile,sep=""),sep="\t",header=FALSE)[,2]
}
rownames(readmatrix)=read.table(paste(path,myfiles[i],"/",truefile,sep=""),sep="\t",header=FALSE)[,1]
colnames(readmatrix)=mynames[myfiles]; htseqcounts=readmatrix
htseqcounts=htseqcounts[,as.character(samples[,"Cells"])]
save(htseqcounts,file=paste(outpath,"STAR-HTSeqTags.-",myname,".RData",sep=""))
write.table(htseqcounts,paste(outpath,"STAR-HTSeqTags.-",myname,".counts.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

#### Start from here if only counts available #####
load(paste(outpath,"STAR-HTSeqTags.-",myname,".RData",sep=""))

htseqinfo=htseqcounts[(c(dim(htseqcounts)[1]-4):dim(htseqcounts)[1]),]
htseqcounts=htseqcounts[-(c(dim(htseqcounts)[1]-4):dim(htseqcounts)[1]),]
htseqinfo=rbind(apply(htseqcounts,2,function(x) sum(x)),htseqinfo)
rownames(htseqinfo)=c("Gene","No_gene","Ambig","lowQ","no_align","Not_unique")

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
#### Overview plots: PCA, heatmaps, clustering, ... ####
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
par(mfrow=c(1,1))
mymatrixOverviewPlots(mymatrix.log,mycols,paste(myname,"-all",sep=""),bootnr=20)
dev.off()

#### output raw for Figure
write.table(genesPerSample.1/1000,paste(path,"FigS4c-raw.BrbseqGenesPerSample.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


###### Process for DE using limma #####
fullCountTable=htseqcounts
colnames(fullCountTable)=gsub("-","m",colnames(fullCountTable))
colnames(fullCountTable)=gsub("[+]","p",colnames(fullCountTable))
fullCountTable = round(fullCountTable[apply(fullCountTable,1,function(x) any(round(x)>0)),] ) # only keep features with counts
dge <- DGEList(counts=fullCountTable, genes=rownames(fullCountTable))
dge <- calcNormFactors(dge); isexpr <- rowSums(cpm(dge) > 2) >= 3
y <- voom(dge[isexpr,],plot=F)
lev <- c(unique(gsub("[+]","p",gsub("-","m",as.character(samples[,"CellTypeTime"]))))); 
lev2 <- c(unique(gsub("[+]","p",gsub("-","m",as.character(samples[,"CellType"])))))
f <- factor(gsub("[+]","p",gsub("-","m",as.character(samples[,"CellTypeTime"]))), levels=lev)
celltype <- factor(gsub("[+]","p",gsub("-","m",as.character(samples[,"CellType"]))), levels=lev2)
batch=samples[,7]; names(batch)=colnames(fullCountTable)
design <- model.matrix(~0+f)
colnames(design) <- levels(f)

mydesign=cbind(rep(1:length(celltype)),as.character(celltype),batch,as.character(samples[,"Time"]))
colnames(mydesign)[1]="sample"; mydesign=data.frame(mydesign)

#### use combat for batch correction ####
library(sva)
mycounts=y$E
#mod = model.matrix(~f, data=mydesign); mod0 = model.matrix(~1,data=mydesign)
modcombat = model.matrix(~1, data=mydesign)
mycounts.combat = ComBat(dat=mycounts, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
fit.combat = lmFit(mycounts.combat,design)
fit <- lmFit(y$E,design)

cm <- makeContrasts(DPvsF3m.0 = DP.0-F3m.0,
	DPvsF3p.0 = DP.0-F3p.0,
	F3mvsF3p.0 = F3m.0-F3p.0,
	F3mvsF3p.5h = F3m.0-F3p.5h,
	F3mvsF3p.24h = F3m.0-F3p.24h,
	DPvsF3m.8 = DP.d8-F3m.d8,
	DPvsF3p.8 = DP.d8-F3p.d8,
	F3mvsF3p.8 = F3m.d8-F3p.d8,
	DPd0vs8.0 = DP.0-DP.d8,
	F3md0vs8.0 = F3m.0-F3m.d8,
	F3pd0vs8.0 = F3p.0-F3p.d8,
	F3md5vs8.0 = F3m.5h-F3m.d8,
	F3pd5vs8.0 = F3p.5h-F3p.d8,
	levels=design) #,
#fit2 <- contrasts.fit(fit, cm); fit2 <- eBayes(fit2)
#fit2.c <- contrasts.fit(fit.combat, cm); fit2.c <- eBayes(fit2.c)

#save(fullCountTable,dge,fit, fit2,mycounts.combat,fit2.c,file=paste(outpath,"Limma-",myname,".RData",sep=""))
load(paste(outpath,"Limma-",myname,".RData",sep=""))

mycounts=mycounts.combat
fit2=fit2.c
#### F3 comparisons
DPvsF3m.0.de=topTable(fit2, adjust="BH",coef="DPvsF3m.0",p.value=0.05,number=10000,lfc=log2(2))
DPvsF3p.0.de=topTable(fit2, adjust="BH",coef="DPvsF3p.0",p.value=0.05,number=10000,lfc=log2(2))
F3mvsF3p.0.de=topTable(fit2, adjust="BH",coef="F3mvsF3p.0",p.value=0.05,number=10000,lfc=log2(2))
DPvsF3m.8.de=topTable(fit2, adjust="BH",coef="DPvsF3m.8",p.value=0.05,number=10000,lfc=log2(2))
DPvsF3p.8.de=topTable(fit2, adjust="BH",coef="DPvsF3p.8",p.value=0.05,number=10000,lfc=log2(2))
F3mvsF3p.8.de=topTable(fit2, adjust="BH",coef="F3mvsF3p.8",p.value=0.05,number=10000,lfc=log2(2))
F3mvsF3p.24h.de=topTable(fit2, adjust="BH",coef="F3mvsF3p.24h",p.value=0.05,number=10000,lfc=log2(2))
F3mvsF3p.5h.de=topTable(fit2, adjust="BH",coef="F3mvsF3p.5h",p.value=0.05,number=10000,lfc=log2(2))

#### diff/undiff
DPd0vs8.0.de=topTable(fit2, adjust="BH",coef="DPd0vs8.0",p.value=0.05,number=10000,lfc=log2(2))
F3md0vs8.0.de=topTable(fit2, adjust="BH",coef="F3md0vs8.0",p.value=0.05,number=10000,lfc=log2(2))
F3pd0vs8.0.de=topTable(fit2, adjust="BH",coef="F3pd0vs8.0",p.value=0.05,number=10000,lfc=log2(2))
F3md5vs8.0.de=topTable(fit2, adjust="BH",coef="F3md5vs8.0",p.value=0.05,number=10000,lfc=log2(2))
F3pd5vs8.0.de=topTable(fit2, adjust="BH",coef="F3pd5vs8.0",p.value=0.05,number=10000,lfc=log2(2))

allDEs=list(DPvsF3m.0.de,DPvsF3p.0.de,F3mvsF3p.0.de,DPvsF3m.8.de,DPvsF3p.8.de,F3mvsF3p.8.de,F3mvsF3p.24h.de,
	F3mvsF3p.5h.de,DPd0vs8.0.de,F3md0vs8.0.de,F3pd0vs8.0.de,F3md5vs8.0.de,F3pd5vs8.0.de)
names(allDEs)=c("DPvsF3m.0","DPvsF3p.0","F3mvsF3p.0","DPvsF3m.8","DPvsF3p.8","F3mvsF3p.8","F3mvsF3p.24h",
		"F3mvsF3p.5h","DPd0vs8","F3md0vs8","F3pd0vs8","F3md5vs8","F3pd5vs8")

results <- decideTests(fit2,adjust.method="BH",p.value=0.05,lfc=log2(2))
##### Write out DE results ####
toprint=cbind(rownames(results[which(!rownames(results)%in%"chrCreERT"),]),
	elementMetadata(spikedgenes84.GR[rownames(results[which(!rownames(results)%in%"chrCreERT"),])])$geneName,
	results[which(!rownames(results)%in%"chrCreERT"),])
colnames(toprint)[1:2]=c("EnsgID","Name")
write.table(toprint,paste(path,"DEResults-Limma-",myname,"-",myset,"-",mydate,".all.txt",sep=""),
	row.names=F,col.names=T,quote=F)

mycols.saved=mycols;mycounts=mycounts.combat
names(mycols)=gsub("[+]","p",gsub("-","m",names(mycols)))

##### vennDiagramms with the results ######### make some overview plots: PCA, heatmaps, clustering, ... ####
pdf(paste(path,"STAR-limma-",myname,".overviewplots.combat.v2.pdf",sep=""),width=7,height=7)
par(mfrow=c(1,1))
###### plots based on y$E values
mymatrixOverviewPlots(mycounts,mycols,paste(myname,"-all",sep=""),bootnr=20)
#### D0 only ####
mymatrixOverviewPlots(mycounts[,grep("[.]0",colnames(mycounts))],mycols,paste(myname,"-D0",sep=""),bootnr=20)
#### D8 only ####
mymatrixOverviewPlots(mycounts[,grep("d8",colnames(mycounts))],mycols,paste(myname,"-D8",sep=""),bootnr=20)
#### without F3+ only ####
mymatrixOverviewPlots(mycounts[,-grep("F3p",colnames(mycounts))],mycols,paste(myname,"-noF3p",sep=""),bootnr=20)
#### without DP only ####
mymatrixOverviewPlots(mycounts[,-grep("DP",colnames(mycounts))],mycols,paste(myname,"-noDP",sep=""),bootnr=20)

######## plots based on differential genes
mymatrix.log=mycounts[unique(unlist(lapply(allDEs,function(x) rownames(x)))),]; par(mfrow=c(1,1))
mymatrixOverviewPlots(mymatrix.log,mycols,paste(myname,"-all-de",sep=""),bootnr=20)
#### D0 only ####
mymatrixOverviewPlots(mymatrix.log[,grep("[.]0",colnames(mymatrix.log))],mycols,paste(myname,"-D0-de",sep=""),bootnr=20)
#### D8 only ####
mymatrixOverviewPlots(mymatrix.log[,grep("d8",colnames(mymatrix.log))],mycols,paste(myname,"-D8-de",sep=""),bootnr=20)
#### without F3+ only ####
mymatrixOverviewPlots(mymatrix.log[,-grep("F3p",colnames(mymatrix.log))],mycols,paste(myname,"-noF3p-de",sep=""),bootnr=20)
#### without DP only ####
mymatrixOverviewPlots(mymatrix.log[,-grep("DP",colnames(mymatrix.log))],mycols,paste(myname,"-noDP-de",sep=""),bootnr=20)

mymatrix.log=mycounts
vennDiagram(results[,c(1:3)]); vennDiagram(results[,c(6:8)])
vennDiagram(results[,c(3:5,8)]); vennDiagram(results[,c(9:13)]); vennDiagram(results[,c(9:11)])
barplot(sapply(allDEs, function(x) dim(x)[1]),main="Nr. of DE regions - FDR 5% FC2",ylab="Nr. genes",las=2,cex.names=0.7)
dev.off()

#### output raw for Figure
write.table(results[,c(1:3,6:8,3:5,8)],paste(path,"FigS4S5j-raw.AregVenns.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(sapply(allDEs, function(x) dim(x)[1]),paste(path,"FigS4e-raw.AregDEs.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


#### set order of the columns
colorder=c("DP10_1.0","DP10_2.0", "DP10_3.0" , "DP10_4.0","DP11_3.0", "DP11_4.0"  ,
	"DP_1.d8", "DP_2.d8","DP_3.d8","DP_4.d8" ,
	"F3m10_1.0","F3m20_1.0","F3m10_2.0" ,"F3m20_2.0","F3m10_3.0","F3m20_3.0","F3m10_4.0", "F3m20_4.0",
	"F3m10_1.5h","F3m10_2.5h","F3m20_2.5h","F3m10_3.5h","F3m20_3.5h","F3m10_4.5h" ,"F3m20_4.5h" ,
	"F3m10_1.24h","F3m10_2.24h" , "F3m20_2.24h","F3m10_3.24h" ,"F3m20_3.24h","F3m10_4.24h","F3m20_4.24h",
	"F3m_1.d8","F3m_2.d8" ,"F3m_3.d8" ,"F3m_4.d8" ,
	"F3p10_1.0","F3p10_2.0" ,"F3p10_3.0","F3p10_4.0" ,"F3p10_1.5h","F3p10_2.5h","F3p10_3.5h" , "F3p10_4.5h",
	 "F3p10_1.24h","F3p10_2.24h",  "F3p10_3.24h", "F3p10_4.24h",
	"F3p_1.d8","F3p_2.d8", "F3p_3.d8","F3p_4.d8"  )

mainname=myname
pdf(paste(path,"STAR-limma-",myname,".heatmaps.v2.pdf",sep=""),width=7,height=7)
#### 1. Validation of SC results & higher adipogenesis in F3- vs. DP? 
# Plot expression of chosen adipogenic genes
plotgenes=unique(c("F3","Mgp","Meox2","Abcg1","Dpp4","Akr1c18","Cd55","Il13ra1","Sparcl1","Adam12","Aoc3","Cd34","Tnfa","Il6","Cd34","Cd31",
	"Pecam1","Mcam","Vcam1","Pparg","Fabp4","Adipoq","Cidec","Cd36","Retn",gsub(" ","",unique(c(apply(mymarkers[c(1:6,9:10),],2,function(x) as.character(x)))))))
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(genes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(genes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Interest.Genes",
	clustering_method="average",scale = "row",fontsize_col=4,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])

# Plot expression of endothelial genes
plotgenes=unique(c(c(gsub(" ","",t(mymarkers["endo",])))))
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(genes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(genes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Interest.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=4,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
#### output raw for Figure
write.table(plotmatrix,paste(path,"FigS4n-raw.EndothelialHeatmap.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


# Plot expression of most F3+ genes before differentiation
#### F3+ low at differentiated stage
plotmatrix=mycounts[rownames(F3mvsF3p.8.de[which(F3mvsF3p.8.de[,1]>0),]),][which(rownames(F3mvsF3p.8.de[which(F3mvsF3p.8.de[,1]>0),])%in%rownames(DPvsF3p.8.de[which(DPvsF3p.8.de[,1]>0),])),colorder]
#DPvsF3m.8.de.len=topTable(fit2, adjust="BH",coef="DPvsF3m.8",p.value=0.1,number=10000,lfc=log2(1.5))
#plotmatrix=mycounts[rownames(DPvsF3m.8.de[which(DPvsF3m.8.de[,1]>0),]),]#[which(rownames(F3mvsF3p.8.de[which(F3mvsF3p.8.de[,1]>0),])%in%rownames(DPvsF3p.8.de[which(DPvsF3p.8.de[,1]>0),]))],colorder]
#plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%unique(unlist(c(GENETYPES84[c(1,2,10,11,16,18,19,20,21,22)])))),]
rownames(plotmatrix)=elementMetadata(genes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F, main="F3+ low D8 selected",
	clustering_method="average",scale = "row",fontsize_col=4,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
#### output raw for Figure
write.table(plotmatrix,paste(path,"FigS4g-raw.F3+lowatdiffHeatmap.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

# Plot expression of most F3+ genes before differentiation
#### F3m vs. F3p D0 & 24h
plotmatrix=mycounts[rownames(F3mvsF3p.0.de[which(F3mvsF3p.0.de[,1]<0),])[which(rownames(F3mvsF3p.0.de[which(F3mvsF3p.0.de[,1]<0),])%in%intersect(rownames(F3mvsF3p.24h.de[which(F3mvsF3p.24h.de[,1]<0),]),
	rownames(F3mvsF3p.5h.de[which(F3mvsF3p.5h.de[,1]<0),])))],colorder]
rownames(plotmatrix)=elementMetadata(genes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="F3+ specific D0, 5h & 24h",
	clustering_method="average",scale = "row",fontsize_col=4,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
#### output raw for Figure
write.table(plotmatrix,paste(path,"Fig2g-raw.F3+highD0D24.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

# Plot expression of most F3+ genes before differentiation
#### F3+ specific TFs, surface markers & secreted factors
plotmatrix=mycounts[rownames(F3mvsF3p.0.de[which(F3mvsF3p.0.de[,1]<0),])[which(rownames(F3mvsF3p.0.de[which(F3mvsF3p.0.de[,1]<0),])%in%rownames(DPvsF3p.0.de[which(DPvsF3p.0.de[,1]<0),]))],]
plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%unique(unlist(c(GENETYPES84[c(1,2,10,11,16,18,19,20,21,22)])))),colorder]
rownames(plotmatrix)=elementMetadata(genes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F, main="F3+ specific D0 selected",
	clustering_method="average",scale = "row",fontsize_col=4,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
#### output raw for Figure
write.table(plotmatrix,paste(path,"FigS3k-raw.F3+highD0TFsecreted.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
dev.off()

pdf(paste(path,"STAR-limma-",myname,".GO.pdf",sep=""),width=7,height=7)
allgenes.background=rownames(mycounts);myp=0.001;showgraph=10; 
write.table(allgenes.background,paste(path,"BackgroundGenes-Limma-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(spikedgenes84.GR[which(names(spikedgenes84.GR)%in%allgenes.background)])$geneName,paste(path,"BackgroundGenes-Limma-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),
	row.names=F,col.names=F,quote=F)

#### F3+ low at differentiated stage
plotmatrix=mycounts[rownames(F3mvsF3p.8.de[which(F3mvsF3p.8.de[,1]>0),]),][which(rownames(F3mvsF3p.8.de[which(F3mvsF3p.8.de[,1]>0),])%in%rownames(DPvsF3p.8.de[which(DPvsF3p.8.de[,1]>0),])),colorder]
myIDs.GO=runGO.84(rownames(plotmatrix),allgenes.background,outpath,ensembl84,"BP","elimCount",myp,paste(mainname,"GO-pcut-",i,sep="_"),ontoGOs.84,showgraph)
write.table(rownames(plotmatrix),paste(path,"F3pLowDiff-Limma-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(spikedgenes84.GR[which(names(spikedgenes84.GR)%in%rownames(plotmatrix))])$geneName,paste(path,"F3pLowDiff-Limma-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
temp=plotmatrix[,colorder]; rownames(temp)=elementMetadata(spikedgenes84.GR[which(names(spikedgenes84.GR)%in%rownames(plotmatrix))])$geneName
write.table(temp,paste(path,"F3pLowDiff-Limma-",myname,"-",myset,"-",mydate,".details.txt",sep=""),
	row.names=T,col.names=T,quote=F)

#### F3m vs. F3p D0 & 24h & 5h
plotmatrix=mycounts[rownames(F3mvsF3p.0.de[which(F3mvsF3p.0.de[,1]<0),])[which(rownames(F3mvsF3p.0.de[which(F3mvsF3p.0.de[,1]<0),])%in%intersect(rownames(F3mvsF3p.24h.de[which(F3mvsF3p.24h.de[,1]<0),]),
	rownames(F3mvsF3p.5h.de[which(F3mvsF3p.5h.de[,1]<0),])))],colorder]
myIDsF3pstr.GO=runGO.84(rownames(plotmatrix),allgenes.background,outpath,ensembl84,"BP","elimCount",myp,paste(mainname,"GO-pcut-F3pstr-",i,sep="_"),ontoGOs.84,showgraph)
write.table(rownames(plotmatrix),paste(path,"F3pStr-Limma-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(spikedgenes84.GR[which(names(spikedgenes84.GR)%in%rownames(plotmatrix))])$geneName,paste(path,"F3pStr-Limma-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
temp=plotmatrix[,colorder]; rownames(temp)=elementMetadata(spikedgenes84.GR[which(names(spikedgenes84.GR)%in%rownames(plotmatrix))])$geneName
write.table(temp,paste(path,"F3pStr-Limma-",myname,"-",myset,"-",mydate,".details.txt",sep=""),
	row.names=T,col.names=T,quote=F)


#### F3+ specific D0
plotmatrix=mycounts[rownames(F3mvsF3p.0.de[which(F3mvsF3p.0.de[,1]<0),])[which(rownames(F3mvsF3p.0.de[which(F3mvsF3p.0.de[,1]<0),])%in%rownames(DPvsF3p.0.de[which(DPvsF3p.0.de[,1]<0),]))],]
myIDsF3pD0.GO=runGO.84(rownames(plotmatrix),allgenes.background,outpath,ensembl84,"BP","elimCount",myp,paste(mainname,"GO-pcut-F3p-",i,sep="_"),ontoGOs.84,showgraph)
write.table(rownames(plotmatrix),paste(path,"F3p-Limma-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(spikedgenes84.GR[which(names(spikedgenes84.GR)%in%rownames(plotmatrix))])$geneName,paste(path,"F3p-Limma-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
temp=plotmatrix[,colorder]; rownames(temp)=elementMetadata(spikedgenes84.GR[which(names(spikedgenes84.GR)%in%rownames(plotmatrix))])$geneName
write.table(temp,paste(path,"F3p-Limma-",myname,"-",myset,"-",mydate,".details.txt",sep=""),
	row.names=T,col.names=T,quote=F)

#D0 F3+ simple
plotmatrix=mycounts[rownames(F3mvsF3p.0.de[which(F3mvsF3p.0.de[,1]<0),]),]
myIDsF3pD0s.GO=runGO.84(rownames(plotmatrix),allgenes.background,outpath,ensembl84,"BP","elimCount",myp,paste(mainname,"GO-pcut-F3p-D0s-",i,sep="_"),ontoGOs.84,showgraph)
write.table(rownames(plotmatrix),paste(path,"F3p-D0s-Limma-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(spikedgenes84.GR[which(names(spikedgenes84.GR)%in%rownames(plotmatrix))])$geneName,paste(path,"F3p-D0s-Limma-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
#5h F3+ simple
plotmatrix=mycounts[rownames(F3mvsF3p.5h.de[which(F3mvsF3p.5h.de[,1]<0),]),]
myIDsF3p5hs.GO=runGO.84(rownames(plotmatrix),allgenes.background,outpath,ensembl84,"BP","elimCount",myp,paste(mainname,"GO-pcut-F3p-5hs-",i,sep="_"),ontoGOs.84,showgraph)
write.table(rownames(plotmatrix),paste(path,"F3p-5hs-Limma-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(spikedgenes84.GR[which(names(spikedgenes84.GR)%in%rownames(plotmatrix))])$geneName,paste(path,"F3p-5hs-Limma-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),
	row.names=F,col.names=F,quote=F)

#24h F3+ simple
plotmatrix=mycounts[rownames(F3mvsF3p.24h.de[which(F3mvsF3p.24h.de[,1]<0),]),]
myIDsF3p24hs.GO=runGO.84(rownames(plotmatrix),allgenes.background,outpath,ensembl84,"BP","elimCount",myp,paste(mainname,"GO-pcut-F3p-24hs-",i,sep="_"),ontoGOs.84,showgraph)
write.table(rownames(plotmatrix),paste(path,"F3p-24hs-Limma-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(spikedgenes84.GR[which(names(spikedgenes84.GR)%in%rownames(plotmatrix))])$geneName,paste(path,"F3p-24hs-Limma-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),
	row.names=F,col.names=F,quote=F)

#5h & 24h F3+ simple
#### F3m vs. F3p D0 & 24h & 5h
plotmatrix=mycounts[rownames(F3mvsF3p.5h.de[which(F3mvsF3p.5h.de[,1]<0),])[which(rownames(F3mvsF3p.5h.de[which(F3mvsF3p.5h.de[,1]<0),])%in%intersect(rownames(F3mvsF3p.24h.de[which(F3mvsF3p.24h.de[,1]<0),]),
	rownames(F3mvsF3p.5h.de[which(F3mvsF3p.5h.de[,1]<0),])))],colorder]
myIDsF3pinvitro.GO=runGO.84(rownames(plotmatrix),allgenes.background,outpath,ensembl84,"BP","elimCount",myp,paste(mainname,"GO-pcut-F3pinvitro-",i,sep="_"),ontoGOs.84,showgraph)
write.table(rownames(plotmatrix),paste(path,"F3pInvitro-Limma-",myname,"-",myset,"-",mydate,".ensgID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
write.table(elementMetadata(spikedgenes84.GR[which(names(spikedgenes84.GR)%in%rownames(plotmatrix))])$geneName,paste(path,"F3pInvitro-Limma-",myname,"-",myset,"-",mydate,".ID.txt",sep=""),
	row.names=F,col.names=F,quote=F)
temp=plotmatrix[,colorder]; rownames(temp)=elementMetadata(spikedgenes84.GR[which(names(spikedgenes84.GR)%in%rownames(plotmatrix))])$geneName
write.table(temp,paste(path,"F3pInvitro-Limma-",myname,"-",myset,"-",mydate,".details.txt",sep=""),
	row.names=T,col.names=T,quote=F)

# Plot expression of interest genes:
#daniel=elementMetadata(hsagenes84.GR[c("ENSG00000173546", "ENSG00000143248", "ENSG00000175084", "ENSG00000113721")])$geneName
daniel=c("Cspg4","Rgs5","Des","Pdgfrb")
tgfbgenes="GO:0007179"; fgfgenes="GO:0008543"; vegfgenes="GO:0038084"
adipogenes=c("GO:0045444", "GO:0060612")
hhgenes=c("GO:0007224")

tgfbgenes=getBM( c("external_gene_name","ensembl_gene_id"),filters=c("go_id"),values=tgfbgenes,mart = ensembl84)[,1]
fgfgenes=getBM( c("external_gene_name","ensembl_gene_id"),filters=c("go_id"),values=fgfgenes,mart = ensembl84)[,1]
vegfgenes=getBM( c("external_gene_name","ensembl_gene_id"),filters=c("go_id"),values=vegfgenes,mart = ensembl84)[,1]
adipogenes=getBM( c("external_gene_name","ensembl_gene_id"),filters=c("go_id"),values=adipogenes,mart = ensembl84)[,1]
hhgenes=getBM( c("external_gene_name","ensembl_gene_id"),filters=c("go_id"),values=hhgenes,mart = ensembl84)[,1]

# Plot expression of chosen TGFB genes
plotgenes=daniel
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(spikedgenes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(spikedgenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Daniel.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])

# Plot expression of tgfb genes
plotgenes=tgfbgenes
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(spikedgenes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(spikedgenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Tgfb.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
#### only F3+ vs. F3- 
plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%elementMetadata(spikedgenes84.GR[unique(unlist(sapply(allDEs[c(3,6,7,8)],function(x) rownames(x))))])$geneName),]
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Tgfb.Genes sign.DE",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
plotmatrix=results[which(rownames(results)%in%getEnsgID(spikedgenes84.GR,rownames(plotmatrix))),]
rownames(plotmatrix)=elementMetadata(spikedgenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Tgfb.Genes sign.DE",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "none",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1)

# Plot expression of fgf genes
plotgenes=fgfgenes
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(spikedgenes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(spikedgenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Fgf.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
#### only F3+ vs. F3- 
plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%elementMetadata(spikedgenes84.GR[unique(unlist(sapply(allDEs[c(3,6,7,8)],function(x) rownames(x))))])$geneName),]
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Fgf.Genes sign.DE",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
plotmatrix=results[which(rownames(results)%in%getEnsgID(spikedgenes84.GR,rownames(plotmatrix))),]
rownames(plotmatrix)=elementMetadata(spikedgenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Fgf.Genes sign.DE",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "none",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1)

# Plot expression of vegf genes
plotgenes=vegfgenes
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(spikedgenes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(spikedgenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Vegf.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
#### only F3+ vs. F3- 
plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%elementMetadata(spikedgenes84.GR[unique(unlist(sapply(allDEs[c(3,6,7,8)],function(x) rownames(x))))])$geneName),]
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Vegf.Genes sign.DE",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
plotmatrix=results[which(rownames(results)%in%getEnsgID(spikedgenes84.GR,rownames(plotmatrix))),]
rownames(plotmatrix)=elementMetadata(spikedgenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Vegf.Genes sign.DE",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "none",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1)


# Plot expression of Hh genes
plotgenes=hhgenes
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(spikedgenes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(spikedgenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Hh.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
#### only F3+ vs. F3- 
plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%elementMetadata(spikedgenes84.GR[unique(unlist(sapply(allDEs[c(3,6,7,8)],function(x) rownames(x))))])$geneName),]
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Hh.Genes sign.DE",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
plotmatrix=results[which(rownames(results)%in%getEnsgID(spikedgenes84.GR,rownames(plotmatrix))),]
rownames(plotmatrix)=elementMetadata(spikedgenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Hh.Genes sign.DE",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "none",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1)


# Plot expression of adipo genes
plotgenes=adipogenes
plotmatrix=mycounts[rownames(mycounts)[which(rownames(mycounts)%in%getEnsgID(spikedgenes84.GR,plotgenes))],colorder]
rownames(plotmatrix)=elementMetadata(spikedgenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Adipo.Genes",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
#### only F3+ vs. F3- 
plotmatrix=plotmatrix[which(rownames(plotmatrix)%in%elementMetadata(spikedgenes84.GR[unique(unlist(sapply(allDEs[c(3,6,7,8)],function(x) rownames(x))))])$geneName),]
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Adipo.Genes sign.DE",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "row",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1,ColSideColors=mycols[colnames(plotmatrix)])
plotmatrix=results[which(rownames(results)%in%getEnsgID(spikedgenes84.GR,rownames(plotmatrix))),]
rownames(plotmatrix)=elementMetadata(spikedgenes84.GR[rownames(plotmatrix)])$geneName
pheatmap(plotmatrix,clustering_distance_rows="correlation",cluster_cols = F,main="Adipo.Genes sign.DE",
	col=rev(colorRampPalette(c("tomato4","sienna2","ivory","steelblue1","steelblue3"))(216)),
	clustering_method="average",scale = "none",fontsize_col=6,fontsize_row=7,show_rownames=T,cexRow=0.1,cexCol=0.1)
dev.off()


#### Save counts & matrix for other processing 
fullCountTable.brbseq=fullCountTable
normcounts.brbseq=mycounts
save(normcounts.brbseq,fullCountTable.brbseq,results,allDEs, file=paste(outpath,"Limma-",myname,".BrbseqCounts.RData",sep=""))
#mapped/brbseq/STAR/plots/Limma-Nov16.BrbseqCounts.RData

#### Keep any cytokine or secreted factor ####
#### uniprot_swissprot
tab84.uniprot = getBM( c("ensembl_gene_id","external_gene_name","uniprot_swissprot"),mart = ensembl84)
secreted.cur=read.table("myarchive/utils/funseckb2_search_results_currated.txt",sep="\t",colClasses=c("V1"="character"))
secreted.highconf=read.table("myarchive/utils/funseckb2_search_results_highconf.txt",sep="\t",colClasses=c("V1"="character"))
secreted.lowconf=read.table("myarchive/utils/funseckb2_search_results_lowconf.txt",sep="\t",colClasses=c("V1"="character"))
secreted.medconf=read.table("myarchive/utils/funseckb2_search_results_medconf.txt",sep="\t",colClasses=c("V1"="character"))
allsecreted=c(secreted.cur[,1],secreted.highconf[,1],secreted.medconf[,1])
allsecreted.IDs=tab84.uniprot[which(tab84.uniprot[,3]%in%allsecreted),1]




