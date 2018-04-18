
##### General overview of samples: MDS, PCA & clustering
mymatrixOverviewPlots = function(mymatrix.log,mycols,myname,bootnr=20) {
	plotMDS(scale(mymatrix.log),col=mycols[colnames(mymatrix.log)],top=500,cex=0.5,
		main=paste("MDS scaled",myname,sep=""))
	mypca=PCA(t(scale(mymatrix.log)),graph=FALSE)
	plot.PCA(mypca,col.hab=mycols[colnames(mymatrix.log)],choix="ind",habillage="ind",title=paste(myname," scaled",sep=""))
	plot.PCA(mypca,choix="var",title=paste(myname," scaled",sep=""),select="contrib 20",col.var="gray20",unselect=0.3)
	par(mfrow=c(1,1)); #bootnr=20
	fit.pvclust <- pvclust(scale(mymatrix.log), method.hclust="ward.D2", method.dist="correlation",nboot=bootnr)
	plot(fit.pvclust ,cex=0.5) # dendogram with p values
	fit.pvclust <- pvclust(mymatrix.log, method.hclust="ward.D2", method.dist="correlation",nboot=bootnr)
	plot(fit.pvclust ,cex=0.5) # dendogram with p values
}


##### read the htseq aigned reads ####
#readSTARHTSeqASC(starpath,"cut",mynames,names(mynames),mydate) 
readSTARHTSeqASC = function(starpath,mypattern,mynames,shortnames,mydate) {
	myfiles=list.files(starpath,pattern=mypattern)
	myfiles=myfiles[which(!myfiles%in%myfiles[grep("Transcriptome",myfiles)])]
	myfiles=sapply(mynames,function(x) myfiles[grep(x,myfiles)])
#	myfiles=myfiles[sapply(mysplit(mynames,".fastq-cut",2,1),function(x) grep(x,myfiles))]
	readmatrix=matrix(data=NA,ncol=length(myfiles),nrow=47824)
	readstats=matrix(data=NA,ncol=length(myfiles),nrow=5)
	kok=NULL
	for(i in 1:length(myfiles)) {
         	truefile=list.files(paste(starpath,myfiles[i],"/",sep=""),pattern="gene")
		info = file.info(paste(starpath,myfiles[i],"/",truefile,sep=""))[1]
		if (info[[1]]>0) {
		kok=c(kok,i)
		temp=read.table(paste(starpath,myfiles[i],"/",truefile,sep=""),sep="\t",header=FALSE)
         	readmatrix[,i]=temp[1:47824,2]
		readstats[,i]=temp[47825:47829,2]}		
		}
	rownames(readmatrix)=temp[1:47824,1]
	rownames(readstats)=temp[47825:47829,1]
	readexons=matrix(data=NA,ncol=length(myfiles),nrow=453487)
	for(i in 1:length(myfiles)) {
         	truefile=list.files(paste(starpath,myfiles[i],"/",sep=""),pattern="exon")
		info = file.info(paste(starpath,myfiles[i],"/",truefile,sep=""))[1]
		if (info[[1]]>0) {
		temp=read.table(paste(starpath,myfiles[i],"/",truefile,sep=""),sep="\t",header=FALSE)
         	readexons[,i]=temp[1:453487,2]
		}
		}
	rownames(readexons)=temp[1:453487,1]
	colnames(readmatrix)=colnames(readexons)=shortnames
	htseqgenes=readmatrix
	htseqexons=readexons
	htseqstats=t(rbind(apply(htseqgenes,2,function(x) sum(x)),readstats))
	colnames(htseqstats)=c("Genes","No_genes","Ambig","Low_qual","Not_align","Not_unique")
	save(htseqgenes,htseqexons,htseqstats,file=paste(outpath,"SSAll-STARHTSeqTags.-",mydate,".RData",sep=""))
}

#### Interesting genes
mymarkers=read.table("utils/Sep16_markers.txt",head=F,sep=",",fill=T)
mymarkers=mymarkers[-9,]; rownames(mymarkers)=mysplit(mymarkers[,1],":",2,1)
mymarkers[,1]=mysplit(mymarkers[,1],": ",2,2)
allmarkers=gsub(" ","",unlist(apply(mymarkers[1:7,],1,function(x) as.vector(x))))

myMarkerAll=unique(c("chrtdRFP","Dlk1","Fabp4", "Sox9","Peg3",
	"Pparg", "Runx2","Cd36","Cxcl9","Mgp","Meox2","C7","F3","Tspan11","Adcyap1r1",
	"Abcc9","Abcg1","Cd34","Cd29", "Cd200","Cd45","Cd146","Cd90","Cd73","Cd105","Cd44",
	"Bmp4","Bmp3","Bmp5","Bmp6","Bmp7","Bmp2","Bmp8a","Bmp8b","Bmp9","Bmp10","Bmpr1a","Bmpr2","Bmpr1b","BmprII",
	"Tgfbi","Akr1c18","Qk","Sparcl1","Sept11","Duoxa1","Col5a3",
	"Zfp423","Myf5","Cd137","Tmem26","Tbx1","Prdm16","Ppargc1a","Fgf21",
	"L6c1","Bgn","Cxcl1","Cxcl12","Cxcl13","C3","Chrdl1","Plin2","Plin4","Cd24a","Mmd2","Tusc5","Hdx",
	"Acta2", "Tagln", "Myh11", "Myl9","Zfp423","Myf5","Cnn1","Pi16","Adipoq", "Pdgfra","Pdgfrb","Il4ra","Il4",
	"Il5","Il13","Il33","Il18","Ccl2","Ccr2","Prex1","Ednrb1","Pdk4","Dio2","Fabp3","Cxcl14",
	"Adam12","Lpl","Cdh11","Tmem100","Vcam1","Cd36","Mmp3","Chrna3","Chrna2","Chrna1","Myoc","Oc","Alpl","Bsp","Dcn","Bglap",
	"Hes1","Egr2", "Meg3","Ryan","Aoc3","Gpm6b","Lum","Ramp1","Epha4","Cd9","Sulf1","Thbd","Islr2",
	"Cd14","Cd38","Mas1","Csf1r","Creb5","Hdac7","Il33","Il13ra1","Gpr133","Tmem158","Tgfbr2","Scara5","Dpp4","Ppap2b",
	"Scara3","Cd55","Car8","Axl","Anxa1","Fosb","Bhlhe40","Cd9","Cd63","Mme","Apoe",
	"Sparcl1","Col15a1","Col4a2","Prnp","Lrpap1","Mdh2","Gdf10","Anxa2","Anxa3","Sept2","Csde1","Nptn","Lgr4",
	"Snap23","Lrpap1","Plxna1","Adam11","Zkscan4","Creb3l1","Zfp28","Dnmt3b","Slc12a9","Nt5e","Wnt16","Slc25a36","Tek",
	"Ap2m1","Vpm1","Lgals9","Heg1","Plxna4","Hoxc5","Zfp30","Bcl6","Sox6","Notch2","Pla1a","Tek","Fndc1","Fosb","Irf1","Ugt1a6a",
	"Chrdl1","Fgf7","Nrep","Ptprc","Acta2","Pecam1","Mcam","Cd44","Cd90","Ly6a","Cd29","Eng","Zfp423","Mki67","Tie2","mt-Rnr2",
	"mt-Nd1","mt-Co1","mt-Co3","Malat1","mt-Atp6","Inmt","Ucp1","Ucp2","Crem","Hoxc5","Nfe2l3","Prrx2","Stat5a","Cenpv","Ncor","Ikzf2","Mxi1",
	"Ybx2","Bcl6","Sox6","Sox18","Nfia","Thrsp","Pcx","Pck1","Impa2","Ces1d", "Grem2","Arx","Slc11a1","Kcns3","Hhip","Ebf2","Edn1","Dbp","Hopx",
	"Postn","March6", "Arx","H2afz","Hmgb1","Hmgn2","Cdh5","Atp5g3","Atp5b","Atp5g1","Atp5g2", "Sox2","Nanog","Pou5f1","Lep","Pparg",
	"Vegfa","Alp1","Osx","Igf1","Abcg1","Abcg2","Abcg3","Abcg4","Abcg5","Abcg8","Aldh1","Aldh1a1","Aldh1b1","Aldh1l2","Aldh1l1","Aldh2","Aldh1a2", 
	"Myh11","Myocd","Acta2","Tagln","Cnn1","Mkl2"))
	#,unlist(genesDiff.TFs[c(2,3,6)]),unlist(genesDiff.Membrane[c(2,3,6)])))

mychoiceMarkerAll=unique(c("chrtdRFP","Dlk1","Cd34","Pdgfra","Pdgfrb",
"Zfp423","Cd24a","Mki67","Pecam1","Cd44",
"Cxcl12","Cdh11","Lpl","Slc25a5","Lum",
"Sparcl1","Col4a2","Col15a1","Egr2","Apoe",
"Bhlhe40","Fabp4","Peg3","Adam12","Aoc3","Wt1",
"Pparg","Cd36","Mme","Gdf10","Adcyap1r1",
"Mgp","Meox2","F3","Tspan11","Alpl",
"C7","C2","Serping1","Igf1","Car8",
"Bmp7","Bmp2","Duoxa1","Il18","Il33",
"Sulf1","Il13ra1","Tek","Creb5","Hdac7",
"Notch2","Cxcl13","Pi16","Scara3","Bmp7",
"Plxna4","Adam11","Hoxc5","Zkscan4","Nt5e"))

#runGO.84(myIDs,allgenes.background,plotpath,ensembl84,"BP","elimCount",myp,paste(mainname,"GO-pcut-",i,sep="_"),ontoGOs.84,showgraph)
runGO.84 <- function(genes,mybackground,path,ensembl,ont,algorithm,myp,set,ontoGOs,showgraph) {
	mygenes <-getBM(attributes=c("ensembl_gene_id","external_gene_name","description"), filters="ensembl_gene_id", values=genes, mart=ensembl)
	mm10GO <- do.call("rbind", ontoGOs)
	mm10GO=mm10GO[which(mm10GO[,1]%in%mybackground),]
	mm10.gene2GO <- with(mm10GO, split(go, ensembl_gene_id))
	hRes <- sigGOTable(mygenes$ensembl_gene_id,names(mm10.gene2GO),mm10.gene2GO,ont,myp,algorithm,set,showgraph)
	write.table(hRes,paste(path,"GO/",set,"_topGO_",myp,".",ont,".",algorithm,".txt",sep=""),append=FALSE,sep="\t",quote=FALSE, eol="\n",na="NA",dec=".",row.names=TRUE,col.names=TRUE)
		
	myGO=hRes[,1]
	mygenes=cbind(mygenes,"")
	mygenes=cbind(mygenes,"")
	colnames(mygenes)[4:5]=c("GO.ID","GO.Description")
	
	mygenes=as.matrix(mygenes)
	for(i in 1:dim(hRes)[1]) {
		x=hRes[i,1]
		x.info=mm10GO[which(mm10GO[,2]%in%x),]
		mygenes[which(mygenes[,1]%in%x.info[,1]),4]=paste(mygenes[which(mygenes[,1]%in%x.info[,1]),4],hRes[i,1],sep=";")
		mygenes[which(mygenes[,1]%in%x.info[,1]),5]=paste(mygenes[which(mygenes[,1]%in%x.info[,1]),5],hRes[i,2],sep=";")
	}	
	return(mygenes)
}

sigGOTable <- function(selGenes, GOgenes, gene2GO,ont, maxP,mymethod,set,showgraph){
	inGenes <- factor(as.integer(GOgenes %in% selGenes))
	names(inGenes) <- GOgenes
	GOdata <- new("topGOdata", ontology=ont, allGenes=inGenes,annot=annFUN.gene2GO, gene2GO=gene2GO)
	myTestStat <- new(mymethod, testStatistic=GOFisherTest,name="Fisher test", cutOff=maxP)
	mySigGroups <- getSigGroups(GOdata, myTestStat)
	sTab <- GenTable(GOdata, mySigGroups, topNodes=length(usedGO(GOdata)))
	names(sTab)[length(sTab)] <- "p.value"
	sTab <- subset(sTab, as.numeric(p.value) < maxP)
	sTab$Term <- sapply(mget(sTab$GO.ID, env=GOTERM), Term)
	if ((showgraph%in%c(4:10))&(dim(sTab)[1]>0)) {
		showSigOfNodes(GOdata, topGO::score(mySigGroups), firstSigNodes = showgraph, useInfo = "def")
		textplot( sTab, valign="top"  )
		title(paste(set,ont,sep="_"))
	}
	return(sTab)}

########## GET THE ENSEMBL ID OF A GENE #########
getEnsgID <- function(genes.GR,geneName) {
	return(names(genes.GR[which(elementMetadata(genes.GR)$geneName%in%geneName)]))
}

valuesToColorsAbs=function (myvector,choosecols,mynames,mymin,mymax) {
		       myvector=c(mymin,mymax,myvector); names(myvector)[1:2]=c("A1","A2"); mynames=c("A1","A2",mynames)
                       myvector[which(is.na(as.numeric(as.vector(myvector))))]=min(myvector[which(!is.na(as.numeric(as.vector(myvector))))])
  		       myvector=round(myvector,digits=2)
		       colormap=seq(min(myvector),max(myvector),by=0.01)
                       colormap=colorRampPalette(choosecols)(length(colormap))
                       names(colormap)=round(seq(min(myvector),max(myvector),by=0.01),digits=2)
                       rpfmycol=as.vector(colormap[match(myvector,names(colormap))])
                       names(rpfmycol)=mynames
                       return(rpfmycol[which(!names(rpfmycol)%in%c("A1","A2"))])
               }