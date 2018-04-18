mysplit <- function (vector,splitsign,numbers,position) {
	return(unlist(strsplit(as.character(vector),splitsign))[seq(position,length(vector)*numbers,by=numbers)])
}

replaceMinus <- function(mytemp) {
	mytemp[which(mytemp < c(-5))]= -5
	return(mytemp)
}
######## function correcting a matrix: -inf to min of row, inf to max of row, nan to NA
replaceValuesCufflinks=function(x) {
	x[which(x<-16.8525)]=-16.8525
	return(x)
}

exploratoryExpress <- function(genes.expr,path,myset,what,mynames,pcacols,interestGenes,genes.expr.withnames,mylimit,mytypeY,myscale,genes.expr.mcse,genes.expr.minu,allin,toscale) {
	temp=genes.expr
	if(toscale %in%"scale") {
	#	myscaling=apply(genes.expr,2,median)/mean(apply(genes.expr,2,median))
	#	temp=sapply(1:dim(genes.expr)[2],function(x) genes.expr[,x]/myscaling[x])
		temp=normalizeBetweenArrays(temp,method="quantile")
		colnames(temp)=mynames
	}
	genes.expr=temp
	
	my.abs = abs(cor(genes.expr))
	write.table(my.abs,paste(path,myset,"-",what,"EstimateCorrelations.txt",sep=""),sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)

	pdf(paste(path,myset,"-exploratory",what,"Plots.pdf",sep=""),width=10,height=10)
	
	boxplot(genes.expr,main=what,col="darkgrey",names=colnames(genes.expr),ylim=c(min(min(apply(tab,2,function(x) summary(x))[2,])/2,min(apply(genes.expr,2,function(x) summary(x))[2,])*2),
		max(apply(genes.expr,2,function(x) summary(x))[5,]+apply(genes.expr,2,function(x) (summary(x)[5]-summary(x)[2])*2))))
	#display all the values as histogramms
	mydim=min(as.integer(dim(genes.expr)[2]/3),3)
	par(mfrow=c(mydim,3))
	for (i in 1:dim(genes.expr)[2]) {
		truehist(genes.expr[,i],main=colnames(genes.expr)[i],col="grey",xlab=c("Expression Estimate"))
	}
	# display the correlations as heatmaps
	heatmap.2(my.abs,scale="none",density.info="none",col=rev(heat.colors(200)),cellnote=apply(my.abs,1,function(x) prettyNum(x,digits=3)),
		notecex=1,notecol="black",breaks=seq(0,1,by=0.005),trace="none",cexRow=0.9,cexCol=0.9)

	#display the correlations per sample
	mylim=dim(genes.expr)[2]
	myplots=list(NA)
	for (i in 1:(mylim-1)) {
		for (j in (i+1):mylim) {
			mycor=cor(genes.expr[,i],genes.expr[,j])
			mytitle=paste(paste(colnames(genes.expr)[c(i,j)],collapse=":"),prettyNum(mycor,digits=3),sep=":",collapse=":")
			smoothScatter(genes.expr[,i],genes.expr[,j],main=mytitle,xlab=colnames(genes.expr)[i],ylab=colnames(genes.expr)[j])
			}
	}
	
	# plot PCA Now
	par(mfrow=c(1,1))
	plotPCA(genes.expr,addtext=mynames,legend=FALSE,col=c(pcacols),pch=rep(19,length(mynames)))
	plotPCA(apply(genes.expr, 2,function(x) replaceMinus(x)),addtext=mynames,legend=FALSE,col=c(pcacols),pch=rep(19,length(mynames)))
	
	#### display clustering for selected sets of genes
	if (mytypeY%in%"long") {
		for (i in 1:length(interestGenes)) {
			select=which(genes.expr.withnames[,dim(genes.expr.withnames)[2]]%in%interestGenes[[i]])
			if (length(select)>mylimit) {
				select=select[1:mylimit]
			}
		#	hmcol = colorRampPalette(brewer.pal(9, "Reds"))(100)
			hmcol=rev(heat.colors(256))
			temp=genes.expr[select,]
			rownames(temp)=genes.expr.withnames[match(rownames(temp),rownames(genes.expr.withnames)),dim(genes.expr.withnames)[2]]
			if (dim(temp)[1]>3) {
			blueredcol=bluered(200)
			heatmap.2(apply(temp, 2,function(x) replaceMinus(x)), cexRow=0.5, col = blueredcol[61:213], scale=myscale,
				breaks=seq(-5,15,by=0.13),dendrogram="row",Colv=FALSE,trace="none", margin=c(10, 6),main=names(interestGenes)[i])
		}
		}
		dev.off()	
					
			
	##### display the expression for selected sets of genes (with errors) ######
		if (allin%in%"allin") {
		pdf(paste(path,myset,"-exploratory",what,"LinePlots.pdf",sep=""),width=12,height=6)
		for (i in 1:length(interestGenes)) {
			select=which(genes.expr.withnames[,dim(genes.expr.withnames)[2]]%in%interestGenes[[i]])
			if (length(grep("mmseq",myset))==1) {
				temp=genes.expr[select,c(1:9)]
				temp.mcse=genes.expr.mcse[select,c(1:9)]
				temp=apply(temp, 2,function(x) replaceMinus(x))
				if (length(select)>mylimit) {
					mymeans=cbind(apply(temp[,1:3],1,mean),apply(temp[,4:6],1,mean),apply(temp[,7:9],1,mean))
					mydiff=apply(mymeans,1,function(x) max(abs(x[1]-x[2]),abs(x[1]-x[3]),abs(x[2]-x[3])))					
#					mydiff=log((temp[,1]+1)/(temp[,3]+1))
					temp=temp[order(-(mydiff))[1:mylimit],]
					temp.mcse=temp.mcse[order(-(mydiff))[1:mylimit],]
				}
				myseq=1:length(temp[,4])
				par(mfrow=c(1,2))
				temp=temp[order(temp[,4]),]
				temp.mcse=temp.mcse[order(temp[,4]),]
				##### first one is brite-brown
				plotCI(myseq,temp[,4],uiw=temp.mcse[,4],liw=temp.mcse[,4],col=pcacols[4],pch=19,main=names(interestGenes)[i],
					xlab=names(interestGenes)[i],barcol=pcacols[4],gap=0,ylab="MMseq expression",xlim=c(0,c(length(myseq)+as.integer(length(myseq)/10))))
				for (k in 5:6) {
					plotCI(myseq,temp[,k],uiw=temp.mcse[,k],liw=temp.mcse[,k],col=pcacols[k],add=TRUE,pch=19,barcol=pcacols[k],gap=0)
				}
				for (k in 1:3) {
					plotCI(myseq,temp[,k],uiw=temp.mcse[,k],liw=temp.mcse[,k],col=pcacols[k],add=TRUE,pch=19,barcol=pcacols[k],gap=0)
				}
				mymeans=cbind(apply(temp[,1:3],1,median),apply(temp[,4:6],1,median),apply(temp[,7:9],1,median))
				tolabel=which((c(c(mymeans[,2]+1)/c(mymeans[,1]+1))<=0.7692308)|(c(c(mymeans[,2]+1)/c(mymeans[,1]+1))>=1.3))
				if (length(tolabel)>0) {
					text(myseq[tolabel],temp[tolabel,2],labels=genes.expr.withnames[match(rownames(temp)[tolabel],rownames(genes.expr.withnames)),dim(genes.expr.withnames)[2]],cex=0.8,adj=0.1,pos=4,offset=0.2)
				}
				legend("topleft",c("brown","brite","white"),col=c("tan4","bisque3","slategray1"),lwd=2)
					
				##### second one is brite-white
				temp=temp[order(temp[,4]),]
				temp.mcse=temp.mcse[order(temp[,4]),]
				plotCI(myseq,temp[,4],uiw=temp.mcse[,4],liw=temp.mcse[,4],col=pcacols[4],pch=19,main=names(interestGenes)[i],
					xlab=names(interestGenes)[i],barcol=pcacols[4],gap=0,ylab="MMseq expression",xlim=c(0,c(length(myseq)+as.integer(length(myseq)/10))))
				for (k in 5:6) {
						plotCI(myseq,temp[,k],uiw=temp.mcse[,k],liw=temp.mcse[,k],col=pcacols[k],add=TRUE,pch=19,barcol=pcacols[k],gap=0)
				}
				for (k in 7:9) {
					plotCI(myseq,temp[,k],uiw=temp.mcse[,k],liw=temp.mcse[,k],col=pcacols[k],add=TRUE,pch=19,barcol=pcacols[k],gap=0)
				}
				mymeans=cbind(apply(temp[,1:3],1,median),apply(temp[,4:6],1,median),apply(temp[,7:9],1,median))
				tolabel=which((c(c(mymeans[,2]+1)/c(mymeans[,3]+1))<=0.7692308)|(c(c(mymeans[,2]+1)/c(mymeans[,1]+1))>=1.3))
				if (length(tolabel)>0) {
					text(myseq[tolabel],temp[tolabel,8],labels=genes.expr.withnames[match(rownames(temp)[tolabel],rownames(genes.expr.withnames)),dim(genes.expr.withnames)[2]],cex=0.8,adj=0.1,pos=4,offset=0.2)
				}
				legend("topleft",c("brown","brite","white"),col=c("tan4","bisque3","slategray1"),lwd=2)
			} else {
					temp=genes.expr[select,c(1:9)]
					temp=apply(temp, 2,function(x) replaceMinus(x))
					if (length(select)>mylimit) {
						mymeans=cbind(apply(temp[,1:3],1,mean),apply(temp[,4:6],1,mean),apply(temp[,7:9],1,mean))
						mydiff=apply(mymeans,1,function(x) max(abs(x[1]-x[2]),abs(x[1]-x[3]),abs(x[2]-x[3])))					
	#					mydiff=log((temp[,1]+1)/(temp[,3]+1))
						temp=temp[order(-(mydiff))[1:mylimit],]
					}
					myseq=1:length(temp[,4])
					par(mfrow=c(1,2))
					temp=temp[order(temp[,4]),]
					##### first one is brite-brown
					plotCI(myseq,temp[,4],col=pcacols[4],pch=19,main=names(interestGenes)[i],
						xlab=names(interestGenes)[i],barcol=pcacols[4],gap=0,ylab="MMseq expression",xlim=c(0,c(length(myseq)+as.integer(length(myseq)/10))))
					for (k in 5:6) {
						plotCI(myseq,temp[,k],col=pcacols[k],add=TRUE,pch=19,barcol=pcacols[k],gap=0)
					}
					for (k in 1:3) {
						plotCI(myseq,temp[,k],col=pcacols[k],add=TRUE,pch=19,barcol=pcacols[k],gap=0)
					}
					mymeans=cbind(apply(temp[,1:3],1,median),apply(temp[,4:6],1,median),apply(temp[,7:9],1,median))
					tolabel=which((c(c(mymeans[,2]+1)/c(mymeans[,1]+1))<=0.7692308)|(c(c(mymeans[,2]+1)/c(mymeans[,1]+1))>=1.3))
					if (length(tolabel)>0) {
						text(myseq[tolabel],temp[tolabel,2],labels=genes.expr.withnames[match(rownames(temp)[tolabel],rownames(genes.expr.withnames)),dim(genes.expr.withnames)[2]],cex=0.8,adj=0.1,pos=4,offset=0.2)
					}
					legend("topleft",c("brown","brite","white"),col=c("tan4","bisque3","slategray1"),lwd=2)

					##### second one is brite-white
					temp=temp[order(temp[,4]),]
					plotCI(myseq,temp[,4],col=pcacols[4],pch=19,main=names(interestGenes)[i],
						xlab=names(interestGenes)[i],barcol=pcacols[4],gap=0,ylab="MMseq expression",xlim=c(0,c(length(myseq)+as.integer(length(myseq)/10))))
					for (k in 5:6) {
							plotCI(myseq,temp[,k],col=pcacols[k],add=TRUE,pch=19,barcol=pcacols[k],gap=0)
					}
					for (k in 7:9) {
						plotCI(myseq,temp[,k],col=pcacols[k],add=TRUE,pch=19,barcol=pcacols[k],gap=0)
					}
					mymeans=cbind(apply(temp[,1:3],1,median),apply(temp[,4:6],1,median),apply(temp[,7:9],1,median))
					tolabel=which((c(c(mymeans[,2]+1)/c(mymeans[,3]+1))<=0.7692308)|(c(c(mymeans[,2]+1)/c(mymeans[,1]+1))>=1.3))
					if (length(tolabel)>0) {
						text(myseq[tolabel],temp[tolabel,8],labels=genes.expr.withnames[match(rownames(temp)[tolabel],rownames(genes.expr.withnames)),dim(genes.expr.withnames)[2]],cex=0.8,adj=0.1,pos=4,offset=0.2)
					}
					legend("topleft",c("brown","brite","white"),col=c("tan4","bisque3","slategray1"),lwd=2)
			}
		} 
		dev.off()
		} 
	} else {
		dev.off()
	}
}	


exploratoryDESeqAll <- function(path,myset,what,cds,nrtests,genes.GR,selectgenes,ugroups, interestGenes, mylimit,sign,sign2,allin,pcacols) {
		
	pdf(paste(path,myset,"-DESeq",what,"Plots.pdf",sep=""),width=10,height=10)
	boxplot(counts(cds),main=what,col="darkgrey",names=colnames(counts(cds)),cex.names=0.6,las=2)
	
	plotDispEsts(cds)
	
	cdsFullBlind = estimateDispersions( cds, method = "blind" )
	vsdFull = varianceStabilizingTransformation( cds )

	#select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:100]
	#take only adipogenic genes#
	for (i in 1:length(interestGenes)) {
		specgenes=names(genes.GR)[which(elementMetadata(genes.GR)$gene%in%interestGenes[[i]])]
		select=which(rownames(counts(cds))%in%specgenes)
		if (length(select)>mylimit) {
			select=select[1:mylimit]
		}
	#	hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
		hmcol=rev(heat.colors(256))
		temp=exprs(vsdFull)[select,]
		rownames(temp)=elementMetadata(genes.GR[match(rownames(exprs(vsdFull)[select,]),names(genes.GR))])$gene
		heatmap.2(temp, col = hmcol, trace="none", margin=c(10, 6),main=names(interestGenes)[i])
	}
	
	dists = dist( t( exprs(vsdFull) ) )
	mat = as.matrix( dists )
	rownames(mat) = colnames(mat) = with(pData(cdsFullBlind), condition)
	heatmap.2(mat, trace="none", col = hmcol, margin=c(13, 13))
#	print(plotPCA2(vsdFull,mylimit),main="vdsFull")
	plotPCA(exprs(vsdFull),addtext=mynames,legend=FALSE,col=c(pcacols[c(1:2,4:5,7:8)]),pch=rep(19,6))
	
#	print(plotPCA2(cdsFullBlind),main="cdsFullBlind")
	
	for(i in 1:nrtests) {
		ncu = counts( cds, normalized=TRUE )[ , conditions(cds)==ugroups[i] ]
		plotMA(data.frame(baseMean = rowMeans(ncu),log2FoldChange = log2( ncu[,2] / ncu[,1] )), col = "black",main=ugroups[i])
	}
	dev.off()
	
	if (allin%in%"allin") {
	pdf(paste(path,myset,"-DESeq",what,"LinePlots.pdf",sep=""),width=12,height=6)
	for (i in 1:length(interestGenes)) {
		specgenes=names(genes.GR)[which(elementMetadata(genes.GR)$gene%in%interestGenes[[i]])]
		select=which(rownames(exprs(vsdFull))%in%specgenes)
		temp=exprs(vsdFull)[select,]
		if (length(select)>mylimit) {
			mydiff=log((temp[,1]+1)/(temp[,3]+1))
			temp=temp[order(-abs(mydiff))[1:mylimit],]
		}
		myseq=1:length(temp[,3])
		par(mfrow=c(1,2))
		temp=temp[order(temp[,4]),]
		plotCI(myseq,temp[,3],uiw=0,liw=0,col=pcacols[4],pch=19,main=paste(names(interestGenes)[i],"BrownBrite",sep=":"),xlab=names(interestGenes)[i],gap=0,ylab="Deseq counts",xlim=c(0,c(length(myseq)+as.integer(length(myseq)/10))))
		plotCI(myseq,temp[,4],uiw=0,liw=0,col=pcacols[4],add=TRUE,pch=19,gap=0)
		plotCI(myseq,temp[,1],uiw=0,liw=0,col=pcacols[1],add=TRUE,pch=19,gap=0)
		plotCI(myseq,temp[,2],uiw=0,liw=0,col=pcacols[1],gap=0,add=TRUE,pch=19)
		
		tolabel=which(rownames(temp)%in% sign[,1])
		if (length(tolabel)>0) {
			text(myseq[tolabel],temp[tolabel,1],labels=elementMetadata(genes.GR)$gene[match(rownames(temp)[tolabel],names(genes.GR))],cex=0.8,adj=0.1,pos=4,offset=0.2)
		}
		legend("topleft",c("brite","brown","white"),col=pcacols[c(1,4,7)],lwd=2)
		temp=temp[order(temp[,4]),]
		plotCI(myseq,temp[,3],uiw=0,liw=0,col=pcacols[4],pch=19,main=paste(names(interestGenes)[i],"WhiteBrite",sep=":"),xlab=names(interestGenes)[i],gap=0,ylab="Deseq counts",xlim=c(0,c(length(myseq)+as.integer(length(myseq)/10))))
		plotCI(myseq,temp[,4],uiw=0,liw=0,col=pcacols[4],add=TRUE,pch=19,gap=0)
		plotCI(myseq,temp[,5],uiw=0,liw=0,col=pcacols[7],add=TRUE,pch=19,gap=0)
		plotCI(myseq,temp[,6],uiw=0,liw=0,col=pcacols[7],gap=0,add=TRUE,pch=19)
		tolabel=which(rownames(temp)%in%sign2[,1])
		if (length(tolabel)>0) {
			text(myseq[tolabel],temp[tolabel,5],labels=elementMetadata(genes.GR)$gene[match(rownames(temp)[tolabel],names(genes.GR))],cex=0.8,adj=0.1,pos=4,offset=0.2)
		}
		legend("topleft",c("brite","brown","white"),col=pcacols[c(1,4,7)],lwd=2)
		} 
	} 
	dev.off()

}


#### Plot for DeSEQ ####
plotDispEsts <- function( cds )
{
   plot(
      rowMeans( counts( cds, normalized=TRUE ) ),
      fitInfo(cds)$perGeneDispEsts,
      pch = '.', log="xy" )
   xg <- 10^seq( -.5, 5, length.out=300 )
   lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
}
varianceStabilizingTransformation = function (cds)
{
  new("ExpressionSet",
      exprs          = getVarianceStabilizedData(cds),
      phenoData      = phenoData(cds),
      featureData    = featureData(cds),
      experimentData = experimentData(cds),
      annotation     = annotation(cds),
      protocolData   = protocolData(cds))
}
plotMA = function(x, ylim,
  col = ifelse(x$padj>=mycut, "gray32", "red3"),
  linecol = "#ff000080",
  xlab = "mean of normalized counts", mycut=0.05,ylab = expression(log[2]~fold~change),
  log = "x", cex=0.45, ...)
{
  if (!(is.data.frame(x) && all(c("baseMean", "log2FoldChange") %in% colnames(x))))
    stop("'x' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")

  x = subset(x, baseMean!=0)
  py = x$log2FoldChange
  if(missing(ylim))
      ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
  plot(x$baseMean, pmax(ylim[1], pmin(ylim[2], py)),
       log=log, pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
       cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
  abline(h=0, lwd=4, col=linecol)
}



plotPCA2 = function(x, intgroup, ntop=mylimit)
{
  require("lattice")
  require("genefilter")
  rv = rowVars(exprs(x))
  select = order(rv, decreasing=TRUE)[seq_len(ntop)]
  pca = prcomp(t(exprs(x)[select,]))

  fac = factor(apply(pData(vsdFull)[, intgroup], 1, paste, collapse=" : "))
  colours = brewer.pal(nlevels(fac), "Paired")

  pcafig = xyplot(PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=2,
    aspect = "iso", col=colours,
    main = draw.key(key = list(
      rect = list(col = colours),
      text = list(levels(fac)),
      rep = FALSE)))
}

##### function that capitalizes #####
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}

###### function that adds the gene names to a deseq table where the first row is ensembl geneIDs #####
addGeneName <- function(de.DESeq.shZeb1D0,genes.GR) {
	temp=cbind(de.DESeq.shZeb1D0,elementMetadata(genes.GR[match(de.DESeq.shZeb1D0[,1],names(genes.GR))])$gene)
	colnames(temp[dim(temp)[2]])="Name"
	return(temp)
}
######## function correcting a matrix: -inf to min of row, inf to max of row, nan to NA
replaceValuesDeseq=function(x) {
	x[which(x%in%c(-Inf))]=min(x[which(!x%in%c(-Inf,Inf,NaN))])
	x[which(x%in%c(Inf))]=max(x[which(!x%in%c(-Inf,Inf,NaN))])
	x[which(x%in%"NaN")]=NA
	return(x)
}

######## function correcting a matrix: -inf to min of row, inf to max of row, nan to NA
replaceValuesCufflinks=function(x) {
	x[which(x < c(-16.8525))]=-16.8525
	return(x)
}

#### make a smooth scatterplot (with the correlation in the title)
#takes as arguments two vectors & thein names 
smoothPlot <- function(a,b,namea,nameb,mymethod) {
	myind=which((!a%in%c(NA,"NaN","-Inf","Inf"))&(!b%in%c(NA,"NaN","-Inf","Inf")))
	mycor=cor(a[myind],b[myind],method=mymethod)
	mytitle=paste(paste(namea,nameb,collapse="|",sep="|"),prettyNum(mycor,digits=4),sep=":",collapse=":")
	smoothScatter(a[myind],b[myind],main=mytitle,xlab=namea,ylab=nameb)
	return(mycor)
}
#### make a normal scatterplot (with the correlation in the title)
#takes as arguments two vectors & thein names
#plots text on top 
normalPlot <- function(a,b,namea,nameb,mymethod) {
	myind=which((!a%in%c(NA,"NaN","-Inf","Inf"))&(!b%in%c(NA,"NaN","-Inf","Inf")))
	mycor=cor(a[myind],b[myind],method=mymethod)
	mytitle=paste(paste(namea,nameb,collapse="|",sep="|"),prettyNum(mycor,digits=4),sep=":",collapse=":")
	f=abs(max(a[myind])-min(a[myind]))/4
	myxlim=c(c(min(a[myind])-abs(f/4)), max(a[myind]+abs(f/4)))
	f=abs(max(b[myind])-min(b[myind]))/4
	myylim=c(c(min(b[myind])-abs(f/4)), max(b[myind]+abs(f/4)))
	
	plot(a[myind],b[myind],main=mytitle,xlab=namea,ylab=nameb,col="darkblue",pch=20,xlim=c(min(myxlim),max(myxlim)),ylim=c(min(myylim),max(myylim)))
	text(a[myind],b[myind],labels=names(a[myind]),pos=4,offset=0.2)
	return(mycor)
}

library(biomaRt); library(RColorBrewer); library(affycoretools)
library(GenomicRanges); library(MASS)
library(DESeq); library(gplots)
library(topGO); library(cummeRbund)

ensembl=useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl",host="www.ensembl.org")
