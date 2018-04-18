#### make a normal scatterplot (with the correlation in the title)
#takes as arguments two vectors & thein names
#plots text on top 
scatterNormalPlot <- function(a,b,namea,nameb,mymethod) {
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

#### make a smooth scatterplot (with the correlation in the title)
#takes as arguments two vectors & their names 
scatterSmoothPlot <- function(a,b,namea,nameb,mymethod) {
	myind=which((!a%in%c(NA,"NaN","-Inf","Inf"))&(!b%in%c(NA,"NaN","-Inf","Inf")))
	mycor=cor(a[myind],b[myind],method=mymethod)
	mytitle=paste(paste(namea,nameb,collapse="|",sep="|"),prettyNum(mycor,digits=4),sep=":",collapse=":")
	smoothScatter(a[myind],b[myind],main=mytitle,xlab=namea,ylab=nameb)
	return(mycor)
}

###### Make the heatmap needed in expranalysis #######
makeprefheatmaprlog = function(genetable,mygenes.GR,colvals="grey",mainname="row.scaled",
	rowvals="grey") {
	mymatrix=genetable
	if (length(grep("ENS",rownames(mymatrix)))>0) {
		rownames(mymatrix)[grep("ENS",rownames(mymatrix))]=elementMetadata(mygenes.GR[rownames(mymatrix)[grep("ENS",rownames(mymatrix))]])$geneName
	}
	rowvals=valuesToColors(apply(mymatrix,1,function(x) median(x)),c("ivory","darkred"),names(apply(mymatrix,1,function(x) median(x))))
	mymatrix=mymatrix[which(!duplicated(rownames(mymatrix))),]
	
	heatmap.2(mymatrix,scale="none",col=rev(colorRampPalette(c("red","orange","ivory","lightblue","darkblue"))(216)),
		main=mainname,density.info="none",notecex=0.8,notecol="black",keysize = 1,trace="none",cexRow=0.6,cexCol=0.4,
		ColSideColors=colvals[colnames(mymatrix)],RowSideColors=rowvals[rownames(mymatrix)])
	heatmap.2(mymatrix,scale="row",col=rev(colorRampPalette(c("red","ivory","lightblue"))(128)),
		main=mainname,density.info="none",notecex=0.8,notecol="black",keysize = 1,trace="none",cexRow=0.6,cexCol=0.4,
		ColSideColors=colvals[colnames(mymatrix)],RowSideColors=rowvals[rownames(mymatrix)])
	#### prettier heatmap #####
	pheatmap(mymatrix,scale="none",clustering_distance_rows="correlation",cluster_cols = F,main=mainname,
		clustering_method="average",fontsize_col=4,fontsize_row=7,show_rownames=T,cexRow=0.1,
		cexCol=0.1,ColSideColors=colvals[colnames(mymatrix)],RowSideColors=rowvals[rownames(mymatrix)])
	pheatmap(mymatrix,clustering_distance_rows="correlation",cluster_cols = F,main=mainname,
		clustering_method="average",scale = "row",fontsize_col=4,fontsize_row=7,show_rownames=T,cexRow=0.1,
		cexCol=0.1,ColSideColors=colvals[colnames(mymatrix)],RowSideColors=rowvals[rownames(mymatrix)])
	
}
