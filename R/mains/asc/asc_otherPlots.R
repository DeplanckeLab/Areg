###### Initialise variables and load functionality #####GeneNr
source("plot-util.R")
source("util-asc.R")

##### set paths and libraries 
mmseqpath="src/mmseq-latest/"
myset="bow-mmseqEnsg84-runCut-20161121"
mmset=myset
mydate="161121"
choosecols=c("ivory3","orange", "red")
choosecolsgenes=c("steelblue3","yellow", "red")
#### Load interesting sets of genes ######
load("mygenomes/annotation/GeneTypes84.R")
myname="ASC"
exonLength84=width(exons84.GR)
names(exonLength84)=names(exons84.GR)

path="adipo/ETH/data/other/"
outpath="myadipo/ETH/data/other/plots/"
myfiles=list.files(path,pattern="csv")

#### plot and output
pdf(paste(outpath,"Otherplots-",myname,"-",myset,"-",mydate,".complete.pdf",sep=""),height=10,width=12)
par(mfrow=c(3,3))

##### qPCR FC from validation
myfile=read.csv(paste(path,"Mouse_validation_nrs.csv",sep=""),head=T)
mylogs=myfile[,2]; names(mylogs)=mysplit(myfile[,1],"_",2,2)
barplot(log2(mylogs[grep("P1",myfile[,1])]),main="P1",col="darkgreen",ylab="log2FC",las=2)
barplot(log2(mylogs[grep("P2",myfile[,1])]),main="P2",col="red",ylab="log2FC",las=2)
barplot(log2(mylogs[grep("P3",myfile[,1])]),main="P3.ABCG1",col="blue",ylab="log2FC",las=2)
barplot(log2(mylogs[grep("P4",myfile[,1])]),main="P3.F3",col="blue",ylab="log2FC",las=2)
#### raw Figure
write.table(myfile,paste(path,"Fig2bS3c-raw.PopulationValidationqPCRFC.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


#### qPCR differentiated data
myfile=read.csv(paste(path,"Mouse_differentiation_qPCR.csv",sep=""),head=T,fill=T)
mylogs=myfile[,5]; names(mylogs)=myfile[,2]
barplot(log2(mylogs[grep("CD55",myfile[,1])]),main="P1",col="darkgreen",ylab="log2FC",las=2)
barplot(log2(mylogs[grep("AOC3",myfile[,1])]),main="P2",col="red",ylab="log2FC",las=2)
barplot(log2(mylogs[grep("F3",myfile[,1])])[1:8],main="P3.F3ABCG1",col="blue",ylab="log2FC",las=2)
barplot(log2(mylogs[grep("F3",myfile[,1])])[9:16],main="P3.F3",col="blue",ylab="log2FC",las=2)
#### raw Figure
write.table(myfile,paste(path,"FigS3f-raw.DifferentiationqPCRFC.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


##### differentiation % from differentiation 
myfile=read.csv(paste(path,"Mouse_differentiation_nrs.csv",sep=""),head=T)
beanplot(lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P1"]),
	main="P1",col="darkgreen",ylab="Fr. diff",ylim=c(0.1,0.55))
	#1 vs. 2 0.03961; 1 vs. 3 0.02843; 2 vs. 3  0.02843
	#t-test 0.02142; 0.001867; 0.000286
temp=lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P1"])
pval1=c(t.test(temp[[1]],temp[[2]])$p.val,t.test(temp[[1]],temp[[3]])$p.val,t.test(temp[[2]],temp[[3]])$p.val)
c(wilcox.test(temp[[1]],temp[[2]])$p.val,wilcox.test(temp[[1]],temp[[3]])$p.val,wilcox.test(temp[[2]],temp[[3]])$p.val)
# ttest 0.0214170191 0.0018667186 0.0002862941
# wtest 0.03960870 0.02842954 0.02842954
beanplot(lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P2"]),
	main="P2",col="red",ylab="Fr. diff",ylim=c(0.1,0.55))
	#1 vs. 2  0.4372; 1 vs. 3 0.1532; 2 vs. 3  0.2795
temp=lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P2"])
pval2=c(t.test(temp[[1]],temp[[2]])$p.val,t.test(temp[[1]],temp[[3]])$p.val,t.test(temp[[2]],temp[[3]])$p.val)
#0.5486532 0.4161466 0.6385386
#0.8743671 0.3064922 0.5589857
beanplot(lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P3"]),
	main="P3",col="blue",ylab="Fr. diff",ylim=c(0.1,0.55))
	# 1 vs. 2 0.01141 1 vs. 3 0.01193; 2 vs. 3 0.01167
	# 0.0005379; 7.668e-06; 8.379e-06
temp=lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P3"])
pval3=c(t.test(temp[[1]],temp[[2]])$p.val,t.test(temp[[1]],temp[[3]])$p.val,t.test(temp[[2]],temp[[3]])$p.val)
# ttest 5.379250e-04 7.668050e-06 8.378689e-06
# wtest 0.02857143 0.02940105 0.02940105
beanplot(lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P3.ABCG1"]),
	main="P3.ABCG1",col="blue",ylab="Fr. diff",ylim=c(0,0.3))
	# 1 vs. 2  0.01429; 1 vs. 3 0.02857; 2 vs. 3 0.01429
	# 0.01125; 0.0298; 0.001484
temp=lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P3.ABCG1"])
pval4=c(t.test(temp[[1]],temp[[2]])$p.val,t.test(temp[[1]],temp[[3]])$p.val,t.test(temp[[2]],temp[[3]])$p.val)
# ttest 0.022492663 0.059609729 0.002968105
# wtest 0.02940105 0.07959409 0.02857143
beanplot(lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P3.F3"]),
	main="P3.F3",col="blue",ylab="Fr. diff",ylim=c(0.1,0.55))
	# 1 vs. 2  0.02857; 1 vs. 3 0.0294; 2 vs. 3 0.0294
	#4.807e-07; 1.644e-06; 5.298e-09
temp=lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P3.F3"])
pval5=c(t.test(temp[[1]],temp[[2]])$p.val,t.test(temp[[1]],temp[[3]])$p.val,t.test(temp[[2]],temp[[3]])$p.val)
# 9.613251e-07 3.287659e-06 1.059507e-08
# 0.01141204 0.01192523 0.01166731

mypvals=rbind(pval1,pval2,pval3,pval4,pval5)
rownames(mypvals)=c("P1.CD55","P2.AOC3","P3.ABCG1F3","P3.ABCG1","P3.F3")
colnames(mypvals)=c("P1-P2","P1-P3","P2-P3")
#### raw Figure
write.table(mypvals,paste(path,"Fig2dS3e-raw.DifferentiationPvals.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"Fig2dS3e-raw.Differentiation.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


##### nuclei from differentiation 
myfile=read.csv(paste(path,"Mouse_differentiation_nuclei.csv",sep=""),head=T)
beanplot(lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P1"]),
	main="P1",col="darkgreen",ylab="Fr. diff")
temp=lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P1"])
pval1=c(t.test(temp[[1]],temp[[2]])$p.val,t.test(temp[[1]],temp[[3]])$p.val,t.test(temp[[2]],temp[[3]])$p.val)
beanplot(lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P2"]),
	main="P2",col="red",ylab="Fr. diff")
temp=lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P2"])
pval2=c(t.test(temp[[1]],temp[[2]])$p.val,t.test(temp[[1]],temp[[3]])$p.val,t.test(temp[[2]],temp[[3]])$p.val)
beanplot(lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P3"]),
	main="P3",col="blue",ylab="Fr. diff")
temp=lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P3"])
pval3=c(t.test(temp[[1]],temp[[2]])$p.val,t.test(temp[[1]],temp[[3]])$p.val,t.test(temp[[2]],temp[[3]])$p.val)
beanplot(lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P3.F3"]),
	main="P3.F3",col="blue",ylab="Fr. diff")
temp=lapply(c("pp","p-","p+"),function(x) myfile[which(myfile[,1]%in%x),"P3.F3"])
pval4=c(t.test(temp[[1]],temp[[2]])$p.val,t.test(temp[[1]],temp[[3]])$p.val,t.test(temp[[2]],temp[[3]])$p.val)

mypvals=rbind(pval1,pval2,pval3,pval4)
rownames(mypvals)=c("P1.CD55","P2.AOC3","P3.ABCG1F3","P3.F3")
colnames(mypvals)=c("P1-P2","P1-P3","P2-P3")

#### raw Figure
write.table(mypvals,paste(path,"FigS3g-raw.DifferentiationNucleiPvals.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"FigS3g-raw.DifferentiationNuclei.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


##### obesity nrs.
myfile=read.csv(paste(path,"Mouse_obesity_nrs.csv",sep=""),head=T)
mylogs=myfile[,2]
beanplot(list(myfile[grep("iW",myfile[,1]),2][1:4],myfile[grep("iW",myfile[,1]),2][5:8]),names=c("Lean","Obese"),main="subq",col="blue",ylab="% P3")
# p-val 0.02914 paired t.test
pval1=t.test(myfile[grep("iW",myfile[,1]),2][1:4],myfile[grep("iW",myfile[,1]),2][5:8],paired=T)$p.val
beanplot(list(myfile[grep("vis",myfile[,1]),2][1:4],myfile[grep("vis",myfile[,1]),2][5:8]),names=c("Lean","Obese"),main="visc",col="blue",ylab="% P3")
#p-val 0.002416 paired t.test
pval2=t.test(myfile[grep("vis",myfile[,1]),2][1:4],myfile[grep("vis",myfile[,1]),2][5:8],paired=T)$p.val
beanplot(list(myfile[grep("iW",myfile[,1]),2][1:4],myfile[grep("iW",myfile[,1]),2][5:8]),
	list(myfile[grep("vis",myfile[,1]),2][1:4],myfile[grep("vis",myfile[,1]),2][5:8]),names=c("Subq.Lean","Subq.Obese","Visc.Lean","Visc.Obese")
	,col="blue",ylab="% P3")
pval3=t.test(myfile[grep("iW",myfile[,1]),2][1:4],myfile[grep("vis",myfile[,1]),2][1:4],paired=T)$p.val
pval4=t.test(myfile[grep("iW",myfile[,1]),2][5:8],myfile[grep("vis",myfile[,1]),2][5:8],paired=T)$p.val
# lean subq vs. visc 0.01406 #obese subq vs. visc 0.0001203

mypvals=c(pval1,pval2,pval3,pval4); names(mypvals)=c("iW-leanvsob","vis-leanvsob","lean-iWvsovis","ob-iWvsob")
#### raw Figure
write.table(mypvals,paste(path,"Fig4e-raw.ObesityPvals.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"Fig4e-raw.Obesity.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


##### titration nrs. 
myfile=read.csv(paste(path,"Mouse_titration_nrs.csv",sep=""),head=T)
myvals=lapply(c("F0_","F5_","F10_","F25_","F50_","F80_","F100_"),function(x) myfile[grep(x,myfile[,1]),2])
beanplot(myvals,names=c("0","5","10","25","50","80","100"),ylab="Fr. diff",xlab="% F3")
avg=sapply(myvals,function(x)mean(x))
sdev=sapply(myvals,function(x)sd(x))
plot(c(0,5,10,25,50,80,100), avg,ylim=range(c(avg-sdev, avg+sdev)),
	    pch=19, xlab="% P3", ylab="Mean Fr. Diff +/- SD",main="Titration")
arrows(c(0,5,10,25,50,80,100), avg-sdev, c(0,5,10,25,50,80,100), avg+sdev, length=0.05, angle=90, code=3)
#### and add line
points(c(0,50,100),c(avg[1],c(avg[1]-avg[7])/2+avg[7],avg[7]),col="blue",pch=19)
abline(lsfit(c(0,50,100),c(avg[1],c(avg[1]-avg[7])/2+avg[7],avg[7])),col="blue",lty=2)
plot(c(0,5,10,25,50,80,100), sapply(myvals,function(x) median(x)), xlab="% P3", ylab="median Fr. Diff",main="Titration")
#### raw Figure
write.table(myfile,paste(path,"Fig2f-raw.Titration.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


##### titration nuclei nrs. 
myfile=read.csv(paste(path,"Mouse_titration_nuclei.csv",sep=""),head=T)
myvals=lapply(c("F0_","F5_","F10_","F25_","F50_","F80_","F100_"),function(x) myfile[grep(x,myfile[,1]),2])
beanplot(myvals,names=c("0","5","10","25","50","80","100"),ylab="Nuclei nr.",xlab="% F3")
avg=sapply(myvals,function(x)mean(x))
sdev=sapply(myvals,function(x)sd(x))
plot(c(0,5,10,25,50,80,100), avg,ylim=range(c(avg-sdev, avg+sdev)),
	    pch=19, xlab="% P3", ylab="Mean Nuclei +/- SD",main="Titration")
arrows(c(0,5,10,25,50,80,100), avg-sdev, c(0,5,10,25,50,80,100), avg+sdev, length=0.05, angle=90, code=3)
plot(c(0,5,10,25,50,80,100), sapply(myvals,function(x)median(x)), xlab="% P3", ylab="median Fr. Diff",main="Titration")
#### raw Figure
write.table(myfile,paste(path,"FigS4a-raw.TitrationNuclei.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

##### in vivo nrs. 
myfile=read.csv(paste(path,"Mouse_injection_nrs.csv",sep=""),head=T)
beanplot(list(apply(myfile[1:3,2:5],1,function(x) mean(x)),apply(myfile[4:6,2:5],1,function(x) mean(x))),
	col="blue",ylab="Fr. diff",names=c("+","-"),main="In vivo")
# p-val 0.05
#### raw Figure
write.table(myfile,paste(path,"Fig4hS7j-raw.InVivo.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

#### Enrichr plot
P1=c(8,2,0,0,0,0); P2=c(0,8,1,1,0,0); P3=c(0,1,0,0,7,2)
names(P1)=names(P2)=names(P3)=c("skin","adipose_tissue","stomach","breast","blood_vessel","nerve")
barplot(P1,col="darkgreen",las=2,main="P1");barplot(P2,col="red",las=2,main="P2")
barplot(P3,col="blue",las=2,main="P3")
myfile=cbind(P1,P2,P3);rownames(myfile)=c("skin","adipose_tissue","stomach","breast","blood_vessel","nerve")
colnames(myfile)=c("P1","P2","P3")
#### raw Figure
write.table(myfile,paste(path,"FigS1g-raw.Enrichr.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

##### human in vivo nrs. 
myfile=read.csv(paste(path,"human_i07_exvivo.csv",sep=""),head=T)
beanplot(lapply(2:4,function(x) as.numeric(myfile[,x])),
	col="yellow",ylab="Fr. diff",names=c("DN","QN","DNDP"),main="Ex vivo i07")
# wilcox 1-2 0.07085 1-3 4.826e-06 2-3 2.757e-10 i07 p-val 0.05
# ttest 1-2 0.03055 1-3  5.256e-06 2-3 4.287e-09
temp=lapply(2:4,function(x) as.numeric(myfile[,x]))
mypval=c(wilcox.test(temp[[1]],temp[[2]])$p.val,wilcox.test(temp[[1]],temp[[3]])$p.val,wilcox.test(temp[[2]],temp[[3]])$p.val)
names(mypval)=c("DN-QN","DN-DNDP","QN-DNDP")
#### raw Figure
write.table(mypval,paste(path,"Fig3b-raw.humanDiffi07Pvals.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"Fig3b-raw.humanDiffi07.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### review data #####
#### transwell 
## exp1
myfile=read.csv(paste(path,"Mouse_transwell_KD.csv",sep=""),head=T)[1:24,]
beanplot(lapply(c("ctrl_siRNA","fgf12_siRNA","rtp3_siRNA","spink2_siRNA"),
	function(x) myfile[which(myfile[,1]%in%x),"Diff"]),
	names=c("Ctrl","Fgf12","Rtp3","Spink2"),
	main="Transwell KD",col="yellow",ylab="Fr. diff",ylim=c(0.05,0.35))
temp=lapply(c("ctrl_siRNA","fgf12_siRNA","rtp3_siRNA","spink2_siRNA"),
		function(x) myfile[which(myfile[,1]%in%x),"Diff"])
myp1=sapply(2:4,function(x) t.test(temp[[1]],temp[[x]])$p.val)
# ttest 0.0056596027 0.0002438112 0.0063106924
# wtest 0.004329004 0.004329004 0.008658009
##### nuclei from differentiation 
beanplot(lapply(c("ctrl_siRNA","fgf12_siRNA","rtp3_siRNA","spink2_siRNA"),
	function(x) myfile[which(myfile[,1]%in%x),"Nuclei"]),
	names=c("Ctrl","Fgf12","Rtp3","Spink2"),
	main="Transwell KD",col="yellow",ylab="Fr. diff",ylim=c(4000,6500))
temp=lapply(c("ctrl_siRNA","fgf12_siRNA","rtp3_siRNA","spink2_siRNA"),
		function(x) myfile[which(myfile[,1]%in%x),"Nuclei"])
myp1.nuc=sapply(2:4,function(x) t.test(temp[[1]],temp[[x]])$p.val)
names(myp1)=names(myp1.nuc)=c("Fgf12","Rtp3","Spink2")
myp=cbind(myp1,myp1.nuc)
#### raw Figure
write.table(myp,paste(path,"FigS5hi-raw.transwellKDR2pval.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"FigS5hi-raw.transwellKDR2l.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

## exp2
myfile=read.csv(paste(path,"Mouse_transwell_KD.csv",sep=""),head=T)[25:54,]
beanplot(lapply(c("ctrl_siRNA","fgf12_siRNA","rtp3_siRNA","spink2_siRNA","vit-siRNA"),
	function(x) myfile[which(myfile[,1]%in%x),"Diff"]),
	names=c("Ctrl","Fgf12","Rtp3","Spink2","Vit"),
	main="Transwell KD",col="yellow",ylab="Fr. diff",ylim=c(0.05,0.3))
temp=lapply(c("ctrl_siRNA","fgf12_siRNA","rtp3_siRNA","spink2_siRNA","vit-siRNA"),
			function(x) myfile[which(myfile[,1]%in%x),"Diff"])
myp2=sapply(2:5,function(x) t.test(temp[[1]],temp[[x]])$p.val)
# ttest 1.906324e-01 3.083478e-05 5.573071e-03
# wtest 0.240259740 0.002164502 0.008658009
##### nuclei from differentiation 
beanplot(lapply(c("ctrl_siRNA","fgf12_siRNA","rtp3_siRNA","spink2_siRNA","vit-siRNA"),
	function(x) myfile[which(myfile[,1]%in%x),"Nuclei"]),
	names=c("Ctrl","Fgf12","Rtp3","Spink2","Vit"),
	main="Transwell KD",col="yellow",ylab="Fr. diff",ylim=c(3500,5000))
temp=lapply(c("ctrl_siRNA","fgf12_siRNA","rtp3_siRNA","spink2_siRNA","vit-siRNA"),
				function(x) myfile[which(myfile[,1]%in%x),"Nuclei"])
myp2.nuc=sapply(2:5,function(x) t.test(temp[[1]],temp[[x]])$p.val)
names(myp2)=names(myp2.nuc)=c("Fgf12","Rtp3","Spink2","Vit")
myp=cbind(myp2,myp2.nuc)
#### raw Figure
write.table(myp,paste(path,"Fig2iS5ef-raw.transwellKDR1pval.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"Fig2iS5ef-raw.transwellKDR1.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

##### qPCR transwell differentiation 
myfile=read.csv(paste(path,"Mouse_transwell_KD_diffqPCR.csv",sep=""),head=T)
myp.wilc=myp.ttest=matrix(data=1,ncol=4,nrow=8);k=1
for (mygene in c("Adipoq", "Fabp4", "Cd36" ,"Fasn", "Cebpa", "Cebpb" ,"Leptn","Pparg2")) {
	toplot=cbind(	sapply(seq(1,12,by=2), function(x) mean(myfile[x:(x+1),mygene])),
			sapply(seq(13,24,by=2), function(x) mean(myfile[x:(x+1),mygene])),
			sapply(seq(25,36,by=2), function(x) mean(myfile[x:(x+1),mygene])),
			sapply(seq(37,48,by=2), function(x) mean(myfile[x:(x+1),mygene])),
			sapply(seq(49,60,by=2), function(x) mean(myfile[x:(x+1),mygene])))
	colnames(toplot)=c("Cntrl","Fgf12","Rtp3","Spink2","Vit"); rownames(toplot)=c(1:6)
	barx=barplot(apply(toplot,2,function(x) mean(x)),main=mygene); error.bar(barx,apply(toplot,2,function(x) mean(x)),apply(toplot,2,function(x) sd(x)))
	points(rep(barx,rep(6,5)),c(toplot), pch=19, col="black")	
	myp.wilc[k,]=p.adjust(sapply(2:5,function(x) wilcox.test(toplot[,1],toplot[,x])$p.val),n=32,method="BH")
	myp.ttest[k,]=p.adjust(sapply(2:5,function(x) t.test(toplot[,1],toplot[,x])$p.val),n=32,method="BH")
	k=k+1
}
colnames(myp.wilc)=colnames(myp.ttest)=c("Fgf12","Rtp3","Spink2","Vit")
rownames(myp.wilc)=rownames(myp.ttest)=c("Adipoq", "Fabp4", "Cd36" ,"Fasn", "Cebpa", "Cebpb" ,"Leptn","Pparg2")
write.table(myp.wilc,paste(path,"Mouse_transwell_KD_diffqPCR.pwilcox.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
write.table(myp.ttest,paste(path,"Mouse_transwell_KD_diffqPCR.ttest.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
#### raw Figure
write.table(myp.ttest,paste(path,"FigS5j-raw.transwellqPCR.ttest.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myp.wilc,paste(path,"FigS5j-raw.transwellqPCR.wilcoxtest.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"FigS5j-raw.transwellqPCR.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


##### qPCR transwell genes 
myfile=read.csv(paste(path,"Mouse_transwell_KD_qPCR.csv",sep=""),head=T)
k=1
myfile.red=myfile
par(mfrow=c(4,6))
myp=NULL;myp.wilc=NULL
while(dim(myfile.red)[1]>1) {
mylen=grep("kd",myfile.red[,1])[1]-1
control.sub=myfile.red[k:c(k+mylen-1),]
factor.sub=myfile.red[c(k+mylen):(k+2*mylen-1),]
toplot=cbind(	sapply(seq(1,mylen,by=2), function(x) mean(control.sub[x:(x+1),3])),
		sapply(seq(1,mylen,by=2), function(x) mean(factor.sub[x:(x+1),3])))
colnames(toplot)=c("C",mysplit(factor.sub[,1][1],"-kd",2,1))
barx=barplot(apply(toplot,2,function(x) mean(x))); error.bar(barx,apply(toplot,2,function(x) mean(x)),apply(toplot,2,function(x) sd(x)))
points(rep(barx,rep(dim(toplot)[1],dim(toplot)[2])),c(toplot), pch=19, col="black")
myfile.red=myfile.red[-(k:c(k+mylen*2-1)),];myp=c(myp,t.test(toplot[,1],toplot[,2])$p.val)
myp.wilc=c(myp.wilc,wilcox.test(toplot[,1],toplot[,2])$p.val)}
names(myp)=names(myp.wilc)=c("Fgf12","Rtp3","Spink2","Vit")
myp=rbind(myp,myp.wilc);rownames(myp)=c("t.test","wilcox.test")
#### raw Figure
write.table(myp,paste(path,"FigS5g-raw.transwellqPCRPval.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"FigS5g-raw.transwellqPCR.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


##### 
myfile=read.csv(paste(path,"Mouse_differentiation_KD.csv",sep=""),head=T)
beanplot(lapply(c("Cntrl","Fgf12","Rtp3","Spink2","Vit"),
	function(x) myfile[1:4,x]),names=c("Cntrl","Fgf12","Rtp3","Spink2","Vit"),
	main="KD.1",col="yellow",ylab="Fr. diff",ylim=c(0.05,0.25),las=2)
beanplot(lapply(c("Cntrl","Fgf12","Rtp3","Spink2","Vit"),
	function(x) myfile[5:7,x]),names=c("Cntrl","Fgf12","Rtp3","Spink2","Vit"),
	main="KD.2",col="yellow",ylab="Fr. diff",ylim=c(0.05,0.25),las=2)
beanplot(lapply(c("Cntrl","Fgf12","Rtp3","Spink2"),
	function(x) myfile[8:11,x]),names=c("Cntrl","Fgf12","Rtp3","Spink2"),
	main="KD.1 Exp2",col="yellow",ylab="Fr. diff",ylim=c(0.05,0.25),las=2)
beanplot(lapply(c("Cntrl","Fgf12","Rtp3","Spink2"),
	function(x) myfile[12:15,x]),names=c("Cntrl","Fgf12","Rtp3","Spink2"),
	main="KD.2 Exp2",col="yellow",ylab="Fr. diff",ylim=c(0.05,0.25),las=2)
### and merged
beanplot(lapply(c("Cntrl","Fgf12","Rtp3","Spink2","Vit"),
	function(x) myfile[1:7,x]),names=c("Cntrl","Fgf12","Rtp3","Spink2","Vit"),
	main="KD",col="yellow",ylab="Fr. diff",ylim=c(0.05,0.25),las=2)
	#t.test 4.400382e-05 4.400382e-05 4.400382e-05 4.400382e-05
	#wilcox.test 0.0005827506 0.0005827506 0.0005827506 0.0005827506
temp=lapply(c("Cntrl","Fgf12","Rtp3","Spink2","Vit"),function(x) myfile[1:7,x])
myp.ttest=sapply(2:5,function(x) t.test(temp[[1]],temp[[x]])$p.val)
myp.wilcox=sapply(2:5,function(x) wilcox.test(temp[[1]],temp[[x]])$p.val)
beanplot(lapply(c("Cntrl","Fgf12","Rtp3","Spink2"),
	function(x) myfile[8:15,x]),names=c("Cntrl","Fgf12","Rtp3","Spink2"),
	main="KD Exp2",col="yellow",ylab="Fr. diff",ylim=c(0.05,0.25),las=2)
temp=lapply(c("Cntrl","Fgf12","Rtp3","Spink2"),function(x) myfile[8:15,x])
myp.ttest2=sapply(2:4,function(x) t.test(temp[[1]],temp[[x]])$p.val)
myp.wilcox2=sapply(2:4,function(x) wilcox.test(temp[[1]],temp[[x]])$p.val)
#t.test 3.986509e-05 3.986509e-05 3.986509e-05
#wilcox.test 0.0001554002 0.0001554002 0.0001554002
#### and vit
beanplot(lapply(c("Cntrl","Vit"),
	function(x) myfile[16:21,x]),names=c("Cntrl","Vit"),
	main="KD.Vit",col="yellow",ylab="Fr. diff",ylim=c(0.1,0.25),las=2)
temp=lapply(c("Cntrl","Vit"),function(x) myfile[16:21,x])
myp.ttest3=sapply(2,function(x) t.test(temp[[1]],temp[[x]])$p.val)
myp.wilcox3=sapply(2,function(x) wilcox.test(temp[[1]],temp[[x]])$p.val)

# ttest 0.0004884
# wilcox.test 0.002165

myp=rbind(myp.ttest,c(myp.ttest2,myp.ttest3),myp.wilcox,c(myp.wilcox2,myp.wilcox3))
colnames(myp)=c("Fgf12","Rtp3","Spink2","Vit")
rownames(myp)=c("t.testSupp","t.testMain","wilcox.Supp","wilcox.Main")
#### raw Figure
write.table(myp,paste(path,"FigS5e-raw.KDDifferentiationPval.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"FigS5e--raw.KDDifferentiation.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


##### Nuclei for the differentiation
myfile=read.csv(paste(path,"Mouse_differentiation_KD_nuclei.csv",sep=""),head=T)
beanplot(lapply(c("Cntrl","Fgf12","Rtp3","Spink2","Vit"),
	function(x) myfile[1:4,x]),names=c("Cntrl","Fgf12","Rtp3","Spink2","Vit"),
	main="KD.1",col="yellow",ylab="Nr. nuclei",ylim=c(2800,9000),las=2)
beanplot(lapply(c("Cntrl","Fgf12","Rtp3","Spink2","Vit"),
	function(x) myfile[5:7,x]),names=c("Cntrl","Fgf12","Rtp3","Spink2","Vit"),
	main="KD.2",col="yellow",ylab="Nr. nuclei",ylim=c(2800,9000),las=2)
beanplot(lapply(c("Cntrl","Fgf12","Rtp3","Spink2"),
	function(x) myfile[8:11,x]),names=c("Cntrl","Fgf12","Rtp3","Spink2"),
	main="KD.1 Exp2",col="yellow",ylab="Nr. nuclei",ylim=c(2800,9000),las=2)
beanplot(lapply(c("Cntrl","Fgf12","Rtp3","Spink2"),
	function(x) myfile[12:15,x]),names=c("Cntrl","Fgf12","Rtp3","Spink2"),
	main="KD.2 Exp2",col="yellow",ylab="Nr. nuclei",ylim=c(2800,9000),las=2)
#### together
beanplot(lapply(c("Cntrl","Fgf12","Rtp3","Spink2","Vit"),
	function(x) myfile[1:7,x]),names=c("Cntrl","Fgf12","Rtp3","Spink2","Vit"),
	main="KD",col="yellow",ylab="Nr. nuclei",ylim=c(2800,9000),las=2)
temp=lapply(c("Cntrl","Fgf12","Rtp3","Spink2","Vit"),function(x) myfile[1:7,x])
myp.ttest=sapply(2:5,function(x) t.test(temp[[1]],temp[[x]])$p.val)
myp.wilcox=sapply(2:5,function(x) wilcox.test(temp[[1]],temp[[x]])$p.val)
beanplot(lapply(c("Cntrl","Fgf12","Rtp3","Spink2"),
	function(x) myfile[8:15,x]),names=c("Cntrl","Fgf12","Rtp3","Spink2"),
	main="KD.1 Exp2",col="yellow",ylab="Nr. nuclei",ylim=c(2800,9000),las=2)
temp=lapply(c("Cntrl","Fgf12","Rtp3","Spink2"),function(x) myfile[8:15,x])
myp.ttest2=sapply(2:4,function(x) t.test(temp[[1]],temp[[x]])$p.val)
myp.wilcox2=sapply(2:4,function(x) wilcox.test(temp[[1]],temp[[x]])$p.val)

#### and vit
beanplot(lapply(c("Cntrl","Vit"),
	function(x) myfile[16:21,x]),names=c("Cntrl","Vit"),
	main="KD.Vit",col="yellow",ylab="Nr. nuclei",ylim=c(4000,6000),las=2)
temp=lapply(c("Cntrl","Vit"),function(x) myfile[16:21,x])
myp.ttest3=sapply(2,function(x) t.test(temp[[1]],temp[[x]])$p.val)
myp.wilcox3=sapply(2,function(x) wilcox.test(temp[[1]],temp[[x]])$p.val)

myp=rbind(myp.ttest,c(myp.ttest2,myp.ttest3),myp.wilcox,c(myp.wilcox2,myp.wilcox3))
colnames(myp)=c("Fgf12","Rtp3","Spink2","Vit")
rownames(myp)=c("t.testSupp","t.testMain","wilcox.Supp","wilcox.Main")

#### raw Figure
write.table(myp,paste(path,"FigS5f-raw.KDDifferentiationNucPval.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"FigS5f-raw.KDDifferentiationNuc.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)



####### Differentiation for all factors
myfile=read.csv(paste(path,"Mouse_differentiation_KD_all.csv",sep=""),head=T)
allexp=as.vector(unique(myfile[,1]))
for(exp in allexp) {
	myfile.sub=myfile[which(myfile[,1]%in%exp),]
	allsamples=unique(myfile.sub[,2])
	beanplot(lapply(allsamples,function(x) as.numeric(myfile.sub[which(myfile.sub[,2]%in%x),8])),
	names=as.vector(sapply(allsamples,function(x) unique(myfile.sub[which(myfile.sub[,2]%in%x),3]))),
	main=unique(sapply(allsamples,function(x) unique(myfile.sub[which(myfile.sub[,2]%in%x),5]))),
	col="yellow",ylab="Fr. diff",las=2)
	beanplot(lapply(allsamples,function(x) as.numeric(myfile.sub[which(myfile.sub[,2]%in%x),7])),
	names=as.vector(sapply(allsamples,function(x) unique(myfile.sub[which(myfile.sub[,2]%in%x),3]))),
	main=unique(sapply(allsamples,function(x) unique(myfile.sub[which(myfile.sub[,2]%in%x),5]))),
	col="yellow",ylab="Nr. nuclei",las=2)
}

##### qPCR transwell differentiation 
myfile=read.csv(paste(path,"Mouse_KD_qPCR.csv",sep=""),head=T)
k=1
myfile.red=myfile
par(mfrow=c(4,6))
myp=NULL;myp.wilc=NULL
while(dim(myfile.red)[1]>1) {
mylen=grep("C-",myfile.red[,1])[1]-1
factor.sub=myfile.red[k:c(k+mylen-1),]
control.sub=myfile.red[c(k+mylen):(k+2*mylen-1),]
toplot=cbind(	sapply(seq(1,mylen,by=2), function(x) mean(control.sub[x:(x+1),2])),
		sapply(seq(1,mylen,by=2), function(x) mean(factor.sub[x:(x+1),2])))
colnames(toplot)=c("C",mysplit(factor.sub[,1][1],"-kd",2,1))
barx=barplot(apply(toplot,2,function(x) mean(x))); error.bar(barx,apply(toplot,2,function(x) mean(x)),apply(toplot,2,function(x) sd(x)))
points(rep(barx,rep(dim(toplot)[1],dim(toplot)[2])),c(toplot), pch=19, col="black")	
myfile.red=myfile.red[-(k:c(k+mylen*2-1)),]
myp=c(myp,t.test(toplot[,1],toplot[,2])$p.val)
myp.wilc=c(myp.wilc,wilcox.test(toplot[,1],toplot[,2])$p.val)
}
names(myp)=names(myp.wilc)=c("Wisp1","Slamf9","Tspan11","Klf5","Wfdc1","Angpt4","Fmod","Vit","Mamdc2","Ndufa4l2","Dner",
"Ptch2","Hhip","Hopx","Rtp3","Il20rb","Slc6a2","Spink2","Fgf12","Dapk2","Kcnj8","Cpxm2","Akrlc12","Vit")

myp=rbind(myp,myp.wilc,p.adjust(myp,method="BH"),p.adjust(myp.wilc,method="BH"))
rownames(myp)=c("t.test","wilcox.test","adj.ttest","adj.wilcox.test")
#### raw Figure
write.table(myp,paste(path,"FigS5c-raw.AregKDqPCRpval.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"FigS5c-raw.AregKDqPCR.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


####### Differentiation from human
##### differentiation % from differentiation 
par(mfrow=c(3,3))
myfile=read.csv(paste(path,"Hsa_exvivo.csv",sep=""),head=T)
myind=unique(myfile[,1])
myp.wcox=myp.ttest=matrix(data=NA,ncol=3,nrow=length(myind)); k=1
for(i in myind) {
	temp=lapply(2:4,function(x) myfile[which(myfile[,1]%in%i),x])
	beanplot(temp,main=i,col="yellow",ylab="Adipogenic diff. [a.u.]")
	myp.wcox[k,]=p.adjust(c(wilcox.test(temp[[1]],temp[[2]])$p.val,wilcox.test(temp[[1]],temp[[3]])$p.val,
		+ wilcox.test(temp[[2]],temp[[3]])$p.val),n=18,method="BH")
	myp.ttest[k,]=p.adjust(c(t.test(temp[[1]],temp[[2]])$p.val,t.test(temp[[1]],temp[[3]])$p.val,t.test(temp[[2]],temp[[3]])$p.val),n=18,method="BH");k=k+1

#	p.adjust(c(wilcox.test(temp[[1]],temp[[2]])$p.val,wilcox.test(temp[[1]],temp[[3]])$p.val,
#	wilcox.test(temp[[2]],temp[[3]])$p.val),n=18,method="BH")
	
}
colnames(myp.wcox)=colnames(myp.ttest)=c("1-2","1-3","2-3"); 
rownames(myp.wcox)=rownames(myp.ttest)=myind
write.table(myp.wcox,paste(path,"Human_diff.pwilcox.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
write.table(myp.ttest,paste(path,"Human_diff.ttest.txt",sep=""),quote=F,row.names=F,col.names=F,sep="\t")

#### raw Figure
write.table(myp.ttest,paste(path,"Fig3bd-raw.humanDiffPval.ttest.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myp.wcox,paste(path,"Fig3bd-raw.humanDiffPval.wilcox.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"Fig3bd-raw.humanDiff.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


###### Human qPCR nrs. 
myfile=read.csv(paste(path,"Hsa_exvivo_qPCR.csv",sep=""),head=T)
myind=unique(myfile[,1])
undiffvals=sapply(myind[1:3],function(i) log2(myfile[which(myfile[,1]%in%i),4]/myfile[which(myfile[,1]%in%i),3]))

par(mfrow=c(3,3))
rownames(undiffvals)=myfile[1:5,2];colnames(undiffvals)=myind[1:3]
diffvals=sapply(myind[4:6],function(i) log2(myfile[which(myfile[,1]%in%i),4]/myfile[which(myfile[,1]%in%i),3]))
rownames(diffvals)=myfile[16:20,2];colnames(diffvals)=myind[1:3]
for (i in 1:3) {
	barplot(undiffvals[,i],main=colnames(undiffvals)[i],names=rownames(undiffvals),ylab=c("log2 +/-"))
	barplot(diffvals[,i],main=colnames(diffvals)[i],names=rownames(diffvals),ylab=c("log2 +/-"))
}	
barx=barplot(apply(undiffvals,1,mean),main="Undiff",names=rownames(undiffvals),ylab=c("log2 +/-"))
error.bar(barx,apply(undiffvals,1,mean),apply(undiffvals,1,sd))
points(rep(barx,rep(dim(undiffvals)[2],dim(undiffvals)[1])),c(t(undiffvals)), pch=19, col="black")
barx=barplot(apply(diffvals,1,mean),main="Diff",names=rownames(diffvals),ylab=c("log2 +/-"))
error.bar(barx,apply(diffvals,1,mean),apply(diffvals,1,sd))
points(rep(barx,rep(dim(diffvals)[2],dim(diffvals)[1])),c(t(diffvals)), pch=19, col="black")

#### raw Figure
colnames(myfile)[3:4]=c("CD142-","CD142+")
write.table(myfile,paste(path,"FigS6b-raw.humanqPCR.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

undiff.p=sapply(c("F3","Mgp","Meox2","Abcg1"),
	function(x) t.test(myfile[grep(x,myfile[,2])[1:3],3],
	myfile[grep(x,myfile[,2])[1:3],4],paired=T,alternative="less")$p.val )
diff.p=c(sapply(c("F3","Mgp"),
	function(x) t.test(myfile[grep(x,myfile[,2])[4:6],3],
	myfile[grep(x,myfile[,2])[4:6],4],paired=T,alternative="greater")$p.val ),
	sapply(c("Pparg","Fabp4"),function(x) t.test(myfile[grep(x,myfile[,2])[1:3],3],
		myfile[grep(x,myfile[,2])[1:3],4],paired=T,alternative="greater")$p.val ))
myp=c(undiff.p,diff.p);names(myp)=c("undiff.F3","undiff.Mgp","undiff.Meox2","undiff.Abcg1","diff.F3","diff.Mgp","diff.Pparg","diff.Fabp4")
write.table(myp,paste(path,"FigS6b-raw.humanqPCRPval.ttest.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

#### Wilcox
undiff.p=sapply(c("F3","Mgp","Meox2","Abcg1"),
	function(x) wilcox.test(myfile[grep(x,myfile[,2])[1:3],3],
	myfile[grep(x,myfile[,2])[1:3],4],paired=T,alternative="less")$p.val )
diff.p=c(sapply(c("F3","Mgp"),
	function(x) wilcox.test(myfile[grep(x,myfile[,2])[4:6],3],
	myfile[grep(x,myfile[,2])[4:6],4],paired=T,alternative="greater")$p.val ),
	sapply(c("Pparg","Fabp4"),function(x) wilcox.test(myfile[grep(x,myfile[,2])[1:3],3],
		myfile[grep(x,myfile[,2])[1:3],4],paired=T,alternative="greater")$p.val ))
myp=c(undiff.p,diff.p);names(myp)=c("undiff.F3","undiff.Mgp","undiff.Meox2","undiff.Abcg1","diff.F3","diff.Mgp","diff.Pparg","diff.Fabp4")
write.table(myp,paste(path,"FigS6b-raw.humanqPCRPval.wilcoxtest.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

#### Mouse human transwells 
myfile=read.csv(paste(path,"Mmu_hsa_transwell.csv",sep=""),head=T)
beanplot(myfile[,2:3],main="Mouse-human transwell",col="yellow",ylab="Adipogenic diff. [a.u.]")
myp=c(wilcox.test(myfile[,2],myfile[,3])$p.val,t.test(myfile[,2],myfile[,3])$p.val)
names(myp)=c("Wilcox","ttest")
write.table(myp,paste(path,"FigS6k-raw.mousehuman.pvals.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"FigS6k-raw.mousehuman.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


#### Mouse simple transwells 
myfile=read.csv(paste(path,"Transwell_Diff_simple.csv",sep=""),head=T)
beanplot(myfile[grep("S1",myfile[,1]),2:4],main="Mouse transwell S1",col="yellow",ylab="Adipogenic diff. [a.u.]")
beanplot(myfile[grep("S2",myfile[,1]),2:4],main="Mouse transwell S2",col="yellow",ylab="Adipogenic diff. [a.u.]")

myp=rbind(c(wilcox.test(myfile[grep("S1",myfile[,1]),2],myfile[grep("S1",myfile[,1]),3])$p.val,
	wilcox.test(myfile[grep("S1",myfile[,1]),2],myfile[grep("S1",myfile[,1]),4])$p.val,
	wilcox.test(myfile[grep("S1",myfile[,1]),3],myfile[grep("S1",myfile[,1]),4])$p.val),
	c(wilcox.test(myfile[grep("S2",myfile[,1]),2],myfile[grep("S2",myfile[,1]),3])$p.val,
		wilcox.test(myfile[grep("S2",myfile[,1]),2],myfile[grep("S2",myfile[,1]),4])$p.val,
		wilcox.test(myfile[grep("S2",myfile[,1]),3],myfile[grep("S2",myfile[,1]),4])$p.val),
	c(t.test(myfile[grep("S1",myfile[,1]),2],myfile[grep("S1",myfile[,1]),3])$p.val,
		t.test(myfile[grep("S1",myfile[,1]),2],myfile[grep("S1",myfile[,1]),4])$p.val,
		t.test(myfile[grep("S1",myfile[,1]),3],myfile[grep("S1",myfile[,1]),4])$p.val),
		c(t.test(myfile[grep("S2",myfile[,1]),2],myfile[grep("S2",myfile[,1]),3])$p.val,
			t.test(myfile[grep("S2",myfile[,1]),2],myfile[grep("S2",myfile[,1]),4])$p.val,
			t.test(myfile[grep("S2",myfile[,1]),3],myfile[grep("S2",myfile[,1]),4])$p.val))
colnames(myp)=c("ASPCs vs. -","ASPCs. vs. +", "- vs. +")
rownames(myp)=c("S1.wilcox","S2.wilcox","S1.ttest","S2.ttest")
write.table(myp,paste(path,"FigS5b-raw.transwellSimple.pvals.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"FigS5b-raw.transwellSimple.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)

#### In vivo mouse plugs
myfile=read.csv(paste(path,"Mouse_invivo.csv",sep=""),head=T)
beanplot(list(myfile[1:5,2],myfile[8:12,2]),names=c("C","CD142-"),col="yellow",main="SVF")
beanplot(list(myfile[1:7,3],myfile[8:14,3]),names=c("C","CD142-"),col="yellow",main="SVF Lin- SCA1+")

temp1=list(myfile[1:5,2],myfile[8:12,2]) # 0.01593
temp2=list(myfile[1:7,3],myfile[8:14,3]) # 0.008441

t.test([temp1[[1]],temp1[[2]],paired=T)
myp=c(t.test(temp1[[1]],temp1[[2]],paired=T)$p.val,t.test(temp2[[1]],temp2[[2]],paired=T)$p.val)
names(myp)=c("SVF","Lin-SVF")
#### raw Figure
write.table(myp,paste(path,"Fig4hS7i-raw.invivoPlugsPval.ttest.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"Fig4hS7i-raw.invivoPlugs.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


#####
myfile=read.csv(paste(path,"Vessels.csv",sep=""),head=T)
beanplot(list(myfile[,1],myfile[,2]),names=c("C","CD142-"),col="yellow",main="IB4")
temp1=list(myfile[,1],myfile[,2]) #0.9815

myp=t.test(myfile[,1],myfile[,2],paired=T)$p.val
#### raw Figure
write.table(myp,paste(path,"FigS7l-raw.VesselsPval.ttest.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(myfile,paste(path,"FigS7l-raw.Vessels.txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)


dev.off()
