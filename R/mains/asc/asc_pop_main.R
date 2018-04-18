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
myname="Pop"
exonLength84=width(exons84.GR)
names(exonLength84)=names(exons84.GR)


##### READ MMSEQ FILES AS REQUIRED ######
####### Read the data from the allmmseq directory #######
path="adipo/ETH/mapped/pop/mmseq/"
outpath="adipo/ETH/mapped/pop/mmseq/plots/"
mynames=list.files(path,pattern="gene.mmseq")

##### Read information about samples (decide which cells not to include) #####
shortnames=mynames=substr(mynames,1,12)
names(mynames)=mynames

#### save the parameters ####
##### Read mmseq cells ####
source(paste("src/mmseq-latest/","/src/R/mmseq.R",sep=""))
mmseqpath=path
scmmseq_files=paste(mmseqpath,list.files(mmseqpath,pattern=".gene.mmseq"),sep="")
scmmseq_files=sapply(mynames,function(x) scmmseq_files[grep(x,scmmseq_files)])
#scfiles = readmmseq(scmmseq_files,names(mynames),normalize=F)
#scfiles.mmseqnorm = readmmseq(scmmseq_files,names(mynames),normalize=T)
#save(scfiles,scfiles.mmseqnorm,scmmseq_files,mynames,shortnames,file=paste(outpath,myname,"-",myset,"-",mydate,".RData",sep=""))
load(paste(outpath,myname,"-",myset,"-",mydate,".RData",sep="")) #scfiles and scfiles.mmseqnorm #SC-bow-mmseqEnsg75-runCut-20140609-140717.RData
mmseqcounts=apply(scfiles[[10]], 2, function(x) round(x))



##### Read STAR cells ####
# path to the gsnap/htseq results
starpath="adipo/ETH/mapped/pop/STAR/"
#readSTARHTSeqASC(starpath,"cut",mynames,names(mynames),mydate) 
load(paste(outpath,"SSAll-STARHTSeqTags.-",mydate,".RData",sep="")) #htseqgenes and htseqexons #SSAll-STARHTSeqTags.-140717.RData
htseqgenes.pop=htseqgenes;htseqstats.pop=htseqstats
save(htseqgenes.pop,htseqstats.pop=htseqstats,file=paste(outpath,"SSAll-STARHTSeqTags.-",mydate,".Counts.RData",sep=""))


