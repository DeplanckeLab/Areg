##### Run HTseq ######
### run as sh submit-htseqcount-fromSTAR.sh $myset $ggenome $ggenfile $inpath $outpath "$myfiles"

myset=$1
ggenome=$2
ggenfile=$3
inpath=$4
outpath=$5
myfiles="$6"

for myfile in $myfiles
	do
	bsub -M 10000000 -R "rusage[mem=10000]" -o htseq-exons-$myset-$myfile-%J "samtools view ${inpath}/$myfile/Aligned.sortedByCoord.out.bam | htseq-count -m intersection-
nonempty -s no -a 10 -t exon -i exon_id - $ggenfile > ${outpath}/${myfile}/Aligned.sortedByCoord-exons.txt"
	bsub -M 10000000 -R "rusage[mem=10000]" -o htseq-genes-$myset-$myfile-%J "samtools view ${inpath}/$myfile/Aligned.sortedByCoord.out.bam | htseq-count -m intersection-
nonempty -s no -a 10 -t exon -i gene_id - $ggenfile > ${outpath}/${myfile}/Aligned.sortedByCoord-genes.txt"
	done