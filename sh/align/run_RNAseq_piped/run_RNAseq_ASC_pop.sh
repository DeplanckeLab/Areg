###### Set the directory ########
toolpath="toolbox/sh/"
rundir="asc/"
cd $rundir


###### Clean up the reads and QC ########																																									                         
#sh ${toolpath}/submit-prinseq.sh $inpath $outpath "$myfiles"
#sh ${toolpath}/submit-cutadapt.sh $inpath $outpath "$myfiles"
#sh ${toolpath}/submit-fastqc.sh $inpath $outpath "$myfiles"
#sh ${toolpath}/get_qc_info.sh ${outpath}/QC/ "$myfiles" ##### Information in ${outpath}/QC/Bigfilter_report.txt FastQC_quality_report.txt CutAdapt_report.txt

###### Choose the files to work with ########
inpath="data/pop/"
outpath="data/pop/"
myfiles="BD_P1_ALL_ALL_L1_R1_001.fastq-cut.fastq BD_P2_ALL_ALL_L1_R1_001.fastq-cut.fastq BD_P3_ALL_ALL_L1_R1_001.fastq-cut.fastq NextSub1Rm.fastq-cut.fastq NextSub1Rp.fastq-cut.fastq NextSub2Rm.fastq-cut.fastq NextSub2Rp.fastq-cut.fastq NextSub3Rm.fastq-cut.fastq NextSub3Rp.fastq-cut.fastq NextSub4Rm.fastq-cut.fastq NextSub4Rp.fastq-cut.fastq NextSub5Rm.fastq-cut.fastq NextSub6Rm.fastq-cut.fastq"
myname="Pop"


###### Choose the files to work with ########
inpath="data/pop/"
outpath="mapped/pop/"

###### Run BOWTIE|MMSEQ ALINGMENT ##############
### params
#myset="bow-mmseqEnsg75-runCut-$(date +%Y%m%d)"
myset="bow-mmseqEnsg84-runCut-20161109"
transcripts="genomes/annotation/Mus_musculus.GRCm38.84.mmseq.wbwt"
transcriptsFa="genomes/annotation/Mus_musculus.GRCm38.84.mmseq.fa"

sh ${toolpath}/submit-bowtie.sh $myset $transcripts $inpath $outpath "$myfiles"
sh ${toolpath}/get_bowtie_info.sh $outpath/logs/ $outpath/logs/Bowtie-$myname-$myset-alninfo-pop.txt #### followed by process-bowtie.report.R for plotting
sh ${toolpath}/submit-bamtohits.sh $myset $transcriptsFa $inpath $outpath "$myfiles"
sh ${toolpath}/submit-mmseq.sh $myset $inpath $outpath "$myfiles"
sh ${toolpath}/get-bam2hits-info.sh $outpath/logs/ $outpath/logs/Bam2Hits-$myname-$myset-final.txt


####### Run STAR ALIGNMENT #####################
inpath="data/pop/"
outpath="mapped/pop/STAR/"
myset="star-Ensg84-runCut-20161109"
myext=""
genomedir="genomes/STAR-Mus_musculus.GRCm38.84_mm10_spiked_rfp/"
transcriptsFa="genomes/annotation/Mus_musculus.GRCm38.84.mmseq.fa"

##### Submit for both transcriptome and genome mapping ######
sh ${toolpath}/submit-star-send.sh $myset $genomedir $inpath $outpath "$myext" "$myfiles"
sh ${toolpath}/get_star_info.sh $outpath $outpath/STAR-align-pop.txt cut ###### Get information about the read mapping --- add to bigtable eventually

outpath="mapped/pop/STAR/"
inpath="mapped/pop/STAR/"
ggenfile="genomes/annotation/Mus_musculus.GRCm38.84.red.spiked-rfp-gfp.withexonid.for-tophat.gtf"
ggenome="genomes/mm10-spiked-rfp.fa"
sh ${toolpath}/submit_htseqcount-fromSTAR.sh $myset $ggenome $ggenfile $inpath $outpath "$myfiles"

ext="fastq.gz"
while read -r old new; do
     mv "${old}.${ext}" "BD_Jan15_mAreg_${new}.${ext}"
done < map.txt
