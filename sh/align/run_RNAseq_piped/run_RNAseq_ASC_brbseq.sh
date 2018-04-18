
##### Start every script with copying over the toolbox utils (sh functions from main) ######
# they are in 

###### Set the directory ########
toolpath="/toolbox/"
rundir="/home/pschwali/asc/"
cd $rundir

###### Choose the files to work with ########
inpath="data/brbseq/"
outpath="data/brbseq/"
myname="Nov16"
mysamples="A10.fastq A11.fastq A12.fastq A1.fastq A2.fastq A3.fastq A4.fastq A5.fastq A6.fastq A7.fastq A8.fastq A9.fastq B10.fastq B11.fastq B12.fastq B1.fastq B2.fastq B3.fastq B4.fastq B5.fastq B6.fastq B7.fastq B8.fastq B9.fastq C10.fastq C11.fastq C12.fastq C1.fastq C2.fastq C3.fastq C4.fastq C5.fastq C6.fastq C7.fastq C8.fastq C9.fastq D10.fastq D11.fastq D12.fastq D1.fastq D2.fastq D3.fastq D4.fastq D5.fastq D6.fastq D7.fastq D8.fastq D9.fastq E1.fastq E2.fastq E3.fastq E4.fastq"

###### Clean up the reads and QC ########																																									                         
sh ${toolpath}/submit-prinseq.sh $inpath $outpath "$mysamples"
for myfile in $mysamples
	do
	myfor="$inpath/$myfile"	
	bsub -M 4000000 -R  "rusage[mem=4000]" -o cutadapt-$myfile.out -e cutadapt-$myfile.err -J cutadapt-$myfile "cutadapt -b GTACTCTGCGTTGATACCAC -m 25 -q 20 $myfor-good.fastq > $myfor-cut.fastq"
done
for mysample in $mysamples
	do
	### se
	myfor="$inpath/${mysample}-cut.fastq"
	trimopt="-q 15 --fastqc --stringency 3 -e 0.05 --length 30 -o $inpath/trimmed/" ### for single end
	bsub -M 4000000 -R  "rusage[mem=4000]" -o fastqc-$mysample.out -e fastqc-$mysample.err -J fastqc-$mysample "trim_galore $trimopt $myfor"
	done
	
###### and check the trimming stats ######
sh $toolpath/get_trimmed_reads.sh $inpath/trimmed/ $inpath/trimmed/${myname}_TrimGalore_trimming_report.txt $mysamples


###### Choose the files to work with ########
inpath="data/brbseq/"
outpath="mapped/brbseq/"
myname="Nov16"

###### Run BOWTIE|MMSEQ ALINGMENT ##############
### params
#myset="bow-mmseqEnsg75-runCut-$(date +%Y%m%d)"
#myset="bow-mmseqEnsg75-runCut-20160626"
#transcripts="genomes/annotation/Mus_musculus.GRCm38.75.mmseq.spiked-tdRFP-GFP"
#transcriptsFa="genomes/annotation/Mus_musculus.GRCm38.75.mmseq.spiked-tdRFP-GFP.fa"
#sh ${toolpath}/submit-bowtie.sh $myset $transcripts $inpath $outpath "$myfiles"
#sh ${toolpath}/get_bowtie_info.sh $outpath/logs/ $outpath/logs/Bowtie-$myname-$myset-alninfo-final.txt #### followed by process-bowtie.report.R for plotting
#sh ${toolpath}/submit-bamtohits.sh $myset $transcriptsFa $inpath $outpath "$myfiles"
#sh ${toolpath}/submit-mmseq.sh $myset $inpath $outpath "$myfiles"
#sh ${toolpath}/get-bam2hits-info.sh $outpath/logs/ $outpath/logs/Bam2Hits-$myname-$myset-final.txt

####### Run STAR ALIGNMENT #####################
inpath="data/brbseq/trimmed/"
outpath="mapped/brbseq/STAR/"
myset="star-Ensg84-runCut-20161111"
genomedir="genomes/STAR-Mus_musculus.GRCm38.84_mm10_spiked_rfp/"
transcriptsFa="genomes/annotation/Mus_musculus.GRCm38.84.mmseq.fa"
myfiles="A10 A11 A12 A1 A2 A3 A4 A5 A6 A7 A8 A9 B10 B11 B12 B1 B2 B3 B4 B5 B6 B7 B8 B9 C10 C11 C12 C1 C2 C3 C4 C5 C6 C7 C8 C9 D10 D11 D12 D1 D2 D3 D4 D5 D6 D7 D8 D9 E10 E11 E12 E1 E2 E3 E4"

##### Submit for both transcriptome and genome mapping ######
sh ${toolpath}/submit-star-send.sh $myset $genomedir $inpath $outpath ".fastq-cut_trimmed.fq" "$myfiles"
sh ${toolpath}/get_star_info.sh $outpath $outpath/STAR-align.txt "" ###### Get information about the read mapping 

outpath="mapped/brbseq/STAR/"
inpath="mapped/brbseq/STAR/"
ggenfile="genomes/annotation/Mus_musculus.GRCm38.84.red.spiked-rfp-gfp.withexonid.for-tophat.gtf"
ggenome="genomes/mm10-spiked-rfp.fa"
sh ${toolpath}/submit-htseqcount-fromSTAR.sh $myset $ggenome $ggenfile $inpath $outpath "$myfiles"

#### 
#mv to /bam
myfiles="A10 A11 A12 A1 A2 A3 A4 A5 A6 A7 A8 A9 B10 B11 B12 B1 B2 B3 B4 B5 B6 B7 B8 B9 C10 C11 C12 C1 C2 C3 C4 C5 C6 C7 C8 C9 D10 D11 D12 D1 D2 D3 D4 D5 D6 D7 D8 D9 E10 E11 E12 E1 E2 E3 E4"
for myfile in $myfiles
do
mv $myfile/Aligned.sortedByCoord.out.bam bam/${myfile}.STAR.sortedByCoord.bam
done

while read -r old new; do
     mv "${old}.STAR.sortedByCoord.bam" "BD_Nov16_mAreg_${new}.STAR.sortedByCoord.bam"
done < map.txt



