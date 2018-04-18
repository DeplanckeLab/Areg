inpath="hsabrbseq/nxid1712/"
outpath="hsabrbseq/Lib10ng/"

##### Start every script with copying over the toolbox utils (sh functions from main) ######
# they are in 
###### Set the directory ########
toolpath="toolbox/sh/"
rundir="hsaasc/"
cd $rundir

###### Choose the files to work with ########
inpath="hsabrbseq/Demultiplexed_Lib6/"
outpath="hsabrbseq/Demultiplexed_Lib6/"
myname="Jan17"
mysamples="i04_T0_F3-.fastq i04_T0_F3+.fastq i04_T0_Total.fastq i05_T0_F3-.fastq i05_T0_F3+.fastq i05_T0_Total.fastq i06_T0_F3-.fastq i06_T0_F3+.fastq i06_T0_Total.fastq i08_T0_F3-.fastq i08_T0_F3+.fastq i08_T0_Total.fastq"

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



myname="Jan17"
###### Run BOWTIE|MMSEQ ALINGMENT ##############
### params
#myset="bow-mmseqEnsg75-runCut-$(date +%Y%m%d)"
#myset="bow-mmseqEnsg75-runCut-20160626"
#transcripts="genomes/annotation/Mus_musculus.GRCm38.75.mmseq.spiked-tdRFP-GFP"
#transcriptsFa="genomes/annotation/Mus_musculus.GRCm38.75.mmseq.spiked-tdRFP-GFP.fa"
#myfiles="A10_S.fastq-cut.fastq A11_S.fastq-cut.fastq A1_S1.fastq-cut.fastq A2_S2.fastq-cut.fastq A3_S3.fastq-cut.fastq A4_S4.fastq-cut.fastq A5_S5.fastq-cut.fastq A6_S6.fastq-cut.fastq A7_S7.fastq-cut.fastq A8_S8.fastq-cut.fastq A9_S9.fastq-cut.fastq B10_S.fastq-cut.fastq B11_S.fastq-cut.fastq B1_S1.fastq-cut.fastq B2_S1.fastq-cut.fastq B3_S1.fastq-cut.fastq B4_S1.fastq-cut.fastq B5_S1.fastq-cut.fastq B6_S1.fastq-cut.fastq B7_S1.fastq-cut.fastq B8_S1.fastq-cut.fastq B9_S2.fastq-cut.fastq C1_S2.fastq-cut.fastq C2_S2.fastq-cut.fastq C3_S2.fastq-cut.fastq C4_S2.fastq-cut.fastq C5_S2.fastq-cut.fastq C6_S2.fastq-cut.fastq C7_S2.fastq-cut.fastq C8_S3.fastq-cut.fastq C9_S3.fastq-cut.fastq D10_S.fastq-cut.fastq D11_S.fastq-cut.fastq D3_S3.fastq-cut.fastq D4_S3.fastq-cut.fastq D5_S3.fastq-cut.fastq D6_S3.fastq-cut.fastq D7_S3.fastq-cut.fastq D8_S3.fastq-cut.fastq D9_S3.fastq-cut.fastq E10_S.fastq-cut.fastq E11_S.fastq-cut.fastq E1_S4.fastq-cut.fastq E2_S4.fastq-cut.fastq E3_S4.fastq-cut.fastq E4_S4.fastq-cut.fastq E5_S4.fastq-cut.fastq E6_S4.fastq-cut.fastq E7_S4.fastq-cut.fastq E8_S4.fastq-cut.fastq E9_S4.fastq-cut.fastq F10_S.fastq-cut.fastq F11_S.fastq-cut.fastq F1_S5.fastq-cut.fastq F2_S5.fastq-cut.fastq F3_S5.fastq-cut.fastq F4_S5.fastq-cut.fastq F5_S5.fastq-cut.fastq F6_S5.fastq-cut.fastq F7_S5.fastq-cut.fastq F8_S5.fastq-cut.fastq F9_S6.fastq-cut.fastq G10_S.fastq-cut.fastq G11_S.fastq-cut.fastq G1_S6.fastq-cut.fastq G2_S6.fastq-cut.fastq G3_S6.fastq-cut.fastq G4_S6.fastq-cut.fastq G5_S6.fastq-cut.fastq G6_S6.fastq-cut.fastq G7_S6.fastq-cut.fastq G8_S7.fastq-cut.fastq G9_S7.fastq-cut.fastq H1_S7.fastq-cut.fastq H2_S7.fastq-cut.fastq H3_S7.fastq-cut.fastq H4_S7.fastq-cut.fastq H5_S7.fastq-cut.fastq H6_S7.fastq-cut.fastq H7_S8.fastq-cut.fastq H8_S8.fastq-cut.fastq"
#sh ${toolpath}/submit-bowtie.sh $myset $transcripts $inpath $outpath "$myfiles"
#sh ${toolpath}/get_bowtie_info.sh $outpath/logs/ $outpath/logs/Bowtie-$myname-$myset-alninfo-final.txt #### followed by process-bowtie.report.R for plotting
#sh ${toolpath}/submit-bamtohits.sh $myset $transcriptsFa $inpath $outpath "$myfiles"
#sh ${toolpath}/submit-mmseq.sh $myset $inpath $outpath "$myfiles"
#sh ${toolpath}/get-bam2hits-info.sh $outpath/logs/ $outpath/logs/Bam2Hits-$myname-$myset-final.txt

####### Run STAR ALIGNMENT #####################
inpath="hsabrbseq/Demultiplexed_Lib6/trimmed/"
outpath="mapped/hsabrbseq/Demultiplexed_Lib6/"

myset="star-Ensg84-runCut-20170106"
genomedir="genomes/STAR-Homo_sapiens.GRCh38.84_mm10_spiked_rfp/"
transcriptsFa="genomes/annotation/Homo_sapiens.GRCh38.84.mmseq.fa"
myfiles="i04_T0_F3- i04_T0_F3+ i04_T0_Total i05_T0_F3- i05_T0_F3+ i05_T0_Total i06_T0_F3- i06_T0_F3+ i06_T0_Total i08_T0_F3- i08_T0_F3+ i08_T0_Total"

##### Submit for both transcriptome and genome mapping ######
sh ${toolpath}/submit-star-send.sh $myset $genomedir $inpath $outpath ".fastq-cut_trimmed.fq" "$myfiles"
sh ${toolpath}/get_star_info.sh $outpath $outpath/STAR-align.txt "" ###### Get information about the read mapping --- add to bigtable eventually

inpath="mapped/hsabrbseq/Demultiplexed_Lib6/"
outpath="mapped/hsabrbseq/Demultiplexed_Lib6/"

####myfiles="BD_P2_A2_CGTACTAG-GCGTAAGA_L003_R1_100-cut.fastq"
ggenfile="genomes/annotation/Homo_sapiens.GRCh38.84.red.spiked-rfp-gfp.withexonid.for-tophat.gtf"
ggenome="genomes/hg19-spiked-rfp.fa"
sh ${toolpath}/submit_htseqcount-fromSTAR.sh $myset $ggenome $ggenfile $inpath $outpath "$myfiles"

myfiles="i04_T0_F3- i04_T0_F3+ i04_T0_Total i05_T0_F3- i05_T0_F3+ i05_T0_Total i06_T0_F3- i06_T0_F3+ i06_T0_Total i08_T0_F3- i08_T0_F3+ i08_T0_Total"
for myfile in $myfiles
do
cp $myfile/Aligned.sortedByCoord.out.bam bam/${myfile}.STAR.sortedByCoord.bam
done

ext="STAR.fastq.gz"
while read -r old new; do
     mv "${old}.${ext}" "BD_Jan17_mAreg_${new}.${ext}"
done < map.txt

