#### run STAR on single end data ###
# run as sh ${toolpath}/submit-star-send.sh $myset $genomedir $inpath $outpath "$myfiles"

#myset="star-Ensg75-runCut-20140610"
#genomedir="STAR-Mus_musculus.GRCm38.75_mm10_spiked_rfp_gfp/"
#inpath="sc/"
#outpath="sc/STAR/"
#myfiles="BD_P3_H9_GCTACGCT-CTAAGCCT_L001_R1_001.fastq-cut.fastq"

myset=$1
genomedir=$2
inpath=$3
outpath=$4
myext=$5
myfiles="$6"

for myfile in $myfiles
        do
        staropt="--runThreadN 4 --runMode alignReads --genomeDir $genomedir --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20  --alignIntronMax 1000000 --alignMatesGapMax 1000000 g --outSAMtype BAM SortedByCoordinate"
        infiles="--readFilesIn $inpath/${myfile}$myext"
        mkdir $outpath/${myfile}/
        bsub -M 30000000 -R  "rusage[mem=30000] span[ptile=4]" -n 4 -o $myset-$myfile.out -e $myset-$myfile.err -J $myset-$myfile "STAR  $staropt $infiles --outFileNamePrefix $outpath/${myfile}/"
        mkdir $outpath/Transcriptome-${myfile}/
        staropt="--runThreadN 4 --runMode alignReads --genomeDir $genomedir --outFilterType BySJout --outFilterMultimapNmax 100 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20  --alignIntronMax 1000000 --alignMatesGapMax 1000000 g --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM"
        bsub -M 30000000 -R  "rusage[mem=30000] span[ptile=4]" -n 4 -o $myset-$myfile.out -e $myset-$myfile.err -J $myset-$myfile "STAR  $staropt $infiles --outFileNamePrefix $outpath/Transcriptome-${myfile}/"
        done
