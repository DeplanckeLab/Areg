###### run as sh submit-prinseq.sh $inpath $outpath test.fastq

inpath=$1
outpath=$2
myfiles="$3"

for myfile in $myfiles
	do
	myfor="$inpath/$myfile"	
	bsub -M 4000000 -R  "rusage[mem=4000]" -q long -o bigfilter-$myfile.out -e bigfilter-$myfile.err -J bigfilter-$myfile "perl /home/pschwali/src/prinseq-lite-0.20.3/prinseq-lite.pl -fastq $myfor -custom_params 'A 70%;T 70%;G 70%;C 70%'  -trim_tail_left 36  -trim_tail_right 36 -lc_method dust -lc_threshold 45 -min_gc 1 -out_format 3 -verbose -out_good $myfor-good -out_bad $myfor-bad"
done
