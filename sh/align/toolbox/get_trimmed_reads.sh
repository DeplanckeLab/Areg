#### goes through the trimmed reads and reports stuff out ####
## run as get_trimmed_reads.sh $inpath/ TrimGalore_trimming_report.txt $myfiles

path=$1
outfile=$2
myfiles="$3"

for filename in ${path}/*trimming_report.txt
	do
	myname=`basename $filename`
	myinfo1=`grep "Processed reads:" ${filename}`
	myinfo2=`grep "Trimmed reads:" ${filename}`
	myinfo3=`grep "Quality-trimmed:" ${filename}`
	myinfo4=`grep "shorter than" ${filename}`
	myinfo=`echo ${myname} -- ${myinfo1} ${myinfo2} ${myinfo3} ${myinfo4}`
	echo $myinfo >> $outfile
	done

##### and also the fastqc result after that #####
outfile="${path}/FastQC_quality_report.txt"
#ls *fq > myfiles
#myfiles=`more myfiles`
for myfile in $myfiles
	do
	temp=`more ${path}/${myfile}-good_fastqc/summary.txt`
	echo $temp >> $outfile
	done

outfile="${path}/FastQC_quality_report.txt"
#ls *fq > myfiles
#myfiles=`more myfiles`
for myfile in $myfiles
	do
	temp=`more ${path}/${myfile}_trimmed_fastqc/summary.txt`
	echo $temp >> $outfile
	done	
	
	
outfile="${path}/FastQC_quality_report_paired.txt"
#ls *fq > myfiles
#myfiles=`more myfiles`
for myfile in $myfiles
	do
	temp=`more ${path}/${myfile}_R1_val_1_fastqc/summary.txt`
	echo $temp >> $outfile
	temp=`more ${path}/${myfile}_R2_val_2_fastqc/summary.txt`
	echo $temp >> $outfile
	done

outfile="${path}/FastQC_quality_report_paired.txt"
#ls *fq > myfiles
#myfiles=`more myfiles`
for myfile in $myfiles
	do
	temp=`more ${path}/${myfile}_1_val_1_fastqc/summary.txt`
	echo $temp >> $outfile
	temp=`more ${path}/${myfile}_2_val_2_fastqc/summary.txt`
	echo $temp >> $outfile
	done		
	
	
