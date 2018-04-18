
#### get info of all the STAR ran alignments in a folder ####
## run as get_star_info.sh /adipo/ adipo/Star-mmm10-run1-alninfo.txt STAR
path=$1
outfile=$2
mypattern="$3"

for filename in ${path}/*${mypattern}*/
        do
        myname=`basename $filename`
        myinfo=`grep "Number of input reads" $filename/Log.final.out`
         echo $myinfo >> $outfile
        myinfo=`grep "Uniquely mapped reads number" $filename/Log.final.out`
         echo $myinfo >> $outfile
        myinfo=`grep "Number of splices: Total" $filename/Log.final.out`
         echo $myinfo >> $outfile
        myinfo=`grep "Number of splices: Annotated (sjdb)" $filename/Log.final.out`
         echo $myinfo >> $outfile
        myinfo=`grep "Number of reads mapped to multiple loci" $filename/Log.final.out`
         echo $myinfo >> $outfile
        myinfo=`grep "Number of reads mapped to too many loci" $filename/Log.final.out`
         echo $myinfo >> $outfile
        myinfo=`grep "reads unmapped: too many mismatches" $filename/Log.final.out`
         echo $myinfo >> $outfile
        myinfo=`grep "of reads unmapped: too short" $filename/Log.final.out`
         echo $myinfo >> $outfile
        myinfo=`grep "reads unmapped: other" $filename/Log.final.out`
         echo $myinfo >> $outfile
        myinfo=`echo $myname | echo ${myname}`
        echo $myinfo >> $outfile
        done