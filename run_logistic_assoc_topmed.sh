#!/bin/bash
# Script will go over all files specified in $vcflist file, split them by the $numlines splits, do association analysis, then recombine the association results in similar order and number of files as in $vcflist, and move the resulting files as specified in $vcflist's 2nd column
# Syntax: bash run_linear_assoc_topmed.sh <vcflist file> <number of max lines per file> <r2 filter> <maf filter>
# Example: bash run_linear_assoc_topmed.sh vcflist.txt 400000 0.4 0.01
# Depends on: decompress_and_filter_vcf.sh, split_and_assoc.sh, run_linear_assoc_topmed.pbs, linear_assoc_topmed.R, combine.sh, simple_move.sh

vcflist=$1
numlines=$2
r2=$3
maf=$4

date
me=`basename "$0"`
command=$(echo "$me $vcflist $numlines $r2 $maf")
echo
echo "Starting $command"

# check inputs ok:
echo "Checking inputs ok"
#[ ! $(wc -l $vcflist |cut -d' ' -f1) -eq $(wc -l $outlist |cut -d' ' -f1) ] && echo "$vcflist and $outlist do not have an equal number of lines"
while IFS= read -r line; do
  vcf=$(echo $line |cut -d' ' -f1)
  [ ! -f $vcf ] && echo "$vcf does not exist" && exit 1
done < $vcflist
echo "Inputs pass"

# conda path (conda activate won't work):
export PATH="$PATH:/hpf/largeprojects/struglis/naim/miniconda3/envs/rbioenv/bin:/hpf/largeprojects/struglis/naim/miniconda3/condabin:/usr/local/bin:/usr/local/sbin:/usr/lib64/qt-3.3/bin:/opt/moab/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/sbin:/usr/sbin:/home/naim/bin:/hpf/largeprojects/struglis/naim/miniconda3/envs/rbioenv/bin:/hpf/largeprojects/struglis/naim/miniconda3/condabin"

export LC_NUMERIC="en_US" # to print thousands comma separator in numeric variables
tempnum=2012243
tempdir="tempfiles${tempnum}"
[ ! -d $tempdir ] && mkdir $tempdir


while IFS= read -r line; do # for each vcf file (likely split per chromosome)
  vcf=$(echo $line |cut -d' ' -f1)
  outfile=$(echo $line |cut -d' ' -f2)
  afile="$tempdir/$(basename $vcf)"
  afile=${afile%.gz}
  echo
  #echo "Sending prepjob to decompress and filter $vcf"
  #bash code/decompress_and_filter_vcf.sh $vcf "R2" 0.4 0.01 $afile # filter for R2>=0.4 and MAF>0.01, and extract the dosage (DS) field
  tmpvar=$(basename $afile)
  tmpdir="`pwd`/job_output/"
  #prepjob=$(qsub -F "$vcf R2 $r2 $maf $afile" -l walltime=23:59:00 -l nodes=1:ppn=1 -l mem=20g -l vmem=20g -o ${tmpdir} -e $tmpdir -d `pwd` -N "decompress_$tmpvar" code/decompress_and_filter_vcf.sh)
  echo "Sending association jobs for $afile"
  #assoc_job=$(qsub -W depend=afterok:$prepjob -F "$afile $numlines $outfile" -l walltime=23:59:00 -l nodes=1:ppn=1 -l mem=10g -l vmem=10g -o ${tmpdir} -e $tmpdir -d `pwd` -N "splitassoc_$tmpvar" code/split_and_assoc.sh)
  assoc_job=$(qsub -F "$afile $numlines $outfile" -l walltime=23:59:00 -l nodes=1:ppn=1 -l mem=10g -l vmem=10g -o ${tmpdir} -e $tmpdir -d `pwd` -N "splitassoc_$tmpvar" code/split_and_assoc.sh)
done < "$vcflist"

echo "Done $command"
