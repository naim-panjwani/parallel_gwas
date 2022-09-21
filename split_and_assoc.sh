#!/bin/bash

vcffile=$1
maxsnplines=$2
outfile=$3
tmpdir="`pwd`/job_output/"

date
me=`basename "$0"`
command=$(echo "$me $vcffile $maxsnplines $outfile")
echo
echo "Starting $command"

export PATH="$PATH:/hpf/largeprojects/struglis/naim/miniconda3/envs/rbioenv/bin:/hpf/largeprojects/struglis/naim/miniconda3/condabin:/usr/local/bin:/usr/local/sbin:/usr/lib64/qt-3.3/bin:/opt/moab/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/sbin:/usr/sbin:/home/naim/bin:/hpf/largeprojects/struglis/naim/miniconda3/envs/rbioenv/bin:/hpf/largeprojects/struglis/naim/miniconda3/condabin"

echo
echo "Splitting $vcffile"
echo
#split -d -l $maxsnplines $vcffile ${vcffile%.vcf}_split_ --verbose
for i in $(ls ${vcffile%.vcf}_split_[0-9][0-9]); do
#    echo "Adding header to $i"
    j=${i}_withHeader
#    cat ${vcffile%.vcf}_header.vcf $i > $j
    echo "Submitting association analysis job for $j"
    ajob=$(qsub -v infile=$j,outfile=${i}_assoc.txt -l walltime=23:59:00 -l nodes=1:ppn=1 -l mem=20g -l vmem=20g -o ./job_output/ -e ./job_output/ -d `pwd` -N photosens_assoc_$i code/run_logistic_assoc_topmed.pbs)
    jobs="${jobs}:${ajob}"
done
jobs=$(echo $jobs |sed 's/^://g')
echo
echo "Association analysis for $vcffile. Job IDs:"
echo "$jobs"
tmpvar=$(basename $vcffile)
bjob=$(qsub -W depend=afterok:$jobs -F "$vcffile" -l walltime=23:59:00 -l nodes=1:ppn=1 -l mem=10g -l vmem=10g -o ${tmpdir} -e $tmpdir -d `pwd` -N "assocbind_$tmpvar" code/combine.sh)
infile=${vcffile%.vcf}_assoc.txt.gz
cjob=$(qsub -W depend=afterok:$bjob -F "$infile $outfile" -l walltime=23:59:00 -l nodes=1:ppn=1 -l mem=8g -l vmem=8g -o ${tmpdir} -e $tmpdir -d `pwd` -N "move_$tmpvar" code/simple_move.sh)
echo
echo "Done $command"
date
