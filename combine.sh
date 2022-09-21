#!/bin/bash

vcffile=$1 # original vcf filename from which association analysis was done
outfile=${vcffile%.vcf}_assoc.txt


export PATH="$PATH:/hpf/largeprojects/struglis/naim/miniconda3/envs/rbioenv/bin:/hpf/largeprojects/struglis/naim/miniconda3/condabin:/usr/local/bin:/usr/local/sbin:/usr/lib64/qt-3.3/bin:/opt/moab/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/sbin:/usr/sbin:/home/naim/bin"

date
me=`basename "$0"`
command=$(echo "$me $vcffile $outfile")
echo
echo "Starting $command"
echo
echo "Combining into $outfile"
head -1 ${vcffile%.vcf}_split_00_assoc.txt > $outfile
for i in $(ls ${vcffile%.vcf}_split_*assoc.txt); do
  sed '1d' $i >> $outfile
done
echo "Compressing $outfile"
bgzip $outfile
echo "Tabix indexing $outfile"
tabix -p vcf "${outfile}.gz"
echo "Done compressing and indexing $outfile"
echo "Done $command"
date
