#!/bin/bash
# Decompresses and filters the vcf file for the given r2 and maf filters and extracts the dosage (DS) field
# Syntax: bash decompress_and_filter_vcf.sh <vcf filename> <r2 info field name> <r2 filter> <maf filter> <output filename>
# Example: bash decompress_and_filter_vcf.sh data/raw_data/imputed_data_topmed_chr1.vcf.gz R2 0.4 0.01 tempfiles2012242/imputed_data_topmed_chr1.vcf


vcf=$1
r2info=$2
r2=$3
maf=$4
afile=$5

date
me=`basename "$0"`

echo "Starting $me"
echo "vcf: $vcf"
echo "r2info: $r2info"
echo "r2: $r2"
echo "maf: $maf"
echo "output: $afile"

export PATH="$PATH:/hpf/largeprojects/struglis/naim/miniconda3/envs/rbioenv/bin:/hpf/largeprojects/struglis/naim/miniconda3/condabin:/usr/local/bin:/usr/local/sbin:/usr/lib64/qt-3.3/bin:/opt/moab/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/sbin:/usr/sbin:/home/naim/bin:/hpf/largeprojects/struglis/naim/miniconda3/envs/rbioenv/bin:/hpf/largeprojects/struglis/naim/miniconda3/condabin"

r2str=$(echo "INFO/${r2info}<${r2}")
mafstr=$(echo "${maf}:minor")
echo
echo "Decompressing ${vcf} for ${r2info}>=${r2} and maf>${maf}"
vcftmp1="${afile%.vcf}_tmp1.vcf"
bcftools filter -Ou -e $r2str $vcf |bcftools view -q $mafstr > $vcftmp1
numheaderlines=$(head -1000 $vcftmp1 |grep ^"#" |wc -l)
head -n $numheaderlines $vcftmp1 > "${afile%.vcf}_header.vcf"
echo "Extracting dosage field"
paste <(sed "1,${numheaderlines}d" $vcftmp1 |cut -f1-9) <(bcftools query -f '[\t%DS]\n' $vcftmp1 |sed 's/^\t//g')  > $afile
rm -f $vcftmp1
echo "Done $me $vcf"
date
