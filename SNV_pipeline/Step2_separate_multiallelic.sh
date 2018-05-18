#!/bin/bash

	# BASH Script for Separating VCF into Multiallelic

INDIR=$1 # Input Directory of Initial Run PLATYPUS Calling
OUTDIR=$2

#find $INDIR -maxdepth 1 -iregex '.+\.vcf' > $OUTDIR/2_filter_individual/vcf_list.txt
ls -1 -d $INDIR/*default.vcf >> $OUTDIR/2_filter_individual/vcf_list.txt
LIST=$OUTDIR/2_filter_individual/vcf_list.txt

while read FILE; do 
            NAME=`basename $FILE`
		if [ -f "$OUTDIR/2_filter_individual/"${NAME%.*}".multiallelic.vcf" ] || [ -f "$OUTDIR/2_filter_individual/"${NAME%.*}".multiallelic.vcf.gz" ]; then 
	   		echo "$NAME  exists\n"
         	else
        		echo "$NAME"
        	
		# 1. Break any Multi-allelic Variants into SNVs and merge it into the existing VCF
				vcfbreakmulti $FILE  >$OUTDIR/2_filter_individual/"${NAME%.*}".multiallelic.vcf
				bgzip $FILE
		fi
done<$LIST
