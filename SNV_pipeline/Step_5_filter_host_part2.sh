#!/bin/bash

# FILTERING USING GENOTYPED HOST TUMOR FILE

VCFLISTDIR=$1
OUTDIR=$2

	ls -1 $VCFLISTDIR/final*multiallelic.vcf > $OUTDIR/vcflist.txt
	
	while read FILE; do
		echo "$FILE"
		Rscript 5_Variants_filtering.R $FILE $OUTDIR

		bgzip $FILE

	done<$OUTDIR/vcflist.txt
