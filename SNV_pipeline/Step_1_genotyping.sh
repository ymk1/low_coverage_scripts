#!/bin/bash

	# BASH Script for Running Initial Platypus Variant Calling Across Samples

OUTDIR=$1
REFERENCE=$2
BAMLIST=$3
        
while read FILE; do 
            NAME=`basename $FILE`
	   if [ -f "$OUTDIR/1_individual_calls/platypusVariants_"${NAME%.*}"_default.vcf" ]; then 
	    	echo "$NAME  exists\n"
	   else
            echo -e "\nCalling on $NAME\n"
            
             Platypus.py callVariants \
            		--logFileName=$OUTDIR/1_individual_calls/"${NAME%.*}"_default.log \
            		--refFile=$REFERENCE \
            		--bamFiles=$FILE \
            		--minPosterior=0 \
	    		--minBaseQual=30 \
	    		--filteredReadsFrac=.8 \
	    		--badReadsThreshold=30 \
	    		--badReadsWindow=15 \
	    		--minFlank=0 \
            		--minReads=1 \
            		--nCPU=8 \
            		-o $OUTDIR/1_individual_calls/platypusVariants_"${NAME%.*}"_default.vcf
	 fi
        done <$BAMLIST

