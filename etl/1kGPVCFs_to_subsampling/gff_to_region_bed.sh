#!/bin/bash
# given a region type, get the coordinates in non overlapping bed format
set -e

region="exon"
sourcedir="/data/1KG_GRCh38/"
gff=$sourcedir"/Homo_sapiens.GRCh38.111.gff3"
bed=${sourcedir}${region}.bed

# select lines in gff where the region is an exon/gene/etc
# get the chromosome name and the coordinates: beware, bed files are 0 based opebed intervals
# while gff are 1 based closed intervals
# sort and unique lines to lighten the burden of the bcftools merge
awk -v r="$region" '$3==r {print $1, $4-1, $5}' $gff | sort -V -k1,1 -k2,2 -k3,3 | uniq > $bed
# make the sepatator of bed a tab
awk '{print "chr"$1,$2,$3}' $bed  | sed -e "s/ \+/\t/g" -e 's/\t*$//' > ${bed}.temp && \
mv ${bed}.temp ${bed}
# merge intervals
# beware, once this is done, the intervals no longer represent finite unique genes/exons/etc
# as overlapping genes/exons/etc are merged
bedtools merge -i ${bed} > ${bed}.temp && mv ${bed}.temp ${bed}