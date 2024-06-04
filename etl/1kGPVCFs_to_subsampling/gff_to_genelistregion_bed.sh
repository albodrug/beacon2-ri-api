#!/bin/bash
# given a region type, get the coordinates in non overlapping bed format
set -e

genelist=$1 # file with list of genes
region="exon"
sourcedir="/data/1KG_GRCh38/"
gff=$sourcedir"Homo_sapiens.GRCh38.111.gff3"
regionbed=${sourcedir}${region}.bed
genelistbed=${sourcedir}$(basename "$genelist" .lst).bed

echo "Working directory: $sourcedir"
echo "Using $genelist"
echo "Using $gff"
echo "Output: $regionbed"
echo "Output: $genelistbed"

# extract gene coordinates of gene with ids in the genelist file
while read geneid; do 
  grep "Name=${geneid}" $gff ;
done < $genelist | awk '{print $1, $4-1, $5}' | sort -V -k1,1 -k2,2 -k3,3 | uniq > $genelistbed
awk '{print "chr"$1,$2,$3}' $genelistbed  | sed -e "s/ \+/\t/g" -e 's/\t*$//' > ${genelistbed}.temp && \
mv ${genelistbed}.temp ${genelistbed}
bedtools merge -i ${genelistbed} > ${genelistbed}.temp && mv ${genelistbed}.temp ${genelistbed}

# select lines in gff where the region is an exon/gene/etc
# get the chromosome name and the coordinates: beware, bed files are 0 based opebed intervals
# while gff are 1 based closed intervals
# sort and unique lines to lighten the burden of the bcftools merge
awk -v r="$region" '$3==r {print $1, $4-1, $5}' $gff | sort -V -k1,1 -k2,2 -k3,3 | uniq > $regionbed
# make the sepatator of bed a tab
awk '{print "chr"$1,$2,$3}' $regionbed  | sed -e "s/ \+/\t/g" -e 's/\t*$//' > ${regionbed}.temp && \
mv ${regionbed}.temp ${regionbed}
# merge intervals
# beware, once this is done, the intervals no longer represent finite unique genes/exons/etc
# as overlapping genes/exons/etc are merged
bedtools merge -i ${regionbed} > ${regionbed}.temp && mv ${regionbed}.temp ${regionbed}

# intersect the region bed and the gene list bed to get only exons/cds/etc of genes in list
bedtools intersect -a $regionbed -b $genelistbed > ${sourcedir}${region}_$(basename "$genelist" .lst).bed
echo "Finished."