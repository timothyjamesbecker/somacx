#!/bin/bash
#usage ./minimap2_align.sh ./somacx_N_bam ref.mmi sample 4

BIO=$1 #has path to minimap2,samtools,sambamba
DIR=$2 #the output directory to store BAM files and tmp folders/files for alignment/sorting/merging/duplicate
REF=$3 #the reference fasta/mmi from minimap2
SM=$4  #the sample name to uses as the bam file identifier
TH=$5  #the number of threads to pass to alignment, sorting, merging and duplicate marking tools

RG1="@RG\tID:"$SM"_lane1_N\tLB:"$SM"_lane1_N\tPL:ILLUMINA\tPU:"$SM"_lane1_N\tSM:"$SM
$BIO/minimap2 -ax sr -Y -t $TH -R $RG1 $REF $DIR/lane1.1.fq.gz $DIR/lane1.2.fq.gz | $BIO/samtools view - -Shb > $DIR/$SM.lane1_N.bam
$BIO/samtools sort -l 9 -@ $TH -T _sort -o $DIR/$SM.lane1_N.sorted.bam $DIR/$SM.lane1_N.bam
rm $DIR/$SM.lane1_N.bam || true

RG2="@RG\tID:"$SM"_lane2_N\tLB:"$SM"_lane2_N\tPL:ILLUMINA\tPU:"$SM"_lane2_N\tSM:"$SM
$BIO/minimap2 -ax sr -Y -t $TH -R $RG2 $REF $DIR/lane2.1.fq.gz $DIR/lane2.2.fq.gz | $BIO/samtools view - -Shb > $DIR/$SM.lane2_N.bam
$BIO/samtools sort -l 9 -@ $TH -T _sort -o $DIR/$SM.lane2_N.sorted.bam $DIR/$SM.lane2_N.bam
rm $DIR/$SM.lane2_N.bam || true

$BIO/samtools merge -l 9 -@ $TH -O BAM $DIR/$SM.N.merged.bam $DIR/$SM.lane1_N.sorted.bam $DIR/$SM.lane2_N.sorted.bam
rm $DIR/$SM.lane*_N.sorted.bam || true

$BIO/sambamba markdup -l 9 -t $TH --tmpdir=_sort --sort-buffer-size=8096 --overflow-list-size=2000000 $DIR/$SM.N.merged.bam $DIR/$SM.N.final.bam
rm $DIR/$SM.N.merged.bam*  || true
