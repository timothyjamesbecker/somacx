#!/bin/bash
#usage ./minimap2_align.sh ./somacx_T_bam ref.mmi sample 4

BIO=$1 #has path to minimap2,samtools,sambamba
DIR=$2 #the output directory to store BAM files and tmp folders/files for alignment/sorting/merging/duplicate
REF=$3 #the reference fasta/mmi from minimap2
SM=$4  #the sample name to uses as the bam file identifier
TH=$5  #the number of threads to pass to alignment, sorting, merging and duplicate marking tools

RG1="@RG\tID:"$SM"_lane1_T\tLB:"$SM"_lane1_T\tPL:ILLUMINA\tPU:"$SM"_lane1_T\tSM:"$SM
$BIO/minimap2 -ax sr -Y -t $TH -R $RG1 $REF $DIR/lane1.1.fq.gz $DIR/lane1.2.fq.gz | $BIO/samtools view - -Shb > $DIR/$SM.lane1_T.bam
$BIO/samtools sort -l 9 -@ $TH -T _sort -o $DIR/$SM.lane1_T.sorted.bam $DIR/$SM.lane1_T.bam
rm $DIR/$SM.lane1_T.bam || true

RG2="@RG\tID:"$SM"_lane2_T\tLB:"$SM"_lane2_T\tPL:ILLUMINA\tPU:"$SM"_lane2_T\tSM:"$SM
$BIO/minimap2 -ax sr -Y -t $TH -R $RG2 $REF $DIR/lane2.1.fq.gz $DIR/lane2.2.fq.gz | $BIO/samtools view - -Shb > $DIR/$SM.lane2_T.bam
$BIO/samtools sort -l 9 -@ $TH -T _sort -o $DIR/$SM.lane2_T.sorted.bam $DIR/$SM.lane2_T.bam
rm $DIR/$SM.lane2_T.bam || true

RG3="@RG\tID:"$SM"_lane3_T\tLB:"$SM"_lane3_T\tPL:ILLUMINA\tPU:"$SM"_lane3_T\tSM:"$SM
$BIO/minimap2 -ax sr -Y -t $TH -R $RG3 $REF $DIR/lane3.1.fq.gz $DIR/lane3.2.fq.gz | $BIO/samtools view - -Shb > $DIR/$SM.lane3_T.bam
$BIO/samtools sort -l 9 -@ $TH -T _sort -o $DIR/$SM.lane3_T.sorted.bam $DIR/$SM.lane3_T.bam
rm $DIR/$SM.lane3_T.bam || true

RG4="@RG\tID:"$SM"_lane4_T\tLB:"$SM"_lane4_T\tPL:ILLUMINA\tPU:"$SM"_lane4_T\tSM:"$SM
$BIO/minimap2 -ax sr -Y -t $TH -R $RG4 $REF $DIR/lane4.1.fq.gz $DIR/lane4.2.fq.gz | $BIO/samtools view - -Shb > $DIR/$SM.lane4_T.bam
$BIO/samtools sort -l 9 -@ $TH -T _sort -o $DIR/$SM.lane4_T.sorted.bam $DIR/$SM.lane4_T.bam
rm $DIR/$SM.lane4_T.bam || true

$BIO/samtools merge -l 9 -@ $TH -O BAM $DIR/$SM.T.merged.bam $DIR/$SM.lane1_T.sorted.bam $DIR/$SM.lane2_T.sorted.bam $DIR/$SM.lane3_T.sorted.bam $DIR/$SM.lane4_T.sorted.bam
rm $DIR/$SM.lane*_T.sorted.bam || true

$BIO/sambamba markdup -l 9 -t $TH --tmpdir=_sort --sort-buffer-size=8096 --overflow-list-size=2000000 $DIR/$SM.T.merged.bam $DIR/$SM.T.final.bam
rm $DIR/$SM.T.merged.bam* || true
