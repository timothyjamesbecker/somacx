#!/bin/bash
#usages pb_T_align.sh tool_path fastq_in_dir ref_path bam_out_dir sample threads

BIO=$1      #has path to minimap2,samtools,sambamba
IN_DIR=$2   #the fastq input directory
REF=$3      #the reference fasta to align to
OUT_DIR=$4  #the output directory to store the resulting BAM files in
SM=$5       #the bam sample ID that will also be the name of the bam file
TH=$6       #the number of threads for use with minimap2,samtools,sambamba

RG1="@RG\tID:"$SM"_lane1_T\tLB:"$SM"_lane1_T\tPL:PACBIO\tPU:"$SM"_lane1_T\tSM:"$SM
$BIO/minimap2 -ax map-pb -Y -t $TH -R $RG1 $REF $IN_DIR/lane1.fq.gz | $BIO/samtools view - -Shb > $OUT_DIR/$SM.lane1_T.bam
$BIO/samtools sort -l 9 -@ $TH -T _sort -o $OUT_DIR/$SM.lane1_T.sorted.bam $OUT_DIR/$SM.lane1_T.bam
rm $OUT_DIR/$SM.lane1_T.bam || true

RG2="@RG\tID:"$SM"_lane2_T\tLB:"$SM"_lane2_T\tPL:PACBIO\tPU:"$SM"_lane2_T\tSM:"$SM
$BIO/minimap2 -ax map-pb -Y -t $TH -R $RG2 $REF $IN_DIR/lane2.fq.gz | $BIO/samtools view - -Shb > $OUT_DIR/$SM.lane2_T.bam
$BIO/samtools sort -l 9 -@ $TH -T _sort -o $OUT_DIR/$SM.lane2_T.sorted.bam $OUT_DIR/$SM.lane2_T.bam
rm $OUT_DIR/$SM.lane2_T.bam || true

RG3="@RG\tID:"$SM"_lane3_T\tLB:"$SM"_lane3_T\tPL:PACBIO\tPU:"$SM"_lane3_T\tSM:"$SM
$BIO/minimap2 -ax map-pb -Y -t $TH -R $RG3 $REF $IN_DIR/lane3.fq.gz | $BIO/samtools view - -Shb > $OUT_DIR/$SM.lane3_T.bam
$BIO/samtools sort -l 9 -@ $TH -T _sort -o $OUT_DIR/$SM.lane3_T.sorted.bam $OUT_DIR/$SM.lane3_T.bam
rm $OUT_DIR/$SM.lane3_T.bam || true

RG4="@RG\tID:"$SM"_lane4_T\tLB:"$SM"_lane4_T\tPL:PACBIO\tPU:"$SM"_lane4_T\tSM:"$SM
$BIO/minimap2 -ax map-pb -Y -t $TH -R $RG4 $REF $IN_DIR/lane4.fq.gz | $BIO/samtools view - -Shb > $OUT_DIR/$SM.lane4_T.bam
$BIO/samtools sort -l 9 -@ $TH -T _sort -o $OUT_DIR/$SM.lane4_T.sorted.bam $OUT_DIR/$SM.lane4_T.bam
rm $OUT_DIR/$SM.lane4_T.bam || true

$BIO/samtools merge -l 9 -@ $TH -O BAM $OUT_DIR/$SM.T.merged.bam $OUT_DIR/$SM.lane1_T.sorted.bam $OUT_DIR/$SM.lane2_T.sorted.bam $OUT_DIR/$SM.lane3_T.sorted.bam $OUT_DIR/$SM.lane4_T.sorted.bam
rm $OUT_DIR/$SM.lane*_T.sorted.bam || true

$BIO/sambamba markdup -l 9 -t $TH --tmpdir=_sort --sort-buffer-size=1024 --overflow-list-size=20000 $OUT_DIR/$SM.T.merged.bam $OUT_DIR/$SM.T.final.bam
rm $OUT_DIR/$SM.T.merged.bam* || true