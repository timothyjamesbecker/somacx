#!/bin/bash
#usage: ./art_sim.sh tool_path genome.fa out_dir 20 1 #this will do depth of 20 in 2 lanes starting at 1

BIO=$1      #has path to art_illumina
FASTA_IN=$2 #input fasta sequence
OUT_DIR=$3  #output directory for FASTQ files
DEP=$4      #average depth of coverage
LN=$5       #name it for a lane (multi-lane bam files need this)

DEP=$(awk 'BEGIN {print '$DEP'/2.0}')
mkdir -p $OUT_DIR/artsim

#MEAN=$(awk 'BEGIN {srand(); printf("%3.0f",175+50*rand())}')
#SD=$(awk 'BEGIN {srand(); printf("%3.0f",75+50*rand())}')
$BIO/art_illumina -ss HS25 -i $FASTA_IN -l 125 -f $DEP -m 200 -s 75 -o $OUT_DIR/lane$LN.
rm $OUT_DIR/*.aln || true
pigz -9 $OUT_DIR/lane$LN.1.fq
pigz -9 $OUT_DIR/lane$LN.2.fq
mv $OUT_DIR/lane$LN.1.fq.gz $OUT_DIR/artsim/lane$LN.1.fq.gz
mv $OUT_DIR/lane$LN.2.fq.gz $OUT_DIR/artsim/lane$LN.2.fq.gz

LN=$(awk 'BEGIN{print '$LN'+1;}')
$BIO/art_illumina -ss HS25 -i $FASTA_IN -l 125 -f $DEP -m 175 -s 110 -o $OUT_DIR/lane$LN.
rm $OUT_DIR/*.aln || true
pigz -9 $OUT_DIR/lane$LN.1.fq
pigz -9 $OUT_DIR/lane$LN.2.fq
mv $OUT_DIR/lane$LN.1.fq.gz $OUT_DIR/artsim/lane$LN.1.fq.gz
mv $OUT_DIR/lane$LN.2.fq.gz $OUT_DIR/artsim/lane$LN.2.fq.gz
