#!/bin/bash
#usage: ./pb_sim.sh tool_path source_path.fa out_dir 20 1 #depth=20 and this will use two lanes starting at lane 1

BIO=$1      #has path to pbsim
FASTA_IN=$2 #input fasta sequence
OUT_DIR=$3  #output directory for FASTQ files
DEP=$4      #average depth of coverage
LN=$5       #name it for a lane (multi-lane bam files need this)

DEP=$(awk 'BEGIN {print '$DEP'/2.0}')
mkdir -p $OUT_DIR/pbsim
rm $OUT_DIR/*.maf.gz || true
rm $OUT_DIR/*.ref || true

$BIO/pbsim --strategy wgs --method qshmm --qshmm $BIO/pbsim_data/QSHMM-RSII.model --depth $DEP --prefix $OUT_DIR/lane$LN --genome $FASTA_IN
rm $OUT_DIR/lane$LN*.ref || true
rm $OUT_DIR/lane$LN*.maf.gz || true
cat $OUT_DIR/lane$LN*.gz > $OUT_DIR/pbsim/lane$LN.fq.gz
rm $OUT_DIR/lane$LN*.gz || true

LN=$(awk 'BEGIN{print '$LN'+1;}')
$BIO/pbsim --strategy wgs --method qshmm --qshmm $BIO/pbsim_data/QSHMM-RSII.model --depth $DEP --prefix $OUT_DIR/lane$LN --genome $FASTA_IN
rm $OUT_DIR/lane$LN*.ref || true
rm $OUT_DIR/lane$LN*.maf.gz || true
cat $OUT_DIR/lane$LN*.gz > $OUT_DIR/pbsim/lane$LN.fq.gz
rm $OUT_DIR/lane$LN*.gz || true
