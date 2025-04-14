![GitHub All Releases](https://img.shields.io/github/downloads/timothyjamesbecker/somacx/total.svg) [![DOI](https://zenodo.org/badge/185091540.svg)](https://zenodo.org/badge/latestdoi/185091540) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Alt text](images/somacx_logo.png?raw=true "somacx")
## A complex generative genome modeling framework
Copyright (C) 2020-2025 Timothy James Becker
![Alt text](images/clone_tree_inv.png?raw=true "somacx") <br>

## Installation

### Pip Installation
Requirements: python 3.10+, cython 0.29+, pysam 0.15+, numpy 1.18+<br>
Instructions: use the pip command line as shown below from the current release
```bash
python -m pip install https://github.com/timothyjamesbecker/somacx/releases/download/0.1.3/somacx-0.1.3.tar.gz
```

### Docker Installation (Dockerfile)
Requirements: docker version 28+
Instructions: Clone the repo, change directory into it and then build the image
```bash
git clone https://github.com/timothyjamesbecker/somacx.git
cd somacx
docker build -t somacx ./
```

### Docker Image Usage (no install method)
Requirments: docker version 28+
Instructions: use the docker command line to download the prebuilt image
```bash
docker pull timothyjamesbecker/somacx
```

## Usage

### Docker Usage
```bash
docker run -v /media/yourdata:/data -it timothyjamesbecker/somacx \
generator.py -r /data/hg38.fa -j hg38 -C chr11,chr22,chrX -o /data/out_dir --model 0.25,0.75 --cov 2,6 
```

### Normal Usage
Pip will add the bin directly of this repo into the environmental PATH so you can directly call generator.py:
```bash
generator.py -r hg38.fa -j hg38 -C chr11,chr22,chrX -o out_dir --model 0.25,0.75 --cov 2,6
```
<b>Argment Details:</b><br>

<b>-r</b> the reference fasta file such as GRCh38/Hg38, hg19/human_g1k_v37decoy, the user will need to download this file and then pass its full path in this argument<br>

<b>-j</b> the full JSON control file which includes SV type and sizes rates (including separate linked versus unlinked SNV/MNVs), 
A gene map such as refGene.txt.gz, encoded weighted gene list or mask files which we absract to the WCU data structure that governs 
the gain and loss probabilities of genomic regions and finally an optional user specified clonal tree topology.
Use presets values: hg38 for Hg38/GRCh38 coordinates (no alternate contigs) or hv37d for the 1000 Genomes Phase 3 reference 
fasta used in the FusorSV paper.
Full json presents can be generated using the included tools.

<b>-C</b> the chromosomes from the given reference fasta that will be output (good for testing or producing a small dataset)

<b>-o</b> the output directory where the germline genome file (~2X the size of the reference fasta) as well as the somatic genome file 
(2/cov) will be stored as well as the matching truth VCF files for each. The somatic VCF file includes the clonal tree topology 
encoded using the PEDIGREE entries in the header.

#### Options
all options can now be given uniform random ranges (these values will override JSON presets)

<b>--model</b> a parametric control that controls the tree topology 0.0 will simulate Cancer Stem Cell (CSC) and 1.0 will 
yeild a balanced subclonal topology

<b>--branch</b> the chance of each node to branch every cycle

<b>--decay</b> the chance of each node to die off every cycle

<b>--cov</b> the desired sequencing coverage so that the user can control the approximate noise level threshold for having 
alleles of the lowest frequency present

## Detailed Examples (step 1 generate FASTA file)
### Download a reference FASTA if you don't have one
```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
unzip hg38.fa.gz
```
### Generate one germline genome (chr1,chr2,...,chrX,chrY) using hg38
```bash
generator.py -r hg38.fa -j hg38 --cov 20,30 --o out_dir --germ
```

### Generate one germline and one somatic sub-clonal genome 
```bash
generator.py -r hg38.fa -j hg38 -o out_dir --model 0.0 --branch 0.5 --cov 30
```

### Generate one germline and one somatic cancer-stem cell genome
```bash
generator.py -r hg38.fa -j hg38 -o out_dir --model 1.0 --branch 0.5 --cov 30
```

### Generate 10 somatic genomes between cancer stem cell and sub clonal
```bash
for i in $(seq 1 10);
do
  generator.py -r hg38.fa -j hg38 -o out_dir --model 0.0,1.0 --branch 0.5 --cov 30  
done
```

## Detailed Examples (step 2 use somaCX FASTA with a read simulator to make FASTQ reads)
Next you take the FASTA genome file that is generated in the somaCX framework and simulate sequence reads. 
FASTA is the universal input for a wide variety of read simulators.  
We have also included [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) 
and [PBSIM3](https://github.com/yukiteruono/pbsim3) examples scripts that simulate two lanes for a 
normal tissue sequence experiment or 4 lanes (at 2X coverage) for tumor tissue sequencing. 
### Simulate Illumina Paired-end normal genome generated from somaCX

## Detail Examples (step3 use FASTQ reads to process with minimap2, samtools, sambamba to make BAM)
These examples will also require a sequence alignment tool such as [minimap2](), read utilities [samtools]() and [sambamba]().
### align FASTQ reads to genome

## Advanced JSON Usage
The full distribution for a somatic genome generation will contain germline distribution parameters. Together the default
human JSON file called full.hg38.json.gz contains the following data slots (keys in JSON nomenclature):<br><br>
#### full.hg38.json.gz (all distributions packed up into one file)
gene_map `` - a JSON object used for gene symbol coordinates (has seq, gene and wcu keys for query)``<br><br>
g_var_map``- a JSON object used for linked and unlinked germline multi-nucleotide and structural variation 
(type/size/frequency)``<br><br>
s_var_map -`` a JSON object used for linke and unlinked somatic germlin multi-nucleotide and structural variation
(type/size/frequency)``<br><br>
clone_tree ``- a JSON object that defines the structure of the simulated clone tree in the somaic genome generation 
phase``<br><br>
somatic_aneuploidy ``- the multinomial distribution for each chromosome to generate aneuploidy (abnormal number
of chromosomes)``<br><br>
wcus (weighted class units below) ``- JSON objects with chromosomes as keys with arrays for each chromosome containing:
[start, end, length, {class_id:genomic_id}]. class_id are the g1kp3, onco, svmask, mitcp, mmej, nhej, apot identifiers and 
genomic_id are a gene SYMBOL or other unique genomic region id``<br><br>
germline_[genelist]_loss_wcu(s)``- g1kp3, onco, svmask, mitcp, mmej, nhej, apot gene list based ranges that
will become the positional distributions for germline loss generation (for types like DEL, INV, TRA, etc)``<br><br>
germline_[genelist]_gain_wcu(s)``- g1kp3, onco, svmask, mitcp, mmej, nhej, apot gene list based ranges that
will become the positional distributions for germline gain generation (for types like DUP, etc)``<br><br>
somatic_[genelist]_loss_wcu(s)``- g1kp3, onco, svmask, mitcp, mmej, nhej, apot gene list based ranges that
will become the positional distributions for somatic loss generation (for types like DEL, INV, TRA, etc)``<br><br>
somatic_[genelist]_gain_wcu(s)``- g1kp3, onco, svmask, mitcp, mmej, nhej, apot gene list based ranges that
will become the positional distributions for somatic gain generation (for types like DUP, etc)``<br><br>

In order to design distribution parameters for new genomes, users would need to build a new JSON file and set values as above. We have some development plans for building JSON presets for mm genomes.

More details on the individual JSON objects: <br>

#### gene_map
``allows query by sequence name or by gene symbol. seq: will be each chromosome followed by an array of all genes and start,end
coordinates so that all genes for a given chromosome can easily be integrated into gain or loss distribtions. gene: are the offciale gene symbols
such as DDX11L1, TP53, etc and will resolve to one or more genomic coordinates.``<br>

#### g_var_map, s_var_map
`` keys denote the layer of the variation where 1 will be a MNV before layer 2 (SV layer) these SNV will be linked to the SVs
that are intersecting, (such as DUP that would also duplicate the SNV) wheras layer 3 MNVs are applied after the SV events which
will make them unlinked (such as a DUP with one copy containing a MNV but the other doesn't). In layer 2 (SV layer) there are
parameters that govern the distributions of each type detailed further:``

DUP ``l:n - distribution of the lengths (chance of a SV of lenght l per bp of genome). CN - number of copies allowed (2,..,8)
default. CNP - probabilty of each copy number. TYPE - DUP subtypes: TANDEM, DISPERSED, INVERTED. TP - probabilty of each 
type (type distribution)``<br>
TRA ``l:n - distribution of the lengths (chance of a SV of lenght l per bp of genome). TYPE - subtype for traslocations:
DEL or INS which are balanced or unbalanced. TP - is the type distribution for TRA``<br>
INV ``l:n - distribution of the lengths (chance of a SV of lenght l per bp of genome). TYPE - is the subtype: PERFECT or
COMPLEX. PERFECT have very simple signatures (perfect DNA repair) and reverse complement. COMPLES on the other hand
have DEL and DUP posisble on the break ends as well as nove INS inside the INV SV call region. s:p - sets the complex INV
subtype DEL/DUP/INS length distribution (50bp,100bp,etc)``<br>
DEL ``l:n - distribution of the lengths (chance of a SV of lenght l per bp of genome).``<br>
INS  ``l:n - distribution of the lengths (chance of a SV of lenght l per bp of genome).``

#### clone_tree
``clone tree keeps the same of the somatic evolution in tree form but is used practically to spool the FASTA files and
generate each clone. The root key keeps track of the genome name (sample name) while the tree key keeps the nested JSON tree 
structure for each clone which will be nodes. Alive is an array that determines if a specific clone will be generated
which is useful for simulation of ancesters that dropped out from selective pressure. freq is an array that keeps track of
the number of ancestors which is used to allocate variation (used for allele frequency generation)``<br>

## Web Browser JSON Editor
``Visual web broswer JSON editor is in development``

## Python API
``gene_effects, read_utils, variant_utils and sim python modules offer much greater flexibility in simulating mammalian genomes
``

























