[![Build Status](https://api.travis-ci.org/timothyjamesbecker/somacx.svg)](https://travis-ci.com/timothyjamesbecker/somacx) ![GitHub All Releases](https://img.shields.io/github/downloads/timothyjamesbecker/somacx/total.svg) [![DOI](https://zenodo.org/badge/185091540.svg)](https://zenodo.org/badge/latestdoi/185091540) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Alt text](images/somacx_logo.png?raw=true "somacx") <br>
## A pathway based somatic genome generation framework
Copyright (C) 2020 Timothy Becker
![Alt text](images/clone_tree.png?raw=true "somacx") <br>


## Requirements
python 2.7.11+, cython 0.29+, pysam 0.15+, numpy 1.18+ or<br>
python 3.6.8+, cython 0.29+, pysam 0.15+, numpy 1.18+

## PIP Installation
```bash
pip install https://github.com/timothyjamesbecker/somacx/releases/download/0.1.1/somacx-0.1.1.tar.gz
```

## Docker Installation
Our repo has a Dockerfile that will build somacx for you or you can easily grab a pre-built image for testing or cloud use:
```bash
docker pull timothyjamesbecker/somacx
```

## Basic Usage
```bash
generator.py -r ref_fasta -j DEFAULT -c 11,22,X -C 11,22,X -o out_dir --cov 2
```
<b>-r</b> the reference fasta file such as hg19, Hg38, human_g1k_v37decoy, etc<br>

<b>-j</b> the full JSON control file which includes SV type and sizes rates (including separate linked versus unlinked SNV/MNVs), A gene map such as refGene.txt.gz, encoded weighted gene list or mask files which we absract to the WCU data structure that governs the gain and loss probabilities of genomic regions and finally an optional user specified clonal tree topology.

<b>-c</b> the chroms from the given reference fasta that should have Variation modeled (default is all)

<b>-C</b> the chroms from the given reference fasta that will be output (good for testing or producing a small dataset)

<b>-o</b> the output directory where the germline genome file (~2X the size of the reference fasta) as well as the somatic genome file (2/cov) will be stored as well as the matching truth VCF files for each. The somatic VCF file includes the clonal tree topology encoded using the PEDIGREE entries in the header.

#### options
<b>--model</b> a parametric control that controls the tree topology 0.0 will simulate Cancer Stem Cell (CSC) and 1.0 will yeild a balanced subclonal topology

<b>--branch</b> the chance of each node to branch every cycle

<b>--decay</b> the chance of each node to die off every cycle

<b>--cov</b> the desired sequencing coverage so that the user can control the approximate noise level threshold for having alleles of the lowest frequency present

## Advanced Usage
``JSON file format API coming soon...``

## Python API
``gene_effects, read_utils, variant_utils and sim python modules offer much greater flexibility in simulating beyond somatic mammalian genomes
``
