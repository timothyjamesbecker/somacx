 ![GitHub All Releases](https://img.shields.io/github/downloads/timothyjamesbecker/somacx/total.svg) [![DOI](https://zenodo.org/badge/185091540.svg)](https://zenodo.org/badge/latestdoi/185091540) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Alt text](images/somacx_logo.png?raw=true "somacx") <br>
## A pathway based somatic genome generation framework
Copyright (C) 2020 Timothy Becker
![Alt text](images/clone_tree.png?raw=true "somacx") <br>


## Requirements
<b>(python 3.6+ is now the prefered platform as of release 0.1.2)</b><br>
python 3.6.8+, cython 0.29+, pysam 0.15+, numpy 1.18+ or<br>
python 2.7.11+, cython 0.29+, pysam 0.15+, numpy 1.18+
## PIP Installation
```bash
python3 -m pip install https://github.com/timothyjamesbecker/somacx/releases/download/0.1.2/somacx-0.1.2.tar.gz
```

## Docker Installation
This repo has a Dockerfile that will build somacx for you or you can easily grab a pre-built image for testing or cloud use:
```bash
docker pull timothyjamesbecker/somacx
```

## Basic Usage
```bash
python3 generator.py -r ref_fasta -j hg38 -C 11,22,X -o out_dir --model 0.25,0.75 --cov 2,6
```
<b>-r</b> the reference fasta file such as GRCh38/Hg38, hg19/human_g1k_v37decoy<br>

<b>-j</b> the full JSON control file which includes SV type and sizes rates (including separate linked versus unlinked SNV/MNVs), A gene map such as refGene.txt.gz, encoded weighted gene list or mask files which we absract to the WCU data structure that governs the gain and loss probabilities of genomic regions and finally an optional user specified clonal tree topology.
Use presets values: hg38 for Hg38/GRCh38 coordinates (no alternate contigs) or hv37d for the 1000 Genomes Phase 3 reference fasta used in the FusorSV paper.
Full json presents can be generated using the included tools.
<b>-C</b> the chroms from the given reference fasta that will be output (good for testing or producing a small dataset)
<b>-o</b> the output directory where the germline genome file (~2X the size of the reference fasta) as well as the somatic genome file (2/cov) will be stored as well as the matching truth VCF files for each. The somatic VCF file includes the clonal tree topology encoded using the PEDIGREE entries in the header.

#### options
<b>all options can now be given uniform random ranges (these values will override JSON presets)</b>

<b>--model</b> a parametric control that controls the tree topology 0.0 will simulate Cancer Stem Cell (CSC) and 1.0 will yeild a balanced subclonal topology

<b>--branch</b> the chance of each node to branch every cycle

<b>--decay</b> the chance of each node to die off every cycle

<b>--cov</b> the desired sequencing coverage so that the user can control the approximate noise level threshold for having alleles of the lowest frequency present

## Web Browser JSON Editor
``Visual web broswer JSON editor coming soon``

## Advanced JSON Usage
``JSON file format API coming soon...``

## Python API
``gene_effects, read_utils, variant_utils and sim python modules offer much greater flexibility in simulating mammalian genomes
``
