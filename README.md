![Alt text](images/somacx_logo.png?raw=true "somacx") <br>
## A pathway based somatic genome generation framework
(c) 2019 Timothy Becker

![Alt text](images/clone_tree.png?raw=true "somacx") <br>


## Requirements
python 2.7.11+, cython 0.28+, pysam 0.15+, numpy 1.16+<br>

## Installation
```bash
pip install https://
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
## Advanced Usage
``JSON file format API coming soon...``

## Python API
``gene_effects, read_utils, variant_utils and sim python modules offer much greater flexibility in simulating beyond somatic mammalian genomes
``
