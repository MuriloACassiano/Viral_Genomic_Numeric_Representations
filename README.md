# Numeric representations and alignment-free techniques for sequence clustering as tools for unsupervised grouping:
## Applications in coronavirus phylogeny and Brazilian lineages of SARS-CoV-2

In this study, we aimed to evaluate the effectiveness of several K-mer-based numerical representations, including Triplet Frequency, K-mer Natural Vector, Fast Vector, and Magnus Genomic Representation, in analyzing coronaviridae sequences (both cured and published) and approximately 86 thousand sequenced genomes in Brazil, obtained from the GISAID EpiCov database. Our analyses included comparing measures that summarize structural characteristics, cofeneic correlations, and distances between trees constructed using Euclidean distances of the numerical representations and the tree built from multiple sequence alignment and subsequent phylogenetic estimation by maximum likelihood. We also assessed the ability of each numerical representation to carry biological information known from the sequences, such as taxonomic group or viral lineage, through dimensionality reduction techniques.

________________________________________________________________________________

## Getting Started

We tested our scripts in a desktop PC (16GB RAM, Intel® Core™ i5-10400F CPU @ 2.90GHz × 12, Ubuntu 22.04.2 LTS) and in a cloud-based machine (128 GB RAM, Intel® Xeon® CPU E7- 2870  @ 2.40GHz x 32, Ubuntu 20.04.5 LTS).

We used python, R and shell scripting to run our analysis.

________________________________________________________________________________

### Prerequisites

#### Python

These are the python libraries required to run the codes:

numpy==1.23.5

pandas==1.4.2

pyfastx==0.8.4

ray==1.12.1

regex==2022.4.24

scikit-learn==1.1.2

umap-learn==0.5.3

It is required Python 3.8+, pip correclty installed in your machine.
See (https://pip.pypa.io/en/stable/installation/)


To install:
```
pip install numpy==1.23.5 pandas==1.4.2 pyfastx==0.8.4 ray==1.12.1 regex==2022.4.24 scikit-learn==1.1.2 umap-learn==0.5.3
#or
pip3 install numpy==1.23.5 pandas==1.4.2 pyfastx==0.8.4 ray==1.12.1 regex==2022.4.24 scikit-learn==1.1.2 umap-learn==0.5.3
```


#### R

We suggest the usage of rpy2, as it simplifies the runnning of python and R at the same notebook.

##### One must properly install rpy2 before continue:
https://rpy2.github.io/doc/latest/html/overview.html#install-installation

At the end of installation, one must be able to run the cell bellow:
```
%load_ext rpy2.ipython
```

These are the R libraries required to run the codes:

###### Install and load the Bioconductor manager
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
```

###### Install the required Bioconductor packages
```
BiocManager::install(c("phyloTop", "treeio", "castor", "ggtree"))
```

###### Install the required CRAN packages
```
install.packages(c("tidyverse", "ape", "tidytree", "ggpubr"))
install.packages(c("ggExtra", "ggpmisc", "reshape2"))
install.packages(c("phylotools", "cowplot", "dendextend"))
```

#### Third-Party Software

### * MAFFT: a multiple sequence alignment program for unix-like operating systems
  * [Homepage](https://mafft.cbrc.jp/alignment/software/)
  * [LINUX Installation](https://mafft.cbrc.jp/alignment/software/linux.html)
  * [MAC-OS Installation](https://mafft.cbrc.jp/alignment/software/macosx.html)
  * [WSL windows Installation](https://mafft.cbrc.jp/alignment/software/windows.html)

### * trimAL: a tool for the automated removal of spurious sequences or poorly aligned regions from a multiple sequence alignment
  * [Homepage](http://trimal.cgenomics.org/trimal)
  * [Installation](http://trimal.cgenomics.org/downloads)

### * IQ-TREE: Efficient software for phylogenomic inference
  * [Homepage](http://www.iqtree.org/)
  * [Documentation](http://www.iqtree.org/doc/)

### * Augur: A bioinformatic toolkit for phylogenetic analysis
  * [Homepage](https://docs.nextstrain.org/projects/augur/en/stable/index.html)
  * [Installation](https://docs.nextstrain.org/projects/augur/en/stable/installation/installation.html)

### * Nextclade: Clade assignment, mutation calling, and sequence quality checks
  * [Homepage](https://clades.nextstrain.org/)
  * [Download](https://github.com/nextstrain/nextclade/releases)

### * phylonium: fast and accurate estimation of evolutionary distances
  * [Homepage](https://github.com/EvolBioInf/phylonium)
  * [Manual](https://github.com/EvolBioInf/phylonium/blob/master/documentation/manual.pdf)

________________________________________________________________________________

### Numeric encoding

Python Script [numeric_encoding.py](https://github.com/MuriloACassiano/Viral_Genomic_Numeric_Representations/blob/main/numeric_encoding.py)

As input, this script takes a fasta file, containing whole genome sequences in ACGTN characters (not aligned sequences).
As output, this software provide a csv file, with the first column being the fasta header from the original fasta file.

usage:

numeric_encoding.py [-h] [--fasta FASTA] [--repr REPRESENTATION] [--thread THREAD]

optional arguments:

  -h, --help            show this help message and exit
  
  --fasta FASTA         location/name of the input fasta file
  
  --repr REPRESENTATION
                        Choosen numeric representation:
                        
                        triplet - Triplet dictionary
                        
                        magnus - Magnus Representation
                        
                        Fast - Fast Vector
                        
                        4-mer - 4-mer natural vector
                        
                        6-mer - 6-mer natural vector
                        
                        c-4-mer - Cumulative 4-mer natural vector
                        
                        c-6-mer - Cumulative 6-mer natural vector
  --thread THREAD       Number of threads to use
  

Example:

```
python3 numeric_encoding.py --fasta my.sequences.fasta --repr "Magnus" --thread 12
python3 numeric_encoding.py --fasta my.sequences.fasta --repr "Fast" --thread 18
```

Or executing uncoupled from terminal (unix based systems):

```
nohup python3 numeric_encoding.py --fasta my.sequences.fasta --repr "4-mer" --thread 32 &
nohup python3 numeric_encoding.py --fasta my.sequences.fasta --repr "6-mer" --thread 32 &
```
________________________________________________________________________________

### Orthocoronaviridae

Notebook with Python, Shell and R Scripts [coronavirinae.ipynb](https://github.com/MuriloACassiano/Viral_Genomic_Numeric_Representations/blob/main/coronavirinae.ipynb)

This notebook constains a detailed description of the analysis of 69 genomes of alpha, beta, delta and gammacoronaviruses.

________________________________________________________________________________

### Brazilian SARS-CoV-2 Genomes

Notebook with Python, Shell and R Scripts [BR_SARS_CoV_2.ipynb](https://github.com/MuriloACassiano/Viral_Genomic_Numeric_Representations/blob/main/BR_SARS_CoV_2.ipynb)

This notebook constains a detailed description of the analysis of SARS-CoV-2 genomes, from the ~86K sequences subsampled to 3K sequences in order to make the comparisons possible.

________________________________________________________________________________

### Brazilian SARS-CoV-2 Genomes

Notebook with Shell and R Scripts [phylonium_output.ipynb](https://github.com/MuriloACassiano/Viral_Genomic_Numeric_Representations/blob/main/phylonium_output.ipynb)

This notebook constains a detailed description of the analysis of ~3K SARS-CoV-2 genomes using a program for estimating the evolutionary distances between closely related genomes.

According it's developers and recend independent benchmarks, phylonium is much faster than alignment based approaches for phylogeny reconstruction and more accurate than similar alignment-free methods.
________________________________________________________________________________

## Authors

* **Murilo H. Anzolini Cassiano** - *Conceptualization and execution* - [Murilo Cassiano](https://github.com/MuriloACassiano)

* **Daniel Macedo de Melo Jorge** - *Co-orientation (Bioinformatics)*

* **Eurico Arruda** - *Orientation (Virology)*

________________________________________________________________________________

## License

This project is licensed under the GNU General Public License (GPL).

________________________________________________________________________________

## Acknowledgments

* The present work was carried out with the support of the Coordination for the Improvement of Higher Education Personnel - Brazil (CAPES).
