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
install.packages(c("phylotools", "dende"))
```

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

Nootbook with Python, Shell and R Scripts [coronavirinae.ipynb](https://github.com/MuriloACassiano/Viral_Genomic_Numeric_Representations/blob/main/coronavirinae.ipynb)

This notebook constains a detailed description of the analysis of 69 genomes of alpha, beta, delta and gammacoronaviruses.

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
