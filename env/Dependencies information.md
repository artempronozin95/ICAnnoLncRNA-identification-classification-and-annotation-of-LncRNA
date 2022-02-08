# programs.yaml 
## [EMBOSS](http://emboss.open-bio.org/html/use/ch02s07.html)
"The European Molecular Biology Open Software Suite". EMBOSS is a free Open Source software analysis package specially developed for the needs of the molecular biology (e.g. EMBnet) user community. The software automatically copes with data in a variety of formats and even allows transparent retrieval of sequence data from the web. Also, as extensive libraries are provided with the package, it is a platform to allow other scientists to develop and release software in true open source spirit. EMBOSS also integrates a range of currently available packages and tools for sequence analysis into a seamless whole. EMBOSS breaks the historical trend towards commercial software packages.

`conda install -c bioconda emboss`
## [GMAP](https://academic.oup.com/bioinformatics/article/21/9/1859/409207)
The program maps and aligns a single sequence with minimal startup time and memory requirements, and provides fast batch processing of large sequence sets. The program generates accurate gene structures, even in the presence of substantial polymorphisms and sequence errors, without using probabilistic splice site models. Methodology underlying the program includes a minimal sampling strategy for genomic mapping, oligomer chaining for approximate alignment, sandwich DP for splice site detection, and microexon identification with statistical significance testing.

`conda install -c bioconda gmap`
## [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)
gffcompare can be used to compare, merge, annotate and estimate accuracy of one or more GFF files (the “query” files), when compared with a reference annotation (also provided as GFF).

`conda install -c bioconda gffcompare`
## [bedops](https://bedops.readthedocs.io/en/latest/)
BEDOPS is an open-source command-line toolkit that performs highly efficient and scalable Boolean and other set operations, statistical calculations, archiving, conversion and other management of genomic data of arbitrary scale. Tasks can be easily split by chromosome for distributing whole-genome analyses across a computational cluster.

`conda install -c bioconda bedops`
## [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
The Snakemake workflow management system is a tool to create reproducible and scalable data analyses. Workflows are described via a human readable, Python based language. They can be seamlessly scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition. Finally, Snakemake workflows can entail a description of required software, which will be automatically deployed to any execution environment.

`conda install -c bioconda snakemake`
### [seqtk](https://docs.csc.fi/apps/seqtk/)
Seqtk is a fast and lightweight tool for processing sequences in the FASTA or FASTQ format. It seamlessly parses both FASTA and FASTQ files which can also be optionally compressed by gzip.

`conda install -c bioconda seqtk`
## Python
### [pandas](https://pandas.pydata.org)
pandas aims to be the fundamental high-level building block for doing practical, real world data analysis in Python. Additionally, it has the broader goal of becoming the most powerful and flexible open source data analysis / manipulation tool available in any language.

`conda install -c anaconda pandas`
### [more-itertools](https://pypi.org/project/more-itertools/)
Python’s itertools library is a gem - you can compose elegant solutions for a variety of problems with the functions it provides. In more-itertools we collect additional building blocks, recipes, and routines for working with Python iterables.

`conda install -c anaconda more-itertools`
### [Biopython](https://biopython.org)
Biopython is a set of freely available tools for biological computation written in Python by an international team of developers.

`conda install -c conda-forge biopython`
### [matplotlib](https://matplotlib.org)
Matplotlib is a comprehensive library for creating static, animated, and interactive visualizations in Python. Matplotlib makes easy things easy and hard things possible.

`conda install -c conda-forge matplotlib`
### [gtfparse](https://pypi.org/project/gtfparse/)
Parsing tools for GTF (gene transfer format) files.

`conda install -c bioconda gtfparse`
### [scikit](https://scikit-learn.org/stable/index.html)
Simple and efficient tools for predictive data analysis

`conda install -c anaconda scikit-learn`
### [seaborn](https://seaborn.pydata.org/)
Seaborn is a Python data visualization library based on matplotlib. It provides a high-level interface for drawing attractive and informative statistical graphics.

`conda install -c anaconda seaborn`
# cpc.yaml
## [CPC2](https://github.com/biocoder/CPC2/blob/master/README.md)
Coding Potential Calculator (CPC), can assess the protein-coding potential of a transcript based on six biologically meaningful sequence features. Tenfold cross-validation on the training dataset and further testing on several large datasets showed that CPC can discriminate coding from noncoding transcripts with high accuracy.
# lncfinder.yaml
## [LncFinder](https://cran.r-project.org/web/packages/LncFinder/index.html)
Long non-coding RNAs identification and analysis. Default models are trained with human, mouse and wheat datasets by employing SVM. Features are based on intrinsic composition of sequence, EIIP value (electron-ion interaction pseudopotential), and secondary structure. 

`conda install -c heyarti r-lncfinder`
## [libgfortran](https://github.com/gcc-mirror/gcc/tree/master/libgfortran)
The libgfortran contains a Fortran shared library, needed to run Fortran dynamically-linked programs, required for the Vertica LAPACK library addition.

`conda install -c conda-forge libgfortran`
# alignment.yaml
## [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
BLAST finds regions of similarity between biological sequences. The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance.

`conda install -c bioconda blast`
## [diamond](https://github.com/bbuchfink/diamond)
DIAMOND is a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data.

`conda install -c bioconda diamond`
