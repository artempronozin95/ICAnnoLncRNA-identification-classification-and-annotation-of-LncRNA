## Table of contents
* [Introduction](#introduction)
* [Installation](#installation)
* [Input](#input)
  * [Configuration](#configuration-file)
  * [Folders](#folders)
* [Running the pipeline](#work-start)
* [Models](#models)
* [Structure information](#lncrna-structure-information)
* [Known LncRNA](#known-lncrna-for-database)
* [Tissue analysis](#tissue-analysis)
* [Output](#output)
## Introduction
Long non-coding RNAs (lncRNAs) typically defined as transcripts of more than 200 nucleotides in length and without any protein coding potential. These RNAs are involved in important plant development processes such as phosphate homeostasis, flowering, photomorphogenesis and stress response. However, their structural and functional properties are not clear. Information about lncRNA sequences and their expression typically obtained from RNA-seq experiments. To process these data relevant tools for automated lncRNA are required. 
Here we propose ICAnnoLncRNA, pipeline for automatic identification, classification and annotation of plant lncRNAs. ICAnnoLncRNA, work with RNA-seq data and it was tested on 877 transcriptome libraries of *Zea mays*.
The pipeline includes the following steps: 
#### 1. Pre-processing.
+ Gmap index building.
+ LncFinder model building.
#### 2. LncRNA filtering.
+ LncFinder predicting.
+ CPC2 predicting.
+ ORF detection.
+ Transmembrane domains detection.
+ GMAP alignment on reference genome.
#### 3. LncRNA annotation.
+ Classification gffcompare.
+ Blastn alignment.

The pipeline is implemented using the workflow management system [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), which provides ability to platform-independent installation and execution of the software.

## Schematic diagram
![Test Image 1](https://github.com/artempronozin95/ICAnnoLncRNA---identification-classification-and-annotation-of-LncRNA/blob/main/image/Pipeline.png)

## Installation 
```
1. wget https://github.com/artempronozin95/ICAnnoLncRNA---identification-classification-and-annotation-of-LncRNA/archive/refs/heads/main.zip
2. unzip main.zip
3. cd ./ICAnnoLncRNA---identification-classification-and-annotation-of-LncRNA-main
4. conda env create --file env/programs.yaml
5. conda activate ICAnnoLncRNA
```
After these steps all necessary packages are installed. If you need update packages (**not recommended**), change the version of  packages after “=” (example - snakemake=4.0.1 -> snakemake=6.0.0), then “conda env update --file ./programs.yaml”. All necessary packages will be updated. 
More detailed description provided [here](https://github.com/artempronozin95/ICAnnoLncRNA---identification-classification-and-annotation-of-LncRNA/blob/main/env/Dependencies%20information.md)

It is necessary to install (download program in pipeline folder): 
+ [CPC2](https://github.com/biocoder/CPC2)
```
1. cd ./ICAnnoLncRNA---identification-classification-and-annotation-of-LncRNA-main
2. wget https://github.com/biocoder/CPC2/archive/refs/heads/master.zip
3. unzip master.zip
4. mv CPC2-master CPC2
```
+ [TMHMM](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0)
```
1. click on Download.
2. register.
3. check up your mail.
4. cd ./ICAnnoLncRNA---identification-classification-and-annotation-of-LncRNA-main
5. wget the link from mail.
6. tar -xf tmhmm-2.0c.Linux.tar.gz
7. mv tmhmm-2.0c tmhmm
8. cd tmhmm/bin
9. open tmhmm.
10. change $opt_plot = 1 on $opt_plot = 0 and $opt_short = 0 on $opt_short = 1.
11. save.
```
## Input
### Genome sequence
1. Known LncRNA transcripts of the species in `FASTA` format. File with lncRNA from the public database of lncRNA. For example [Ensembl](https://www.ensembl.org/index.html), [EVLncRNAs v2](https://www.sdklab-biophysics-dzu.net/EVLncRNAs2/).
2. Known protein coding sequences of the species in `FASTA` format. File with protein coding sequences from a public database. For example [Ensembl](https://www.ensembl.org/index.html).
3. RNA-seq transcripts of the species in `FASTA` format. This file contains transcripts obtained after de novo assembly of transcriptome. For example file that we provide, [30k_zey.fasta](https://github.com/artempronozin95/ICAnnoLncRNA---identification-classification-and-annotation-of-LncRNA/blob/main/data/input/30k_zey.zip), obtained by [Trinity](https://github.com/trinityrnaseq/trinityrnaseq) assembly.
All three sets are necessary for build 
### Reference genome
1. Reference genome of the species in `FASTA` format.
2. Annotation of the species in `GFF/GTF` format.

### Configuration file
Input all necessary files into configuration file “config.yaml”:
+ `lnc:` - known LncRNA sequences of studied organisms in `FASTA` format. (Example: `lnc: "data/input/lincrna.fasta"`)
+ `cds:` - known CDS (coding) sequences of studied organisms in `FASTA` format. (Example: `cds: "data/input/7000_len_mRNA.fasta"`)
+ `model:` - pipeline check if model (lncFinder) for this organism already exist and use it, if it's not, pipeline will create a new one (in output folder). Model that exist represented in [Models](#models). (Example: `model: "Zea_mays"`)
+ `sequence:` - sequences that need to study in `FASTA` format (Example: `sequence: "data/input/Zea.fasta"`)
+ `structure:` - need to choose whether to use secondary structure in model building or not. Choose between `DNA` or `SS` (secondary structure). (Example: `structure: "DNA"`)
+ gmap:
  + `gff_reference:` - reference genome annotation in `GFF` format. (Example: `gff_reference: "data/reference/Zea_mays.AGPv4.40.gff"`)
  + `reference:` - reference genome in `FASTA` format. (Example: `reference: "data/reference/Zea_mays.AGPv4.dna.toplevel.fa"`)
  + `out:` - output file. (Example: `out: "data/output/alignm.gff"`)
+ `gff:` - gff file of reference genome. (Example: `gff: "data/reference/Zea_mays.AGPv4.40.gff"`)
+ diamond: - alignment of lncRNA on protein database (recommended to run after main pipeline)
  + `option:` - "on" (activate alignment step) or "off" (deactivate alignment step) (Example:  `option: "off"`)
  + `database:` - protein database in `FASTA` format. (Example: `database: "data/reference/data_base/Arabidopsis_thaliana.TAIR10.pep.all.fa"`)
  + `query:` - lncRNA transcripts in `FASTA` format. (Example: `query: "data/output/unknown.fasta"`)
  + `out:` - output file in `outfmt6` format. (Example: `out: "data/output/diamond.outfmt6"`)
+ tissue:
  + `organism:` - choose between organisms in [Tissue analysis](#tissue-analysis) or input your organism. (Example: `organism: "ZM"`)
### Folders
#### data/input
Contain: 
1. Known LncRNA transcripts of the species in `FASTA` format.
2. Known protein coding sequences of the species in `FASTA` format.
3. RNA-seq transcripts of the species in `FASTA` format.
#### data/reference
Contain:
+ reference genome of studied organism in `FASTA` format.
+ annotation of studied organism in `GFF/GTF` format.
+ `data_index` folder with library of known lncRNA from databases.
+ `intron_data` folder with structure information of organisms.
+ `models` folder with model for lncFinder.
## Work start
  #### 1. `snakemake -j 2`
  `-j` or  --cores -  Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores.
  #### 2. `snakemake -nr` 
  `-n` - Do not execute anything, and display what would be done. If you have a very large workflow, use –dry-run –quiet to just print a summary of the DAG of jobs.
  `-r` - Print the reason for each executed rule.
  #### 3. `snakemake --use-conda`
  `--use-conda` - Use additional conda environment.
  #### 4. Recommended run: 
  `snakemake -j 2 --use-conda`
## Models
+ *Zea_mays*
+ *Arabidopsis_thaliana*

## LncRNA structure information
+ *Zea_mays*
+ *Arabidopsis_thaliana*
+ *Oryza_sativa*

## Known LncRNA for database
Folder `data/reference/data_index` contain lncRNA library that build on base of 5 lncRNA databases (PNRD, CANTATAdb, GREENC, PlncDB, EVLncRNA)

## Tissue analysis
File `tissue/SRX_all_org.tsv`, contains information about transcript experiment libraries with respect to specific tissue. In config.yaml need to choose between organisms that presented in SRX file:
+ HV - *Hordeum vulgare*
+ OS - *Oryza sativa*
+ SL - *Solanum lycopersicum*
+ ST - *Solanum tuberosum*
+ ZM - *Zea mays*
Or input your organism in the same format file.
## Output
A typical structure of `Output` is follows:

    output/
        ├── alignm.bed
        ├── alignm_filter.gff
        ├── anti.png                                                                                                                         
        ├── blast.outfmt6                                                                                                                    
        ├── class                                                                                                                                                
        │   ├── class_org.png                                                                                                                                    
        ├── classes.png                                                                                                                      
        ├── Coding.fasta                                                                                                                     
        ├── cpc.txt
        ├── exon_size.png
        ├── filter_alignm.bed
        ├── gffcmp.alignm_filter.gff.tmap
        ├── gffcmp.loci
        ├── intron_size.png
        ├── itron_coordin.tsv
        ├── lncFinder.csv
        ├── new_lncrna.fasta
        ├── Noncoding.fasta
        ├── Noncoding_trans_out.fasta
        ├── not_trans.csv
        ├── number_of_exon.png
        ├── number_of_lncRNA.png
        ├── ORF.orf
        ├── pipe
        │   ├── pipe.png
        ├── reference.bed
        ├── statistic_bed.tsv
        ├── statistic.csv
        ├── tissue
        │   ├── tissue_org.csv
        │   ├── tissue_org.png
        │   ├── transc_cod.csv
        │   ├── transc.csv
        │   └── transc_non.csv
        ├── tmhmm.csv
        └── trans.csv
        
* Folder "class", contains lncRNA class distribution.
* Folder "pipe", contains lncRNA and protein coding genes distribution after prediction.
* Folder "tissue", contains lncRNA distribution in tissue.

Groups of output files

**Image:** 
+ `intron_size.png` - lncRNA intron size distribution.
+ `classes.png` - distribution of lncRNA classes.    
+ `anti.png` - allocation of antisense lncRNA alignment to target gene structure.       
+ `exon_size.png` - lncRNA exon size distribution. 
+ `number_of_exon.png` - the ratio of the number of exons per lncRNA.
+ `number_of_lncRNA.png` - number of lncRNA per chromosome.                                                                                                 

**Prediction:**
+ `lncFinder.csv` - Transcripts predicted as lncRNA by lncFinder, `CSV` format.
+ `cpc.txt` - Transcripts predicted as lncRNA by CPC2, `TXT` format.
+ `Coding.fasta` - Transcripts predicted as protein coding by lncFinder, `FASTA` format. 
+ `Noncoding.fasta` - Transcripts predicted as lncRNA by lncFinder, `FASTA` format. 

**Transmembrane segment:**
+ `not_trans.csv` - id of transcripts without transmembrane domains, `CSV` format.
+ `ORF.orf` - reading frames of predicted lncRNA, `ORF` format.
+ `trans.csv` - id of transcripts with transmembrane domains, `CSV` format.
+ `tmhmm.csv` - TMHMM results, `CSV` format.
+ `Noncoding_trans_out.fasta` -set of lncRNA transcripts without transcripts with transmembrane domains, `FASTA` format.

**Alignment:**                                                                                         
+ `alignm.bed` - results of lncRNA transcripts alignment on reference genome by GMAP, 'BED' format.
+ `alignm_filter.gff` - results of lncRNA transcripts alignment on reference genome without long intron transcripts, `GFF` format.
+ `itron_coordin.tsv` - coordinates lncRNA introns, `TSV` format.
+ `statistic_bed.tsv` - lncRNA intron statistic, `TSV` format.
+ `statistic.csv` - lncRNA structure statistic, `CSV` format.
+ `blast.outfmt6` - BLASTn results in `outfmt6` format.

**Classification:**
+ `gffcmp.alignm_filter.gff.tmap` - lncRNA classification by gffcompare, `TMAP` format.
+ `gffcmp.loci` - lncRNA loci.
+ `reference.bed` - reference genome annotation, `BED` format.
+ `new_lncrna.fasta` - classified lncRNA transcripts,`FASTA` format.
