## Table of contents
* [Introduction](#identification-classification-and-annotation-of-lncrna)
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
## ICAnnoLncRNA - identification, classification and annotation of LncRNA
Long non-coding RNAs (lncRNAs) typically defined as transcripts of more than 200 nucleotides in length and without any protein coding potential. These RNAs are involved in important plant development processes such as phosphate homeostasis, flowering, photomorphogenesis and stress response. However, their structural and functional properties are not clear. Information about lncRNA sequences and their expression typically obtained from RNA-seq experiments. To process these data relevant tools for automated lncRNA are required. 
Here we propose ICAnnoLncRNA, pipeline for automatic identification, classification and annotation of plant lncRNAs based on their localization in the genome. ICAnnoLncRNA, work with RNA-seq data and it was tested on 877 transcriptome libraries of *Zea mays*.
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

![Test Image 1](https://github.com/artempronozin95/ICAnnoLncRNA---identification-classification-and-annotation-of-LncRNA/blob/main/image/Pipeline.png)
## Installation 
    1. Open folder env
    2. conda env create --file ./programs.yaml
    3. conda activate ICAnnoLncRNA
After these steps all necessary packages are installed. If you need update packages (**not recommended**), change the version of  packages after “=” (example - snakemake=4.0.1 -> snakemake=6.0.0), then “conda env update --file ./programs.yaml”. All necessary packages will be updated. 
More detailed description provided [here](https://github.com/artempronozin95/ICAnnoLncRNA---identification-classification-and-annotation-of-LncRNA/blob/main/env/Dependencies%20information.md)

It is nessesary to install (download program in pipeline folder): 
+ [CPC2](https://github.com/biocoder/CPC2)
+ [TMHMM](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0)
## Input
### Configuration file
Input all necessary files into configuration file “config.yaml”:
+ lnc: - known LncRNA sequences of studied organisms in FASTA format.
+ cds: - known CDS (coding) sequences of studied organisms in FASTA format.
+ model: - pipeline check if model (lncFinder) for this organism already exist and use it, if it's not, pipeline will create new one (in output folder). Model that exist represented in [Models](#models).
+ sequence: - sequences that need to study in FASTA format
+ structure: - need to choose use secondary structure in model building or not. Choose between DNA or SS (secondary structure).
+ gmap:
  + gff_reference: - reference genome annotation in GFF format.
  + reference: - reference genome in FASTA format.
  + out: - output file.
+ gff: - gff file of reference genome.
+ diamond: - alignment of lncRNA on proteine database (recomendent to run after main pipeline)
  + option: - "on" (activate alignment step) or "off" (deactivate alignment step)
  + database: - proteine database in FASTA format.
  + query: - lncRNA transcripts in FASTA format.
  + out: - output file in outfmt6 format.
+ tissue:
  + organism: - choose between organisms in [Tissue analysis](#tissue-analysis) or input your organism.
### Folders
#### data/input
Contain: 
+ known LncRNA sequences of studied organisms.
+ known coding sequences of studied organisms.
+ sequences that need to study.
#### data/reference
Contain:
+ reference genome of studied organism in 'fasta' format.
+ annotation of studied organism in 'gff/gtf' format.
+ `data_index` folder with library of known lncRNA from databases.
+ `intron_data` folder with structure information of organisms.
+ `models` folder with model for lncFinder.
## Work start
  #### 1. `snakemake -j 2`
  j or  --cores -  Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores.
  #### 2. `snakemake -nr` 
  n - Do not execute anything, and display what would be done. If you have a very large workflow, use –dry-run –quiet to just print a summary of the DAG of jobs.
  r - Print the reason for each executed rule.
  #### 3. `snakemake --use-conda`
  use-conda - Use additional conda environment.
  #### 4. Recomendent run: 
  `snakemake -j 2 --use-conda`
## Models
+ *Zea_mays*
+ *Arabidopsis_thaliana*

## LncRNA structure information
+ *Zea_mays*
+ *Arabidopsis_thaliana*
+ *Oryza_sativa*

## Known LncRNA for database
Folder `data/reference/data_index` containe lncRNA library that build on base of 5 lncRNA databases (PNRD, CANTATAdb, GREENC, PlncDB, EVLncRNA)

## Tissue analysis
File `tissue/SRX_all_org.tsv`, contain information about transcript experiment libraries respect to specific tissue. In config.yaml need to choose between organisms that presented in SRX file:
+ HV - *Hordeum vulgare*
+ OS - *Oryza sativa*
+ SL - *Solanum lycopersicum*
+ ST - *Solanum tuberosum*
+ ZM - *Zea mays*
Or input your organism in same format file.
## Output
A typical structure of `Output` is follows:

    output/
        ├── alignm.bed
        ├── alignm_filter.gff
        ├── anti.png                                                                                                                         
        ├── blast.outfmt6                                                                                                                    
        ├── class                                                                                                                                                
        │   ├── class_org.png                                                                                                                                    
        ├── classes.png                                                                                                                      
        ├── Coding.fasta                                                                                                                     
        ├── cpc.txt
        ├── exon_size.png
        ├── filter_alignm.bed
        ├── gffcmp.alignm_filter.gff.tmap
        ├── gffcmp.loci
        ├── intron_size.png
        ├── itron_coordin.tsv
        ├── lncFinder_train.csv
        ├── new_lncrna.fasta
        ├── Noncoding.fasta
        ├── Noncoding_trans_out.fasta
        ├── not_trans.csv
        ├── number_of_exon.png
        ├── number_of_lncRNA.png
        ├── ORF.orf
        ├── pipe
        │   ├── pipe.png
        ├── reference.bed
        ├── statistic_bed.tsv
        ├── statistic.csv
        ├── tissue
        │   ├── tissue_org.csv
        │   ├── tissue_org.png
        │   ├── transc_cod.csv
        │   ├── transc.csv
        │   └── transc_non.csv
        ├── tmhmm.csv
        └── trans.csv
        
* Folder "class", containe lncRNA class destribution.
* Folder "pipe", containe lncRNA and proteine coding genes destribution after prediction.
* Folder "tissue", containe lncRNA destribution in tissue.

Groups of output files

**Image:** 
+ `intron_size.png` - lncRNA intor size distribution.
+ `classes.png` - lncRNA class distribution.    
+ `anti.png` - antisense lncRNA distribution.       
+ `exon_size.png` - lncRNA exon size distribution. 
+ `number_of_exon.png` - lncRNA exon distribution across one lncRNA.
+ `number_of_lncRNA.png` - lncRNA distribution across chromosome.                                                                                                 

**Prediction:**
+ `lncFinder.csv` - Transcripts predicted as lncRNA by lncFinder. 
+ `cpc.txt` - Transcripts predicted as lncRNA by CPC2.
+ `Coding.fasta` - Transcripts predicted as protein coding by lncFinder, `FASTA` format. 
+ `Noncoding.fasta` - Transcripts predicted as lncRNA by lncFinder, `FASTA` format. 

**Transmembrane segment:**
+ `not_trans.csv` - id of transcripts without transmembrane domains.
+ `ORF.orf` - reading frames of predicted lncRNA, `ORF` format.
+ `trans.csv` - id of transcripts with transmembrane domains.
+ `tmhmm.csv` - TMHMM results, `CSV` format.
+ `Noncoding_trans_out.fasta` -set of lncRNA transcripts without transcripts with transmembrane domains, `FASTA` format.

**Alignment:**                                                                                         
+ `alignm.bed` - results of lncRNA transcripts alignment on reference genome by GMAP, 'BED' format.
+ `alignm_filter.gff` - results of lncRNA transcripts alignment on reference genome witout long intron transcripts, `GFF` format.
+ `itron_coordin.tsv` - coordinates lncRNA introns, `TSV` format.
+ `statistic_bed.tsv` - lncRNA intron statistic.
+ `statistic.csv` - lncRNA structure statistic, `CSV` format.
+ `blast.outfmt6` - BLASTn results in `outfmt6` format.

**Classification:**
+ `gffcmp.alignm_filter.gff.tmap` - lncRNA classification by gffcompare, 'tmap' format.
+ `gffcmp.loci` - lncRNA loci.
+ `reference.bed` - reference genome annotation, 'BED' format.
+ `new_lncrna.fasta` - classified lncRNA transcripts,`FASTA` format.
