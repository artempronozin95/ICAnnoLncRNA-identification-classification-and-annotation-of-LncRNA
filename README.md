
## Table of contents
* [Introduction](#identification-classification-and-annotation-of-lncrna)
* [Environment](#before-work-start-create-environment)
* [Configuration](#input-all-necessary-files-into-configuration)
* [Running the pipeline](#work-start)
* [Models](#models)
* [Structure information](#lncrna-structure-information)
* [Known LncRNA](#known-lncrna-for-database)
* [Tissue analysis](#tissue-analysis)
* [Output](#output)
## ICAnnoLncRNA - identification, classification and annotation of LncRNA
Pipeline for automatic identification, classification and annotation of plant lncRNAs based on their localization in the genome.
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
## Before work start, create environment:
    1. Open folder env
    2. conda env create --file env/programs.yaml
    3. conda activate pipeline
After these steps all necessary packages are installed. If you need update packages (**not recommended**), change the version of  packages after “=” (example - snakemake=4.0.1 -> snakemake=6.0.0), then “conda env update --file env/programs.yaml”. All necessary packages will be updated. 
It is nessesary install download [CPC2](https://github.com/biocoder/CPC2) program in pipeline folder.
## Input all necessary files into configuration file “config.yaml”:
+ lnc: - known LncRNA sequences of studied organisms.
+ cds: - known CDS (coding) sequences of studied organisms.
+ model: - pipeline check if model (lncFinder) for this organism already exist and use it, if it's not, pipeline will create new one (in output folder)
+ sequence: - sequences that need to study
+ structure: - need to choose use secondary structure in model building or not. Choose between DNA or SS (secondary structure).
+ gmap:
  + gff_reference: - gff file of reference genome.
  + reference: - reference genome in fasta format.
  + out: - output file.
+ gff: - gff file of reference genome.
+ diamond: - alignment of lncRNA on proteine database (recomendent to run after main pipeline)
  + option: - "on" (activate alignment step) or "off" (deactivate alignment step)
  + database: - proteine database.
  + query: - lncRNA
  + out: - output file in outfmt6 format.
+tissue:
  + organism: - choose between organisms in "Tissue analysis" or input your organism.
## Work start:
  #### 1. `snakemake -j 2`
  j or  --cores -  Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores.
  #### 2. `snakemake -nr` 
  n - Do not execute anything, and display what would be done. If you have a very large workflow, use –dry-run –quiet to just print a summary of the DAG of jobs.
  r - Print the reason for each executed rule.
  #### 3. `snakemake --use-conda`
  use-conda - Use additional conda environment.
  #### 4. Recomendent run: 
  `snakemake -j 2 --use-conda`
## Models:
+ *Zea_mays*
+ *Arabidopsis_thaliana*

## LncRNA structure information:
+ *Zea_mays*
+ *Arabidopsis_thaliana*
+ *Oryza_sativa*

## Known LncRNA for database:
Folder data/referencedata_index containe lncRNA library that build on base of 5 lncRNA databases (PNRD, CANTATAdb, GREENC, PlncDB, EVLncRNA)

## Tissue analysis:
File tissue/SRX_all_org.tsv, contain information about transcript experiment libraries respect to specific tissue. In config.yaml need to choose between organisms that presented in SRX file:
+ HV - *Hordeum vulgare*
+ OS - *Oryza sativa*
+ SL - *Solanum lycopersicum*
+ ST - *Solanum tuberosum*
+ ZM - *Zea mays*
Or input your organism in same format file.
## Output:
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
        ├── piep
        │   ├── piep.png
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
* Folder "piep", containe lncRNA and proteine coding genes destribution after prediction.
* Folder "tissue", containe lncRNA destribution in tissue.

**output**:
+ alignm.bed - results of GMAP alignment in 'bed' format.
+ alignm_filter.gff - filter of long intron sequences.
+ anti.png - antisense lncRNA distribution.                                                                                                                       
+ blast.outfmt6 - BLASTn results in 'outfmt6' format.
+ classes.png - containe lncRNA class distribution.                                                                                                               
+ Coding.fasta - All protein coding sequences in `fasta` format.                                                                                                   
+ cpc.txt - lncRNA predicted by CPC2.
+ exon_size.png - lncRNA exon size distribution.
+ gffcmp.alignm_filter.gff.tmap - lncRNA classification by gffcompare in 'tmap' format.
+ gffcmp.loci - lncRNA locus.
+ intron_size.png - lncRNA intor size distribution.
+ itron_coordin.tsv - introns of lncRNA coordinates on reference genome.
+ lncFinder_train.csv -lncRNA predicted by lncFinder.
+ new_lncrna.fasta - Classified lncRNA sequences in `fasta` format.
+ Noncoding.fasta - All lncRNA sequences in `fasta` format.
+ Noncoding_trans_out.fasta -lncRNA without sequences with transmembrane domains in `fasta` format.
+ not_trans.csv - id sequences without transmembrane domains.
+ number_of_exon.png - number of lncRNA exon distribution across one lncRNA.
+ number_of_lncRNA.png - lncRNA distribution across chromosome.
+ ORF.orf - reading frames of predicted lncRNA.
+ reference.bed - reference genome annotation in 'bed' format.
+ statistic_bed.tsv - lncRNA intron statistic.
+ statistic.csv - lncRNA structure statistic.
+ tmhmm.csv - TMHMM results.
+ trans.csv - id sequences with transmembrane domains.

