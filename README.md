## Table of contents
* [Introduction](#introduction)
* [Installation](#installation)
* [Input](#input)
* [Data](#data)
  * [Configuration](#configuration-file)
  * [Folders](#folders)
* [Running the pipeline](#work-start)
* [Models](#models)
* [Structure information](#lncrna-structure-information)
* [Known LncRNA](#known-lncrna-for-database)
* [Species](#species)
* [Tissue analysis](#tissue-analysis)
* [Expression and entropy](#expression-and-entropy)
* [Output](#output)
* [Errors](#errors)
## Introduction
Long non-coding RNAs (lncRNAs) are RNA molecules longer than 200 nucleotides that do not encode proteins. Experimental studies have shown the diversity and importance of lncRNA functions in plants. To expand knowledge about lncRNAs in other species, computational pipelines that allow for standardised data-processing steps in a mode that does not require user control up to the final result were actively developed recently. These advancements enable wider functionality for lncRNA data identification and analysis. In the present work, we propose the ICAnnoLncRNA pipeline for automatic identification, classification and annotation of plant lncRNAs in transcriptomic sequences assembled from high-throughput RNA sequencing (RNA-seq) data.

#### This pipeline is only applicable to the Linux operating system.

The pipeline includes the following steps: 
#### 1. Pre-processing.
+ Gmap index building.
+ LncFinder model building.
#### 2. LncRNA filtering.
+ Length filtering
+ LncFinder predicting.
+ CPC2 predicting.
+ GMAP alignment on reference genome.
+ Filtering erros/noise.
+ Filtering possible transposable elements (TEs).
#### 3. LncRNA annotation.
+ Classification gffcompare.
+ Blastn alignment.
+ Analysis of lncRNA expression.

The pipeline is implemented using the workflow management system [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), which provides ability to platform-independent installation and execution of the software.

## Schematic diagram
![Test Image 1](https://github.com/artempronozin95/ICAnnoLncRNA-identification-classification-and-annotation-of-LncRNA/blob/main/image/Pipeline.png)

## Installation 
# Automatic
**recommended for clusters/servers**

Install only **programs.yaml** environment. Other environments will install automatically when ICAnnoLncRNA start work.
```
1. wget https://github.com/artempronozin95/ICAnnoLncRNA-identification-classification-and-annotation-of-LncRNA/archive/refs/heads/main.zip
2. unzip main.zip
3. cd ./ICAnnoLncRNA-identification-classification-and-annotation-of-LncRNA-main
4. conda env create --file env/programs.yaml
5. conda activate ICAnnoLncRNA
```
After these steps all necessary packages are installed. If you need update packages (**not recommended**), change the version of  packages after “=” (example - `snakemake=4.0.1 -> snakemake=6.0.0`), then `conda env update --file ./programs.yaml`. All necessary packages will be updated. Recomended on clusters, requires a lot of  processing power.
# Step method
**recommended for personal computer**
```
1. conda update conda.
2. conda create -n ICAnnoLncRNA python=3.6
3. conda activate ICAnnoLncRNA
4. conda install -c bioconda emboss
5. install next packeges from file below
```
Install all packeges represented **[here](https://github.com/artempronozin95/ICAnnoLncRNA---identification-classification-and-annotation-of-LncRNA/blob/main/env/Dependencies%20information.md)**

# It is necessary to install and download (download program in pipeline folder): 
+ [CPC2](https://github.com/biocoder/CPC2)
```
1. cd ./ICAnnoLncRNA---identification-classification-and-annotation-of-LncRNA-main
2. wget https://github.com/biocoder/CPC2/archive/refs/heads/master.zip
3. unzip master.zip
4. mv CPC2-master CPC2
```
+ [LncAPDB.fasta](https://data.mendeley.com/datasets/fnk8pmp2yz/3)

download LncAPDB.fasta into `data/reference/data_index` folder.

## Input
### Genome sequence
1. Known LncRNA transcripts of the species in `FASTA` format. File with lncRNA from the public database of lncRNA. For example [Ensembl](https://www.ensembl.org/index.html), [EVLncRNAs v2](https://www.sdklab-biophysics-dzu.net/EVLncRNAs2/).
2. Known protein coding sequences of the species in `FASTA` format. File with protein coding sequences from a public database. For example [Ensembl](https://www.ensembl.org/index.html).
3. RNA-seq transcripts of the species in `FASTA` format. This file contains transcripts obtained after assembly of transcriptome. For example file that we provide, [30k_zey.fasta](https://github.com/artempronozin95/ICAnnoLncRNA---identification-classification-and-annotation-of-LncRNA/blob/main/data/input/30k_zey.zip), obtained by [Trinity](https://github.com/trinityrnaseq/trinityrnaseq) assembly.
All three sets are necessary for build 
### Reference genome
1. Reference genome of the species in `FASTA` format.
2. Annotation of the species in `GFF/GTF` format.
### Additional data
1. TE coordinats on the reference genome, example [Zea_N_merged.bed](https://github.com/artempronozin95/ICAnnoLncRNA-identification-classification-and-annotation-of-LncRNA/blob/main/data/reference/N_coords/Zea_N_merged.bed)
2. File with sequence ID from public lncRNA databases connected to IDs in library of known lncRNAs. example [index_and_newindex.csv](https://github.com/artempronozin95/ICAnnoLncRNA-identification-classification-and-annotation-of-LncRNA/blob/main/tissue/index_and_newindex.csv)
3. Annotation of the organism transcriptome libraries by tissue type, example [SRX_all_org.tsv](https://github.com/artempronozin95/ICAnnoLncRNA-identification-classification-and-annotation-of-LncRNA/blob/main/tissue/SRX_all_org.tsv)

## Data
Additional data is presented here: https://data.mendeley.com/datasets/fnk8pmp2yz/3

It includes: 
+ blast.outfmt6 - Blastn results. Contain homologs with known lncRNAs sequences from the LncAPDB library.
+ index_and_newindex.fasta - index of PNRD, CANTATAdb, GREENC, PlncDB, EVLncRNA databases compared with new index for LncAPDB library. 
+ LncAPDB.fasta - lncRNA sequences of LncAPDB library in fasta format.
+ lncrna_class.tmap - novel lncRNAs divided into gffcompare classes.
+ lncrna_coordinats.bed - coordinates of novel lncRNAs on chromosomes.
+ lncrna_intron_size.tsv - intron size of novel lncRNAs and their coordinates on the genome.
+ new_lncrna.fasta - novel lncRNA sequences in fasta format.
+ new_lncrna_locus.loci - locus of lncRNA sequences on genome.
+ transcriptome_lib.txt - Maize transcriptome libraries. 

### Configuration file
Input all necessary files into configuration file “config.yaml”:
+ `lnc:` - known LncRNA sequences of studied organisms in `FASTA` format. 
  + (Example: `lnc: "data/input/lincrna.fasta"`)
+ `cds:` - known CDS (coding) sequences of studied organisms in `FASTA` format. 
  + (Example: `cds: "data/input/7000_len_mRNA.fasta"`)
+ `model:` - pipeline check if model (lncFinder) for this organism already exist and use it, if it's not, pipeline will create a new one (in output folder). Model that exist represented in [Models](#models). 
  + (Example: `model: "Zea_mays"`)
+ `sequence:` - sequences that need to study in `FASTA` format.
  + (Example: `sequence: "data/input/Zea.fasta"`)
+ `structure:` - need to choose whether to use secondary structure in model building or not. Choose between `DNA` or `SS` (secondary structure). 
  + (Example: `structure: "DNA"`)
+ `gmap:`
  + `gff_reference:` - reference genome annotation in `GFF` format. 
    + (Example: `gff_reference: "data/reference/Zea_mays.AGPv4.40.gff"`)
  + `reference:` - reference genome in `FASTA` format.
    +  (Example: `reference: "data/reference/Zea_mays.AGPv4.dna.toplevel.fa"`)
  + `out:` - output file. 
    + (Example: `out: "data/output/alignm.gff"`)
+ `gff:` - gff file of reference genome. 
  + (Example: `gff: "data/reference/Zea_mays.AGPv4.40.gff"`)
+ `blast:` - standard BLAST parameters
  + `evalue:` - evalue
    + (Example: `evalue: 1e-50`)
  + `max_target:` - max_target_seqs
    + (Example: `max_target: 1`)
  + `identity:` - perc_identity
    + (Example: `identity: 30`)
  + `threads:` - num_threads
    + (Example: `threads: 1`)
+ `TE_coords:` - if there is no TE coordinats information - recomended to turn off this step.
   + `option:` "on" - turn on or turn off this step (on/off) 
   + `coords:` - TE coordinat
     + (Example: `coords: "data/reference/N_coords/Zea_N_merged.bed"`)
+ `tissue:`
  + `organism:` - choose species that you study form species [list](#species), or input your own species ID in same format
    + (Example: `organism: "ZMAY"`)
  + `exp:` - choose between organisms tissue experiments in [Tissue analysis](#tissue-analysis) or input your organism. 
    + (Example: `exp: "ZM"`)
  + `LncAPDB:` - File with sequence ID from public lncRNA databases connected to IDs in library of known lncRNAs.
    + (Example: `LncAPDB: "tissue/index_and_newindex.csv"`)
  + `SRX_exp:` - Annotation of the transcriptome libraries by tissue type and gene expression levels.
    + (Example: `SRX_exp: "tissue/SRX_all_org.tsv"`)
+ `expression`:
  + `option`: "on" - turn on or turn off this step (on/off) 
  + `path`: - path to expression data.
    + (Example: `path: "expressin/Kallisto"`)
  + `calculation_type`: - type of expression data. There three types: 1.Kallisto 2.Htseq 3.Own type of data.
    + (Example: `path: "Kallisto"`)
    
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
+ `N_coords` folder with TE coordinats
## Work start
  #### 1. `snakemake -j 2`
  `-j` or  `--cores` -  Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores.
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
Folder `data/reference/data_index` contain lncRNA library (LncAPDB.fasta) that build on base of 5 lncRNA databases (PNRD, CANTATAdb, GREENC, PlncDB, EVLncRNA). The library you can find here: https://data.mendeley.com/v1/datasets/fnk8pmp2yz/draft?a=13de1dfb-b631-42f3-8108-f3de2f50fd90

## Species
+ *Medicago truncatula* - `MTRU`
+ *Glycine max* - `GMAX`
+ *Populis trichocarpa* - `PTRI`
+ *Arabidopsis thaliana* - `ATHA`
+ *Vitis vivifera* - `VVIN`
+ *Oryza sativa japonica* - `OSAT`
+ *Brachipodiumdistachion* - `BDIS`
+ *Sorghum bicolor* - `SBIC`
+ *Zea mays* - `ZMAY`
+ *Selaginella moellendorfii* - `SMOE`
+ *Physcomitrella patens* - `PPAT`
+ *Ostreococcus tauri*- `OTAU`
+ *Volvox* - `VOLV`
+ *Amborell trichopoda* - `ATRI`

## Tissue analysis
File `tissue/SRX_all_org.tsv`, contains information about 1241 transcript experiment libraries for 5 organisms with respect to specific tissue. In config.yaml need to choose between organisms that presented in SRX file:
+ HV - *Hordeum vulgare*
+ OS - *Oryza sativa*
+ SL - *Solanum lycopersicum*
+ ST - *Solanum tuberosum*
+ ZM - *Zea mays*

If the libraries you are investigating (for other organism for exemple) are not in this file. You can add them to the data file in the same format 
```
ST	SRX1478098	leaf
ZM	SRX339769	stomatal division zone
ZM	SRX339787	ear
ZM	SRX711129	leaf
OS	SRX2582231	unknown
ZM	SRX711024	ear
ZM	SRX711053	leaf

```
Where 1 column - id, 2 - organism that you study, 3 - SRX library, 4 - tissue. 

## Expression and entropy
For this step there is necessary to provide expression data for every transcript. Thus recommended `turn off` this step at first run of the pipeline. After receiving expression data `turn on` this step and pipeline will calculete entropy and espression statistic. In config.yaml need to choose between three types of table building:
+ `Kallisto` - result are received by means of Kallisto program. Results should be presented as on example:
```
├── SRR957466
│   ├── abundance.h5
│   ├── abundance.tsv
│   └── run_info.json
├── SRR957467
│   ├── abundance.h5
│   ├── abundance.tsv
│   └── run_info.json
└── SRR957468
    ├── abundance.h5
    ├── abundance.tsv
    └── run_info.json
```  
+ `Htseq` - result are received by means of Htseq program. Results should be presented as on example:
```
├── SRR1586686.txt
└── SRR504468.txt
```
+ `own` - In this case user should build table by himself. Table should looks like as [table](https://github.com/artempronozin95/ICAnnoLncRNA-identification-classification-and-annotation-of-LncRNA/blob/main/example/expr_table.tsv) presented in example folder. Where:
+ `type` - `mRNA` - is protein coding genes. `noncons` - nonconserved lncRNAs. `consv` - conserved lncRNAs.
+ all tissue that needed. In example table, value is a median of expression for libraries that belong to tissue in column name.
+ `target_id` - ID of the transcript. it is necessary that the transcripts correspond to the results of the pipeline.

## Output
A typical structure of `Output` is follows:
```
├── gffcompare_first
    ├── gffcmp.annotated.gtf
    ├── gffcmp.filter_alignm.bed.refmap
    ├── gffcmp.filter_alignm.bed.tmap
    ├── gffcmp.loci
    ├── gffcmp.stats
    └── gffcmp.tracking
├── gffcompare_second
    ├── gff.annotated.gtf
    ├── gff.filter_alignm.bed.refmap
    ├── gff.filter_alignm.bed.tmap
    ├── gff.loci
    ├── gff.stats
    └── gff.tracking
├── GMAP
    ├── alignm.bed
    ├── alignm_filter.gff
    ├── filter_alignm.bed
    ├── gmap_build.err.log
    ├── gmap_build.out.log
    ├── gmap.out.log
    ├── reference.bed
    └── statistic.csv
├── lncRNA_prediction
    ├── Coding.fasta
    ├── cpc.txt
    ├── lncFinder.csv
    └── Noncoding.fasta
├── lncRNA_structure
    ├── anti.png
    ├── exon_size.png
    ├── intron_size.png
    ├── itron_coordin.tsv
    ├── long_transcripts.bed
    ├── number_of_exon.png
    ├── number_of_lncRNA.png
    └── statistic_bed.tsv
├── loci
    ├── lncRNA_before_loci.bed
    └── lncRNA_loci.bed
├── new_lncRNA
    ├── blast.outfmt6
    ├── classes.png
    ├── LncAPDB_vs_blast.csv
    ├── new_lncrna.fasta
    └── True_lncRNA.bed
├── TE_finder
    ├── Lnc_aling_with_TE.tsv
    └── lncRNA_after_loci.bed
├── tissue
    ├── all.txt
    ├── cod.txt
    ├── cons.txt
    ├── id.txt
    ├── non.txt
    ├── tissue_org.csv
    ├── tissue_org.png
    ├── transc_cod.csv
    ├── transc.csv
    └── transc_non.csv
├── expression
    ├── box_plot.png
    ├── entropy.png
    ├── expr_table.tsv
    ├── hist_expr.png
    └── table_type.txt
```
* Folder "tissue", contains lncRNA distribution in tissue.

Groups of output files

**1. Pre-processing.**

Pre-processing performed in input directory:
+ filtered_seq.fasta - filter transcripts with length < 200 nt.
+ new_vs_old_id.tsv - new ID and old ID of transcripts.
+ test_train - folder of lncFinder model studie. 

**2. LncRNA filtering.**

***lncRNA_prediction***

`lncRNA_prediction`
+ Coding.fasta - Transcripts predicted as protein coding by lncFinder, `FASTA` format. 
+ cpc.txt - Transcripts predicted as lncRNA by CPC2, `TXT` format.
+ lncFinder.csv - Transcripts predicted as lncRNA by lncFinder, `CSV` format.
+ Noncoding.fasta - Transcripts predicted as lncRNA by lncFinder, `FASTA` format.

***GMAP alignment on reference genome.***

`GMAP`
+ alignm.bed - results of lncRNA transcripts alignment on reference genome by GMAP, 'BED' format.
+ alignm_filter.gff - results of lncRNA transcripts alignment on reference genome without long intron transcripts, `GFF` format.
+ filter_alignm.bed - results of lncRNA transcripts alignment on reference genome without long intron transcripts, `BED` format.
+ gmap_build.err.log
+ gmap_build.out.log
+ gmap.out.log
+ reference.bed - reference genome annotation, `BED` format.
+ statistic.csv - lncRNA structure statistic, `CSV` format.

***Filtering erros/noise.***

`gffcompare_first` - classification of the candidate lncRNA sequences by their location in the genome relative to a protein-coding genes, more information can be found on [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)
+ gffcmp.annotated.gtf - Gffcompare reports a GTF file containing the “union” of all transcript IDs in each sample.
+ gffcmp.filter_alignm.bed.refmap - This tab-delimited file lists, for each reference transcript, which query transcript either fully or partially matches that reference transcript.
+ gffcmp.filter_alignm.bed.tmap - This tab delimited file lists the most closely matching reference transcript for each query transcript. 
+ gffcmp.loci - Represent loci of transcripts. 
+ gffcmp.stats
+ gffcmp.tracking - This file matches transcripts up between samples.

`loci`
+ lncRNA_before_loci.bed - the candidate lncRNA sequences before merging into loci, `BED` format.
+ lncRNA_loci.bed - loci of the candidate lncRNA, `BED` format.

`gffcompare_second`- classification of the candidate lncRNA sequences by their location in the genome relative to lncRNA loci.
+ gff.annotated.gtf
+ gff.filter_alignm.bed.refmap
+ gff.filter_alignm.bed.tmap
+ gff.loci
+ gff.stats
+ gff.tracking

***Filtering possible transposable elements (TEs).***

`TE_finder`
+ Lnc_aling_with_TE.tsv - the candidate lncRNA sequences that mapped on TE, `TSV` format.
+ lncRNA_after_loci.bed - the candidate lncRNA sequences that are completely matched (`=` - class of gff compare) with lncRNA loci, `BED` format.

**3. LncRNA annotation.**

`lncRNA_structure`
+ anti.png - allocation of antisense lncRNA alignment to target gene structure.
+ exon_size.png - lncRNA exon size distribution. 
+ intron_size.png - lncRNA intron size distribution.
+ itron_coordin.tsv - coordinates lncRNA introns, `TSV` format.
+ long_transcripts.bed - transcripts with length > 5000 nt, `BED` format.
+ number_of_exon.png - the ratio of the number of exons per lncRNA.
+ number_of_lncRNA.png - number of lncRNA per chromosome.
+ statistic_bed.tsv - lncRNA intron statistic, `TSV` format.

`new_lncRNA`
+ blast.outfmt6 - BLASTn results in `outfmt6` format.
+ classes.png - distribution of lncRNA classes. 
+ LncAPDB_vs_blast.csv - known lncRNAs that were predicted.
+ new_lncrna.fasta - identified lncRNA transcripts,`FASTA` format. (**final result**)
+ True_lncRNA.bed - identified lncRNA transcripts,`BED` format.

`tissue`
+ all.txt - lncRNA that have homologs with known lncRNAs.
+ cod.txt - ID of proteine coding transcripts.
+ cons.txt - ID of the lncRNAs similar to lncRNAs from other species (conserved sequences)
+ id.txt - ID of all input transcripts
+ non.txt - ID of the lncRNAs similar only to known studied organism sequences (non-conserved sequences)
+ tissue_org.csv - all in one (proteine coding, conserved lncRNAs, non-conserved lncRNAs) 
+ tissue_org.png - heatmap of transcript libraries tissue type.
+ transc_cod.csv - proteine coding transcript libraries tissue type.
+ transc.csv - conserved lncRNAs transcript libraries tissue type.
+ transc_non.csv - non-conserved lncRNAs transcript libraries tissue type.
+ Blastn alignment.
+ Analysis of lncRNA expression.

`expression`
+ box_plot.png - box plot diagrams of gene expression.
+ entropy.png - density plot of Shannon entropy.
+ expr_table.tsv - table with expression data for every transcript in every tissue.
+ hist_expr.png - density plot of maximum expression levels of transcripts.
+ table_type.txt - type of expression results.

## Errors
```
rule model_lncFinder:
    input: data/input/test_train/sets.txt                                    
    output: data/input/test_train/dirs.txt                                                                                                            
    jobid: 2    

Activating conda environment /home/pronozin/ICAnnoLncRNA---identification-classification-and-annotation-of-LncRNA-main/.snakemake/conda/767d2b45.
/bin/bash: activate:
```
**Solution**: type `conda install conda` in ICAnnoLncRNA environment.

*************************************************************************************
```
rule Gmap_build:
    input: data/reference/Zea_mays.AGPv4.dna.toplevel.fa, data/reference/Zea_mays.AGPv4.40.gff
    output: data/output/statistic.csv
    jobid: 5
...

Error in rule Gmap_build:
    jobid: 5
    output: data/output/statistic.csv

RuleException:
CalledProcessError in line 123 of ...:
Command ' set -euo pipefail; ...'
```
**Solution**: update gtfparse package. `Example: conda install -c bioconda gtfparse`

*************************************************************************************
```
error while loading shared libraries: libgfortran.so.3
```
This error connected with updating of libraries of libgfortran at Conda database.

**Solution**: Try these two commands:
```
locate libgfortran.so.3
locate libgfortran
```
If you manage to locate libgfortran.so.3, you can make a symbolic link to it inside a known library directory. If you find this library but it ends in a higher number than 3, say libgfortran.so.4, go to that directory and make a symbolic link :
```
ln -s libgfortran.so.4 libgfortran.so.3
```
If error still occurs, then download this files [libgfortran.so.3](https://github.com/artempronozin95/ICAnnoLncRNA-identification-classification-and-annotation-of-LncRNA/blob/main/example/libgfortran.so.3), [libgfortran.so.3.0.0](https://github.com/artempronozin95/ICAnnoLncRNA-identification-classification-and-annotation-of-LncRNA/blob/main/example/libgfortran.so.3.0.0) and paste them into **lib** folder of your conda environment where R situated.
For example ```conda/60023710/lib``` or ```miniconda3/envs/lnc/lib```


