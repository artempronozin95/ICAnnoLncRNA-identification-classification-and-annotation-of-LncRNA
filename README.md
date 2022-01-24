# ICAnnoLncRNA - identification, classification and annotation of LncRNA
Pipeline for automatic identification, classification and annotation of plant lncRNAs based on their localization in the genome.
## Before work start, create environment:
    1. Open folder env
    2. conda env create --file env/programs.yaml
    3. conda activate pipeline
After these steps all necessary packages are installed. If you need update packages (not recommended), change the version of  packages after “=” (example - snakemake=4.0.1 -> snakemake=6.0.0), then “conda env update --file env/programs.yaml”. All necessary packages will be updated. 
## Input all necessary files into configuration file “config.yaml”:
+ lnc: - known LncRNA sequences of studied organisms.
+ cds: - known CDS (coding) sequences of studied organisms.
+ model: - pipeline check if model (lncFinder) for this organism already exist and use it, if it's not, pipeline will create new one (in output folder)
+ sequence: - sequences that need to study
+ structure: - need to choose use secondary structure in model building or not. Choose between DNA or SS (secondary structure).
+ gmap:
  + gff_reference: - gff file of reference genome.
  + reference: - reference genome in fasta format.
        out: - output file.
    gff: - gff file of reference genome.
    diamond: - alignment of lncRNA on proteine database (recomendent to run after main pipeline)
        option: - "on" (activate alignment step) or "off" (deactivate alignment step)
        database: - proteine database.
        query: - lncRNA
        out: - output file in outfmt6 format.
    tissue:
        organism: - choose between organisms in "Tissue analysis" or input your organism.
