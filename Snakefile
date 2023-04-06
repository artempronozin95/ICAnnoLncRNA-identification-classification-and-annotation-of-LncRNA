# read config info into this namespace
configfile: "config.yaml"
print (config['lnc'], config['cds'], config['sequence'])



rule all:
    input:
        "data/input/test_train/sets.txt",
        "data/output/lncRNA_prediction/cpc.txt",
        "data/input/test_train/best_model.txt",
        "data/input/test_train/dirs.txt",
        "data/input/filtered_seq.fasta",
        "data/output/lncRNA_prediction/Noncoding.fasta",
        "data/output/lncRNA_prediction/lncFinder.csv",
        "data/output/GMAP/statistic.csv",
        "data/output/GMAP/alignm_filter.gff",
        "data/output/GMAP/alignm.bed",
        "data/output/GMAP/filter_alignm.bed",
        "data/output/gffcompare_first/gffcmp.filter_alignm.bed.tmap",
        "data/output/loci/lncRNA_loci.bed",
        "data/output/gffcompare_second/gff.filter_alignm.bed.tmap",
        "data/output/new_lncRNA/True_lncRNA.bed",
        "data/output/lncRNA_structure/itron_coordin.tsv",
        "data/output/lncRNA_structure/statistic_bed.tsv",
        "data/output/new_lncRNA/blast.outfmt6",
        "data/output/tissue/tissue_org.png"

rule test_train:
    input:
       expand("{lnc}", lnc=config['lnc']),
       expand("{cds}", cds=config['cds'])
    output:
       "data/input/test_train/sets.txt"
    run:
       shell("python scripts/test_train.py {input} > {output}")
       
rule model_lncFinder:
    input:
       set="data/input/test_train/sets.txt"
    params:
       path=ancient("data/input/test_train"),
       str=expand("{strc}", strc=config['structure']),
       lnc="data/input/test_train/train_lnc.fasta"
    output:
        "data/input/test_train/dirs.txt"
    conda:
        "env/lncfinder.yaml"
    shell:
        "Rscript scripts/model_lncFind.r {params.path} {params.lnc} {params.str} >> {output}"
rule F1:
    input:
      set="data/input/test_train/dirs.txt"
    params:
      "data/input/test_train/"
    output:
      "data/input/test_train/best_model.txt"
    shell:
      "python scripts/TP_FN.py {params}"

rule filter_short_seq:
    input:
      expand("{sequence}", sequence=config['sequence'])
    params:
      "data/input/new_vs_old_id.tsv"
    output:
      "data/input/filtered_seq.fasta"
    shell:
      "python scripts/rename.py {input} {output} {params}"

rule lncFinder:
    input:
       "data/input/test_train/best_model.txt",
       "data/input/filtered_seq.fasta",
    params:
       ancient("data/output/lncRNA_prediction"),
       expand("{strc}", strc=config['structure'])
    output:
        "data/output/lncRNA_prediction/lncFinder.csv"
    conda:
        "env/lncfinder.yaml"
    shell:
        "scripts/Chunk_dataframe.sh {input} {params}"

rule cpc:
    input:
        "data/input/filtered_seq.fasta"
    params:
        "data/output/lncRNA_prediction/cpc"
    output:
        "data/output/lncRNA_prediction/cpc.txt"
    conda:
        "env/cpc.yaml"
    shell:
        "CPC2/bin/CPC2.py -i {input} -r -o {params}"
        
rule take_noncoding:
    input:
        ancient("data/output/lncRNA_prediction/lncFinder.csv"),
        "data/input/filtered_seq.fasta"
    output:
        "data/output/lncRNA_prediction/Noncoding.fasta"
    shell:
        "python scripts/take_noncoding.py {input} {output}"
  
rule Gmap_build:
    input:
        expand("{reference}", reference=config['gmap']['reference']),
        expand("{gff_reference}", gff_reference=config['gmap']['gff_reference'])
    output:
       "data/output/GMAP/statistic.csv"
    shell:
       "python scripts/GMAP_build.py {input} {output}"

rule Gmap_aling: 
    input: 
        expand("{reference}", reference=config['gmap']['reference']),
        "data/output/lncRNA_prediction/Noncoding.fasta",
        "data/output/GMAP/statistic.csv"
    output:
        "data/output/GMAP/alignm_filter.gff"
    shell:
        "python scripts/GMAP.py {input} {config[model]} {output}"
        
rule to_bed:
    input:
        query="data/output/GMAP/alignm_filter.gff",
        reference=expand("{gff_reference}", gff_reference=config['gmap']['gff_reference'])
    output:
        query_out="data/output/GMAP/alignm.bed",
        reference_out="data/output/GMAP/reference.bed"
    run:
        shell("gff2bed < {input.query} > {output.query_out}")
        shell("gff2bed < {input.reference} > {output.reference_out}")

rule filter_long_introns:
    input:
       query_out="data/output/GMAP/alignm.bed",
       reference_out="data/output/GMAP/reference.bed"
    output:
       "data/output/GMAP/filter_alignm.bed",
       "data/output/lncRNA_structure/itron_coordin.tsv"
    shell:
       "python scripts/gmap_filter.py {input} {output}"

rule gffcomare:
    input:
        expand("{annot}", annot=config['gff']),
        "data/output/GMAP/filter_alignm.bed"
    params:
        "data/output/gffcompare_first/gffcmp"
    output:
        "data/output/gffcompare_first/gffcmp.filter_alignm.bed.tmap"
    run:
        shell("gffcompare -r {input} -o {params}")
        shell("mv data/output/GMAP/gffcmp.filter_alignm.bed.tmap data/output/gffcompare_first/")
        shell("mv data/output/GMAP/gffcmp.filter_alignm.bed.refmap data/output/gffcompare_first/")

if config['TE_coords']['option'] == 'on':
  rule lncRNA_loci:
    input:
        "data/output/gffcompare_first/gffcmp.filter_alignm.bed.tmap",
        "data/output/GMAP/filter_alignm.bed"
    params:
        "data/output/loci/lncRNA_before_loci.bed"
    output:
        "data/output/loci/lncRNA_loci.bed"
    run:
        shell("python scripts/lncRNA_class.py {input} {params} {output}")

  rule gffcomare_with_loci:
    input:
        "data/output/loci/lncRNA_loci.bed",
        "data/output/GMAP/filter_alignm.bed"
    params:
        "data/output/gffcompare_second/gff"
    output:
        "data/output/gffcompare_second/gff.filter_alignm.bed.tmap"
    run:
        shell("gffcompare -r {input} -o {params}")
        shell("mv data/output/GMAP/gff.filter_alignm.bed.tmap data/output/gffcompare_second/")
        shell("mv data/output/GMAP/gff.filter_alignm.bed.refmap data/output/gffcompare_second/")

  rule TE_finder:
    input:
       gffcomp="data/output/gffcompare_second/gff.filter_alignm.bed.tmap",
       aling="data/output/GMAP/filter_alignm.bed",
       coords=expand("{TE_coords}", TE_coords=config['TE_coords']['coords'])
    params:
        option=expand("{option}", option=config['TE_coords']['option'])
    output:
       loci="data/output/TE_finder/lncRNA_after_loci.bed",
       TE="data/output/TE_finder/Lnc_aling_with_TE.tsv",
       true_lnc="data/output/new_lncRNA/True_lncRNA.bed"
    shell:
       "python scripts/bed_intersect.py {input.gffcomp} {input.aling} {params.option} {output.loci} {input.coords} {output.TE} {output.true_lnc}"

elif config['TE_coords']['option'] == 'off':
  rule lncRNA_loci:
    input:
        "data/output/gffcompare_first/gffcmp.filter_alignm.bed.tmap",
        "data/output/GMAP/filter_alignm.bed"
    params:
        "data/output/loci/lncRNA_before_loci.bed"
    output:
        "data/output/loci/lncRNA_loci.bed"
    run:
        shell("python scripts/lncRNA_class.py {input} {params} {output}")

  rule gffcomare_with_loci:
    input:
        "data/output/loci/lncRNA_loci.bed",
        "data/output/GMAP/filter_alignm.bed"
    params:
        "data/output/gffcompare_second/gff"
    output:
        "data/output/gffcompare_second/gff.filter_alignm.bed.tmap"
    run:
        shell("gffcompare -r {input} -o {params}")
        shell("mv data/output/GMAP/gff.filter_alignm.bed.tmap data/output/gffcompare_second/")
        shell("mv data/output/GMAP/gff.filter_alignm.bed.refmap data/output/gffcompare_second/")
   
  rule true_lncRNA:
        input:
           "data/output/gffcompare_second/gff.filter_alignm.bed.tmap",
           "data/output/GMAP/filter_alignm.bed"
        params:
           expand("{option}", option=config['TE_coords']['option'])
        output:  
           "data/output/new_lncRNA/True_lncRNA.bed"
        shell:
           "python scripts/bed_intersect.py {input} {params} {output}"
           
rule statistic:
    input:
        "data/output/new_lncRNA/True_lncRNA.bed",
        "data/output/GMAP/reference.bed",
        "data/output/gffcompare_first/gffcmp.filter_alignm.bed.tmap"
    output:
        "data/output/lncRNA_structure/statistic_bed.tsv"
    shell:
        "python scripts/gff.py {input}"

rule blast:
    input:
       "data/output/new_lncRNA/True_lncRNA.bed",
       "data/output/gffcompare_first/gffcmp.filter_alignm.bed.tmap",
       "data/reference/data_index/LncAPDB.fasta",
       expand("{sequence}", sequence=config['sequence'])
    params:
       evalue=expand("{evalue}", evalue=config['blast']['evalue']),
       max_target=expand("{max_target}", max_target=config['blast']['max_target']),
       identity=expand("{identity}", identity=config['blast']['identity']),
       threads=expand("{threads}", threads=config['blast']['threads']),
       old_new="data/input/new_vs_old_id.tsv"
    output:
       fasta="data/output/new_lncRNA/new_lncrna.fasta",
       blast="data/output/new_lncRNA/blast.outfmt6"
    conda:
       "env/alignment.yaml"
    shell:
       "python scripts/blast.py {input} {output.fasta} {params.evalue} {params.max_target} {params.identity} {params.threads} {params.old_new}"

rule tissue:
    input:
       blast="data/output/new_lncRNA/blast.outfmt6",
       lnc="data/output/lncRNA_prediction/lncFinder.csv",
       new="data/output/new_lncRNA/new_lncrna.fasta"
    params:
       org=expand("{org}", org=config['tissue']['organism']),
       exp=expand("{exp}", exp=config['tissue']['exp']),
       LncAPDB=expand("{LncAPDB}", LncAPDB=config['tissue']['LncAPDB']),
       old_new="data/input/new_vs_old_id.tsv"
    output:
       "data/output/tissue/tissue_org.png"
    run:
       shell("""grep -v NonCoding {input.lnc} | awk -F ',' '{{print $1}}' | tr -d '"' > data/output/tissue/cod.txt""")
       shell("grep -v  {params.org} {input.blast} | awk -F '\t' '{{print $1}}' | sort | uniq > data/output/tissue/cons.txt")
       shell("grep {params.org} {input.blast} | awk -F '\t' '{{print $1}}' | sort | uniq > data/output/tissue/non.txt")
       shell("grep '>' {input.new} | awk -F '>' '{{print $2}}' > data/output/tissue/id.txt")
       shell("awk -F '\t' '{{print $1}}' {input.blast} | sort | uniq > data/output/tissue/all.txt")
       shell("grep -v -w -f data/output/tissue/all.txt data/output/tissue/id.txt > data/output/tissue/non.txt")
       shell("python scripts/tissue.py {params.exp} {input.blast} {params.LncAPDB} {params.old_new}")
     #  shell("rm data/output/tissue/cod.txt data/output/tissue/cons.txt data/output/tissue/non.txt")
     
     

