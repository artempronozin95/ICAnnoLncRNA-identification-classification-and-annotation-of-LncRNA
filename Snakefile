# read config info into this namespace
configfile: "config.yaml"
print (config['lnc'], config['cds'], config['sequence'])


if config['diamond']['option'] == 'on':
 rule all_1:
    input:
        "data/input/test_train/sets.txt",
        "data/output/cpc.txt",
        "data/input/test_train/best_model.txt",
        "data/output/Noncoding.fasta",
        "data/output/lncFinder.csv",
        "data/output/statistic.csv",
        "data/output/Noncoding_trans_out.fasta",
        "data/output/ORF.orf",
        "data/input/test_train/dirs.txt",
        "data/output/tmhmm.csv",
        "data/output/alignm_filter.gff",
        "data/output/gffcmp.alignm_filter.gff.tmap",
        "data/output/filter_alignm.bed",
        "data/output/itron_coordin.tsv",
        "data/output/statistic_bed.tsv",
        "data/output/tissue/tissue_org.png",
        expand("{out}", out=config['diamond']['out']),
        blast="data/output/blast.outfmt6",
        query="data/output/alignm.bed"

if config['diamond']['option'] == 'off':
 rule all_2:
    input:
        "data/input/test_train/sets.txt",
        "data/output/cpc.txt",
        "data/output/Noncoding.fasta",
        "data/output/lncFinder.csv",
        "data/output/statistic.csv",
        "data/input/test_train/best_model.txt",
        "data/output/Noncoding_trans_out.fasta",
        "data/output/ORF.orf",
        "data/output/tmhmm.csv",
        "data/input/test_train/dirs.txt",
        "data/output/alignm_filter.gff",
        "data/output/gffcmp.alignm_filter.gff.tmap",
        "data/output/filter_alignm.bed",
        "data/output/itron_coordin.tsv",
        "data/output/statistic_bed.tsv",
        "data/output/tissue/tissue_org.png",
        blast="data/output/blast.outfmt6",
        query="data/output/alignm.bed"

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

rule lncFinder:
    input:
       "data/input/test_train/best_model.txt",
       expand("{sequence}", sequence=config['sequence']),
       ancient("data/output/")
    params:
       expand("{strc}", strc=config['structure'])
    output:
        "data/output/lncFinder.csv"
    conda:
        "env/lncfinder.yaml"
    shell:
        "scripts/Chunk_dataframe.sh {input} {params}"

rule cpc:
    input:
        expand("{sequence}", sequence=config['sequence'])
    output:
        "data/output/cpc.txt"
    conda:
        "env/cpc.yaml"
    shell:
        "CPC2/bin/CPC2.py -i {input} -r -o data/output/cpc"
        
rule take_noncoding:
    input:
        ancient("data/output/lncFinder.csv"),
        expand("{sequence}", sequence=config['sequence'])
    output:
        "data/output/Noncoding.fasta"
    shell:
        "python scripts/take_noncoding.py {input} {output}"
  
rule Gmap_build:
    input:
        expand("{reference}", reference=config['gmap']['reference']),
        expand("{gff_reference}", gff_reference=config['gmap']['gff_reference'])
    output:
       "data/output/statistic.csv"
    shell:
       "python scripts/GMAP_build.py {input}"
    
rule ORF:
    input:
        "data/output/Noncoding.fasta"
    output:
        "data/output/ORF.orf"
    shell:
        "transeq {input} -auto -stdout -frame R -table 0 > {output}"

rule tmhmm:
    input:
        "data/output/ORF.orf",
        ancient("data/output/")
    output:
        "data/output/tmhmm.csv"
    shell:
        "scripts/tmhmm.sh {input}"

rule filter:
    input:
        tmh="data/output/tmhmm.csv",
        fasta="data/output/Noncoding.fasta"
    output:
        "data/output/Noncoding_trans_out.fasta",
        trans="data/output/trans.csv",
        not_trans="data/output/not_trans.csv"
    run:
        shell("python scripts/frames_and_tmtmm.py {input.tmh} {output.trans} {output.not_trans}")
        shell("seqtk subseq {input.fasta} {output.not_trans} > {output}")

rule Gmap_aling: 
    input: 
        expand("{reference}", reference=config['gmap']['reference']),
        "data/output/Noncoding_trans_out.fasta",
        "data/output/statistic.csv"
    output:
        "data/output/alignm_filter.gff"
    shell:
        "python scripts/GMAP.py {input} {config[model]} {output}"
        
rule to_bed:
    input:
        query="data/output/alignm_filter.gff",
        reference=expand("{gff_reference}", gff_reference=config['gmap']['gff_reference'])
    output:
        query_out="data/output/alignm.bed",
        reference_out="data/output/reference.bed"
    run:
        shell("gff2bed < {input.query} > {output.query_out}")
        shell("gff2bed < {input.reference} > {output.reference_out}")

rule filter_long_introns:
    input:
       query_out="data/output/alignm.bed",
       reference_out="data/output/reference.bed"
    output:
       "data/output/filter_alignm.bed",
       "data/output/itron_coordin.tsv"
    shell:
       "python scripts/gmap_filter.py {input} {output}"

rule gffcomare:
    input:
        expand("{annot}", annot=config['gff']),
        "data/output/alignm_filter.gff"
    params:
        "data/output/gffcmp"
    output:
        "data/output/gffcmp.alignm_filter.gff.tmap"
    shell:
        "gffcompare -r {input} -o {params}"

rule statistic:
    input:
        "data/output/filter_alignm.bed",
        "data/output/reference.bed",
        "data/output/gffcmp.alignm_filter.gff.tmap"
    output:
        "data/output/statistic_bed.tsv"
    shell:
        "python scripts/gff.py {input}"

rule blast:
    input:
       "data/output/gffcmp.alignm_filter.gff.tmap",
       "data/reference/data_index/LncAPDB.fasta",
       expand("{sequence}", sequence=config['sequence'])
    params:
       evalue=expand("{evalue}", evalue=config['blast']['evalue']),
       max_target=expand("{max_target}", max_target=config['blast']['max_target']),
       identity=expand("{identity}", identity=config['blast']['identity']),
       threads=expand("{threads}", threads=config['blast']['threads'])
    output:
       fasta="data/output/new_lncrna.fasta",
       blast="data/output/blast.outfmt6"
    conda:
       "env/alignment.yaml"
    shell:
       "python scripts/blast.py {input} {output.fasta} {params.evalue} {params.max_target} {params.identity} {params.threads}"

rule tissue:
    input:
       blast="data/output/blast.outfmt6",
       lnc="data/output/lncFinder.csv",
       new="data/output/new_lncrna.fasta"
    params:
       org=expand("{org}", org=config['tissue']['organism']),
       exp=expand("{exp}", exp=config['tissue']['exp']),
       LncAPDB=expand("{LncAPDB}", LncAPDB=config['tissue']['LncAPDB'])
    output:
       "data/output/tissue/tissue_org.png"
    run:
       shell("""grep -v NonCoding {input.lnc} | awk -F ',' '{{print $1}}' | tr -d '"' > data/output/tissue/cod.txt""")
       shell("grep -v  {params.org} {input.blast} | awk -F '\t' '{{print $1}}' | sort | uniq > data/output/tissue/cons.txt")
       shell("grep {params.org} {input.blast} | awk -F '\t' '{{print $1}}' | sort | uniq > data/output/tissue/non.txt")
       shell("grep '>' {input.new} | awk -F '>' '{{print $2}}' > data/output/tissue/id.txt")
       shell("awk -F '\t' '{{print $1}}' {input.blast} | sort | uniq > data/output/tissue/all.txt")
       shell("grep -v -w -f data/output/tissue/all.txt data/output/tissue/id.txt > data/output/tissue/non.txt")
       shell("python scripts/tissue.py {params.exp} {input.blast} {params.LncAPDB}")
     #  shell("rm data/output/tissue/cod.txt data/output/tissue/cons.txt data/output/tissue/non.txt")

if config['diamond']['option'] == 'on':
  rule diamond:
     input:
        expand("{database}", database=config['diamond']['database']),
        expand("{query}", query=config['diamond']['query'])
     output:
        expand("{out}", out=config['diamond']['out'])
     conda:
       "env/alignment.yaml"
     shell:
        "python scripts/diamond.py {input}"
elif config['diamond']['option'] == 'off':
     pass
