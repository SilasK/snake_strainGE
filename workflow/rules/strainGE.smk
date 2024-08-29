

rule download_genomes:
    output:
        directory("ref_genomes/{genus}")
    conda:
        "../envs/ncbi_download.yaml"
    log:
        "log/{genus}/download_ncbi.log"
    threads:
        4
    shell:
        " mkdir -p {output} ; "
        " ncbi-genome-download bacteria -l complete -g {wildcards.genus} -H -F all -o {output} -P -p {threads} &> {log}"


## Build DB


checkpoint preprare_strainge_db:
    input:
        rules.download_genomes.output
    output:
        dir=directory("strainge_db/{genus}/fasta"),
        meta = "strainge_db/{genus}/fasta/references_meta.tsv"
    params:
        split_plasmids= True
    log:
        "log/{genus}/preprare_strainge_db.log"
    script:
        "../scripts/prepare_strainge_db.py"
        
rule kmerize_genome:
    input:
        "strainge_db/{genus}/fasta/{genome}.fa.gz"
    output:
        "strainge_db/{genus}/hdf5/{genome}.hdf5"
    log:
        "log/{genus}/kmerize/{genome}.log"
    resources:
        mem_mb=1000,
        time_min=20
    threads:
        1
    conda:
        "../envs/strainge.yaml"
    shell:
        "straingst kmerize -o {output} {input} &> {log}"

def kmersim_input(wildcards):

    db_dir = Path(checkpoints.preprare_strainge_db.get(**wildcards).output.dir)
    genomes = glob_wildcards( db_dir /"{genome}.fa.gz").genome
    assert len(genomes) >0, f"No fa.gz files found in {db_dir}"

    return expand(rules.kmerize_genome.output,genome=genomes,**wildcards)    

rule kmersim:
    input:
        kmersim_input
    output:
        "strainge_db/{genus}/similarities.tsv"
    log:
        "log/{genus}/kmersim.log"
    threads:
        8
    conda:
        "../envs/strainge.yaml"
    shell:
        "straingst kmersim "
        " --all-vs-all -t {threads} "
        " -S jaccard -S subset "
        " {input} > {output} 2> {log}"

rule cluster:
    input:
        similarities= rules.kmersim.output,
        hdf= kmersim_input
    output:
        "strainge_db/{genus}/references_to_keep.txt"
    log:
        "log/{genus}/cluster.log"
    conda:
        "../envs/strainge.yaml"
    shell:
        "straingst cluster -i {input.similarities} "
        " -d -C 0.99 -c 0.90 "
        " --clusters-out clusters.tsv "
        " {input.hdf} > {output} 2> {log}"

rule createdb:
    input:
        rules.cluster.output
    output:
        protected("{genus}-db.hdf5")
    log:
        "log/{genus}/createdb.log"
    conda:
        "../envs/strainge.yaml"
    shell:
        "straingst createdb -f {input} -o {output} 2> {log}"



### Run Strain 

rule kmerize_sample:
    input:
        get_reads
    output:
        "samples/{sample}.hdf5"
    log:
        "log/kmerize_sample/{sample}.log"
    conda:
        "../envs/strainge.yaml"
    shell:
        "straingst kmerize -k 23 "
        " -o {output} "
        " {input} "
        " &> {log} "

rule run_straingst:
    input:
        db = "{genus}-db.hdf5",
        sample = rules.kmerize_sample.output
    output:
        "Intermediate_files/{genus}/StrainGST/{sample}.stats.tsv",
        "Intermediate_files/{genus}/StrainGST/{sample}.strains.tsv"
    params:
        out_prefix= lambda wc,output: output[0].replace(".stats.tsv","")
    log:
        "log/{genus}/strainGST/{sample}.log"
    conda:
        "../envs/strainge.yaml"
    shell:
        "straingst run "
        " --separate-output -o {params.out_prefix} "
        " {input.db} {input.sample} "
        " 2> {log}"
    
rule combine_strainGST:
    input:
        expand("Intermediate_files/{{genus}}/StrainGST/{sample}.strains.tsv",sample=get_all_samples())
    output:
        "StrainGST_{genus}.tsv"
    params:
        sample_names = get_all_samples()
    run:

        import pandas as pd


        df_list=[]
        for f in input:
            df = pd.read_csv(f, sep='\t', index_col=1)
            df_list.append(df)



        # Combine all StrainGST results from each sample into a single DataFrame.
        straingst_df = pd.concat(df_list, keys=params.sample_names, names=["Sample"])

        straingst_df.to_csv(output[0],sep='\t')

