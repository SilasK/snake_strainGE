

## Build DB


rule preprare_strainge_db:
    input:
        "ref_genomes/human_readable"
    output:
        dir=directory("strainge_db/fasta"),
        meta = "strainge_db/fasta/references_meta.tsv"
    params:
        split_plasmids= True
    log:
        "log/preprare_strainge_db.log"
    script:
        "scripts/prepare_strainge_db.py"
        
rule kmerize_genome:
    input:
        "strainge_db/fasta/{genome}.fa.gz"
    output:
        "strainge_db/hdf5/{genome}.hdf5"
    log:
        "log/kmerize/{genome}.log"
    resources:
        mem_mb=1000,
        time_min=20
    threads:
        1
    shell:
        "straingst kmerize -o {output} {input} &> {log}"

def kmersim_input(wildcards):

    #db_dir = Path(checkpoints.preprare_strainge_db.get().output.dir)
    db_dir = Path(rules.preprare_strainge_db.output.dir)
    genomes = glob_wildcards( db_dir /"{genome}.fa.gz").genome
    assert len(genomes) >0, f"No fa.gz files found in {db_dir}"

    return expand(rules.kmerize_genome.output,genome=genomes)    

rule kmersim:
    input:
        kmersim_input
    output:
        "strainge_db/similarities.tsv"
    log:
        "log/kmersim.log"
    threads:
        8
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
        "strainge_db/references_to_keep.txt"
    log:
        "log/cluster.log"
    shell:
        "straingst cluster -i {input.similarities} "
        " -d -C 0.99 -c 0.90 "
        " --clusters-out clusters.tsv "
        " {input.hdf} > {output} 2> {log}"

rule createdb:
    input:
        rules.cluster.output
    output:
        protected("pan-genome-db.hdf5")
    log:
        "log/createdb.log"
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
    shell:
        "straingst kmerize -k 23 "
        " -o {output} "
        " {input} "
        " &> {log} "

rule run_straingst:
    input:
        db = "pan-genome-db.hdf5",
        sample = rules.kmerize_sample.output
    output:
        "Intermediate_files/StrainGST/{sample}.stats.tsv",
        "Intermediate_files/StrainGST/{sample}.strains.tsv"
    log:
        "log/strainGST/{sample}.log"
    shell:
        "straingst run "
        " --separate-output -o Intermediate_files/StrainGST/{wildcards.sample}"
        " {input.db} {input.sample} "
        " 2> {log}"
    
rule combine_strainGST:
    input:
        expand("Intermediate_files/StrainGST/{sample}.strains.tsv",sample=get_all_samples())
    output:
        "StrainGST.tsv"
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

