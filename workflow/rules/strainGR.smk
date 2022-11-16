
### Quantify


rule prepare_fasta:
    input:
        strain_list = expand("Intermediate_files/StrainGST/{sample}.strains.tsv",sample=get_all_samples()),
        similarities= rules.kmersim.output,
        fna_dir= rules.preprare_strainge_db.output.dir,
        metadata_file = rules.preprare_strainge_db.output.meta,
    output:
        "Intermediate_files/concatenated_detected_strains.fasta"
    log:
        "log/strainGR/prepare_fasta_for_mapping.log"
    shell:
        "straingr prepare-ref "
        " -s {input.strain_list} "
        ' -p "{input.fna_dir}/{{ref}}.fa.gz" ' 
        " -S {input.similarities} "
        " -o {output} "
        " &> {log} "


"""
straingr prepare-ref 

  Options to change clustering behaviour.

  -S SIMILARITIES, --similarities SIMILARITIES
                        Enable clustering of closely related reference genomes
                        by specifying the path to the k-mer similarity scores
                        as created at the StrainGST database construction
                        step.
  -t THRESHOLD, --threshold THRESHOLD
                        K-mer clustering threshold, the default (0.7) is a bit
                        more lenient than the clustering step for database
                        construction, because for a concatenated reference
                        you'll want the included references not too closely
                        related, due to increased shared content.
"""


rule minimap_index:
    input:
        target=rules.prepare_fasta.output,
    output:
        "Intermediate_files/StrainGR/detected_strains.mmi",
    log:
        "log/alignment/index.log",
    #params:
        #index_size="33G",
    threads: 6
    wrapper:
        "v1.18.3/bio/minimap2/index"



rule map:
    input:
        target=rules.minimap_index.output,
        query=get_reads,
    output:
        "Intermediate_files/StrainGR/bams/{sample}.bam",
    log:
        "log/alignment/minimap/{sample}.log",
    params:
        extra="-x sr",
        sort="coordinate",
    threads: config["threads"]
    # resources:
    #     mem=config["mem"],
    #     mem_mb=config["mem"] * 1000,
    wrapper:
        "v1.18.3/bio/minimap2/aligner"



rule samtools_index:
    input:
        "Intermediate_files/StrainGR/bams/{sample}.bam",
    output:
        "Intermediate_files/StrainGR/bams/{sample}.bam.bai",
    log:
        "log/alignment/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v1.18.3/bio/samtools/index"



rule straingr_call:
    input:
        fasta= rules.prepare_fasta.output,
        bam = "Intermediate_files/StrainGR/bams/{sample}.bam", #rules.map.output,
        bai = "Intermediate_files/StrainGR/bams/{sample}.bam.bai" #rules.samtools_index.output
    output:
        hdf="Intermediate_files/StrainGR/hdf/{sample}.hdf5",
        summary="Intermediate_files/StrainGR/summary/{sample}.tsv"
    log:
        "log/StrainGR/call/{sample}.log"
    shell:
        "straingr call "
        " {input.fasta} {input.bam} "
        " --hdf5-out {output.hdf} "
        " --summary {output.summary} "
        " --tracks all "
        " &> {log}"


checkpoint all_strainGR:
    input:
        expand("Intermediate_files/StrainGR/summary/{sample}.tsv", sample= get_all_samples()),
    output:
        touch("Intermediate_files/StrainGR/all_finished")





# compare two samples
rule straingr_compare:
    input:
        sampleA="Intermediate_files/StrainGR/hdf/{sampleA}.hdf5",
        sampleB="Intermediate_files/StrainGR/hdf/{sampleB}.hdf5"
    output:
        summary="Intermediate_files/StrainGR/comparisons/summary/{sampleA}_{sampleB}_summary.tsv",
        details = "Intermediate_files/StrainGR/comparisons/details/{sampleA}_{sampleB}_details.tsv"
    log:
        "log/StrainGR/compare/{sampleA}_{sampleB}.log"
    threads:
        1
    resources:
        mem_mb=40 *1000
    shell:
        "straingr compare "
        " {input} "
        " -o {output.summary} -d {output.details} "
        " &> {log}"




def get_all_comparisons(wildcards):
    
    Samples = get_all_samples()
    n_samples= len(Samples)

    output_files = []
    for i in range(n_samples):
        for j in range(i+1,n_samples):

            output_files.append( rules.straingr_compare.output.summary.format(sampleA= Samples[i], sampleB=Samples[j]) )

    return output_files


def get_flag(wildcards):
    return checkpoints.all_strainGR.get().output[0]

rule combine_strainGR:
    input:
        flag=get_flag ,
        comparisons= get_all_comparisons
    output:
        touch("StrainGR_compare_finished")

"""

rule straingr_crate_compare_command:
    input:
        expand("Intermediate_files/StrainGR/hdf/{sample}.hdf5", sample= get_all_samples())
    params:
        dir="Intermediate_files/StrainGR/comparisons"
    output:
        "Intermediate_files/StrainGR/compare.sh"
    threads: 1
    log:
        "log/StrainGR/create_command.log"
    shell:
        "mkdir {params.dir} ; "
        " straingr compare "
        " {input} "
        " --all-vs-all"
        " -D {params.dir} "
        " > {output} "
        " 2> {log}"


rule straingr_command:
    input:
        expand("Intermediate_files/StrainGR/hdf/{sample}.hdf5", sample= get_all_samples()),
        command= rules.straingr_crate_compare_command.output[0]
    output:
        touch("Intermediate_files/StrainGR_compare_finished")
    params:
        dir = rules.straingr_crate_compare_command.params.dir
    threads: 1
    log:
        "log/StrainGR/compare.log"
    shell:
        "mkdir {params.dir} ; "
        " bash {input.command} &> {log}"

"""