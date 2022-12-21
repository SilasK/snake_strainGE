
### Quantify


rule prepare_fasta:
    input:
        strain_list = expand("Intermediate_files/{{genus}}/StrainGST/{sample}.strains.tsv",sample=get_all_samples()),
        similarities= rules.kmersim.output,
        fna_dir= rules.preprare_strainge_db.output.dir,
        metadata_file = rules.preprare_strainge_db.output.meta,
    output:
        "Intermediate_files/{genus}/concatenated_detected_strains.fasta"
    log:
        "log/{genus}/strainGR/prepare_fasta_for_mapping.log"
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



rule bwa_mem2_index:
    input:
        rules.prepare_fasta.output
    output:
        multiext("Intermediate_files/{genus}/StrainGR/detected_strains",".0123",
        ".amb",
        ".ann",
        ".bwt.2bit.64",
        ".pac",)
    log:
        "log/{genus}/alignment/bwa_index.log",
    wrapper:
        "v1.19.1/bio/bwa-mem2/index"


rule map:
    input:
        reads=get_reads,
        idx=rules.bwa_mem2_index.output,
    output:
        "Intermediate_files/{genus}/StrainGR/bams/{sample}.bam",
    log:
        "log/{genus}/alignment/minimap/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}' -I 300 ", # define insert size
        sort="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",  # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "v1.19.1/bio/bwa-mem2/mem"



rule samtools_index:
    input:
        "Intermediate_files/{genus}/StrainGR/bams/{sample}.bam",
    output:
        "Intermediate_files/{genus}/StrainGR/bams/{sample}.bam.bai",
    log:
        "log/{genus}/alignment/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v1.18.3/bio/samtools/index"



rule straingr_call:
    input:
        fasta= rules.prepare_fasta.output,
        bam = "Intermediate_files/{genus}/StrainGR/bams/{sample}.bam", #rules.map.output,
        bai = "Intermediate_files/{genus}/StrainGR/bams/{sample}.bam.bai" #rules.samtools_index.output
    output:
        hdf="Intermediate_files/{genus}/StrainGR/hdf/{sample}.hdf5",
        summary="Intermediate_files/{genus}/StrainGR/summary/{sample}.tsv"
    log:
        "log/{genus}/StrainGR/call/{sample}.log"
    shell:
        "straingr call "
        " {input.fasta} {input.bam} "
        " --hdf5-out {output.hdf} "
        " --summary {output.summary} "
        " --tracks all "
        " &> {log}"


checkpoint all_strainGR:
    input:
        expand("Intermediate_files/{{genus}}/StrainGR/summary/{sample}.tsv", sample= get_all_samples()),
    output:
        touch("Intermediate_files/{genus}/StrainGR/all_finished")





# compare two samples
rule straingr_compare:
    input:
        sampleA="Intermediate_files/{genus}/StrainGR/hdf/{sampleA}.hdf5",
        sampleB="Intermediate_files/{genus}/StrainGR/hdf/{sampleB}.hdf5"
    output:
        summary="Intermediate_files/{genus}/StrainGR/comparisons/summary/{sampleA}_{sampleB}_summary.tsv",
        details = "Intermediate_files/{genus}/StrainGR/comparisons/details/{sampleA}_{sampleB}_details.tsv"
    log:
        "log/{genus}/StrainGR/compare/{sampleA}_{sampleB}.log"
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
        touch("StrainGR_compare_{genus}_finished")

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