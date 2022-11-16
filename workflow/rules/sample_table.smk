

def get_sample_table_file(file ="samples.tsv"):

    sample_table_file = Path(file)
    if not sample_table_file.exists():

        logger.error(f"Sample table {sample_table_file} doesn't exist")
        raise Exception(f"Sample table {sample_table_file} doesn't exist")
    
    return sample_table_file


SampleTable = None

def get_sample_table():
    """Return the sample table as a pandas dataframe"""

    
    import pandas as pd

    global SampleTable
    if SampleTable is None:
        sample_table_file = get_sample_table_file()
        SampleTable = pd.read_csv(sample_table_file, sep="\t", index_col=0)

    return SampleTable


def get_all_samples():
    """Return a list of all samples"""
    SampleTable = get_sample_table()
    return SampleTable.index.tolist()


def get_reads(wildcards):
    """Return a list of all fastq files for a given sample"""
    SampleTable = get_sample_table()

    if wildcards.sample not in SampleTable.index:
        raise ValueError(f"Sample {wildcards.sample} not found in sample table")
    return SampleTable.loc[wildcards.sample, ["Reads_QC_R1", "Reads_QC_R2"]].tolist()
