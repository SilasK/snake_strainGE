
# Prepare databases for StrainGe

## Create 2 conda environments
```
mamba create -n strainge strainge mummer pysam=0.19

mamba create -n download_genome  ncbi-genome-download
```

## Download genomes of all **complete** genomes of your genus or genera

Here is an example of all _Staphylococcus_ and _Streptococcus_

Use this in the `download_genome` environment.

``` 
mkdir ref_genomes
ncbi-genome-download bacteria -l complete -g Staphylococcus,Streptococcus -H -F all -o ref_genomes -P -p 4

```

## Convert to format that can be used by strain GE

```
wget https://github.com/broadinstitute/StrainGE/blob/master/bin/prepare_strainge_db.py

```

```
mkdir strainge_db
python3 prepare_strainge_db.py ref_genomes/human_readable -s \
    -o strainge_db > strainge_db/references_meta.tsv
```
## Kmerize with strainge

In the `strainge`environment. 
```
for f in strainge_db/*.fa.gz; do straingst kmerize -o ${f%.fa.gz}.hdf5 $f; done;
```

