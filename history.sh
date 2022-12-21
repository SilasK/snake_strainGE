
mamba install strainge
mamba create -n download_genome  ncbi-genome-download



mkdir ref_genomes
ncbi-genome-download bacteria -l complete -g Staphylococcus,Streptococcus -H -F all -o ref_genomes -P -p 4

wget https://github.com/broadinstitute/StrainGE/blob/master/bin/prepare_strainge_db.py


mkdir strainge_db
python3 prepare_strainge_db.py ref_genomes/human_readable -s \
    -o strainge_db > strainge_db/references_meta.tsv

for f in strainge_db/*.fa.gz; do straingst kmerize -o ${f%.fa.gz}.hdf5 $f; done;


mamba install -y mummer 


pysam error 
mamba install pysam=0.19