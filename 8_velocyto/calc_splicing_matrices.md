# Compute splicing matrices
This is directly run on the 10x cellranger output folder of each sample.

Example for single sample:
```
path=~/users/dk12
sampleID=HAM13007

velocyto run10x --samtools-threads 20 \
                --samtools-memory 2000 \
                -t uint32 \
                -m $path/references/velocyto/mm10_rmsk.gtf \
                $path/brain-regeneration/$sampleID \
                $path/references/refdata-cellranger-mm10-3.0.0/genes/genes.gtf
```
