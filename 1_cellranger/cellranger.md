# Cellranger

Run per sample in directory containing folder with fastq files.

```
cellranger=~/software/cellranger-3.1.0/cellranger
reference=~/references/refdata-cellranger-mm10-3.0.0/

$cellranger count --id=${PWD##*/} --fastqs=fastq --transcriptome=$reference
```