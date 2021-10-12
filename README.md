# NanoEM

Codes for analyzing nanoEM data.

## Reference
https://academic.oup.com/nar/article/49/14/e81/6279847?login=true

## Requirements

- Python3
- pysam
- minimap2
- sambamba
- samtools

## Usage

### Base conversion of reference genome

You should change "ref" to your reference genome. 

```
python src/convert_ref.py ref > output.fa 
```

### Base conversion of nanoEM data

OUTPUT:
- *_CT.fq.gz: C-to-T converted fastq
- *_GA.fq.gz: G-to-A converted fastq

```
python src/convert_reads.py fastq
```

### Mapping converted data to converted reference genome


```
minimap2 -t 8 --split-prefix temp_sam1 -ax map-ont output.fa  *_CT.fq.gz --eqx | samtools view -b | samtools sort -@ 8 -o 1.sorted.bam
samtools index 1.sorted.bam

minimap2 -t 8 --split-prefix temp_sam2 -ax map-ont output.fa  *_GA.fq.gz --eqx | samtools view -b | samtools sort -@ 8 -o 2.sorted.bam
samtools index 2.sorted.bam
```

### Choosing best alignments

```
python src/best_align.py --bam1 1.sorted.bam --bam2 2.sorted.bam  --fastq fastq
samtools view -b output_CT.sam | samtools sort -o output_CT.sorted.bam
samtools view -b output_GA.sam | samtools sort -o output_GA.sorted.bam
samtools index output_CT.sorted.bam
samtools index output_GA.sorted.bam
```

### Methylation calling

```
sambamba mpileup output_CT.sorted.bam -L cpg_sites.bed -o pileup_CT.tsv -t 8 --samtools -f ref
sambamba mpileup output_GA.sorted.bam -L cpg_sites.bed -o pileup_GA.tsv -t 8 --samtools -f ref
python src/call_methylation.py pileup_CT.tsv pileup_GA.tsv > frequency_methylation.tsv
```

### Visualization for bisulfite mode of IGV

```
python script/vis_GA_utilities.py -b output_GA.sorted.bam | samtools view -b | samtools sort -@ 4 -o output_GA_vis.sorted.bam
samtools index output_GA_vis.sorted.bam

samtools merge output_merge.bam output_CT.sorted.bam output_GA_vis.sorted.bam
samtools sort -@ 4 -o output_merge.sorted.bam output_merge.bam
samtools index output_merge.sorted.bam
```
