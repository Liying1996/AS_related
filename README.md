# AS_related

---

There is a specific explanation at the top of each code file.


#### Introduction

get_AS.sh: The main pipeline;

remove_bias_SE.sh/remove_bias_PE.sh: WASP to remove mapping bias;

betabinomial.R: Use Chen's method to obtain ASB sites;

draw_cCRE.py: Draw pie charts of cCREs distribution;

draw_motif.R: Draw histograms of motif enrichment;

best_match.py: Select best Fimo outputs;

match_bins.py: Select Controls from non-ASB SNPs;

trans_alt.py: Convert the fasta with reference allele to fasta with alternative allele;


---

#### Usage and parameters

-s: Single-end reads, cannot specify option -s after specifying option -p;

-p: Paired-end reads, the 2 files are separated by a comma (','), cannot specify option -p after specifying option -s;

-w: The directory of WASP;

-i: The BWA index (Not provided in the example_data/);

-c: The ChromInfo file;

-h: The HDF5 files containing SNPs information converted from VCF, cannot specify option -h after specifying option -v;

-v: The VCF file containing SNPs, cannot specify option -v after specifying option -h;

-f: The peak file of ChIP-seq/DNase-seq data;

-j: The path of snpEff.jar;

-a: Species, include human and mouse. Human: hg38; mouse: mm10;

-o: Output folder;

-m: The motifs of a transcript factor (you can download from [here](https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.21.tgz));

-r: Reference fasta file.


For example, single-end reads:

```shell
path=~/ASBfinder/
bash get_AS.sh \
	-s ${path}/example_data/single_end/ENCFF000OCP.fastq.gz \
	-w ~/WASP/ \
	-i ~/hg38_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
	-c  ${path}/example_data/chromeinfo/GRCh38_EBV.chrom.sizes.tsv \
	-h ${path}/example_data/h5 \
	-f ${path}/example_data/ENCFF801BDJ.bed \
	-j ${path}/data/snpEFF/snpEff/snpEff.jar \
	-a human \
	-o ${path}/example_data/output/single/
	-m ${path}/example_data/CTCF.meme \
	-r ~/hg38_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
```

Paired-end reads:

```shell
path=~/ASBfinder/
bash get_AS.sh \
	-p ${path}/example_data/paired_end/ENCFF340SQP.fastq.gz,${path}/example_data/paired_end/ENCFF587OVW.fastq.gz \
	-w ~/WASP/ \
	-i ~/hg38_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
	-c  ${path}/example_data/chromeinfo/GRCh38_EBV.chrom.sizes.tsv \
	-h ${path}/example_data/h5 \
	-f ${path}/example_data/ENCFF801BDJ.bed \
	-j ${path}/data/snpEFF/snpEff/snpEff.jar \
	-a human \
	-o ${path}/example_data/output/paired/
	-m ${path}/example_data/CTCF.meme \
	-r ~/hg38_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
```

---

#### Dependencies

The pipeline requires WASP to remove the mapping bias, so it must meet the WASP requirements.

Python3 (Highly recommend Anaconda):

- numpy
- scipy
- pysam
- pytables
- matplotlib
- pandas
- argparse
- All [WASP](https://github.com/bmvdgeijn/WASP) needs

R:

- VGAM
- Unicode
- ggplot2


Softwares:

- fastp
- Picard
- tabix
- bedtools intersect (aka. intersectBed)
- snpEFF
- fimo

Conda can install these softare and packages directly:

```shell
conda install picard, fastp,tabix, bedtools
```
