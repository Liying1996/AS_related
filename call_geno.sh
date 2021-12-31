# Calling Gennotype
# Input file: Merged bam files of one donor - input.bam
# Out file: heterozygous SNPs - input_snp_filter.vcf

samtools mpileup -go chr1_samtools.bcf -f ~/GRCh38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta  -t AD -t DP -t SP -Q 20 -I -d 10000 -E  input.bam
bcftools call -vmO z -o input_bcftools_raw.vcf.gz input_samtools.bcf

# Select the reads: QUAL > 50；DP > 15； minor allele > 3 ; minor/major allele > 0.05
bcftools filter -O v -o  input_bcftools_filter.vcf -i 'QUAL> 50 & FMT/DP > 15 & MIN(FMT/AD[*]) > 3 & (MIN(FMT/AD[*])/MAX(FMT/AD[*]) > 0.05)' input_bcftools_raw.vcf.gz

bcftools view -v snps -g het -O v input_bcftools_filter.vcf >  input_snp_filter.vcf # only heterozygous site