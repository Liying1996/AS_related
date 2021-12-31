while getopts 'w:i:1:2:h:o:' arg 
do
        case $arg in
             w)
                wasp_path=$OPTARG 
                ;;
             i)
                bwa_index=$OPTARG
                ;;
             1)
                fastq1=$OPTARG
                ;;
             2)
                fastq2=$OPTARG
                ;;
             h)
				h5_dir=$OPTARG
                ;;
             o)
				output_dir=$OPTARG
                ;;
             ?)  
	            echo "Unkonw argument!"
	        	exit 1
	        	;;
        esac
done

# BWA mapping
echo 'Start: BWA!'

name=`echo $fastq1 | awk -F '/' '{print $NF}' | sed "s/.\(fq\|fastq\)[0-9]*\(\.gz\)*//g"`

mkdir ${output_dir}/map/
bwa aln -t 8 $bwa_index $fastq1  > ${output_dir}/map/${name}.aln1.sai
bwa aln -t 8 $bwa_index $fastq2  > ${output_dir}/map/${name}.aln2.sai
bwa sampe $bwa_index ${output_dir}/map/${name}.aln1.sai ${output_dir}/map/${name}.aln2.sai ${fastq1} ${fastq2}> ${output_dir}/map/${name}.sam

samtools view -@ 8 -bS ${output_dir}/map/${name}.sam > ${output_dir}/map/${name}.tmp.bam
samtools sort -@ 8 ${output_dir}/map/${name}.tmp.bam -o ${output_dir}/map/${name}.bam
samtools index ${output_dir}/map/${name}.bam
rm ${output_dir}/map/${name}.tmp.bam

# WASP remapping
echo 'WASP remapping...'
mkdir  ${output_dir}/find_intersecting_snps/  ${output_dir}/remap/  ${output_dir}/filter_remapped_reads/

python ${wasp_path}/mapping/find_intersecting_snps.py \
     --is_sorted \
     --is_paired_end \
     --output_dir ${output_dir}/find_intersecting_snps \
     --snp_tab ${h5_dir}/snp_tab.h5 \
     --snp_index ${h5_dir}/snp_index.h5 \
     --haplotype ${h5_dir}/haplotypes.h5\
     ${output_dir}/map/${name}.bam

bwa aln -t 8 $bwa_index ${output_dir}/find_intersecting_snps/${name}.remap.fq1.gz  > ${output_dir}/remap/${name}.aln1.sai

bwa aln -t 8 $bwa_index ${output_dir}/find_intersecting_snps/${name}.remap.fq2.gz  > ${output_dir}/remap/${name}.aln2.sai

bwa sampe ${bwa_index} ${output_dir}/remap/${name}.aln1.sai ${output_dir}/remap/${name}.aln2.sai  ${output_dir}/find_intersecting_snps/${name}.remap.fq1.gz  ${output_dir}/find_intersecting_snps/${name}.remap.fq2.gz > ${output_dir}/remap/${name}.remap.sam

samtools view -@ 8 -bS ${output_dir}/remap/${name}.remap.sam > ${output_dir}/remap/${name}.remap.bam

samtools sort -@ 8 ${output_dir}/remap/${name}.remap.bam -o ${output_dir}/remap/${name}.remap.bam

python ${wasp_path}/mapping/filter_remapped_reads.py \
    ${output_dir}/find_intersecting_snps/${name}.to.remap.bam \
    ${output_dir}/remap/${name}.remap.bam \
    ${output_dir}/filter_remapped_reads/${name}.keep.bam

samtools sort -@ 8 ${output_dir}/filter_remapped_reads/${name}.keep.bam -o ${output_dir}/filter_remapped_reads/${name}.keep.sort.bam

samtools index ${output_dir}/filter_remapped_reads/${name}.keep.sort.bam


# Picard remove duplicates
echo 'Picard ...'
mkdir ${output_dir}/dup/
picard MarkDuplicates I=${output_dir}/filter_remapped_reads/${name}.keep.sort.bam O=${output_dir}/dup/${name}.bam M=${output_dir}/dup/dup_metrics_${name}.txt REMOVE_DUPLICATES=true


# remove reads Q < 15
echo 'Remove Q < 15 reads ...'
samtools view -@ 8 -q 15 -bS ${output_dir}/dup/${name}.bam > ${output_dir}/dup/${name}_qual.bam

samtools sort -@ 8 ${output_dir}/dup/${name}_qual.bam > ${output_dir}/dup/${name}_qual_sort.bam

samtools index ${output_dir}/dup/${name}_qual_sort.bam



