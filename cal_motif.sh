# This file can get the number of overlapping motifs between ASB and control non-ASB and p.value.
# Input: [-input] The output file of beta-binomial.R; [-name] The name of sample; [-fasta] The path of reference fasta file; [-motif_file] .meme motif file of one TF; [-output] outoput dir
# Example:
# sh ./cal_motif.sh -i ./WASP/inpeak/ENCFF803YAB_counts_AS_0.1.txt -n ENCFF000VUW -f ~/GRCh38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -m ./CTCF.meme -o ./
# Dependence: best_match.py(or match_bins.py); draw_motif.R
# The first method ("best_match.py") to take controls is to randomly take the same number of nonASB sites as ASB SNPs
# The second method ("match_bins.py"): we randomly selected from the non-allelic specific sites (nonAS)  the same number of SNPs as the control SNPs, and required the quartiles of global MAF of nonAS SNPs,  distance from peak center to the nearest Transcription start sites (TSS) and read coverage of SNPs should be consistent with AS SNPs. MAF information of SNPs was collected from 1000 Genomes Project. The read coverage has a significant impact on the selection of control SNPs. Because the read coverage of control SNPs is significantly different from AS SNPs, so it needs to be considered when selecting control SNPs. Considering directly matching read depth is too stringent, we selected matched control SNPs in bins: first, we looked for exact bin match with bins of 1 for 1-19, bins of 5 from 20-49, bins of 10 from 50-199 and bins of 20 over 200 (eg. if AS is 18, the matched control is 18, if AS is 21, the matched control can be in 21-25...); If there is no control match in the bin, check the next highest bin. Iteratively do this until we hit a match; If there are still no matches, check the next lowest bin. Iteratively do this until we hit a match. 



while getopts 'i:n:f:m:o:' arg 
do
        case $arg in
             i)
                input=$OPTARG 
                ;;
             n)
                name=$OPTARG
                ;;
             f)
                fasta=$OPTARG
                ;;
             m)
                motif_file=$OPTARG
                ;;
             o)
                output=$OPTARG
                ;;
             ?)  
                echo "Unknown argument!"
                exit 1
                ;;
        esac
done

mkdir ${output}/motif/
awk -F '\t' '$NF=="significant"' ${input} > ${output}/motif/${name}_AS.txt 
awk -F '\t' '$NF!="significant"' ${input} > ${output}/motif/${name}_nonAS.txt 

awk -F '\t' '{print $1"\t"$2-21"\t"$2+20"\t"$2"\t"$3"\t"$4}' ${output}/motif/${name}_AS.txt | sed '1d' > ${output}/motif/${name}_AS_snps.bed

awk -F '\t' '{print $1"\t"$2-21"\t"$2+20"\t"$2"\t"$3"\t"$4}' ${output}/motif/${name}_nonAS.txt | sed '1d' > ${output}/motif/${name}_nonAS_snps.bed

cat ${output}/motif/${name}_AS_snps.bed ${output}/motif/${name}_nonAS_snps.bed > ${output}/motif/${name}_snps.bed

bedtools getfasta -fi ${fasta} -bed ${output}/motif/${name}_snps.bed > ${output}/motif/${name}_ref.fasta

curr_path=$(cd `dirname $0`; pwd)
python ${curr_path}/trans_alt.py -name ${name} -input_file ${output}/motif/${name}_snps.bed -fa ${output}/motif/${name}_ref.fasta -output ${output}/motif/

# bg and fg have different results 
fimo --thresh 1 --text ${motif_file} ${output}/motif/${name}_ref.fasta > ${output}/motif/${name}_ref_out.txt

fimo --thresh 1 --text ${motif_file} ${output}/motif/${name}_alt.fasta > ${output}/motif/${name}_alt_out.txt

# only select SNPs seq; 4th is start, 5th is end.
awk -F '[\t]+' '$3<=21 && $4>=21' ${output}/motif/${name}_ref_out.txt > ${output}/motif/${name}_out.txt
awk -F '[\t]+' '$3<=21 && $4>=21' ${output}/motif/${name}_alt_out.txt >> ${output}/motif/${name}_out.txt

python ${curr_path}/best_match.py -n ${name} -input_file ${output}/motif/${name}_out.txt -output ${output}/motif/

# select p.value <= 1e-4
awk -F '\t' '$7<0.0001' ${output}/motif/${name}_best.txt > ${output}/motif/${name}_thresh.txt

awk -F '[\t]+' 'FNR==NR {x[$2];next} ($1":"$2"-"$3 in x)' ${output}/motif/${name}_thresh.txt ${output}/motif/${name}_AS_snps.bed  > ${output}/motif/${name}_AS_inmotif.txt

AS_inmotif=`wc -l ${output}/motif/${name}_AS_inmotif.txt | cut -d ' ' -f1`

awk -F '[\t]+' 'FNR==NR {x[$2];next} ($1":"$2"-"$3 in x)' ${output}/motif/${name}_thresh.txt ${output}/motif/${name}_nonAS_snps.bed  > ${output}/motif/${name}_nonAS_inmotif.txt
nonAS_inmotif=`wc -l ${output}/motif/${name}_nonAS_inmotif.txt | cut -d ' ' -f1`

AS_total=`wc -l ${output}/motif/${name}_AS.txt | cut -d ' ' -f1`
nonAS_total=`wc -l ${output}/motif/${name}_nonAS.txt | cut -d ' ' -f1`

randLoci=0
for num in `seq 0 9999`;do
    shuf -n ${AS_total} ${output}/motif/${name}_nonAS.txt > ${output}/motif/${name}_control.txt  
    awk -F '[\t+]' 'FNR==NR {x[$1"\t"$4];next} ($1"\t"$2 in x)' ${output}/motif/${name}_nonAS_inmotif.txt ${output}/motif/${name}_control.txt > ${output}/motif/${name}_inmotif_control.txt
    control_inmotif=`wc -l ${output}/motif/${name}_inmotif_control.txt | cut -d ' ' -f1`
    echo $control_inmotif >> ${output}/motif/${name}_inmotif_control_count.txt

    if [ "$control_inmotif" -ge "$AS_inmotif" ];then
        randLoci=$(( $randLoci + 1 ))
    fi
done

rm ${output}/motif/${name}_control.txt  
rm ${output}/motif/${name}_inmotif_control.txt

mean=`awk '{sum+=$1} END {print sum/NR}' ${output}/motif/${name}_inmotif_control_count.txt`
sd=`awk '{x[NR]=$0; s+=$0; n++} END{a=s/n; for (i in x){ss += (x[i]-a)^2} sd = sqrt(ss/(n-1)); print sd}' ${output}/motif/${name}_inmotif_control_count.txt`
# pVal=$(echo "$randLoci 10000" | awk '{print $1/$2}')
Rscript ${curr_path}/draw_motif.R ${output}/motif/${name}_inmotif_control_count.txt ${AS_inmotif} ${name} ${output}/motif/${name}_dist.pdf   #generate a tmp file contains p.val

pVal=`cat tmp`
rm tmp 

echo 'name AS_inmotif AS_total nonAS_inmotif nonAS_total control_mean control_sd pVal' > ${output}/motif/${name}_results_inmotif.txt
echo $name $AS_inmotif $AS_total $nonAS_inmotif $nonAS_total $mean $sd $pVal >> ${output}/motif/${name}_results_inmotif.txt

mkdir ${output}/motif/temp_files/
for i in `ls ${output}/motif/${name}* | grep -v 'inmotif'`;do
    mv $i ${output}/motif/temp_files/
done
