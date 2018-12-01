#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe make 1
####### bam to bed
#
bedtools bamtobed -i ear_H3K27ac_1_dedun.bam > ear_H3K27ac_1_dedun.bed
#bedtools bamtobed -i ear_H3K27ac_2_dedun.bam > ear_H3K27ac_2_dedun.bed
perl /NAS4/lien/Chia-pet/ear_H3K27ac/delete_length_less_10kb_PETs.pl ear_H3K27ac_intra_PETs_dedun ear_H3K27ac_1_dedun.bed ear_H3K27ac_2_dedun.bed ear_H3K27ac_1_rmdup_deSL.bed ear_H3K27ac_2_rmdup_deSL.bed
perl /NAS4/lien/Chia-pet/ear_H3K27ac/get_PETs.pl ear_H3K27ac_1_rmdup_deSL.bed ear_H3K27ac_2_rmdup_deSL.bed ear_H3K27ac_intra_PETs_dedun_deSL ear_H3K27ac_inter_PETs_dedun_deSL
#######
awk '$1 != "chr0"' ear_H3K27ac_PETs_peak_peaks.bed > ear_H3K27ac_PETs_peak_peaks_rmchr0.bed
bedtools intersect -a ear_H3K27ac_1_rmdup_deSL.bed -b ear_H3K27ac_PETs_peak_peaks_rmchr0.bed -wo > ear_H3K27ac_1_bowtie_uniq_in_peak
bedtools intersect -a ear_H3K27ac_2_rmdup_deSL.bed -b ear_H3K27ac_PETs_peak_peaks_rmchr0.bed -wo > ear_H3K27ac_2_bowtie_uniq_in_peak
cat ear_H3K27ac_1_rmdup_deSL.bed ear_H3K27ac_2_rmdup_deSL.bed > ear_H3K27ac_rmdup_deSL.bed
#bedtools intersect -b ear_H3K27ac_rmdup_deSL.bed -a ear_H3K27ac_PETs_peak_peaks_rmchr0.bed -wo |awk '{print $1"\t"$2"\t"$3}' |sort -n |uniq -c > ear_H3K27ac_peak_with_PETs
#
#
######
perl /NAS4/lien/Chia-pet/ear_H3K27ac/get_interaction_peak.pl ear_H3K27ac_1_bowtie_uniq_in_peak ear_H3K27ac_2_bowtie_uniq_in_peak ear_H3K27ac_interaction_peak_with_PETs
perl sort_bedpe_based_pos.pl ear_H3K27ac_interaction_peak_with_PETs  - |sort  -k1,1 -k2,2n -k4,4 -k5,5n - |uniq -c > ear_H3K27ac_interaction_peak_with_PETs_uniq
#awk '{print $1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6}' ear_H3K27ac_interaction_peak_with_PETs |sort -n |uniq -c > ear_H3K27ac_PETs_in_peak_both_end
#
perl /NAS4/lien/Chia-pet/ear_H3K27ac/final_summary.pl ear_H3K27ac_PETs_in_peak_both_end ear_H3K27ac_interaction_peak_with_PETs_uniq ear_H3K27ac_interaction_peak_with_PETs_summary
perl /NAS4/lien/Chia-pet/ear_H3K27ac/combine_same_interacted_cluster.pl ear_H3K27ac_interaction_peak_with_PETs_summary ear_H3K27ac_interaction_peak_with_PETs_summary_v1
#R --vanilla --slave < cal_the_FDR_for_interaction_peaks_with_PETs_correct_renew_N.R
#awk '$1~/^[0-9]/' ear_H3K27ac_PETs_peak_peaks.bed  > ear_H3K27ac_PETs_peak_peaks_chr.bed
############get interaction with FDR<0.01 & pets>=3
awk '{print $3"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9"\t"$2"\t"$12}' ../ear_H3K27ac_interaction_peak_with_PETs_summary_with_FDR | awk '$7>=3 && $8<0.01' - | perl -e '$hear=<>;while(<>){print "$_";}' - > ear_H3K27ac_interaction_peak_3pets.bedpe
