###########step1--merged bam from the same repeat
#samtools merge ear_H3K27ac_1.bam EK-2_L1_I601804_1_cutadapt_fil_bowtie_uniq.srt.bam EK-2_L3_I601804_1_cutadapt_fil_bowtie_uniq.srt.bam EK-2_L5_601804_1_cutadapt_fil_bowtie_uniq.srt.bam 0724EK1_L1_I603806_1_cutadapt_fil_bowtie_uniq.srt.bam 0724EK1_L3_I603806_1_cutadapt_fil_bowtie_uniq.srt.bam 0724EK1_L6_603806_1_cutadapt_fil_bowtie_uniq.srt.bam
#samtools merge ear_H3K27ac_2.bam EK-2_L1_I601804_2_cutadapt_fil_bowtie_uniq.srt.bam EK-2_L3_I601804_2_cutadapt_fil_bowtie_uniq.srt.bam EK-2_L5_601804_2_cutadapt_fil_bowtie_uniq.srt.bam 0724EK1_L1_I603806_2_cutadapt_fil_bowtie_uniq.srt.bam 0724EK1_L3_I603806_2_cutadapt_fil_bowtie_uniq.srt.bam 0724EK1_L6_603806_2_cutadapt_fil_bowtie_uniq.srt.bam
##########step2--get pets- bamtobed
#samtools index ear_H3K27ac_1.bam
#samtools index ear_H3K27ac_2.bam
#bedtools bamtobed -i ear_H3K27ac_1.bam > ear_H3K27ac_1_bed
#bedtools bamtobed -i ear_H3K27ac_2.bam > ear_H3K27ac_2_bed
########step2-- get pets-call PETs
#perl /NAS4/lien/Chia-pet/ear_H3K27ac/get_PETs.pl ear_H3K27ac_1_bed ear_H3K27ac_2_bed ear_H3K27ac_intra_PETs ear_H3K27ac_inter_PETs
##########step2--rm redundency PETs
#perl /NAS4/lien/Chia-pet/ear_H3K27ac/de_redundency_PETs.pl ear_H3K27ac_intra_PETs ear_H3K27ac_intra_PETs_dedun
#perl /NAS4/lien/Chia-pet/ear_H3K27ac/de_redundency_PETs.pl ear_H3K27ac_inter_PETs ear_H3K27ac_inter_PETs_dedun
####### step3--get PET peaks--get input bam file
#samtools view -h ear_H3K27ac_1.bam |perl /NAS4/lien/Chia-pet/ear_H3K27ac/get_bam_for_call_peaks.pl ear_H3K27ac_intra_PETs_dedun ear_H3K27ac_inter_PETs_dedun - |samtools view -bS - > ear_H3K27ac_1_dedun.bam
#samtools view -h ear_H3K27ac_2.bam |perl /NAS4/lien/Chia-pet/ear_H3K27ac/get_bam_for_call_peaks.pl ear_H3K27ac_intra_PETs_dedun ear_H3K27ac_inter_PETs_dedun - |samtools view -bS - > ear_H3K27ac_2_dedun.bam
#samtools merge -f ear_H3K27ac_dedun_merge.bam ear_H3K27ac_1_dedun.bam ear_H3K27ac_2_dedun.bam
#samtools sort ear_H3K27ac_dedun_merge.bam -o ear_H3K27ac_dedun_merge.srt.bam 
samtools sort -n ear_H3K27ac_dedun_merge.bam -o ear_H3K27ac_dedun_merge.srt_name.bam
######step3--get PET peaks--MACS call peak
#python /NAS1/lien/software/MACS-1.4.2/bin/macs14 -t ear_H3K27ac_dedun_merge.srt.bam -f BAM -g 2.06e+9 -n ear_H3K27ac_PETs_peak --pvalue 1e-09 --nolambda --nomodel --keep-dup=2
