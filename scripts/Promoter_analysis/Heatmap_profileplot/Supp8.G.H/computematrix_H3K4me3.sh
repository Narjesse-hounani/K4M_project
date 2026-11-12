#!/bin/sh
sbatch -n 6 \
computeMatrix scale-regions \
-R /home/nhounani/DMRs_down_intersect_CpGprom_mm10_TSS_1kb.bed \
/home/nhounani/DMRs_down_intersect_nonCpGprom_mm10_TSS_1kb.bed \
-S /scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_2iD0_TKO_H3_doxPos_CnT_H3K4me3_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_2iD0_TKO_K4M_doxPos_CnT_H3K4me3_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_2iD0_WT_H3_doxPos_CnT_H3K4me3_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_2iD0_WT_K4M_doxPos_CnT_H3K4me3_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_EpiD2_TKO_H3_doxPos_CnT_H3K4me3_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_EpiD2_TKO_K4M_doxPos_CnT_H3K4me3_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_EpiD2_WT_H3_doxPos_CnT_H3K4me3_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_EpiD2_WT_K4M_doxPos_CnT_H3K4me3_q10_merged_norm.bw \
--samplesLabel D0_TKO_H3 D0_TKO_K4M D0_WT_H3 D0_WT_K4M D2_TKO_H3 D2_TKO_K4M D2_WT_H3 D2_WT_K4M \
-m 2000 -b 1000 -a 1000 \
--missingDataAsZero \
--skipZeros \
-bl /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
-p 6 \
--outFileName computematrix_DMRs_down_cginoncgi_allgeno_H3K4me3_v2.gz