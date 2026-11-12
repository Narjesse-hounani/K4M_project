#!/bin/sh
computeMatrix scale-regions \
-R mm10_TSS_TES_plus.bed3  \
-S /scratch/Narjesse/scripts/H3K27me3/heatmap_promoters_merged_replicates/convertedBW/sorted_E14_EpiD2_TKO_H3_doxPos_CnT_H3K27me3_KpDup_q10_b1_s1_LFAP_spikeIn_norm.bw \
/scratch/Narjesse/scripts/H3K27me3/heatmap_promoters_merged_replicates/convertedBW/sorted_E14_EpiD2_TKO_K4M_doxPos_CnT_H3K27me3_KpDup_q10_b1_s1_LFAP_spikeIn_norm.bw \
/scratch/Narjesse/scripts/H3K27me3/heatmap_promoters_merged_replicates/convertedBW/sorted_E14_EpiD2_WT_H3_doxPos_CnT_H3K27me3_KpDup_q10_b1_s1_LFAP_spikeIn_norm.bw \
/scratch/Narjesse/scripts/H3K27me3/heatmap_promoters_merged_replicates/convertedBW/sorted_E14_EpiD2_WT_K4M_doxPos_CnT_H3K27me3_KpDup_q10_b1_s1_LFAP_spikeIn_norm.bw \
--samplesLabel D2_TKO_H3 D2_TKO_K4M D2_WT_H3 D2_WT_K4M \
-bs 100 \
-m 2000 -b 1000 -a 1000 \
--missingDataAsZero \
--skipZeros \
-bl /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
-p 6 \
--outFileName computematrix_heatmap_TSS_TES_allgeno_D2_strand+_H3K27me3.gz