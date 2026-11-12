#!/bin/sh
computeMatrix scale-regions \
-R /scratch/Narjesse/scripts/H3K27me3/profileplot_all_enhancer/H3_vs_K4M_actv_enhancer_all.bed \
/scratch/Narjesse/scripts/H3K27me3/profileplot_all_enhancer/H3_vs_K4M_prim_enhancer_all.bed \
/scratch/Narjesse/scripts/H3K27me3/profileplot_all_enhancer/H3_vs_K4M_inac_enhancer_all.bed \
-S /scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_2iD0_TKO_H3_doxPos_CnT_H3K4me1_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_2iD0_TKO_K4M_doxPos_CnT_H3K4me1_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_2iD0_WT_H3_doxPos_CnT_H3K4me1_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_2iD0_WT_K4M_doxPos_CnT_H3K4me1_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_EpiD2_TKO_H3_doxPos_CnT_H3K4me1_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_EpiD2_TKO_K4M_doxPos_CnT_H3K4me1_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_EpiD2_WT_H3_doxPos_CnT_H3K4me1_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_EpiD2_WT_K4M_doxPos_CnT_H3K4me1_q10_merged_norm.bw \
--samplesLabel D0_TKO_H3 D0_TKO_K4M D0_WT_H3 D0_WT_K4M D2_TKO_H3 D2_TKO_K4M D2_WT_H3 D2_WT_K4M \
-bs 100 \
-m 2000 -b 1000 -a 1000 \
--missingDataAsZero \
--skipZeros \
-bl /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
-p 6 \
--outFileName computematrix_heatmap_all_enhancers_allgeno_H3K4me1.gz