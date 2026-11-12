#!/bin/sh
computeMatrix scale-regions \
-R /scratch/Narjesse/scripts/H3K4me1/profileplot_DMR_enhancer_overlapCORRECTED_binsizeCORRECTED/H3_vs_K4M_up_actv_enhancer_dmrs_corrected_overlap.bed \
/scratch/Narjesse/scripts/H3K4me1/profileplot_DMR_enhancer_overlapCORRECTED_binsizeCORRECTED/H3_vs_K4M_up_prim_enhancer_dmrs_corrected_overlap.bed \
/scratch/Narjesse/scripts/H3K4me1/profileplot_DMR_enhancer_overlapCORRECTED_binsizeCORRECTED/H3_vs_K4M_up_inac_enhancer_dmrs_corrected_overlap.bed \
-S /scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_EpiD2_TKO_H3_doxPos_CnT_H3K4me1_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_EpiD2_TKO_K4M_doxPos_CnT_H3K4me1_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_EpiD2_WT_H3_doxPos_CnT_H3K4me1_q10_merged_norm.bw \
/scratch/Narjesse/scripts/UCSC_convertedBigWig/sorted_E14_EpiD2_WT_K4M_doxPos_CnT_H3K4me1_q10_merged_norm.bw \
--samplesLabel D2_TKO_H3 D2_TKO_K4M D2_WT_H3 D2_WT_K4M \
-m 2000 -b 1000 -a 1000 \
--missingDataAsZero \
--skipZeros \
-bl /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
-p 6 \
--outFileName /scratch/Narjesse/scripts/H3K4me1/profileplot_DMR_enhancer_overlapCORRECTED_binsizeCORRECTED/computematrix_DMR_enhancer_overlapCORRECTED_binsizeCORRECTED_up.gz