#!/bin/sh
computeMatrix scale-regions \
-R  /scratch/Narjesse/scripts/H3K4me1/profileplot_DMR_enhancer_overlapCORRECTED_binsizeCORRECTED/H3_vs_K4M_up_actv_enhancer_dmrs_corrected_overlap.bed \
/scratch/Narjesse/scripts/H3K4me1/profileplot_DMR_enhancer_overlapCORRECTED_binsizeCORRECTED/H3_vs_K4M_up_prim_enhancer_dmrs_corrected_overlap.bed \
/scratch/Narjesse/scripts/H3K4me1/profileplot_DMR_enhancer_overlapCORRECTED_binsizeCORRECTED/H3_vs_K4M_up_inac_enhancer_dmrs_corrected_overlap.bed \
-S /data/bigWigs/E14_EpiD2_WT_H3_doxPos_WGBS_AD_067_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/bigWigs/E14_EpiD2_WT_K4M_doxPos_WGBS_AD_068_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
--samplesLabel D2_WT_H3 D2_WT_K4M \
-bs 100 \
-m 2000 -b 1000 -a 1000 \
-bl /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
-p 6 \
--outFileName computematrix_DMR_enhancer_overlapCORRECTED_binsizeCORRECTED_5mc_up.gz