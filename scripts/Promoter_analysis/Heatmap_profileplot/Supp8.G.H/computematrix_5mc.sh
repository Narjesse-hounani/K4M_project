#!/bin/sh
sbatch -n 6 \
computeMatrix scale-regions \
-R /home/nhounani/DMRs_down_intersect_CpGprom_mm10_TSS_1kb.bed \
/home/nhounani/DMRs_down_intersect_nonCpGprom_mm10_TSS_1kb.bed \
-S /data/bigWigs/E14_EpiD2_WT_H3_doxPos_WGBS_AD_067_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
/data/bigWigs/E14_EpiD2_WT_K4M_doxPos_WGBS_AD_068_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5.bw \
--samplesLabel D2_WT_H3 D2_WT_K4M \
-bs 100 \
-m 2000 -b 1000 -a 1000 \
-bl /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
-p 6 \
--outFileName computematrix_DMRs_down_cginoncgi_allgeno_5mc_v2.gz