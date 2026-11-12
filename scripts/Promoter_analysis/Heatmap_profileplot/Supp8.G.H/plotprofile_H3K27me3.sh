#!/bin/sh
plotProfile \
-m D2_computematrix_DMRs_down_cginoncgi_allgeno_H3K27me3_v2.gz \
--samplesLabel D2_TKO_H3 D2_TKO_K4M D2_WT_H3 D2_WT_K4M \
--colors "#56a3fd" "#00d987" "#ff5d5e" "#ffab28" \
--perGroup \
--yAxisLabel "Norm.H3K27me3" \
--yMin 0 \
--yMax 130 \
--plotFileFormat pdf \
--regionsLabel cgi non_cgi \
--outFileName Plotprofile_D2_DMRs_down_cginoncgi_allgeno_H3K27me3_v2.pdf