#!/bin/sh
plotProfile \
-m D0_computematrix_heatmap_all_cginoncgi_allgeno_H3K27me3_v2.gz \
--samplesLabel D0_TKO_H3 D0_TKO_K4M D0_WT_H3 D0_WT_K4M \
--colors "#56a3fd" "#00d987" "#ff5d5e" "#ffab28" \
--perGroup \
--yAxisLabel "Norm.H3K27me3" \
--yMin 0 \
--plotFileFormat pdf \
--regionsLabel cgi non_cgi \
--outFileName Plotprofile_D0_cginoncgi_allgeno_H3K27me3_v2.pdf