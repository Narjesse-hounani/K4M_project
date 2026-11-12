#!/bin/sh
plotProfile \
-m D2_computematrix_heatmap_all_cginoncgi_allgeno_H3K4me3_v2.gz \
--samplesLabel D2_WT_H3 D2_WT_K4M D2_TKO_H3 D2_TKO_K4M \
--colors "#ff5d5e" "#ffab28" "#56a3fd" "#00d987" \
--perGroup \
--yAxisLabel "Norm.H3K4me3" \
--plotFileFormat pdf \
--regionsLabel cgi non_cgi \
--outFileName Plotprofile_D2_cginoncgi_allgeno_H3K4me3_v2.pdf