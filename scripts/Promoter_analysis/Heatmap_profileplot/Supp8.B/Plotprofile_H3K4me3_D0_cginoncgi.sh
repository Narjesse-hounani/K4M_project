#!/bin/sh
plotProfile \
-m D0_computematrix_heatmap_all_cginoncgi_allgeno_H3K4me3_v2.gz \
--samplesLabel D0_WT_H3 D0_WT_K4M D0_TKO_H3 D0_TKO_K4M \
--colors "#ff5d5e" "#ffab28" "#56a3fd" "#00d987" \
--perGroup \
--yAxisLabel "Norm.H3K4me3" \
--plotFileFormat pdf \
--regionsLabel cgi non_cgi \
--outFileName Plotprofile_D0_cginoncgi_allgeno_H3K4me3_v2.pdf