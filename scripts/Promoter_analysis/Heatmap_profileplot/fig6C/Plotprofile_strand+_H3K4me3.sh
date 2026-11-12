#!/bin/sh
plotProfile \
-m computematrix_heatmap_TSS_TES_allgeno_D2_strand-_H3K4me3.gz \
--samplesLabel D2_TKO_H3 D2_TKO_K4M D2_WT_H3 D2_WT_K4M \
--colors "#56a3fd" "#00d987" "#ff5d5e" "#ffab28" \
--perGroup \
--yAxisLabel "Norm.H3K4me3" \
--plotFileFormat pdf \
--plotFileFormat pdf \
--yMin 0 \
--yMax 55 \
--outFileName plotProfile_TSS_TES_allgeno_D2_strand-_H3K4me3.pdf