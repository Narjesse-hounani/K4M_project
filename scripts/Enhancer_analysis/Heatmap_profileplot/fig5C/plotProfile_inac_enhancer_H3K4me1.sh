#!/bin/sh
plotProfile \
-m computematrix_inac_enhancer_allgeno_H3K4me1_all.gz \
--samplesLabel D2_TKO_H3 D2_TKO_K4M D2_WT_H3 D2_WT_K4M \
--colors "#56a3fd" "#00d987" "#ff5d5e" "#ffab28" \
--perGroup \
--yAxisLabel "Norm.H3K4me1" \
--yMin 0 \
--yMax 140 \
--plotFileFormat pdf \
--outFileName computematrix_inac_enhancer_allgeno_H3K4me1_all.pdf