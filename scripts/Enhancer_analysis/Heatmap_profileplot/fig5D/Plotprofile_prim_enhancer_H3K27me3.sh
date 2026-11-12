#!/bin/sh
plotProfile \
-m computematrix_prim_enhancer_allgeno_H3K27me3_all_bsCORRECTED.gz \
--samplesLabel D2_TKO_H3 D2_TKO_K4M D2_WT_H3 D2_WT_K4M \
--colors "#56a3fd" "#00d987" "#ff5d5e" "#ffab28" \
--perGroup \
--yAxisLabel "Norm.H3K27me3" \
--yMin 0 \
--yMax 40 \
--plotFileFormat pdf \
--outFileName plotProfile_prim_enhancer_allgeno_H3K27me3_all_bsCORRECTED.pdf