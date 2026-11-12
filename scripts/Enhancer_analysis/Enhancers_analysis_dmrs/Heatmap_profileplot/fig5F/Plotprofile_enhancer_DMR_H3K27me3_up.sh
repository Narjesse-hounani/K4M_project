#!/bin/sh
plotProfile \
-m computematrix_DMR_enhancer_overlapCORRECTED_binsizeCORRECTED_up_H3K27me3.gz \
--samplesLabel D2_TKO_H3 D2_TKO_K4M D2_WT_H3 D2_WT_K4M \
--colors "#56a3fd" "#00d987" "#ff5d5e" "#ffab28" \
--perGroup \
--yAxisLabel "Norm.H3K27me3" \
--yMin 0 \
--yMax 50 \
--plotFileFormat pdf \
--regionsLabel active prim inac \
--outFileName Plotprofile_DMR_enhancer_overlapCORRECTED_binsizeCORRECTED_up_H3K27me3.pdf