#!/bin/sh
plotProfile \
-m computematrix_DMR_enhancer_overlapCORRECTED_binsizeCORRECTED_5mc_up.gz \
--samplesLabel D2_WT_H3 D2_WT_K4M \
--colors "#ff5d5e" "#ffab28" \
--perGroup \
--yAxisLabel "5mc" \
--yMin 0 \
--yMax 100 \
--plotFileFormat pdf \
--regionsLabel active prim inac \
--outFileName Plotprofile_DMR_enhancer_overlapCORRECTED_binsizeCORRECTED_5mc_up.pdf