#!/bin/sh
plotProfile \
-m computematrix_DMRs_down_cginoncgi_allgeno_5mc_v2.gz \
--samplesLabel D2_WT_H3 D2_WT_K4M \
--colors "#ff5d5e" "#ffab28" \
--perGroup \
--yAxisLabel "5mc" \
--yMin 0 \
--yMax 100 \
--plotFileFormat pdf \
--regionsLabel cgi noncgi \
--outFileName Plotprofile_DMRs_down_cginoncgi_allgeno_5mc_v2.pdf