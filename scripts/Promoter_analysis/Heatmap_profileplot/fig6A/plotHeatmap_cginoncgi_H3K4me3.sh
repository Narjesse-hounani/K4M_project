#!/bin/sh
plotHeatmap -m computematrix_heatmap_all_cginoncgi_allgeno_H3K4me3_v2.gz \
--colorMap YlOrBr \
--sortRegions descend \
--xAxisLabel "Promoter" \
--yAxisLabel "Norm. K4me3" \
--regionsLabel cgi non_cgi \
--plotFileFormat pdf \
--outFileName heatmap_all_cginoncgi_allgeno_H3K4me3_v2.pdf