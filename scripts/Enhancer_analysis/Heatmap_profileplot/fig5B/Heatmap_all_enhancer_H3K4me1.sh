#!/bin/sh
plotHeatmap -m computematrix_heatmap_all_enhancers_allgeno_H3K4me1.gz \
--colorMap YlOrBr \
--sortRegions descend \
--sortUsingSamples 1 \
--xAxisLabel "Genomic bin" \
--yAxisLabel "Norm. K4me1" \
--regionsLabel actv_enhancers, prim_enhancers, inac_enhancers \
--zMin 0 \
--zMax 100 \
--yMin 0  \
--yMax 170  \
--plotFileFormat pdf \
--outFileName computematrix_heatmap_all_enhancers_allgeno_H3K4me1.pdf
