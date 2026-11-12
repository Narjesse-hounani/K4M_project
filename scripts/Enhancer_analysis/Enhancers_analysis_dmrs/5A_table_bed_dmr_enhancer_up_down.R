## 2026.08.31 Marlet's script
## Getting the enhancer numbers (and intersect with DMR) for Figure 5A
library(dplyr)
library(GenomicRanges)
library(ggplot2)
setwd("/Users/maxim/Documents/MethyLasso_Angelique_v2/")

## -- Read Enhancers
Enhancers <- read.table(file = "WT_Master_Enhancer_Peakome_v6_annotated_TargetGene.txt", header = TRUE, sep = "\t")
## Factorize enhancer classx
Enhancers$D0_class <- factor(x = Enhancers$D0_class, levels = c("actv", "prim", "inac", "othr"))
Enhancers$D2_class <- factor(x = Enhancers$D2_class, levels = c("actv", "prim", "inac", "othr"))
Enhancers$D4_class <- factor(x = Enhancers$D4_class, levels = c("actv", "prim", "inac", "othr"))
## Remove associated gene (remove duplicated enhancers)
Enhancers <- Enhancers |> select(chr:peak_id) |> distinct()
table(Enhancers$D2_class)
## Convert to GRanges
gr_Enhancers <- GRanges(
  seqnames = Enhancers$chr,
  ranges   = IRanges(start = Enhancers$start, end = Enhancers$end))

## -- Read DMRs
DMRs <- read.table(file = "H3_vs_K4M_dmrs.tsv", header = TRUE, sep = "\t")
DMRs_bed <- DMRs %>% select(1:3)
write.table(DMRs_bed, file = "All_DMRs.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## Clean DMRs to simplify operations
DMRs <- DMRs |> select(chr, start, end, diff)
DMRs <- DMRs |> mutate(direction = case_when(
  diff > 0 ~ "down", # meth2 smaller than meth1
  diff < 0 ~ "up"))  # meth2 bigger than meth1
DMRs$direction <- factor(DMRs$direction, levels = c("up", "down"))
table(DMRs$direction)
## Convert to GRanges
gr_DMRs <- GRanges(
  seqnames = DMRs$chr,
  ranges   = IRanges(start = DMRs$start, end = DMRs$end))


DMR_down_bed <- DMRs %>% filter(direction == "down") %>% select(1:3)
DMR_up_bed <- DMRs %>% filter(direction == "up") %>% select(1:3)
write.table(DMR_down_bed, file = "DMRs_down.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(DMR_up_bed, file = "DMRs_up.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## ----------------------
## Plot size distribution

## Enhancers
ggplot(Enhancers, aes(x = end-start+1)) +
  geom_histogram(binwidth = 100, color = "gray30", fill = "gray80") +
  labs(x = "width (bp)", title = "Enhancers") +
  theme_bw()

## DMRs
ggplot(DMRs, aes(x = end-start+1)) +
  geom_histogram(binwidth = 200, color = "gray30", fill = "gray80") +
  labs(x = "width (bp)", title = "DMRs") +
  theme_bw()


## --------------------------------------
## Helper function: intersect two GRanges, keeping exactly the single largest overlap per query feature.
## label_query/label_subject name the two inputs; label_overlap names the actual intersect region.
intersect_max_overlap <- function(gr_query, gr_subject, label_query, label_subject, label_overlap) {
  hits <- findOverlaps(gr_query, gr_subject)
  overlap <- pintersect(gr_query[queryHits(hits)], gr_subject[subjectHits(hits)])
  frac <- width(overlap) / width(gr_query)[queryHits(hits)]
  frac_col <- paste0("frac_", label_query)

  labeled_coords <- function(gr, idx, label) {
    as.data.frame(gr)[idx, c("seqnames", "start", "end", "width")] |>
      setNames(paste0(c("seqnames", "start", "end", "width"), ".", label))
  }

  result <- cbind(
    labeled_coords(overlap, seq_along(overlap), label_overlap),
    labeled_coords(gr_query, queryHits(hits), label_query),
    labeled_coords(gr_subject, subjectHits(hits), label_subject))
  result[[frac_col]] <- frac

  group_cols <- paste0(c("seqnames", "start", "end"), ".", label_query)
  result |>
    group_by(across(all_of(group_cols))) |>
    slice_max(order_by = .data[[frac_col]], n = 1, with_ties = FALSE) |>
    ungroup()
}

## -------------------------------
## -- Intersect DMRs and Enhancers
result_DE <- intersect_max_overlap(gr_query = gr_DMRs, gr_subject = gr_Enhancers,
                                   label_query = "DMR", label_subject = "Enh", label_overlap = "DE")

## Retrieve information back (enhancer coordinates and DMR direction)
result_DE <- left_join(x = result_DE, y = Enhancers, by = c("seqnames.Enh" = "chr",
                                                            "start.Enh" = "start",
                                                            "end.Enh" = "end"))
result_DE <- left_join(x = result_DE, y = DMRs, by = c("seqnames.DMR" = "chr",
                                                       "start.DMR" = "start",
                                                       "end.DMR" = "end")) |> dplyr::rename(direction.DMR = direction)

## Send intersect coordinates to the back
result_DE <- result_DE |> relocate(seqnames.DE:width.DE, .after = last_col())

## Intersect of DMRs with enhancers, classified by enhancer type and DMR
table(result_DE[, c("D2_class", "direction.DMR")], exclude = "othr")

## Number of unique enhancers represented
nrow(result_DE |> filter(D2_class == "actv" & direction.DMR == "up") |> select(peak_id) |> distinct())


## -------------------------------
## -- Intersect Enhancers and DMRs
result_ED <- intersect_max_overlap(gr_query = gr_Enhancers, gr_subject = gr_DMRs,
                                   label_query = "Enh", label_subject = "DMR", label_overlap = "ED")

## Retrieve information back (enhancer coordinates and DMR direction)
result_ED <- left_join(x = result_ED, y = Enhancers, by = c("seqnames.Enh" = "chr",
                                                            "start.Enh" = "start",
                                                            "end.Enh" = "end"))
result_ED <- left_join(x = result_ED, y = DMRs, by = c("seqnames.DMR" = "chr",
                                                       "start.DMR" = "start",
                                                       "end.DMR" = "end")) |> dplyr::rename(direction.DMR = direction)

## Send intersect coordinates to the back
result_ED <- result_ED |> relocate(seqnames.ED:width.ED, .after = last_col())

## Intersect of DMRs with enhancers, classified by enhancer type and DMR
table(result_ED[, c("D2_class", "direction.DMR")], exclude = "othr")

## Number of unique DMRs represented
nrow(result_ED |> filter(D2_class == "actv" & direction.DMR == "up") |> select(seqnames.DMR:end.DMR) |> distinct())

result_ED_up <- result_ED %>% filter(direction.DMR == "up")

up_actv_enhancer <- result_ED_up %>% filter(D2_class == "actv") |> select(seqnames.Enh:end.Enh) |> distinct()

up_inac_enhancer <- result_ED_up %>% filter(D2_class == "inac") |> select(seqnames.Enh:end.Enh) |> distinct()

up_prim_enhancer <- result_ED_up %>% filter(D2_class == "prim") |> select(seqnames.Enh:end.Enh) |> distinct()

write.table(up_actv_enhancer, file = "H3_vs_K4M_up_actv_enhancer_dmrs_corrected_overlap.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(up_inac_enhancer, file = "H3_vs_K4M_up_inac_enhancer_dmrs_corrected_overlap.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(up_prim_enhancer, file = "H3_vs_K4M_up_prim_enhancer_dmrs_corrected_overlap.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

result_ED_down <- result_ED %>% filter(direction.DMR == "down")

down_actv_enhancer <- result_ED_down %>% filter(D2_class == "actv") |> select(seqnames.Enh:end.Enh) |> distinct()

down_inac_enhancer <- result_ED_down %>% filter(D2_class == "inac") |> select(seqnames.Enh:end.Enh) |> distinct()

down_prim_enhancer <- result_ED_down %>% filter(D2_class == "prim") |> select(seqnames.Enh:end.Enh) |> distinct()

write.table(down_actv_enhancer, file = "H3_vs_K4M_down_actv_enhancer_dmrs_corrected_overlap.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(down_inac_enhancer, file = "H3_vs_K4M_down_inac_enhancer_dmrs_corrected_overlap.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(down_prim_enhancer, file = "H3_vs_K4M_down_prim_enhancer_dmrs_corrected_overlap.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
