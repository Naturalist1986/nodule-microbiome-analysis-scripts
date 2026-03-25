## T6SS_operon_diagram.R
## Draws a publication-quality operon diagram of the T6SS locus
## for the 6 Bradyrhizobium symbiont bins that carry vasJ (K11910).
##
## Outputs (saved to the script directory):
##   T6SS_operon_diagram_data.tsv            – raw plot data for all bins
##   T6SS_operon_diagram_combined.pdf/.png/.svg
##   T6SS_operon_diagram_{bin}.pdf/.png      (per-bin)
##
## Requires: RSQLite, readxl, dplyr, ggplot2, gggenes, svglite

suppressPackageStartupMessages({
  library(RSQLite)
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(gggenes)
  library(ggtext)
  library(svglite)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
script_dir <- dirname(normalizePath(if (interactive()) {
  "/mnt/c/Users/owner/My Drive (moshe.alon@mail.huji.ac.il)/Moshe/Efrat_Guy_Project/Boruta_New/bin_diversity_boruta/T6SS_operon_diagram.R"
} else {
  sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])
}))

project_dir  <- file.path(script_dir, "../..")
boruta_dir   <- file.path(script_dir, "..")
pipeline_dir <- file.path(project_dir, "New_Binning_Pipeline_Coassembly")
db_dir       <- file.path(pipeline_dir, "consolidated_bins_anvio")
ko_table     <- file.path(boruta_dir,   "ko_presence_absence_table.xlsx")
abund_tsv    <- file.path(pipeline_dir, "combined_abundance_long_renormalized.tsv")

# ── T6SS KO universe ──────────────────────────────────────────────────────────
t6ss_kos <- c(
  K11891 = "tssM",  K11892 = "dotU",
  K11893 = "tssK",  K11895 = "tssG",  K11896 = "tssF",  K11905 = "tssE",
  K11900 = "tssC",  K11901 = "impB",
  K11903 = "hcp",   K11904 = "vgrG",
  K11906 = "vasD",  K11907 = "clpV",  K11910 = "vasJ"
)

t6ss_class <- c(
  K11891 = "Membrane anchor",    K11892 = "Membrane anchor",
  K11893 = "Baseplate",          K11895 = "Baseplate",
  K11896 = "Baseplate",          K11905 = "Baseplate",
  K11900 = "Contractile sheath", K11901 = "Contractile sheath",
  K11903 = "Secreted / tip",     K11904 = "Secreted / tip",
  K11906 = "Accessory / reg.",   K11907 = "Accessory / reg.",
  K11910 = "Accessory / reg."
)

t6ss_palette <- c(
  "Membrane anchor"    = "#1f78b4",
  "Baseplate"          = "#33a02c",
  "Contractile sheath" = "#e31a1c",
  "Secreted / tip"     = "#ff7f00",
  "Accessory / reg."   = "#6a3d9a",
  "Other"              = "grey72"
)

# Per-gene colour palette (reference-figure style)
gene_palette <- c(
  vasJ  = "#6A3D9A",   # purple
  impB  = "#FF7F00",   # orange
  tssC  = "#2AB07F",   # teal-green
  tssE  = "#B2DF8A",   # light green
  tssF  = "#2078B4",   # blue
  tssG  = "#74C3E8",   # light blue
  vasD  = "#FB9A99",   # salmon
  tssK  = "#E31A1C",   # red
  dotU  = "#FDBF6F",   # peach
  tssM  = "#F768A1",   # pink
  vgrG  = "#9E9AC8",   # lavender
  hcp   = "#FECC5C",   # yellow
  clpV  = "#8B4513",   # brown
  Other = "#BABABA"    # grey
)

# ── Step 1: Identify vasJ bins and species labels ─────────────────────────────
message("Step 1 – Loading KO presence/absence table …")
ko_pa  <- read_excel(ko_table)
ko_col <- names(ko_pa)[1]
vasj_row <- ko_pa[ko_pa[[ko_col]] == "K11910", ]
if (nrow(vasj_row) == 0) stop("K11910 (vasJ) not found in ko_presence_absence_table.xlsx")

bin_counts <- as.numeric(vasj_row[-1])
bin_names  <- names(vasj_row[-1])
vasj_bins  <- bin_names[!is.na(bin_counts) & bin_counts >= 1]

target_bins <- c("hok_binette_bin1", "RH_binette_bin1", "RH_binette_bin2",
                 "carR_binette_bin2", "mtz_binette_bin10", "carK_binette_bin2")
vasj_bins   <- intersect(vasj_bins, target_bins)
message(sprintf("  vasJ bins found: %s", paste(vasj_bins, collapse = ", ")))

message("Step 1 – Loading taxonomy …")
abund <- read.table(abund_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
abund$species_label <- sub(".*s__", "", abund$classification)

species_map <- abund %>%
  filter(Bin_Prefixed %in% vasj_bins) %>%
  select(bin = Bin_Prefixed, species_label) %>%
  distinct(bin, .keep_all = TRUE)

# ── Step 2: Query anvio SQLite databases ─────────────────────────────────────
message("Step 2 – Querying anvio databases …")

query_bin <- function(bin_name) {
  db_path <- file.path(db_dir, paste0(bin_name, ".db"))
  if (!file.exists(db_path)) { warning("DB not found: ", db_path); return(NULL) }
  con <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con))

  gf <- dbGetQuery(con,
    "SELECT gene_callers_id, accession AS KO, [function] AS function_desc
     FROM gene_functions WHERE source = 'KOfam'")
  gc <- dbGetQuery(con,
    "SELECT gene_callers_id, contig, start, stop, direction FROM genes_in_contigs")

  gf_t6 <- gf[gf$KO %in% names(t6ss_kos), ]
  if (nrow(gf_t6) == 0) return(NULL)

  merged      <- merge(gf_t6, gc, by = "gene_callers_id")
  merged$bin  <- bin_name
  merged
}

all_t6ss <- do.call(rbind, Filter(Negate(is.null), lapply(vasj_bins, query_bin)))
if (nrow(all_t6ss) == 0) stop("No T6SS genes found in any bin.")

# ── Steps 3+4: Multi-contig operon — all T6SS contigs ordered by operon position ──
message("Step 3+4 – Building multi-segment operon coordinates …")

# Canonical gene order along the T6SS operon
ko_operon_order <- c(
  K11910 = 1, K11901 = 2, K11900 = 3, K11905 = 4,
  K11896 = 5, K11895 = 6, K11906 = 7, K11893 = 8,
  K11892 = 9, K11891 = 10, K11904 = 11, K11903 = 12
)

BREAK_GAP <- 2500   # display-space gap (bp-equivalent) drawn for each // break

# Order contigs per bin by the median operon rank of their T6SS genes
contig_order <- all_t6ss %>%
  mutate(operon_rank = unname(ko_operon_order[KO])) %>%
  group_by(bin, contig) %>%
  summarise(
    median_rank  = median(operon_rank, na.rm = TRUE),
    dominant_dir = ifelse(mean(direction == "r") > 0.5, "r", "f"),
    contig_g_min = min(start),
    contig_g_max = max(stop),
    .groups = "drop"
  ) %>%
  arrange(bin, median_rank)

# Build concatenated display coordinates for each bin
build_segments <- function(bin_name) {
  segs  <- contig_order[contig_order$bin == bin_name, ]
  genes <- all_t6ss[all_t6ss$bin == bin_name, ]

  x_offset   <- 0
  gene_list  <- vector("list", nrow(segs))
  break_list <- vector("list", max(0L, nrow(segs) - 1L))

  for (s in seq_len(nrow(segs))) {
    seg <- segs[s, ]
    sg  <- genes[genes$contig == seg$contig, ]
    sg  <- sg[order(sg$start), ]

    g_min   <- min(sg$start)
    seg_len <- max(sg$stop) - g_min

    sg$start_norm <- sg$start - g_min
    sg$end_norm   <- sg$stop  - g_min
    sg$forward    <- sg$direction == "f"

    if (seg$dominant_dir == "r") {
      os <- sg$start_norm; oe <- sg$end_norm
      sg$start_norm <- seg_len - oe
      sg$end_norm   <- seg_len - os
      sg$forward    <- TRUE
    }

    sg$start_norm     <- sg$start_norm + x_offset
    sg$end_norm       <- sg$end_norm   + x_offset
    sg$segment_idx    <- s
    sg$segment_contig <- seg$contig
    sg$cluster_start  <- g_min              # genomic min of this segment
    sg$seg_g_max      <- max(sg$stop)       # genomic max of this segment (original stop)
    sg$seg_x_offset   <- x_offset           # display x where this segment starts
    sg$seg_was_flipped <- seg$dominant_dir == "r"

    gene_list[[s]] <- sg

    if (s < nrow(segs)) {
      x_break <- max(sg$end_norm) + BREAK_GAP / 2
      break_list[[s]] <- data.frame(
        bin           = bin_name,
        x_break       = x_break,
        seg_before    = s,
        contig_before = seg$contig,
        contig_after  = segs$contig[s + 1L],
        stringsAsFactors = FALSE
      )
      x_offset <- max(sg$end_norm) + BREAK_GAP
    }
  }

  list(
    genes  = do.call(rbind, gene_list),
    breaks = if (any(lengths(break_list) > 0))
               do.call(rbind, Filter(Negate(is.null), break_list))
             else NULL
  )
}

seg_results <- lapply(unique(all_t6ss$bin), build_segments)
names(seg_results) <- unique(all_t6ss$bin)

# Assemble plot_df and breaks
plot_df <- do.call(rbind, lapply(seg_results, `[[`, "genes"))
breaks_df_raw <- do.call(rbind,
  Filter(Negate(is.null), lapply(seg_results, `[[`, "breaks")))

# Add functional labels and species
plot_df <- plot_df %>%
  left_join(species_map, by = "bin") %>%
  mutate(
    species_label = ifelse(is.na(species_label), bin, species_label),
    gene_label    = unname(t6ss_kos[KO]),
    gene_class    = unname(t6ss_class[KO]),
    y_label       = paste0(bin, "\n", species_label)
  )

# Fix factor level order (hok first → top of plot)
bin_order  <- intersect(target_bins, unique(plot_df$bin))
ylvl_order <- rev(sapply(bin_order,
                  function(b) unique(plot_df$y_label[plot_df$bin == b])))
plot_df$y_label <- factor(plot_df$y_label, levels = ylvl_order)

# ── Step 4b: Extend last segment of each track with downstream flanking genes ─
message("Step 4b – Fetching flanking genes …")

global_max <- max(plot_df$end_norm)

# Parameters taken from the LAST segment of each bin
bin_params <- plot_df %>%
  group_by(bin) %>%
  filter(segment_idx == max(segment_idx)) %>%
  summarise(
    contig        = unique(segment_contig),
    cluster_start = unique(cluster_start),   # genomic min of last segment
    seg_g_max     = unique(seg_g_max),        # genomic max of last segment
    locus_len     = max(end_norm),
    was_flipped   = unique(seg_was_flipped),
    seg_x_offset  = unique(seg_x_offset),
    t6ss_ids      = list(unique(all_t6ss$gene_callers_id[all_t6ss$bin == bin[1]])),
    .groups = "drop"
  )

# Query contig lengths for contig-end marker
bin_params$contig_end_norm <- mapply(function(b, contig, cluster_start, seg_g_max,
                                              locus_len, flipped, seg_x_offset) {
  db_path <- file.path(db_dir, paste0(b, ".db"))
  con     <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con))
  clen <- dbGetQuery(con,
    sprintf("SELECT length FROM contigs_basic_info WHERE contig = '%s'", contig))$length
  if (!length(clen) || is.na(clen)) return(global_max)
  if (!flipped) {
    # display end of contig = seg_x_offset + (contig_length - cluster_start)
    seg_x_offset + (clen - cluster_start)
  } else {
    # after flip, contig starts at genomic 0; maps to seg_x_offset + cluster_start
    seg_x_offset + cluster_start
  }
}, bin_params$bin, bin_params$contig, bin_params$cluster_start,
   bin_params$seg_g_max, bin_params$locus_len,
   bin_params$was_flipped, bin_params$seg_x_offset)

bin_params$contig_end_norm <- pmin(bin_params$contig_end_norm, global_max)

# Helper: extract short gene name from a COG24_FUNCTION description string.
# COG format: "Full description (ShortName) (PDB:...) (PUBMED:...)"
# Take the last parenthetical that looks like a gene name (2-10 alphanum chars,
# not PDB/PUBMED/PMID).
extract_gene_name <- function(desc) {
  vapply(desc, function(d) {
    if (is.na(d) || d == "") return(NA_character_)
    d     <- sub("!!!.*", "", d)                          # keep first if multi
    hits  <- regmatches(d, gregexpr("\\([A-Za-z0-9_-]{2,10}\\)", d))[[1]]
    hits  <- hits[!grepl("PDB|PUBMED|PMID|EC|GO", hits)]
    if (length(hits) == 0) return(NA_character_)
    gsub("[()]", "", hits[length(hits)])                  # last match = gene name
  }, character(1))
}

flank_list <- lapply(seq_len(nrow(bin_params)), function(i) {
  bp  <- bin_params[i, ]
  gap <- global_max - bp$locus_len
  if (gap < 100) return(NULL)

  db_path <- file.path(db_dir, paste0(bp$bin, ".db"))
  con     <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con))

  # Genomic window downstream of the T6SS cluster (transcription direction).
  if (!bp$was_flipped) {
    win_start <- as.integer(bp$seg_g_max)
    win_end   <- as.integer(bp$seg_g_max + gap)
  } else {
    win_start <- as.integer(max(0L, bp$cluster_start - gap))
    win_end   <- as.integer(bp$cluster_start)
  }

  gc <- dbGetQuery(con, sprintf(
    "SELECT gene_callers_id, contig, start, stop, direction
     FROM genes_in_contigs
     WHERE contig = '%s' AND start >= %d AND stop <= %d AND call_type = 1",
    bp$contig, win_start, win_end
  ))
  if (nrow(gc) == 0) return(NULL)
  gc <- gc[!gc$gene_callers_id %in% bp$t6ss_ids[[1]], ]
  if (nrow(gc) == 0) return(NULL)

  # ── Fetch annotations for these flanking genes ─────────────────────────────
  id_list <- paste(gc$gene_callers_id, collapse = ",")

  # COG24_FUNCTION → preferred source for short gene names
  cog <- dbGetQuery(con, sprintf(
    "SELECT gene_callers_id, function AS cog_func
     FROM gene_functions
     WHERE source = 'COG24_FUNCTION' AND gene_callers_id IN (%s)", id_list
  ))
  # One row per gene (take first if multiple COG hits)
  cog <- cog[!duplicated(cog$gene_callers_id), ]

  # KOfam → fallback: use accession (KO number)
  ko_ann <- dbGetQuery(con, sprintf(
    "SELECT gene_callers_id, accession AS ko_acc
     FROM gene_functions
     WHERE source = 'KOfam' AND gene_callers_id IN (%s)", id_list
  ))
  ko_ann <- ko_ann[!duplicated(ko_ann$gene_callers_id), ]

  # Merge annotations
  gc <- merge(gc, cog,   by = "gene_callers_id", all.x = TRUE)
  gc <- merge(gc, ko_ann, by = "gene_callers_id", all.x = TRUE)

  # Build gene_label: COG short name > KO accession > blank
  gc$cog_name <- extract_gene_name(gc$cog_func)
  gc$gene_label <- ifelse(!is.na(gc$cog_name), gc$cog_name,
                   ifelse(!is.na(gc$ko_acc),   gc$ko_acc, ""))

  # ── Normalise coordinates ──────────────────────────────────────────────────
  gc$orig_start_norm <- gc$start - bp$cluster_start
  gc$orig_end_norm   <- gc$stop  - bp$cluster_start
  # Add segment x-offset so flanking genes continue from the last T6SS gene
  gc$orig_start_norm <- gc$orig_start_norm + bp$seg_x_offset
  gc$orig_end_norm   <- gc$orig_end_norm   + bp$seg_x_offset

  if (bp$was_flipped) {
    # Flip relative to the segment's local origin, then add offset
    loc_s <- gc$orig_start_norm - bp$seg_x_offset
    loc_e <- gc$orig_end_norm   - bp$seg_x_offset
    seg_display_len <- bp$locus_len - bp$seg_x_offset
    gc$start_norm <- bp$seg_x_offset + (seg_display_len - loc_e)
    gc$end_norm   <- bp$seg_x_offset + (seg_display_len - loc_s)
    gc$forward    <- gc$direction == "r"
  } else {
    gc$start_norm <- gc$orig_start_norm
    gc$end_norm   <- gc$orig_end_norm
    gc$forward    <- gc$direction == "f"
  }

  gc$start_norm <- pmax(gc$start_norm, 0)
  gc$end_norm   <- pmin(gc$end_norm,   global_max)
  gc <- gc[gc$end_norm > gc$start_norm, ]
  if (nrow(gc) == 0) return(NULL)

  sp <- species_map$species_label[species_map$bin == bp$bin]
  gc$bin           <- bp$bin
  gc$species_label <- if (length(sp)) sp[1] else bp$bin
  gc$contig        <- bp$contig
  gc$cluster_start <- bp$cluster_start
  gc$KO            <- NA_character_
  gc$function_desc <- gc$cog_func
  gc$gene_class    <- "Other"
  gc$segment_idx    <- max(plot_df$segment_idx[plot_df$bin == bp$bin])
  gc$segment_contig <- bp$contig
  gc$cluster_start  <- bp$cluster_start
  gc$seg_g_max      <- bp$seg_g_max
  gc$seg_x_offset   <- bp$seg_x_offset
  gc$seg_was_flipped <- bp$was_flipped
  gc$y_label       <- as.character(unique(plot_df$y_label[plot_df$bin == bp$bin]))
  gc
})

flank_df <- do.call(rbind, Filter(Negate(is.null), flank_list))

# Combine T6SS genes + flanking genes into one data frame for plotting
# (plot_df is kept intact for the summary table and TSV export)
if (!is.null(flank_df) && nrow(flank_df) > 0) {
  # Keep only columns present in plot_df; fill missing flanking columns with NA
  for (col in setdiff(names(plot_df), names(flank_df))) flank_df[[col]] <- NA
  plot_all <- bind_rows(plot_df, flank_df[, names(plot_df)])
  plot_all$y_label <- factor(plot_all$y_label, levels = levels(plot_df$y_label))
} else {
  plot_all <- plot_df
}

# Per-gene fill: T6SS gene label, or "Other" for all flanking genes
plot_all <- plot_all %>%
  mutate(fill_var = ifelse(!is.na(KO), gene_label, "Other"))

message(sprintf("  Added %d flanking genes across %d bins",
                if (!is.null(flank_df)) nrow(flank_df) else 0L,
                length(unique(bin_params$bin[bin_params$locus_len < global_max - 100]))))

# ── Console summary ───────────────────────────────────────────────────────────
message("\n── T6SS Locus Summary ──────────────────────────────────────────────────")
summary_tbl <- plot_df %>%
  group_by(bin, species_label, contig) %>%
  summarise(
    n_T6SS_genes    = n(),
    cluster_span_kb = round((max(stop) - min(start)) / 1000, 1),
    KOs_found       = paste(sort(unique(KO)), collapse = ", "),
    .groups = "drop"
  )
print(as.data.frame(summary_tbl))
message("────────────────────────────────────────────────────────────────────────\n")

# ── Export raw plot data ───────────────────────────────────────────────────────
message("Exporting plot data …")
export_df <- plot_df %>%
  select(
    bin, species_label, contig,
    gene_callers_id, KO, gene_label, gene_class, function_desc,
    start_genomic  = start,
    stop_genomic   = stop,
    direction,
    start_norm, end_norm,
    forward,
    cluster_start
  ) %>%
  arrange(bin, start_norm)

write.table(export_df,
            file      = file.path(script_dir, "T6SS_operon_diagram_data.tsv"),
            sep       = "\t",
            row.names = FALSE,
            quote     = FALSE)
message("  Saved T6SS_operon_diagram_data.tsv")

# ── Step 5: Build y-axis labels (markdown: bold bin, italic species) ─────────
message("Step 5 – Building ggplot (gggenes) …")

# Add markdown y labels to both data frames
add_y_md <- function(df) {
  df %>% mutate(y_md = paste0("**", bin, "**<br>*", species_label, "*"))
}
plot_df  <- add_y_md(plot_df)
plot_all <- add_y_md(plot_all)

# Ordered factor — hok at top (last level → highest y in discrete scale)
md_levels <- sapply(bin_order, function(b) unique(plot_df$y_md[plot_df$bin == b]))
plot_df$y_md  <- factor(plot_df$y_md,  levels = rev(md_levels))
plot_all$y_md <- factor(plot_all$y_md, levels = rev(md_levels))

make_combined_plot <- function(df) {

  n_rows <- length(levels(df$y_md))

  # Build per-bin contig-end info in y_md space
  contig_end_df <- plot_all %>%
    select(bin, y_md) %>%
    distinct() %>%
    left_join(
      bin_params %>% select(bin, contig_end_norm),
      by = "bin"
    ) %>%
    mutate(
      contig_end_norm = pmin(contig_end_norm, global_max),
      truncated       = contig_end_norm < global_max * 0.99
    )
  contig_end_df$y_md <- factor(contig_end_df$y_md, levels = levels(df$y_md))

  # Backbone stops at the contig end (not global_max) for short contigs
  backbone_df <- contig_end_df %>%
    transmute(y_md, xmin = 0, xmax = contig_end_norm)

  # Vertical tick marks where contigs end
  tick_df <- contig_end_df[contig_end_df$truncated, ]

  # Add y_md to breaks
  breaks_plot_df <- if (!is.null(breaks_df_raw) && nrow(breaks_df_raw) > 0) {
    bdf <- merge(breaks_df_raw,
                 unique(plot_all[, c("bin", "y_md")]),
                 by = "bin")
    bdf$y_md <- factor(bdf$y_md, levels = levels(df$y_md))
    bdf
  } else {
    data.frame(x_break = numeric(0), y_md = factor(character(0)),
               contig_before = character(0), contig_after = character(0))
  }

  p <- ggplot(df,
       aes(xmin    = start_norm,
           xmax    = end_norm,
           y       = y_md,
           fill    = fill_var,
           forward = forward,
           label   = gene_label)) +
    # Backbone — stops at contig end
    geom_segment(data = backbone_df,
                 aes(x = xmin, xend = xmax, y = y_md, yend = y_md),
                 colour = "grey70", linewidth = 0.5,
                 inherit.aes = FALSE) +
    # Vertical tick at contig end for truncated bins
    geom_segment(data = tick_df,
                 aes(x = contig_end_norm, xend = contig_end_norm,
                     y = as.numeric(y_md) - 0.42,
                     yend = as.numeric(y_md) + 0.42),
                 colour = "grey30", linewidth = 1.0,
                 inherit.aes = FALSE) +
    # "contig end" text label above the tick
    geom_text(data = tick_df,
              aes(x = contig_end_norm, y = as.numeric(y_md) + 0.55,
                  label = "contig end"),
              colour = "grey30", size = 3.2, fontface = "italic",
              hjust = 0.5, inherit.aes = FALSE) +
    geom_gene_arrow(
      arrowhead_height  = unit(7,   "mm"),
      arrowhead_width   = unit(4,   "mm"),
      arrow_body_height = unit(7,   "mm"),
      colour = "white", linewidth = 0.3
    ) +
    # Gene labels — dark text inside arrows
    geom_gene_label(
      fontface  = "plain",
      size      = 5.5,
      min.size  = 3.0,
      grow      = FALSE,
      colour    = "grey10",
      padding.x = grid::unit(0.3, "mm"),
      padding.y = grid::unit(0.1, "mm")
    ) +
    # // break markers between contig segments
    geom_segment(data = breaks_plot_df,
                 aes(x = x_break, xend = x_break,
                     y = as.numeric(y_md) - 0.48,
                     yend = as.numeric(y_md) + 0.48),
                 colour = "grey20", linewidth = 0.8, linetype = "dashed",
                 inherit.aes = FALSE) +
    geom_text(data = breaks_plot_df,
              aes(x = x_break, y = as.numeric(y_md),
                  label = "//"),
              colour = "grey20", size = 5, fontface = "bold",
              hjust = 0.5, inherit.aes = FALSE) +
    # 5 kb scale bar below lowest track
    annotate("segment",
             x = 0, xend = 5000,
             y = 0.45, yend = 0.45,
             linewidth = 1.2, colour = "grey30") +
    annotate("text",
             x = 2500, y = 0.22,
             label = "5 kb", size = 4.0, colour = "grey30") +
    scale_fill_manual(values = gene_palette, name = "gene",
                      na.value = "#BABABA") +
    scale_x_continuous(
      labels = function(x) paste0(round(x / 1000, 0), " kb"),
      expand = expansion(mult = c(0.01, 0.04)),
      limits = c(-200, global_max * 1.02)
    ) +
    scale_y_discrete(expand = expansion(add = c(0.70, 1.00))) +
    labs(
      title = "T6SS Locus in *Bradyrhizobium* Symbiont Bins",
      x     = "Genomic position (kb, zero-normalised to first T6SS gene)",
      y     = NULL
    ) +
    theme_genes() +
    theme(
      axis.text.y      = element_markdown(size = 13, hjust = 1, lineheight = 1.35),
      axis.text.x      = element_text(size = 12),
      axis.title.x     = element_text(size = 12, margin = margin(t = 7)),
      legend.position  = "right",
      legend.direction = "vertical",
      legend.title     = element_text(face = "bold", size = 11),
      legend.text      = element_text(size = 10),
      legend.key.size  = unit(4.5, "mm"),
      plot.title       = element_markdown(face = "bold", size = 16,
                                          hjust = 0.5, margin = margin(b = 14)),
      plot.margin      = margin(22, 22, 10, 12),
      panel.grid       = element_blank(),
      panel.border     = element_blank()
    )
  p
}

# ── Per-bin figure: single row, faceted for consistent strip label style ─────
make_perbin_plot <- function(df_bin) {
  b  <- unique(df_bin$bin)
  sp <- unique(df_bin$species_label)

  b_name  <- unique(df_bin$bin)
  bdf_bin <- if (!is.null(breaks_df_raw) && nrow(breaks_df_raw) > 0 &&
                  b_name %in% breaks_df_raw$bin)
               breaks_df_raw[breaks_df_raw$bin == b_name, ]
             else data.frame(x_break = numeric(0))

  ggplot(df_bin,
         aes(xmin    = start_norm,
             xmax    = end_norm,
             y       = bin,
             fill    = gene_class,
             forward = forward,
             label   = gene_label)) +
    geom_gene_arrow(
      arrowhead_height  = unit(6,   "mm"),
      arrowhead_width   = unit(4,   "mm"),
      arrow_body_height = unit(6,   "mm"),
      colour = "white", linewidth = 0.25
    ) +
    geom_gene_label(
      fontface  = "bold",
      size      = 3,
      min.size  = 0,
      grow      = FALSE,
      padding.x = grid::unit(1.5, "mm"),
      padding.y = grid::unit(0.5, "mm")
    ) +
    geom_vline(data = bdf_bin,
               aes(xintercept = x_break),
               colour = "grey20", linewidth = 0.8, linetype = "dashed",
               inherit.aes = FALSE) +
    geom_text(data = bdf_bin,
              aes(x = x_break, label = "//"),
              y = Inf, colour = "grey20", size = 5, fontface = "bold",
              vjust = 1.5, hjust = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(values = t6ss_palette, name = "T6SS function") +
    scale_x_continuous(
      labels = function(x) paste0(round(x / 1000, 1), " kb"),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    labs(
      title    = paste0("**", b, "** — *", sp, "*"),
      subtitle = paste0(
        unique(df_bin$contig), "  |  ",
        nrow(df_bin), " T6SS genes  |  ",
        round((max(df_bin$stop) - min(df_bin$start)) / 1000, 1), " kb locus"
      ),
      x = "Genomic position (kb)",
      y = NULL
    ) +
    theme_genes() +
    theme(
      axis.text.y    = element_blank(),
      axis.ticks.y   = element_blank(),
      axis.text.x    = element_text(size = 9),
      legend.position = "right",
      legend.title   = element_text(face = "bold", size = 8),
      legend.text    = element_text(size = 8),
      plot.title     = element_markdown(size = 12, margin = margin(b = 3)),
      plot.subtitle  = element_text(size = 8, colour = "grey40", margin = margin(b = 8)),
      plot.margin    = margin(10, 15, 10, 10)
    )
}

# ── Step 9: Save outputs ───────────────────────────────────────────────────────
message("Step 9 – Saving outputs …")
out_dir <- script_dir

save_plot <- function(p, base_name, width, height) {
  ggsave(file.path(out_dir, paste0(base_name, ".png")),
         plot = p, width = width, height = height, units = "in", dpi = 300)
  ggsave(file.path(out_dir, paste0(base_name, ".svg")),
         plot = p, width = width, height = height, units = "in")
  tryCatch(
    ggsave(file.path(out_dir, paste0(base_name, ".pdf")),
           plot = p, width = width, height = height, units = "in"),
    error = function(e) message("  (PDF skipped – file may be open in viewer: ", conditionMessage(e), ")")
  )
  message(sprintf("  Saved %s (.png/.svg)", base_name))
}

# Combined — single panel, 6 rows (uses plot_all: T6SS + flanking genes)
n_bins     <- length(bin_order)
p_combined <- make_combined_plot(plot_all)
save_plot(p_combined, "T6SS_operon_diagram_combined",
          width = 15, height = 2.0 + n_bins * 0.95)

# Per-bin (also uses flanking genes)
for (b in bin_order) {
  sub <- plot_all[plot_all$bin == b, ]
  p_b <- make_perbin_plot(sub)
  save_plot(p_b, paste0("T6SS_operon_diagram_", b), width = 10, height = 2.8)
}

message("\nDone. All outputs written to: ", out_dir)
