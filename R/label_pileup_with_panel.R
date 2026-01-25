#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
})

opt_list <- list(
  make_option(c("-i","--input"),  type="character", help="Input raw pileup (e.g., pileup/<sample>_raw_pileup.txt) [required]"),
  make_option(c("-p","--panel"),  type="character", help="Panel design TSV with columns: chr, feature_start_position, feature_stop_position, probe_strand, mip_name [required]"),
  make_option(c("-o","--output"), type="character", default=NULL, help="Output path; default = <dirname(input)>/<basename>_labeled.txt")
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$input) || is.null(opt$panel)) stop("Provide -i <pileup> and -p <panel>")

# Read files
pile <- fread(opt$input, sep="\t", header=TRUE, data.table=TRUE)
pan  <- fread(opt$panel, sep="\t", header=TRUE, data.table=TRUE)

# Basic column checks
req_panel <- c("chr","feature_start_position","feature_stop_position","probe_strand","mip_name")
miss <- setdiff(req_panel, names(pan))
if (length(miss)) stop("Panel missing columns: ", paste(miss, collapse=", "))
req_pile <- c("chr","pos","strand","smMIP")
missp <- setdiff(req_pile, names(pile))
if (length(missp)) stop("Pileup missing columns: ", paste(missp, collapse=", "))

# Harmonize chr styles
has_chr <- function(x) any(grepl("^chr", x))
pile_chr_has <- has_chr(pile$chr)
pan_chr_has  <- has_chr(pan$chr)
fix_chr <- function(x, want_chr) {
  x <- as.character(x)
  if (want_chr) {
    ifelse(grepl("^chr", x), x, paste0("chr", x))
  } else {
    sub("^chr", "", x)
  }
}
# Make both match pile's style
pan$chr <- fix_chr(pan$chr, pile_chr_has)

# Make ranges
pile.gr <- GRanges(seqnames = pile$chr,
                   ranges   = IRanges(pile$pos, width=1),
                   strand   = pile$strand)

pan.gr  <- GRanges(seqnames = pan$chr,
                   ranges   = IRanges(pan$feature_start_position, pan$feature_stop_position),
                   strand   = ifelse(pan$probe_strand %in% c("+","-"), pan$probe_strand, "*"),
                   mip_name = pan$mip_name)

# Overlap and choose first hit (positions should belong to at most one smMIP)
hits <- findOverlaps(pile.gr, pan.gr)
best <- rep(NA_integer_, length(pile.gr))
if (length(hits) > 0) {
  # For positions overlapping multiple features, pick the longest overlap
  library(IRanges)
  ov <- pintersect(pile.gr[queryHits(hits)], pan.gr[subjectHits(hits)])
  w  <- width(ov)
  ord <- order(queryHits(hits), -w)
  qh  <- queryHits(hits)[ord]; sh <- subjectHits(hits)[ord]
  first <- !duplicated(qh)
  best[qh[first]] <- sh[first]
}

new_names <- rep(NA_character_, nrow(pile))
ok <- !is.na(best)
if (any(ok)) new_names[ok] <- as.character(mcols(pan.gr)$mip_name[best[ok]])

# Report mapping stats
cat(sprintf("Labeled %d / %d pileup rows with panel mip_name (%.1f%%)\n",
            sum(ok), length(ok), 100*mean(ok)))

# Keep old smMIP column, write a labeled version
setnames(pile, "smMIP", "smMIP_instrument")
pile[, smMIP := new_names]
# Optional gene column parsed from mip_name (before first '/')
pile[, gene := ifelse(!is.na(smMIP), sub("/.*$", "", smMIP), NA_character_) ]

# Default output
if (is.null(opt$output)) {
  dn <- dirname(opt$input); bn <- basename(opt$input)
  out <- file.path(dn, sub("(?:\\.txt)?$", "_labeled.txt", bn, perl=TRUE))
} else {
  out <- opt$output
}
fwrite(pile, out, sep="\t")
cat("Wrote labeled pileup to:", out, "\n")
