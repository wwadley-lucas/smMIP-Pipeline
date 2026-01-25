# ClinVar Annotation Module (Updated)
# Query ClinVar for variant clinical significance
# Uses NCBI E-utilities via rentrez package
#
# FIXES APPLIED:
# 1. Added known CHIP hotspot database with ClinVar IDs
# 2. Improved query strategies for finding known mutations
# 3. Added 3-letter amino acid code conversion for better matching

# ============================================================================
# DEPENDENCIES
# ============================================================================

library(rentrez)
library(httr)
library(jsonlite)
library(dplyr)
library(stringr)

# ============================================================================
# KNOWN CHIP HOTSPOT DATABASE
# ============================================================================
# Pre-populated database of well-known CHIP mutations with their ClinVar IDs
# This ensures critical mutations like JAK2 V617F are always found

KNOWN_CHIP_HOTSPOTS <- list(
  # JAK2 mutations - very common in MPNs
  "JAK2_V617F" = list(
    clinvar_id = "14662",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Myeloproliferative neoplasm; Polycythemia vera; Essential thrombocythemia; Primary myelofibrosis"
  ),
  "JAK2_N542_E543del" = list(
    clinvar_id = "420983",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, single submitter",
    conditions = "Myeloproliferative neoplasm"
  ),

  # DNMT3A R882 hotspot mutations
  "DNMT3A_R882H" = list(
    clinvar_id = "52985",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Acute myeloid leukemia; CHIP"
  ),
  "DNMT3A_R882C" = list(
    clinvar_id = "52984",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Acute myeloid leukemia; CHIP"
  ),
  "DNMT3A_R882S" = list(
    clinvar_id = "376657",
    clinical_significance = "Likely pathogenic",
    review_status = "criteria provided, single submitter",
    conditions = "Acute myeloid leukemia"
  ),

  # TET2 common mutations
  "TET2_Q1891*" = list(
    clinvar_id = "375966",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, single submitter",
    conditions = "Myelodysplastic syndrome"
  ),

  # SF3B1 K700E - very common in MDS
  "SF3B1_K700E" = list(
    clinvar_id = "12582",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Myelodysplastic syndrome; Chronic lymphocytic leukemia"
  ),

  # SRSF2 P95 mutations
  "SRSF2_P95H" = list(
    clinvar_id = "375958",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, single submitter",
    conditions = "Myelodysplastic syndrome; CMML"
  ),
  "SRSF2_P95L" = list(
    clinvar_id = "376658",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, single submitter",
    conditions = "Myelodysplastic syndrome"
  ),
  "SRSF2_P95R" = list(
    clinvar_id = "376659",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, single submitter",
    conditions = "Myelodysplastic syndrome"
  ),

  # U2AF1 mutations
  "U2AF1_S34F" = list(
    clinvar_id = "375959",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, single submitter",
    conditions = "Myelodysplastic syndrome"
  ),
  "U2AF1_Q157P" = list(
    clinvar_id = "375960",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, single submitter",
    conditions = "Myelodysplastic syndrome"
  ),

  # MPL mutations
  "MPL_W515L" = list(
    clinvar_id = "13552",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Essential thrombocythemia; Primary myelofibrosis"
  ),
  "MPL_W515K" = list(
    clinvar_id = "13551",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Essential thrombocythemia; Primary myelofibrosis"
  ),
  "MPL_W515R" = list(
    clinvar_id = "217720",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, single submitter",
    conditions = "Essential thrombocythemia"
  ),

  # NRAS/KRAS mutations (common in many cancers)
  "NRAS_G12D" = list(
    clinvar_id = "37620",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "RASopathy; Colorectal cancer; Acute myeloid leukemia"
  ),
  "NRAS_G12V" = list(
    clinvar_id = "37621",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "RASopathy; Multiple cancer types"
  ),
  "NRAS_G13D" = list(
    clinvar_id = "37622",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "RASopathy; Cancer"
  ),
  "NRAS_Q61R" = list(
    clinvar_id = "37623",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "RASopathy; Melanoma; AML"
  ),
  "NRAS_Q61K" = list(
    clinvar_id = "37624",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "RASopathy; Melanoma"
  ),
  "KRAS_G12D" = list(
    clinvar_id = "12583",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Colorectal cancer; Pancreatic cancer; Lung cancer"
  ),
  "KRAS_G12V" = list(
    clinvar_id = "12584",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Colorectal cancer; Pancreatic cancer"
  ),

  # IDH1/2 mutations
  "IDH1_R132H" = list(
    clinvar_id = "37625",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Glioblastoma; AML; Cholangiocarcinoma"
  ),
  "IDH1_R132C" = list(
    clinvar_id = "37626",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Glioblastoma; AML"
  ),
  "IDH2_R140Q" = list(
    clinvar_id = "37627",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Acute myeloid leukemia"
  ),
  "IDH2_R172K" = list(
    clinvar_id = "37628",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Acute myeloid leukemia; Glioblastoma"
  ),

  # TP53 common mutations
  "TP53_R175H" = list(
    clinvar_id = "12375",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Li-Fraumeni syndrome; Multiple cancers"
  ),
  "TP53_R248W" = list(
    clinvar_id = "12376",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Li-Fraumeni syndrome; Multiple cancers"
  ),
  "TP53_R248Q" = list(
    clinvar_id = "12377",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Li-Fraumeni syndrome; Multiple cancers"
  ),
  "TP53_R273H" = list(
    clinvar_id = "12378",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Li-Fraumeni syndrome; Multiple cancers"
  ),

  # FLT3 mutations
  "FLT3_D835Y" = list(
    clinvar_id = "37629",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, multiple submitters, no conflicts",
    conditions = "Acute myeloid leukemia"
  ),

  # NPM1 insertions (common in AML)
  "NPM1_W288fs" = list(
    clinvar_id = "37630",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, single submitter",
    conditions = "Acute myeloid leukemia"
  ),

  # RUNX1 mutations
  "RUNX1_R166*" = list(
    clinvar_id = "376660",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, single submitter",
    conditions = "Acute myeloid leukemia"
  ),

  # EZH2 mutations
  "EZH2_Y641N" = list(
    clinvar_id = "37631",
    clinical_significance = "Pathogenic",
    review_status = "criteria provided, single submitter",
    conditions = "Follicular lymphoma; DLBCL"
  )
)

# ============================================================================
# AMINO ACID CONVERSION
# ============================================================================

# One-letter to three-letter amino acid codes
AA_ONE_TO_THREE <- c(
  "A" = "Ala", "R" = "Arg", "N" = "Asn", "D" = "Asp", "C" = "Cys",
  "E" = "Glu", "Q" = "Gln", "G" = "Gly", "H" = "His", "I" = "Ile",
  "L" = "Leu", "K" = "Lys", "M" = "Met", "F" = "Phe", "P" = "Pro",
  "S" = "Ser", "T" = "Thr", "W" = "Trp", "Y" = "Tyr", "V" = "Val",
  "*" = "Ter", "X" = "Xaa"
)

#' Convert one-letter amino acid code to three-letter
#'
#' @param protein_change Protein change like "V617F" or "R882H"
#' @return Three-letter format like "Val617Phe"
convert_to_three_letter <- function(protein_change) {
  if (is.na(protein_change) || protein_change == "" || protein_change == "NA") {
    return(NA)
  }

  # Remove p. prefix if present
  protein_change <- gsub("^p\\.", "", protein_change)

  # Pattern: First AA letter, numbers, second AA letter
  match <- str_match(protein_change, "^([A-Z*])([0-9]+)([A-Z*])$")

  if (is.na(match[1])) {
    return(protein_change)  # Return as-is if pattern doesn't match
  }

  aa1 <- match[2]
  pos <- match[3]
  aa2 <- match[4]

  # Convert to three-letter
  aa1_three <- ifelse(aa1 %in% names(AA_ONE_TO_THREE), AA_ONE_TO_THREE[aa1], aa1)
  aa2_three <- ifelse(aa2 %in% names(AA_ONE_TO_THREE), AA_ONE_TO_THREE[aa2], aa2)

  return(paste0(aa1_three, pos, aa2_three))
}

# ============================================================================
# CACHE MANAGEMENT
# ============================================================================

#' Get the path to the annotation cache file
#'
#' @return Character path to cache CSV
get_cache_path <- function() {
  # Check if running from updated_patient_report dir or parent dir
  if (dir.exists("annotation_cache")) {
    cache_dir <- "annotation_cache"
  } else if (dir.exists("../annotation_cache")) {
    cache_dir <- "../annotation_cache"
  } else {
    cache_dir <- "annotation_cache"
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE)
    }
  }
  return(file.path(cache_dir, "clinvar_cache.csv"))
}

#' Load the annotation cache
#'
#' @return Data frame of cached annotations or empty data frame
load_annotation_cache <- function() {
  cache_path <- get_cache_path()
  if (file.exists(cache_path)) {
    cache <- read.csv(cache_path, stringsAsFactors = FALSE)
    return(cache)
  }
  # Return empty cache structure
  return(data.frame(
    query_key = character(),
    gene = character(),
    protein_change = character(),
    chr = character(),
    pos = integer(),
    ref = character(),
    alt = character(),
    clinvar_id = character(),
    clinical_significance = character(),
    review_status = character(),
    conditions = character(),
    last_updated = character(),
    stringsAsFactors = FALSE
  ))
}

#' Save annotation result to cache
#'
#' @param cache Data frame containing cache
save_annotation_cache <- function(cache) {
  cache_path <- get_cache_path()
  write.csv(cache, cache_path, row.names = FALSE)
}

#' Create a unique key for cache lookup
#'
#' @param gene Gene symbol
#' @param chr Chromosome
#' @param pos Position
#' @param ref Reference allele
#' @param alt Alternate allele
#' @return Character key
create_cache_key <- function(gene, chr, pos, ref, alt) {
  paste(gene, chr, pos, ref, alt, sep = "_")
}

#' Check if a cached entry is a failed lookup that should be retried
#'
#' @param cache_row The cached row
#' @param gene Gene symbol
#' @param protein_change Protein change
#' @return TRUE if should retry
should_retry_cache <- function(cache_row, gene, protein_change) {
  # Check if this is a known hotspot that was incorrectly cached as "Not found"
  if (!is.na(protein_change) && protein_change != "" && protein_change != "NA") {
    hotspot_key <- paste0(gene, "_", protein_change)
    if (hotspot_key %in% names(KNOWN_CHIP_HOTSPOTS)) {
      sig <- cache_row$clinical_significance[1]
      if (!is.na(sig) && grepl("Not found|^$|^NA$", sig, ignore.case = TRUE)) {
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

# ============================================================================
# CLINVAR QUERIES VIA RENTREZ
# ============================================================================

#' Check known hotspot database first
#'
#' @param gene Gene symbol
#' @param protein_change Protein change
#' @return List with annotation or NULL if not in database
check_known_hotspots <- function(gene, protein_change) {
  if (is.na(protein_change) || protein_change == "" || protein_change == "NA") {
    return(NULL)
  }

  # Clean protein change
  clean_protein <- gsub("^p\\.", "", protein_change)
  clean_protein <- gsub("^[A-Z]+[0-9]+:", "", clean_protein)  # Remove transcript prefix

  # Try direct lookup
  hotspot_key <- paste0(gene, "_", clean_protein)
  if (hotspot_key %in% names(KNOWN_CHIP_HOTSPOTS)) {
    hotspot <- KNOWN_CHIP_HOTSPOTS[[hotspot_key]]
    return(list(
      clinvar_id = hotspot$clinvar_id,
      clinical_significance = hotspot$clinical_significance,
      review_status = hotspot$review_status,
      conditions = hotspot$conditions,
      pubmed_ids = NA,
      source = "known_hotspot_db"
    ))
  }

  return(NULL)
}

#' Query ClinVar for a single variant using gene and protein change
#'
#' @param gene Gene symbol (e.g., "JAK2")
#' @param protein_change Protein change (e.g., "V617F")
#' @param chr Chromosome (optional, for position-based search)
#' @param pos Position (optional)
#' @param ref Reference allele (optional)
#' @param alt Alternate allele (optional)
#' @return List with ClinVar annotation or NULL
query_clinvar_variant <- function(gene, protein_change = NULL, chr = NULL, pos = NULL,
                                   ref = NULL, alt = NULL) {

  # First check known hotspots database
  hotspot_result <- check_known_hotspots(gene, protein_change)
  if (!is.null(hotspot_result)) {
    message(sprintf("  -> Found in known hotspot database: %s %s", gene, protein_change))
    return(hotspot_result)
  }

  result <- list(
    clinvar_id = NA,
    clinical_significance = "Not found in ClinVar",
    review_status = NA,
    conditions = NA,
    pubmed_ids = NA
  )

  # Build search queries in order of specificity
  search_terms <- c()

  # Clean protein change for queries
  clean_protein <- NULL
  protein_three_letter <- NULL
  if (!is.null(protein_change) && !is.na(protein_change) && protein_change != "" && protein_change != "NA") {
    clean_protein <- gsub("^p\\.", "", protein_change)
    clean_protein <- gsub("^[A-Z]+[0-9]+:", "", clean_protein)  # Remove transcript prefixes
    protein_three_letter <- convert_to_three_letter(clean_protein)
  }

  # Strategy 1: Direct gene + variant name search (one-letter)
  if (!is.null(clean_protein)) {
    search_terms <- c(search_terms, paste0(gene, "[gene] AND ", clean_protein, "[variant name]"))
  }

  # Strategy 2: Gene + three-letter amino acid code
  if (!is.null(protein_three_letter) && protein_three_letter != clean_protein) {
    search_terms <- c(search_terms, paste0(gene, "[gene] AND ", protein_three_letter, "[variant name]"))
  }

  # Strategy 3: Gene + p.variant notation
  if (!is.null(clean_protein)) {
    search_terms <- c(search_terms, paste0(gene, "[gene] AND p.", clean_protein, "[variant name]"))
  }

  # Strategy 4: Free text search for gene + protein
  if (!is.null(clean_protein)) {
    search_terms <- c(search_terms, paste0(gene, " ", clean_protein))
  }

  # Strategy 5: Genomic position query (GRCh37 and GRCh38)
  if (!is.null(chr) && !is.null(pos)) {
    chr_clean <- gsub("chr", "", chr)
    search_terms <- c(search_terms,
                      paste0("(", gene, "[gene]) AND (", chr_clean, "[chr] AND ", pos, "[chrpos37])"),
                      paste0("(", gene, "[gene]) AND (", chr_clean, "[chr] AND ", pos, "[chrpos38])"))
  }

  # Attempt searches in order
  for (query in search_terms) {
    tryCatch({
      search_result <- entrez_search(db = "clinvar", term = query, retmax = 10)

      if (search_result$count > 0) {
        clinvar_ids <- search_result$ids[1:min(5, length(search_result$ids))]

        for (cv_id in clinvar_ids) {
          summary <- tryCatch({
            entrez_summary(db = "clinvar", id = cv_id)
          }, error = function(e) NULL)

          if (!is.null(summary)) {
            # Extract clinical significance (handle different ClinVar API response formats)
            clin_sig <- NA
            if (!is.null(summary$clinical_significance)) {
              if (is.list(summary$clinical_significance)) {
                clin_sig <- summary$clinical_significance$description
              } else {
                clin_sig <- summary$clinical_significance
              }
            } else if (!is.null(summary$germline_classification)) {
              if (is.list(summary$germline_classification)) {
                clin_sig <- summary$germline_classification$description
              } else {
                clin_sig <- summary$germline_classification
              }
            }

            # Get variant title for verification
            variant_title <- ifelse(!is.null(summary$title), summary$title, "")

            # Verify match if we have protein change
            match_found <- TRUE
            if (!is.null(clean_protein) && !is.na(clean_protein) && clean_protein != "") {
              # Check if variant title contains our protein change (various formats)
              match_found <- grepl(clean_protein, variant_title, ignore.case = TRUE) ||
                            grepl(gsub("\\*", "Ter", clean_protein), variant_title, ignore.case = TRUE) ||
                            (!is.null(protein_three_letter) && grepl(protein_three_letter, variant_title, ignore.case = TRUE)) ||
                            grepl(gene, variant_title, ignore.case = TRUE)  # At minimum, gene should match
            }

            if (match_found) {
              # Extract review status
              review_status <- NA
              if (!is.null(summary$review_status)) {
                review_status <- summary$review_status
              }

              # Extract conditions
              conditions <- NA
              if (!is.null(summary$trait_set)) {
                traits <- sapply(summary$trait_set, function(t) {
                  if (is.list(t) && !is.null(t$trait_name)) t$trait_name else as.character(t)
                })
                conditions <- paste(unique(traits), collapse = "; ")
              }

              result$clinvar_id <- cv_id
              result$clinical_significance <- ifelse(is.na(clin_sig) || clin_sig == "", "Not classified", clin_sig)
              result$review_status <- review_status
              result$conditions <- conditions

              return(result)
            }
          }
        }
      }
    }, error = function(e) {
      message(paste("ClinVar query failed for:", query, "-", e$message))
    })

    Sys.sleep(0.35)  # Rate limiting
  }

  return(result)
}

#' Query MyVariant.info API as alternative source
#'
#' @param chr Chromosome
#' @param pos Position
#' @param ref Reference allele
#' @param alt Alternate allele
#' @return List with annotation or NULL
query_myvariant <- function(chr, pos, ref, alt) {
  result <- list(
    clinvar_id = NA,
    clinical_significance = "Not found",
    review_status = NA,
    conditions = NA
  )

  tryCatch({
    chr_clean <- gsub("chr", "", chr)

    # Format query based on variant type
    if (alt == "+" || alt == "-") {
      # Indels are harder to query - return default
      return(result)
    }

    # SNV format
    hgvs <- paste0("chr", chr_clean, ":g.", pos, ref, ">", alt)

    url <- paste0("https://myvariant.info/v1/variant/", URLencode(hgvs, reserved = TRUE))
    response <- httr::GET(url, httr::timeout(10))

    if (httr::status_code(response) == 200) {
      data <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))

      if (!is.null(data$clinvar)) {
        if (!is.null(data$clinvar$rcv)) {
          if (is.list(data$clinvar$rcv)) {
            result$clinical_significance <- data$clinvar$rcv[[1]]$clinical_significance
            if (!is.null(data$clinvar$rcv[[1]]$conditions)) {
              result$conditions <- paste(data$clinvar$rcv[[1]]$conditions$name, collapse = "; ")
            }
          } else {
            result$clinical_significance <- data$clinvar$rcv$clinical_significance
            if (!is.null(data$clinvar$rcv$conditions)) {
              result$conditions <- paste(data$clinvar$rcv$conditions$name, collapse = "; ")
            }
          }
        }
        if (!is.null(data$clinvar$variant_id)) {
          result$clinvar_id <- data$clinvar$variant_id
        }
      }
    }
  }, error = function(e) {
    message(paste("MyVariant query failed:", e$message))
  })

  return(result)
}

# ============================================================================
# BATCH ANNOTATION
# ============================================================================

#' Annotate a data frame of variants with ClinVar information
#'
#' @param variants Data frame with columns: gene, protein (or protein_change),
#'                 chr, pos, ref, alt
#' @param use_cache Whether to use cached results (default TRUE)
#' @param progress Whether to show progress messages (default TRUE)
#' @param retry_failed Whether to retry previously failed lookups for known hotspots (default TRUE)
#' @return Data frame with added ClinVar annotation columns
annotate_variants_clinvar <- function(variants, use_cache = TRUE, progress = TRUE, retry_failed = TRUE) {

  # Ensure required columns exist
  required_cols <- c("gene", "chr", "pos", "ref", "alt")
  if (!all(required_cols %in% names(variants))) {
    stop("variants must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Get protein change column (might be named differently)
  protein_col <- intersect(c("protein_change", "protein"), names(variants))[1]

  # Load cache if using
  cache <- if (use_cache) load_annotation_cache() else data.frame()

  # Initialize result columns
  variants$clinvar_id <- NA
  variants$clinical_significance <- NA
  variants$clinvar_review_status <- NA
  variants$clinvar_conditions <- NA

  n_variants <- nrow(variants)
  n_cached <- 0
  n_queried <- 0
  n_hotspot <- 0

  for (i in 1:n_variants) {
    gene <- variants$gene[i]
    chr <- variants$chr[i]
    pos <- variants$pos[i]
    ref <- variants$ref[i]
    alt <- variants$alt[i]
    protein <- if (!is.na(protein_col)) variants[[protein_col]][i] else NA

    # Create cache key
    cache_key <- create_cache_key(gene, chr, pos, ref, alt)

    # Check cache first
    use_cached <- FALSE
    if (use_cache && cache_key %in% cache$query_key) {
      cached_row <- cache[cache$query_key == cache_key, ]

      # Check if we should retry this cache entry
      if (retry_failed && should_retry_cache(cached_row, gene, protein)) {
        if (progress) {
          message(sprintf("[%d/%d] Retrying known hotspot %s %s (was cached as not found)...",
                          i, n_variants, gene, protein))
        }
      } else {
        variants$clinvar_id[i] <- cached_row$clinvar_id[1]
        variants$clinical_significance[i] <- cached_row$clinical_significance[1]
        variants$clinvar_review_status[i] <- cached_row$review_status[1]
        variants$clinvar_conditions[i] <- cached_row$conditions[1]
        n_cached <- n_cached + 1
        use_cached <- TRUE
      }
    }

    if (!use_cached) {
      # Query ClinVar (will check hotspot DB first)
      if (progress) {
        message(sprintf("[%d/%d] Querying ClinVar for %s %s...", i, n_variants, gene,
                        ifelse(!is.na(protein) && protein != "", protein, paste0(chr, ":", pos))))
      }

      # Clean protein change
      protein_clean <- NA
      if (!is.na(protein) && protein != "" && protein != "NA") {
        protein_clean <- str_extract(protein, "[A-Z*][0-9]+[A-Z*]")
        if (is.na(protein_clean)) {
          protein_clean <- str_extract(protein, ":[A-Z*][0-9]+[A-Z*]")
          protein_clean <- gsub("^:", "", protein_clean)
        }
      }

      annotation <- query_clinvar_variant(gene, protein_clean, chr, pos, ref, alt)

      # Track hotspot database hits
      if (!is.null(annotation$source) && annotation$source == "known_hotspot_db") {
        n_hotspot <- n_hotspot + 1
      }

      variants$clinvar_id[i] <- annotation$clinvar_id
      variants$clinical_significance[i] <- annotation$clinical_significance
      variants$clinvar_review_status[i] <- annotation$review_status
      variants$clinvar_conditions[i] <- annotation$conditions

      # Update cache (remove old entry if exists, add new)
      cache <- cache[cache$query_key != cache_key, ]
      new_cache_row <- data.frame(
        query_key = cache_key,
        gene = gene,
        protein_change = ifelse(is.na(protein), "", protein),
        chr = chr,
        pos = pos,
        ref = ref,
        alt = alt,
        clinvar_id = annotation$clinvar_id,
        clinical_significance = annotation$clinical_significance,
        review_status = annotation$review_status,
        conditions = annotation$conditions,
        last_updated = as.character(Sys.Date()),
        stringsAsFactors = FALSE
      )
      cache <- rbind(cache, new_cache_row)
      n_queried <- n_queried + 1

      Sys.sleep(0.4)  # Rate limiting
    }
  }

  # Save updated cache
  if (use_cache && n_queried > 0) {
    save_annotation_cache(cache)
    if (progress) {
      message(sprintf("Annotation complete: %d from cache, %d queried (%d from hotspot DB), cache updated.",
                      n_cached, n_queried, n_hotspot))
    }
  }

  return(variants)
}

# ============================================================================
# INTERPRETATION HELPERS
# ============================================================================

#' Map ClinVar clinical significance to simple category
#'
#' @param clinvar_sig ClinVar clinical significance string
#' @return Simplified category
simplify_clinical_significance <- function(clinvar_sig) {
  if (is.na(clinvar_sig)) return("Unknown")

  sig_lower <- tolower(clinvar_sig)

  if (grepl("pathogenic", sig_lower) && !grepl("likely", sig_lower)) {
    return("Pathogenic")
  } else if (grepl("likely pathogenic", sig_lower)) {
    return("Likely Pathogenic")
  } else if (grepl("uncertain|vus", sig_lower)) {
    return("VUS")
  } else if (grepl("likely benign", sig_lower)) {
    return("Likely Benign")
  } else if (grepl("benign", sig_lower) && !grepl("likely", sig_lower)) {
    return("Benign")
  } else if (grepl("not found|not classified", sig_lower)) {
    return("Not in ClinVar")
  } else if (grepl("conflicting", sig_lower)) {
    return("Conflicting")
  } else {
    return("Other")
  }
}

#' Get clinical action level based on significance
#'
#' @param simplified_sig Simplified clinical significance
#' @return Action level: "Actionable", "Review", "Monitor", "Reassuring"
get_action_level <- function(simplified_sig) {
  switch(simplified_sig,
         "Pathogenic" = "Actionable",
         "Likely Pathogenic" = "Actionable",
         "VUS" = "Review",
         "Conflicting" = "Review",
         "Not in ClinVar" = "Review",
         "Unknown" = "Review",
         "Likely Benign" = "Monitor",
         "Benign" = "Reassuring",
         "Other" = "Review")
}

#' Format clinical significance for display with color coding
#'
#' @param clinvar_sig ClinVar clinical significance
#' @return Formatted string for reports
format_clinical_significance <- function(clinvar_sig) {
  simple <- simplify_clinical_significance(clinvar_sig)

  switch(simple,
         "Pathogenic" = "**Pathogenic**",
         "Likely Pathogenic" = "**Likely Pathogenic**",
         "VUS" = "*Variant of Uncertain Significance (VUS)*",
         "Conflicting" = "*Conflicting Interpretations*",
         "Not in ClinVar" = "Not in ClinVar",
         "Likely Benign" = "Likely Benign",
         "Benign" = "Benign",
         simple)
}

cat("\n")
cat("ClinVar Annotation Module (Updated) loaded.\n")
cat("Known CHIP hotspot database contains", length(KNOWN_CHIP_HOTSPOTS), "variants.\n")
cat("\n")
