#!/usr/bin/env Rscript
#
# create_genomic_loci.R
#
# Simple R script to create GenomicRiskLoci.txt using a PLINK reference.
# Requirements: R packages data.table, optparse, igraph
#              PLINK accessible as 'plink' (or change plink_cmd)
#
# Usage example:
# Rscript create_genomic_loci.R --bfile /path/to/ref --gwas input.snps \
#   --outdir out --leadP 5e-8 --gwasP 1e-5 --r2 0.6 --r2_2 0.1 \
#   --mergeDist 250 --maf 0.01 --refSNPs 1 --windowKb 500
#
suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(igraph)
})

option_list <- list(
  make_option(c("--bfile"), type="character", help="PLINK prefix for reference panel (bed/bim/fam)"),
  make_option(c("--gwas"), type="character", help="GWAS summary file (tab-delimited). Must have header and columns: chr, bp, p, rsID, a1, a2 (or similar)."),
  make_option(c("--outdir"), type="character", default=".", help="Output directory"),
  make_option(c("--leadP"), type="double", default=5e-8, help="P-value threshold to pick independent significant SNPs"),
  make_option(c("--gwasP"), type="double", default=1e-5, help="p-value cutoff to include GWAS SNPs inside LD window"),
  make_option(c("--r2"), type="double", default=0.6, help="LD r2 threshold for building candidate window"),
  make_option(c("--r2_2"), type="double", default=0.1, help="LD r2 threshold to group independent SNPs into same lead"),
  make_option(c("--mergeDist"), type="integer", default=250, help="Merge distance (kb) to merge nearby lead windows"),
  make_option(c("--maf"), type="double", default=0.01, help="Minimum MAF in reference panel"),
  make_option(c("--refSNPs"), type="integer", default=1, help="1 to include reference-only SNPs, 0 to exclude"),
  make_option(c("--windowKb"), type="integer", default=500, help="Max window (kb) to ask PLINK for LD around a SNP (practical cap)"),
  make_option(c("--threads"), type="integer", default=8, help="Number of threads to pass to PLINK"),
  make_option(c("--plink"), type="character", default="plink", help="PLINK command (default 'plink')")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Create outdir
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

plink_threads_arg <- function() {
  if (!is.null(opt$threads) && opt$threads > 1) {
    sprintf(" --threads %d", opt$threads)
  } else {
    ""
  }
}

# Helper: run plink --freq to make frequency file (if not present)
freq_prefix <- file.path(opt$outdir, "ref_freq")
freq_file <- paste0(freq_prefix, ".frq")
if (!file.exists(freq_file)) {
  message("Computing reference MAF with PLINK ...")
  cmd <- sprintf("%s --bfile %s%s --freq --out %s", opt$plink, opt$bfile, plink_threads_arg(), freq_prefix)
  system(cmd, intern = TRUE)
  if (!file.exists(freq_file)) stop("PLINK frequency file not created; check PLINK and bfile")
}
ref_freq <- fread(freq_file) # columns: CHR SNP A1 A2 MAF ...
setnames(ref_freq, old = c("CHR","SNP","A1","A2","MAF"), new = c("CHR","SNP","A1","A2","MAF"), skip_absent = TRUE)

# Read BIM for mapping (bim columns: chr, rsid, cm, bp, a1, a2)
bim_path <- paste0(opt$bfile, ".bim")
if (!file.exists(bim_path)) stop("Reference .bim not found")
bim <- fread(bim_path, header = FALSE)
setnames(bim, c("chr","rsid","cm","bp","a1","a2"))

# Identify duplicate rsIDs in the reference; PLINK cannot process them in LD mode
dup_rsids <- unique(bim$rsid[duplicated(bim$rsid)])
if (length(dup_rsids)) {
  message(sprintf("Reference panel has %d duplicated rsIDs; excluding them from LD queries", length(dup_rsids)))
}
valid_bim <- bim[!duplicated(rsid)]
if (length(dup_rsids)) {
  valid_bim <- valid_bim[!(rsid %in% dup_rsids)]
}

# Build a quick index: map rsid -> chr, pos, alleles (only unique IDs)
bim_index <- data.table(rsid = valid_bim$rsid, chr = valid_bim$chr, bp = valid_bim$bp, a1 = valid_bim$a1, a2 = valid_bim$a2)
setkey(bim_index, rsid)
valid_rsids <- bim_index$rsid

# Read GWAS summary file
gwas <- fread(opt$gwas)
# Normalize column names: try to find required columns
# Common columns: chr, bp (pos), p, rsID, a1, a2
colnames(gwas) <- gsub("\ufeff", "", colnames(gwas), fixed = TRUE)
colnames(gwas) <- trimws(tolower(colnames(gwas)))

normalize_col <- function(dt, candidates, target) {
  hits <- intersect(trimws(tolower(candidates)), colnames(dt))
  if (!length(hits)) return(FALSE)
  hit <- hits[1]
  if (hit != target) setnames(dt, hit, target)
  TRUE
}

if (!normalize_col(gwas, c("chr","chrom","chromosome"), "chr")) {
  stop("GWAS file must contain a chromosome column (chr/chrom/chromosome)")
}
if (!normalize_col(gwas, c("bp","pos","position","bp_hg19"), "bp")) {
  stop("GWAS file must contain a base-pair position column (bp/pos/position)")
}
if (!normalize_col(gwas, c("p","pval","pvalue"), "p")) {
  stop("GWAS file must contain a P-value column (p/pval/pvalue)")
}
normalize_col(gwas, c("a1","allele1","effect_allele","ea"), "a1")
normalize_col(gwas, c("a2","allele2","other_allele","oa"), "a2")

# attempt to find rsid column
rs_col <- intersect(c("rsid","snp","rs_id","id","markername"), colnames(gwas))
if (!length(rs_col)) stop("GWAS file must contain an rsID column (rsid/snp/id/markername)")
rs_col <- rs_col[1]

message("Detected GWAS columns: ", paste(colnames(gwas), collapse = ", "))

original_rs_col <- rs_col
rsid_cols <- unique(c(original_rs_col, "markername", "snp", "rsid", "id", "variant_id", "variant"))
rsid_cols <- rsid_cols[rsid_cols %in% colnames(gwas)]
gwas[, rsid_match := NA_character_]
for (col in rsid_cols) {
  vals <- trimws(as.character(gwas[[col]]))
  vals[vals == ""] <- NA_character_
  idx <- which(is.na(gwas$rsid_match) & !is.na(vals))
  if (length(idx)) set(gwas, i = idx, j = "rsid_match", value = vals[idx])
}
gwas[rsid_match == "", rsid_match := NA_character_]
if (!any(!is.na(gwas$rsid_match))) {
  stop("Unable to determine rsIDs from GWAS columns; please supply rsid/snp/markername data")
}
rs_col <- "rsid_match"
rsid_candidate_cols <- unique(c("rsid_match", rsid_cols))

# sanitize values
chr_vec <- gwas[["chr"]]
chr_vec <- gsub("^chr", "", chr_vec, ignore.case = TRUE)
special_idx <- chr_vec %in% c("x","y","mt","m")
chr_vec[special_idx] <- toupper(chr_vec[special_idx])
supp_chr <- suppressWarnings(as.integer(chr_vec))
chr_vec <- ifelse(is.na(supp_chr), chr_vec, as.character(supp_chr))
gwas[, chr := chr_vec]
gwas[, bp := suppressWarnings(as.integer(bp))]
if (anyNA(gwas$bp)) {
  stop("Some GWAS bp values are missing or non-integer after coercion")
}
gwas[, p := as.numeric(p)]
if (anyNA(gwas$p)) {
  stop("Some GWAS P-values are missing or non-numeric after coercion")
}

gwas[, resolved_rs := ifelse(rsid_match %chin% valid_rsids, rsid_match, NA_character_)]
matched_count <- sum(!is.na(gwas$resolved_rs))
total_rows <- nrow(gwas)
if (matched_count == 0) {
  stop("None of the GWAS SNP IDs were found in the reference BIM file")
}
if (matched_count < nrow(gwas)) {
  message(sprintf("Retaining %d GWAS rows that match the reference (dropping %d unmatched)", matched_count, total_rows - matched_count))
}
gwas <- gwas[!is.na(resolved_rs)]

# ensure sorted by chr and pos
setkeyv(gwas, c("chr","bp"))

# Quick helper: lookup MAF in ref (by rsid), fallback to bim lookup (rare)
lookup_maf <- function(rsid, pos, chr) {
  rf <- ref_freq[SNP == rsid, MAF]
  if (length(rf) == 1) return(as.numeric(rf))
  # try bim (no MAF in bim) -> NA
  return(NA_real_)
}

# Containers for outputs
ld_rows <- list()
snps_rows <- list()
indsig_rows <- list()
lead_rows <- list()
loci_rows <- list()

# Keep track of assigned uniqIDs
assigned_uids <- new.env(hash = TRUE)

# Helper to make uniqID as "chr:pos:alleleA:alleleB" with alleles sorted
make_uid <- function(chr, pos, a1, a2) {
  alleles <- sort(c(as.character(a1), as.character(a2)))
  paste0(chr, ":", pos, ":", alleles[1], ":", alleles[2])
}

get_candidate_ids <- function(row_dt) {
  ids <- character()
  if ("resolved_rs" %in% colnames(row_dt)) {
    ids <- c(ids, as.character(row_dt[["resolved_rs"]]))
  }
  for (col in rsid_candidate_cols) {
    if (!col %in% colnames(row_dt)) next
    val <- as.character(row_dt[[col]])
    ids <- c(ids, val)
  }
  ids <- trimws(ids)
  ids[ids == ""] <- NA_character_
  ids <- ids[!is.na(ids)]
  unique(ids)
}

# Helper to call plink --ld-snp-list for a set of SNPs (batch LD)
# Returns data.table with columns: SNP_A, BP_A, SNP_B, BP_B, R2
plink_ld_multi <- function(rsids, window_kb, r2cut) {
  rsids <- unique(rsids)
  rsids <- rsids[!is.na(rsids)]
  if (!length(rsids)) return(data.table())
  tmp_pref <- tempfile(pattern = "plink_ld_batch_")
  snp_list <- paste0(tmp_pref, ".snplist")
  fwrite(data.table(rsids), snp_list, col.names = FALSE)
  cmd <- sprintf("%s --bfile %s%s --ld-snp-list %s --ld-window-kb %d --ld-window-r2 %g --r2 --out %s",
                 opt$plink, opt$bfile, plink_threads_arg(), snp_list, window_kb, r2cut, tmp_pref)
  plink_out <- tryCatch(system(cmd, intern = TRUE), warning = function(w) w, error = function(e) e)
  if (inherits(plink_out, c("warning","error"))) {
    warning(sprintf("PLINK batch LD command failed: %s", conditionMessage(plink_out)))
    unlink(Sys.glob(paste0(tmp_pref, "*")), recursive = TRUE, force = TRUE)
    return(data.table())
  }
  ldfile <- paste0(tmp_pref, ".ld")
  if (!file.exists(ldfile)) {
    unlink(Sys.glob(paste0(tmp_pref, "*")), recursive = TRUE, force = TRUE)
    return(data.table())
  }
  ld <- tryCatch(fread(ldfile), error = function(e) data.table())
  keep <- intersect(c("SNP_A","BP_A","SNP_B","BP_B","R2"), colnames(ld))
  if (nrow(ld) == 0 || length(keep) < 5) {
    unlink(Sys.glob(paste0(tmp_pref, "*")), recursive = TRUE, force = TRUE)
    return(data.table())
  }
  ld <- ld[, .(SNP_A, BP_A, SNP_B, BP_B, R2)]
  unlink(Sys.glob(paste0(tmp_pref, "*")), recursive = TRUE, force = TRUE)
  ld
}

plink_ld_snps <- function(rsid, window_kb, r2cut) {
  ld <- plink_ld_multi(rsid, window_kb, r2cut)
  if (!nrow(ld)) return(ld)
  ld[SNP_A == rsid | SNP_B == rsid]
}

# For performance, define a max window in kb (opt$windowKb). Python code used dynamic window size per LD region,
# we will use this practical cap to limit PLINK queries.

# MAIN per-chromosome processing
chroms <- unique(gwas$chr)
chroms <- chroms[order(suppressWarnings(as.integer(chroms)), chroms)]
gidx <- 0
IndSigIdx <- 0
leadIdx <- 0

for (chrom in chroms) {
  message(sprintf("Processing chromosome %s ...", chrom))
  gwas_chr <- gwas[chr == chrom]
  if (nrow(gwas_chr) == 0) next

  # Optionally exclude MHC (user can pass extMHC via change if needed) -- not implemented extra here

  # Identify candidate independent significant SNPs by scanning ordered by P
  # Sort by p ascending
  gwas_ord <- gwas_chr[order(p, bp)]
  lead_rsids_chr <- unique(gwas_ord[p < opt$leadP, resolved_rs])
  lead_rsids_chr <- lead_rsids_chr[!is.na(lead_rsids_chr)]
  ld_pool <- plink_ld_multi(lead_rsids_chr, opt$windowKb, opt$r2)
  ld_expanded <- data.table(lead = character(), lead_bp = integer(), partner = character(), partner_bp = integer(), R2 = numeric())
  if (nrow(ld_pool) > 0) {
    ld_expanded <- rbind(
      ld_pool[, .(lead = SNP_A, lead_bp = BP_A, partner = SNP_B, partner_bp = BP_B, R2)],
      ld_pool[, .(lead = SNP_B, lead_bp = BP_B, partner = SNP_A, partner_bp = BP_A, R2)]
    )
    ld_expanded <- unique(ld_expanded, by = c("lead","partner"))
  }
  setkey(ld_expanded, lead)
  # iterate SNPs with p < leadP
  candidate_inds <- list()
  for (i in seq_len(nrow(gwas_ord))) {
    row <- gwas_ord[i]
    if (is.na(row$p) || row$p >= opt$leadP) break
    lead_ids <- get_candidate_ids(row)
    lead_hits <- lead_ids[lead_ids %chin% bim_index$rsid]
    if (!length(lead_hits)) {
      message(sprintf("Skipping chr%s:%s (no matching SNP in reference)", row$chr, row$bp))
      next
    }
    lead_rs <- lead_hits[1]
    row_copy <- copy(row)
    row_copy[[rs_col]] <- lead_rs
    row_copy$resolved_rs <- lead_rs
    # make uid
    a1 <- if ("a1" %in% colnames(gwas_ord)) row$a1 else NA
    a2 <- if ("a2" %in% colnames(gwas_ord)) row$a2 else NA
    uid <- make_uid(row$chr, row$bp, ifelse(is.na(a1), "NA", a1), ifelse(is.na(a2), "NA", a2))
    if (exists(uid, envir = assigned_uids, inherits = FALSE)) next

    # check MAF in reference
    maf_ref <- lookup_maf(lead_rs, row$bp, row$chr)
    if (!is.na(maf_ref) && maf_ref < opt$maf) {
      next
    }
    # Use plink to get LD partners within windowKb and r2
    # We pass windowKb cap; PLINK will return pairs in that window
    ld_tab <- ld_expanded[.(lead_rs)]
    if (nrow(ld_tab) == 0) {
      # no LD partners above threshold; still record singleton if desired
      # treat the SNP itself as candidate (if p < gwasP or we treat lead itself)
      candidate_list <- data.table(
        uid = uid,
        rsID = lead_rs,
        chr = row$chr,
        bp = row$bp,
        a1 = ifelse("a1" %in% colnames(row), as.character(row$a1), NA_character_),
        a2 = ifelse("a2" %in% colnames(row), as.character(row$a2), NA_character_),
        MAF = ifelse(is.na(maf_ref), NA_real_, maf_ref),
        gwasP = row$p,
        topLead = lead_rs
      )
      candidate_inds[[length(candidate_inds)+1]] <- list(ind_uid = uid, ind_row = row_copy, resolved_rs = lead_rs, candidates = candidate_list, ld_pairs = data.table())
      # mark assigned so we don't re-use
      assign(uid, TRUE, envir = assigned_uids)
      next
    }

    # Build partner list: from ld_tab, collect SNP_B (other SNPs) and SNP_A as needed
    # Combine SNP_A and SNP_B rows that involve our lead
    # rows where SNP_A == lead.rsID or SNP_B == lead.rsID
    partners <- ld_tab
    partners_vec <- unique(partners$partner)
    # Now we want to include GWAS-tagged partners with gwasP < gwasP, and optionally ref-only ones
    candidates_dt <- data.table()
    for (p in partners_vec) {
      # find partner position in bim_index or in GWAS
      # first check GWAS
      gwas_match <- gwas_ord[get(rs_col) == p]
      if (nrow(gwas_match) == 0) {
        for (alt_col in setdiff(rsid_candidate_cols, rs_col)) {
          if (!alt_col %in% colnames(gwas_ord)) next
          alt_vals <- trimws(as.character(gwas_ord[[alt_col]]))
          idx_alt <- which(alt_vals == p)
          if (length(idx_alt)) {
            gwas_match <- gwas_ord[idx_alt]
            break
          }
        }
      }
      if (nrow(gwas_match) > 0) {
        # choose the GWAS row - if multiple at same position, pick the one matching alleles ideally - skip detailed allele matching here
        prow <- gwas_match[1]
        if (!is.na(prow$p) && prow$p < opt$gwasP) {
          # include
          partner_uid <- make_uid(prow$chr, prow$bp, ifelse("a1" %in% colnames(prow), prow$a1, NA), ifelse("a2" %in% colnames(prow), prow$a2, NA))
          candidates_dt <- rbind(candidates_dt, data.table(
            uid = partner_uid,
            rsID = p,
            chr = prow$chr,
            bp = prow$bp,
            a1 = ifelse("a1" %in% colnames(prow), as.character(prow$a1), NA_character_),
            a2 = ifelse("a2" %in% colnames(prow), as.character(prow$a2), NA_character_),
            MAF = lookup_maf(p, prow$bp, prow$chr),
            gwasP = prow$p
          ), fill = TRUE)
        } else {
          # not significant in GWAS for inclusion
          if (opt$refSNPs == 1) {
            # include as ref-only? but it's present in GWAS, so skip
          }
        }
      } else {
        # not in GWAS; if refSNPs==1 include from bim/ref_freq
        if (opt$refSNPs == 1) {
          # try to find in bim_index and freq
          if (p %in% bim_index$rsid) {
            brow <- bim_index[p]
            partner_uid <- make_uid(brow$chr, brow$bp, brow$a1, brow$a2)
            candidates_dt <- rbind(candidates_dt, data.table(
              uid = partner_uid,
              rsID = p,
              chr = brow$chr,
              bp = brow$bp,
              a1 = brow$a1,
              a2 = brow$a2,
              MAF = lookup_maf(p, brow$bp, brow$chr),
              gwasP = NA_real_
            ), fill = TRUE)
          } else {
            # try to find in ref_freq by SNP name
            if (p %in% ref_freq$SNP) {
              rfrow <- ref_freq[SNP == p]
              # no bp/chrom in frq? frq has CHR; but should have, use CHR and SNP
              partner_uid <- paste0(rfrow$CHR, ":", NA, ":", NA, ":", NA)
              candidates_dt <- rbind(candidates_dt, data.table(
                uid = partner_uid,
                rsID = p,
                chr = rfrow$CHR,
                bp = NA_integer_,
                a1 = rfrow$A1,
                a2 = rfrow$A2,
                MAF = rfrow$MAF,
                gwasP = NA_real_
              ), fill = TRUE)
            }
          }
        }
      }
    }

    # also include the lead itself as candidate (topLead)
    lead_candidate <- data.table(
      uid = uid,
      rsID = lead_rs,
      chr = row$chr,
      bp = row$bp,
      a1 = ifelse("a1" %in% colnames(row), as.character(row$a1), NA_character_),
      a2 = ifelse("a2" %in% colnames(row), as.character(row$a2), NA_character_),
      MAF = ifelse(is.na(maf_ref), lookup_maf(lead_rs, row$bp, row$chr), maf_ref),
      gwasP = row$p
    )
    candidates_dt <- rbind(lead_candidate, unique(candidates_dt), fill = TRUE)

    # Build LD pair table (leader -> partner with r2)
    # Use partners dt: map partner to r2
    # ensure partner entries exist; create pair rows with uid format
    ld_pairs_dt <- data.table()
    for (j in seq_len(nrow(partners))) {
      prow <- partners[j]
      partner_rs <- as.character(prow$partner)
      bp_partner <- ifelse(!is.null(prow$partner_bp) && !is.na(prow$partner_bp), as.integer(prow$partner_bp), NA_integer_)
      partner_chr <- chrom
      partner_a1 <- NA
      partner_a2 <- NA
      if (partner_rs %in% bim_index$rsid) {
        brow <- bim_index[partner_rs]
        if (!is.null(brow$bp) && !is.na(brow$bp)) bp_partner <- as.integer(brow$bp)
        if (!is.null(brow$chr) && !is.na(brow$chr)) partner_chr <- brow$chr
        if (!is.null(brow$a1)) partner_a1 <- brow$a1
        if (!is.null(brow$a2)) partner_a2 <- brow$a2
      }
      partner_uid <- make_uid(
        partner_chr,
        ifelse(is.na(bp_partner), NA, bp_partner),
        partner_a1,
        partner_a2
      )
      ld_pairs_dt <- rbind(ld_pairs_dt, data.table(SNP1 = uid, SNP2 = partner_uid, R2 = as.numeric(prow$R2)), fill = TRUE)
    }

    # Save candidate index entry
    candidate_inds[[length(candidate_inds) + 1]] <- list(ind_uid = uid, ind_row = row_copy, resolved_rs = lead_rs, candidates = unique(candidates_dt, by = "uid"), ld_pairs = unique(ld_pairs_dt))

    # mark the candidate SNPs as assigned (so future independent picks skip)
    for (cuid in unique(candidates_dt$uid)) assign(cuid, TRUE, envir = assigned_uids)
  } # end scanning GWAS p-order

  # Convert candidate_inds into data.tables aggregated
  if (length(candidate_inds) == 0) next

  # Build overall ld_pairs, candidate_snps and IndSigSNPs tables for this chromosome
  chrom_ld_pairs <- rbindlist(lapply(candidate_inds, function(x) x$ld_pairs), fill = TRUE)
  chrom_candidates <- rbindlist(lapply(candidate_inds, function(x) {
    # add a column topLead to indicate source ind
    x$candidates[, topLead := x$resolved_rs]
    x$candidates
  }), use.names = TRUE, fill = TRUE)
  chrom_indsig <- rbindlist(lapply(seq_along(candidate_inds), function(ii) {
    x <- candidate_inds[[ii]]
    nSNPs <- nrow(x$candidates)
    nGWASSNPs <- sum(!is.na(x$candidates$gwasP))
    data.table(
      ind_uid = x$ind_uid,
      rsID = as.character(x$resolved_rs),
      chr = x$ind_row$chr,
      pos = x$ind_row$bp,
      p = x$ind_row$p,
      nSNPs = nSNPs,
      nGWASSNPs = nGWASSNPs
    )
  }), fill = TRUE)

  # append to overall outputs
  if (nrow(chrom_ld_pairs) > 0) ld_rows[[length(ld_rows)+1]] <- chrom_ld_pairs
  if (nrow(chrom_candidates) > 0) snps_rows[[length(snps_rows)+1]] <- chrom_candidates
  if (nrow(chrom_indsig) > 0) indsig_rows[[length(indsig_rows)+1]] <- chrom_indsig

  # Build adjacency graph of IndSigSNPs where edges exist if r2 >= r2_2
  # Use chrom_ld_pairs: map SNP2 uid -> see if SNP2 corresponds to another ind_uid
  # Create mapping from ind_uid to index
  ind_uids <- chrom_indsig$ind_uid
  ind_idx <- seq_along(ind_uids)
  uid_to_idx <- as.list(ind_idx)
  names(uid_to_idx) <- ind_uids
  edges <- list()
  if (nrow(chrom_ld_pairs) > 0) {
    for (r in seq_len(nrow(chrom_ld_pairs))) {
      rowr <- chrom_ld_pairs[r]
      # check if SNP2 corresponds to an ind_uid
      # here rowr$SNP1 is lead UID (uid) we made earlier; SNP2 is partner uid; check both directions
      s1 <- rowr$SNP1; s2 <- rowr$SNP2; r2val <- as.numeric(rowr$R2)
      # if both are ind_uids and r2 >= r2_2 -> edge
      if (!is.null(uid_to_idx[[s1]]) && !is.null(uid_to_idx[[s2]])) {
        if (!is.na(r2val) && r2val >= opt$r2_2) {
          edges[[length(edges)+1]] <- c(uid_to_idx[[s1]], uid_to_idx[[s2]])
        }
      }
    }
  }
  # Build graph
  if (length(edges) == 0) {
    # Each ind becomes its own lead
    comps <- as.list(seq_along(ind_uids))
  } else {
    g <- graph_from_edgelist(do.call(rbind, edges), directed = FALSE)
    compsf <- components(g)
    # get membership groups
    member_map <- split(seq_along(ind_uids), compsf$membership)
    comps <- member_map
  }

  # For each component produce a leadSNP record: choose top SNP (smallest p)
  chrom_leads <- list()
  for (comp in comps) {
    comp_inds <- chrom_indsig[comp, ]
    # pick top SNP by p
    top_i <- which.min(comp_inds$p)
    top_uid <- comp_inds$ind_uid[top_i]
    top_rs <- comp_inds$rsID[top_i]
    top_p <- comp_inds$p[top_i]
    top_pos <- comp_inds$pos[top_i]
    # collect all IndSig rsIDs in this lead
    inds_rs <- paste(comp_inds$rsID, collapse = ";")
    nInd <- nrow(comp_inds)
    # collect candidate SNPs that have topLead == any of these ind rsIDs
    cands <- chrom_candidates[topLead %in% comp_inds$rsID]
    # compute start/end
    start_pos <- suppressWarnings(min(cands$bp, na.rm = TRUE))
    end_pos <- suppressWarnings(max(cands$bp, na.rm = TRUE))
    if (!is.finite(start_pos)) {
      start_pos <- NA_integer_
    } else {
      start_pos <- as.integer(start_pos)
    }
    if (!is.finite(end_pos)) {
      end_pos <- NA_integer_
    } else {
      end_pos <- as.integer(end_pos)
    }
    nSNPs <- length(unique(cands$uid))
    nGWASSNPs <- sum(!is.na(cands$gwasP))
    lead_entry <- data.table(
      lead_uid = top_uid,
      lead_rsID = top_rs,
      chr = chrom,
      pos = top_pos,
      p = top_p,
      nIndSigSNPs = nInd,
      IndSigSNPs = inds_rs,
      start = start_pos,
      end = end_pos,
      nSNPs = nSNPs,
      nGWASSNPs = nGWASSNPs,
      candidates = list(unique(cands$uid))
    )
    chrom_leads[[length(chrom_leads)+1]] <- lead_entry
  }
  chrom_leads_dt <- rbindlist(chrom_leads, fill = TRUE)

  # MERGE leads into genomic loci on this chromosome using mergeDist
  if (nrow(chrom_leads_dt) == 0) next
  # Convert mergeDist kb to bp
  merge_bp <- opt$mergeDist * 1000
  # Sort by candidate start
  setorder(chrom_leads_dt, start)
  loci_list <- list()
  cur <- chrom_leads_dt[1]
  cur_candidates <- unlist(cur$candidates)
  cur_inds <- cur$IndSigSNPs
  cur_leads <- cur$lead_rsID
  cur_start <- cur$start
  cur_end <- cur$end
  cur_nSNPs <- cur$nSNPs
  cur_nGWASSNPs <- cur$nGWASSNPs
  cur_top_p <- cur$p
  cur_top_uid <- cur$lead_uid
  cur_top_rs <- cur$lead_rsID

  if (nrow(chrom_leads_dt) > 1) {
    for (li in 2:nrow(chrom_leads_dt)) {
      nxt <- chrom_leads_dt[li]
      # if overlap or gap < merge_bp, merge
      if ((!is.na(nxt$start) && !is.na(cur_end) && nxt$start <= cur_end) || ( (!is.na(nxt$start) && !is.na(cur_end)) && (nxt$start - cur_end < merge_bp) )) {
        # merge
        cur_candidates <- unique(c(cur_candidates, unlist(nxt$candidates)))
        cur_inds <- paste(unique(c(unlist(strsplit(cur_inds, ";")), unlist(strsplit(nxt$IndSigSNPs, ";")))), collapse = ";")
        cur_leads <- paste(unique(c(unlist(strsplit(cur_leads, ";")), nxt$lead_rsID)), collapse = ";")
        cur_start <- min(cur_start, nxt$start, na.rm = TRUE)
        cur_end <- max(cur_end, nxt$end, na.rm = TRUE)
        cur_nSNPs <- length(cur_candidates)
        cur_nGWASSNPs <- NA_integer_ # recompute if needed from chrom_candidates
        # update top lead if nxt has smaller p
        if (!is.na(nxt$p) && nxt$p < cur_top_p) {
          cur_top_p <- nxt$p
          cur_top_uid <- nxt$lead_uid
          cur_top_rs <- nxt$lead_rsID
        }
      } else {
        # finalize current locus
        loci_list[[length(loci_list) + 1]] <- list(
          GenomicLocus = NA_integer_,
          uniqID = cur_top_uid,
          rsID = cur_top_rs,
          chr = chrom,
          pos = NA_integer_,
          p = cur_top_p,
          start = cur_start,
          end = cur_end,
          nSNPs = cur_nSNPs,
          nGWASSNPs = NA_integer_,
          nIndSigSNPs = length(unlist(strsplit(cur_inds, ";"))),
          IndSigSNPs = cur_inds,
          nLeadSNPs = length(unlist(strsplit(cur_leads, ";"))),
          LeadSNPs = cur_leads,
          candidate_uids = list(cur_candidates)
        )
        # start new current
        cur <- nxt
        cur_candidates <- unlist(cur$candidates)
        cur_inds <- cur$IndSigSNPs
        cur_leads <- cur$lead_rsID
        cur_start <- cur$start
        cur_end <- cur$end
        cur_nSNPs <- cur$nSNPs
        cur_nGWASSNPs <- cur$nGWASSNPs
        cur_top_p <- cur$p
        cur_top_uid <- cur$lead_uid
        cur_top_rs <- cur$lead_rsID
      }
    }
  }
  # Add final current locus
  loci_list[[length(loci_list) + 1]] <- list(
    GenomicLocus = NA_integer_,
    uniqID = cur_top_uid,
    rsID = cur_top_rs,
    chr = chrom,
    pos = NA_integer_,
    p = cur_top_p,
    start = cur_start,
    end = cur_end,
    nSNPs = cur_nSNPs,
    nGWASSNPs = NA_integer_,
    nIndSigSNPs = length(unlist(strsplit(cur_inds, ";"))),
    IndSigSNPs = cur_inds,
    nLeadSNPs = length(unlist(strsplit(cur_leads, ";"))),
    LeadSNPs = cur_leads,
    candidate_uids = list(cur_candidates)
  )

  # Append loci to overall list; assign indices sequentially
  for (ll in seq_along(loci_list)) {
    gidx <- gidx + 1
    loci_entry <- loci_list[[ll]]
    loci_entry$GenomicLocus <- gidx
    loci_entry$pos <- NA_integer_ # top lead position could be resolved if needed
    # Compute nGWASSNPs by checking candidate_uids vs chrom_candidates
    cand_uids <- unique(unlist(loci_entry$candidate_uids))
    nGWAS <- sum(chrom_candidates$uid %in% cand_uids & !is.na(chrom_candidates$gwasP))
    loci_entry$nGWASSNPs <- nGWAS
    loci_rows[[length(loci_rows)+1]] <- loci_entry
  }

} # end chromosomes loop

# Combine and write GenomicRiskLoci.txt
if (length(loci_rows) == 0) {
  stop("No loci identified")
}
loci_dt <- rbindlist(lapply(loci_rows, function(x) {
  data.table(
    GenomicLocus = x$GenomicLocus,
    uniqID = x$uniqID,
    rsID = x$rsID,
    chr = x$chr,
    pos = x$pos,
    p = x$p,
    start = x$start,
    end = x$end,
    nSNPs = x$nSNPs,
    nGWASSNPs = x$nGWASSNPs,
    nIndSigSNPs = x$nIndSigSNPs,
    IndSigSNPs = x$IndSigSNPs,
    nLeadSNPs = x$nLeadSNPs,
    LeadSNPs = x$LeadSNPs
  )
}), fill = TRUE)

out_gl <- file.path(opt$outdir, "GenomicRiskLoci.txt")
fwrite(loci_dt, out_gl, sep = "\t", quote = FALSE, na = "NA")
message("Wrote GenomicRiskLoci to: ", out_gl)

# Optionally write other outputs if built
# (This script focuses on GenomicRiskLoci.txt; additional outputs can be written similarly.)

message("Done.")