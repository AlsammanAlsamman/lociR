Target: a clear, actionable description an AI agent can use to reason about, test, or modify the R script create_genomic_loci.R (the PLINK-based locus-identification script you asked for).

Short purpose
- The script identifies genomic risk loci from a GWAS summary file using LD information from a PLINK-format reference panel. It outputs GenomicRiskLoci.txt describing each locus (top lead SNP, locus bounds, counts, assigned independent SNPs, etc.). It implements the same logical steps as FUMA’s getLD.py but in a simplified, PLINK-driven R implementation.

High-level intent (why run it)
- Turn GWAS summary statistics + a reference panel into locus-level regions that aggregate nearby, LD-related association signals. These loci are suitable inputs for downstream annotation (ANNOVAR), gene-mapping and plotting pipelines.

Inputs (what the AI should check exist & are valid)
- --bfile: PLINK reference prefix (.bed/.bim/.fam). Used to compute LD and MAF via PLINK.
- --gwas: GWAS summary statistics (tab-delimited). Must contain columns: chr, bp (position), p (P-value) and an rsID column. Optionally a1/a2 columns for allele matching.
- Parameters (CLI options):
  - leadP (default 5e-8) — P cutoff to pick independent significant SNPs.
  - gwasP (default 1e-5) — P cutoff for including GWAS SNPs into a lead’s LD window.
  - r2 (default 0.6) — r2 threshold used to collect LD partners for each independent SNP.
  - r2_2 (default 0.1) — r2 threshold used to group independent significant SNPs into the same lead.
  - mergeDist (default 250 kb) — distance to merge adjacent lead windows into a single locus.
  - maf (default 0.01) — minimum MAF in reference to include a SNP.
  - refSNPs (0/1) — whether to include reference-only SNPs (not present in GWAS) in candidate lists.
  - windowKb (default 500) — maximum PLINK LD query window (practical cap).
  - plink command (default "plink").

Outputs (what the AI must expect)
- GenomicRiskLoci.txt (tab-delimited): each locus with columns such as GenomicLocus, uniqID (top lead), rsID, chr, pos, p, start, end, nSNPs, nGWASSNPs, nIndSigSNPs, IndSigSNPs, nLeadSNPs, LeadSNPs.
- The script builds intermediate in-memory tables for candidate SNPs, LD pairs and independent significant SNPs; optionally you can extend it to write ld.txt, snps.txt, IndSigSNPs.txt, leadSNPs.txt.

High-level algorithm (step-by-step)
1. Index/reference prep:
   - Generate reference allele frequencies using PLINK --freq (if not available).
   - Read the reference .bim to map rsIDs → positions and alleles.

2. Read and normalize GWAS:
   - Load GWAS summary stats and sort by chr, position and ascending P.

3. Find independent significant SNPs:
   - Iterate GWAS SNPs in ascending p; for each SNP with p < leadP that is not already assigned:
     - Check reference MAF (>= maf). Skip if below threshold.
     - Query PLINK (--ld-snp) around this SNP with windowKb and r2 to get LD partners (r2 >= r2 threshold).
     - Build the candidate set for this independent SNP: include (a) GWAS SNPs in LD with gwasP < gwasP, (b) optionally reference-only partners (refSNPs==1) from the reference.
     - Record IndSigSNP entry and its candidate SNPs and LD pairs. Mark candidates as assigned so later independent picks don’t reassign the same SNPs.

4. Group independent SNPs into lead SNPs:
   - Build a graph where nodes are IndSigSNPs and edges exist if the two IndSigSNPs are in LD ≥ r2_2 (based on LD pairs produced above).
   - Connected components of this graph become leadSNP groups. For each group choose the top SNP (smallest p) as the lead and collect its IndSigSNPs and candidate SNPs.

5. Merge leads into genomic loci:
   - For leadSNPs sorted by genomic location, compute each lead’s window start/end from its candidate SNPs.
   - If a lead’s start ≤ current locus end or the gap < mergeDist, merge it into the current locus: union candidate SNPs, union IndSigSNPs/LeadSNPs, update start/end and counts, and choose the locus top lead as the one with smallest p among merged leads.
   - If merging creates overlap with the previous locus within mergeDist, perform iterative merging.

6. Write GenomicRiskLoci.txt:
   - Output one line per final locus containing the fields described above.

Key implementation details and assumptions the AI should know
- LD queries are made via PLINK’s --ld-snp (one PLINK call per independent candidate). This is simple but may be slow when many independent SNPs exist. Alternative: precompute per-chromosome LD or use bigsnpr to compute correlations in R for speed.
- UID format: chr:pos:alleleA:alleleB with alleles alphabetically ordered. The script constructs UIDs for comparing and mapping candidate SNPs.
- Allele matching: current script uses rsID matching and BIM lookups; it only partially implements allele-disambiguation. For strict correctness you must compare alleles between GWAS and reference (and resolve strand flips).
- MHC handling: the simple script does not automatically exclude MHC unless you add that logic. In full pipelines MHC is typically excluded due to complex LD.
- Input GWAS should be sorted by chr and bp; the script re-sorts but efficient indexing assumes sorted order when doing binary-style lookups.

Failures, edge cases and what the AI should guard against
- PLINK availability / bfile correctness: ensure PLINK is installed and bfile exists with matching sample/population.
- Many PLINK calls: if many independent SNPs exist you may spawn many PLINK processes — monitor runtime and consider batching or precomputing LD.
- Missing rsIDs: if GWAS rsIDs don’t match reference rsIDs, candidate discovery may fail. The script attempts some fallback using BIM, but robust allele/ID harmonization is important.
- Multi-allelic positions or ambiguous alleles (A/T, C/G) require strand-aware matching.
- Positions missing in reference: ref-only SNP inclusion requires reliable mapping of rsID→pos in BIM or reference frequency files.

How an AI agent should interact with / use the script
- Pre-run checks (automatable):
  - Verify bfile (.bed/.bim/.fam) and plink command exist and are readable/executable.
  - Verify GWAS columns present (chr, bp, p, rsID) and that P values are numeric and within (0,1].
  - Validate parameters are in expected ranges (0 < r2 <= 1, 0 < leadP < 1, mergeDist > 0).
- Dry-run mode: run script with smaller windowKb and fewer chromosomes to validate behavior on a test dataset before full run.
- Logging: capture PLINK stdout/stderr and write per-candidate LD query logs for debugging.
- Testing: the AI should run toy unit tests:
  - Known small GWAS with a small reference (e.g., 1000 Genomes chunk) and check expected loci vs. gold-standard.
  - Test refSNPs toggles, mergeDist extremes, and r2/r2_2 thresholds.
- Optimization suggestions to implement automatically:
  - Parallelize PLINK calls when independent SNPs are many (careful with I/O).
  - Replace PLINK per-SNP calls by precomputing LD matrices for windows (bigsnpr) if memory allows.
  - Add robust allele harmonization (match alleles, detect strand flips, flip when needed).
  - Add explicit MHC exclusion option.

What the AI can do next (actionable items)
- Validate inputs and report any mismatches (rsIDs not found in reference; missing columns; invalid parameter values).
- Run a small smoke test with 1–2 chromosomes to confirm outputs.
- If slowness occurs, recommend converting to bigsnpr-based LD extraction or running PLINK in batched mode.
- If needed, extend the script to output the other pipeline files (snps.txt, IndSigSNPs.txt, leadSNPs.txt, ld.txt) in FUMA-compatible format.

Short summary the AI can ingest as metadata
- Name: create_genomic_loci.R
- Purpose: identify independent significant SNPs, group them into leads using LD thresholds, then merge leads into genomic risk loci; produce GenomicRiskLoci.txt for downstream annotation.
- External dependency: PLINK (for LD and MAF queries).
- Key parameters controlling behavior: leadP, gwasP, r2, r2_2, mergeDist, maf, refSNPs.

If you want, I can now:
- Produce a concise JSON “task spec” the AI can use to run tests and verify outputs (input paths, parameter values, expected checks).
- Or produce an extended version of the script that writes all FUMA-style outputs (snps.txt, IndSigSNPs.txt, leadSNPs.txt, ld.txt) and includes allele harmonization and logging. Which do you want?