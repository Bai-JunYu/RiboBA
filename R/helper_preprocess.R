#' Prepare annotation resources for RiboBC
#'
#' Builds a custom `BSgenome` package from the supplied genome FASTA/2bit file
#' and generates transcript-centric objects used by downstream analysis.
#'
#' @param output_dir Directory where intermediate and final artefacts are
#'   written.
#' @param scientific_name Optional scientific name of the organism in the form
#'   `"Genus.species"`. If `NULL`, inferred from file names when possible.
#' @param annotation_version Optional identifier for genome build/annotation
#'   release. If `NULL`, inferred from file names when possible.
#' @param genome_file Path to the reference genome in FASTA or 2bit format.
#' @param gff_file Path to the GFF/GTF annotation file.
#' @param orf_min_aa Minimum ORF length (in amino acids) considered during
#'   candidate generation.
#' @param orf_start_code Start codon sequence used when scanning for ORFs.
#'
#' @return Invisibly returns `NULL`; files are written to `output_dir`.
#'
#' @export
build_annotation_index <- function(
    output_dir,
    scientific_name = NULL,
    annotation_version = NULL,
    genome_file,
    gff_file,
    orf_min_aa = 7L,
    orf_start_code = "ATG") {
  infer_from_filename <- function(genome_file, gff_file) {
    candidates <- c(basename(genome_file), basename(gff_file))
    scientific <- NULL
    anno <- NULL

    for (nm in candidates) {
      no_ext <- sub("\\.(fa|fasta|2bit|gtf|gff3?|gz)$", "", nm, ignore.case = TRUE)
      if (is.null(scientific)) {
        m <- regmatches(no_ext, regexpr("^[A-Z][a-z]+[_\\.][a-z]+", no_ext))
        if (length(m) == 1 && nzchar(m)) {
          scientific <- gsub("_", ".", m)
        }
      }
      if (is.null(anno)) {
        m_anno <- regmatches(no_ext, regexpr("(GRCh[0-9]+|GRCm[0-9]+|Rnor_[0-9.]+|TAIR[0-9]+|WBcel[0-9]+|dm[0-9]+|ce[0-9]+|hg[0-9]+|mm[0-9]+)", no_ext))
        if (length(m_anno) == 1 && nzchar(m_anno)) {
          anno <- m_anno
        }
      }
    }
    if (is.null(scientific) || !grepl("^[A-Z][a-z]+\\.[a-z]+$", scientific)) {
      scientific <- "Unknown.species"
    }
    if (is.null(anno) || !nzchar(anno)) {
      anno <- "v1"
    }
    list(scientific_name = scientific, annotation_version = anno)
  }

  inferred <- infer_from_filename(genome_file = genome_file, gff_file = gff_file)
  if (is.null(scientific_name) || !nzchar(scientific_name)) {
    scientific_name <- inferred$scientific_name
    message("Inferred scientific_name: ", scientific_name)
  }
  if (is.null(annotation_version) || !nzchar(annotation_version)) {
    annotation_version <- inferred$annotation_version
    message("Inferred annotation_version: ", annotation_version)
  }

  # Create BSgenome package
  create_bs <- function(
      out_dir = ".",
      seq_file,
      scientific_name,
      anno_name) {
    seq_format <- base::rev(
      stringr::str_split(
        seq_file, "\\.",
        simplify = TRUE
      )
    )[1]
    if (seq_format %in% base::c("fa", "fasta")) {
      dna_set <- Biostrings::readDNAStringSet(seq_file)
      rtracklayer::export.2bit(
        object = Biostrings::replaceAmbiguities(dna_set),
        con = base::sub(
          x = seq_file,
          pattern = "\\.fa$|\\.fasta", "\\.2bit"
        )
      )
      seq_file <- base::sub(
        x = seq_file,
        pattern = "\\.fa$|\\.fasta", "\\.2bit"
      )
    }
    # Check parameter
    genome_name <- base::basename(seq_file)
    seq_format <- base::rev(
      stringr::str_split(genome_name, pattern = "\\.")[[1]]
    )[1]
    short_name <- base::paste0(
      base::substr(scientific_name, 1, 1),
      base::rev(
        stringr::str_split(
          scientific_name,
          pattern = "\\."
        )[[1]]
      )[1]
    )

    if (!(seq_format %in% base::c("2bit", "fa", "fasta"))) {
      base::stop("Genome file should be *.2bit, *.fa, or *.fasta file!")
    }

    # Create BSgenome package
    seed_text <- base::paste0(
      "Package: BSgenome.", short_name, ".", anno_name, ".Ribo", "\n",
      "Title: Full genome sequences for ", scientific_name, ",", anno_name, "\n",
      "Description: Full genome sequences for ", scientific_name, ",", anno_name, "\n",
      "Version: 0.0.1", "\n",
      "organism: ", scientific_name, "\n",
      "common_name: ", scientific_name, "\n",
      "provider: NA", "\n",
      "genome: ", anno_name, "\n",
      "provider_version: ", anno_name, "\n",
      "release_date: NA", "\n",
      "source_url: NA", "\n",
      "organism_biocview: ", scientific_name, "\n",
      "BSgenomeObjname: ", scientific_name, "\n",
      "seqs_srcdir: ", base::dirname(seq_file), "\n",
      "seqfile_name: ", base::basename(seq_file)
    )

    # Check circle sequence
    df_seqinfo <- GenomeInfoDb::seqinfo(rtracklayer::TwoBitFile(seq_file))
    id_cir <- base::which(
      df_seqinfo@seqnames %in%
        base::c(GenomeInfoDb::DEFAULT_CIRC_SEQS, base::c("Mt", "Pt"))
    )
    seq_cir <- df_seqinfo@seqnames[id_cir]
    df_seqinfo@is_circular[id_cir] <- TRUE
    if (base::length(seq_cir) == 1) {
      seed_text <- base::paste0(
        seed_text, "\n", "circ_seqs: \"", seq_cir, "\""
      )
    } else {
      if (base::length(seq_cir) > 1) {
        seed_text <- base::paste0(
          seed_text, "\n", "circ_seqs: ",
          base::paste0('c("', base::paste0(seq_cir, collapse = '","'), '")')
        )
      } else {
        seed_text <- base::paste0(
          seed_text, "\n", "circ_seqs: character(0)"
        )
      }
    }

    seed_dir <- base::paste0(
      out_dir, "/", scientific_name, ".", anno_name, "_seed"
    )
    base::writeLines(text = seed_text, con = seed_dir)
    # Building BSgenome
    message("Creating BSgenome package source ...")
    BSgenomeForge::forgeBSgenomeDataPkg(
      x = seed_dir,
      destdir = out_dir,
      seqs_srcdir = base::dirname(seq_file)
    )
    devtools::install(
      base::paste(
        out_dir,
        base::paste0("BSgenome.", short_name, ".", anno_name, ".Ribo"),
        sep = "/"
      ),
      upgrade = FALSE
    )
  }
  # get longest tx from txdb
  get_tx_info <- function(
      bs,
      txdb,
      start_atg = TRUE,
      common_stop = TRUE,
      no_frameshift = TRUE,
      min_cds = 30L) {
    # get transcript information
    tx_lens <- GenomicFeatures::transcriptLengths(
      txdb,
      with.utr5_len = TRUE,
      with.utr3_len = TRUE,
      with.cds_len = TRUE
    )

    cds_seqs <- GenomicFeatures::extractTranscriptSeqs(
      bs,
      GenomicFeatures::cdsBy(txdb,
        by = "tx",
        use.names = TRUE
      )
    )
    tx_seqs <- GenomicFeatures::extractTranscriptSeqs(
      bs,
      GenomicFeatures::exonsBy(txdb,
        by = c("tx", "gene"),
        use.names = TRUE
      )
    )

    # remove CDS which not ATG start,
    # not common stop codon, cds length not multiple of 3, less than 30 bp
    # in-frame stop codon(not implemented)
    idx_1 <- idx_2 <- idx_3 <- idx_4 <- vector(length = length(cds_seqs))
    idx_1 <- idx_2 <- idx_3 <- idx_4 <- TRUE
    if (start_atg) {
      idx_1 <- XVector::subseq(cds_seqs, 1, 3) == "ATG"
    }
    if (common_stop) {
      idx_2 <- as.character(
        Biostrings::reverseComplement(
          XVector::subseq(
            Biostrings::reverseComplement(cds_seqs), 1, 3
          )
        )
      ) %in% c("TAA", "TAG", "TGA")
    }
    if (no_frameshift) {
      idx_3 <- IRanges::width(cds_seqs) %% 3L == 0L
    }
    idx_4 <- IRanges::width(cds_seqs) > min_cds
    flt_cds_seqs <- cds_seqs[idx_1 + idx_2 + idx_3 + idx_4 == 4L]

    # choose transcript in a gene (priority longest cds, utr5, utr3)
    flt_tx_lens <- tx_lens[tx_lens$tx_name %in% names(flt_cds_seqs), ]
    tx_len <- gene_id <- tx_name <- cds_len <- utr5_len <- utr3_len <- NULL
    longest_tx <- flt_tx_lens |>
      dplyr::group_by(gene_id) |>
      dplyr::summarise(
        tx_name[
          base::which.max(
            10000L * base::rank(cds_len) +
              100L * base::rank(utr5_len) +
              base::rank(utr3_len)
          )
        ]
      )
    flt_tx_lens <- flt_tx_lens[
      flt_tx_lens$tx_name %in% longest_tx[, 2, drop = TRUE],
    ]

    flt_tx_seqs <- tx_seqs[flt_tx_lens$tx_name]

    # keep transcript with ATCG base
    idx <- rowSums(
      Biostrings::letterFrequency(
        flt_tx_seqs,
        letters = Biostrings::DNA_ALPHABET
      )[, -1:-4]
    ) == 0L
    flt_tx_seqs <- flt_tx_seqs[idx]
    mix_tx <- Biostrings::xscat(unlist(flt_tx_seqs))

    # add structure information
    tx_5pos <- c(
      0L,
      cumsum(
        IRanges::width(flt_tx_seqs)
      )[-length(flt_tx_seqs)]
    )
    tx_str <- data.frame(
      utr5_p5 = tx_5pos + 1L,
      utr5_p3 = tx_5pos + flt_tx_lens$utr5_len,
      utr3_p5 = tx_5pos + flt_tx_lens$utr5_len + flt_tx_lens$cds_len + 1L,
      utr3_p3 = tx_5pos +
        flt_tx_lens$utr5_len + flt_tx_lens$cds_len + flt_tx_lens$utr3_len
    )

    # get lncRNA information
    cds_gene <- unique(
      tx_lens$gene_id[
        tx_lens$tx_name %in% names(cds_seqs)
      ]
    )
    lnc_rna_gene <- setdiff(tx_lens$gene_id, cds_gene)
    # each lncRNA gene, choose the longest transcript
    lnc_tx_lens <- tx_lens[tx_lens$gene_id %in% lnc_rna_gene, ]

    longest_lncrna_tx <- lnc_tx_lens |>
      dplyr::group_by(gene_id) |>
      dplyr::summarise(tx_name[base::which.max(base::rank(tx_len))])

    lnc_tx_lens <- lnc_tx_lens[
      lnc_tx_lens$tx_name %in% longest_lncrna_tx[, 2, drop = TRUE],
    ]
    lnc_tx_seqs <- tx_seqs[lnc_tx_lens$tx_name]

    # keep transcript with ATCG base
    idx <- rowSums(
      Biostrings::letterFrequency(
        lnc_tx_seqs,
        letters = Biostrings::DNA_ALPHABET
      )[, -1:-4]
    ) == 0L
    lnc_tx_seqs <- lnc_tx_seqs[idx]
    lnc_tx_lens <- lnc_tx_lens[idx, ]
    mix_lnc_tx <- Biostrings::xscat(unlist(lnc_tx_seqs))

    tx_5pos <- c(
      0L,
      cumsum(
        IRanges::width(lnc_tx_seqs)
      )[-length(lnc_tx_seqs)]
    )
    lnc_tx_str <- data.frame(
      tx_p5 = tx_5pos + 1L,
      tx_p3 = tx_5pos + IRanges::width(lnc_tx_seqs)
    )

    tx_info <- list(
      flt_tx_seqs = flt_tx_seqs,
      tx_lens = flt_tx_lens,
      mix_tx = mix_tx,
      mix_tx_pos = tx_str
    )

    lncrna_info <- list(
      tx_seqs = lnc_tx_seqs,
      mix_tx = mix_lnc_tx,
      mix_tx_pos = lnc_tx_str
    )

    return(list(
      tx_info = tx_info,
      lncrna_info = lncrna_info
    ))
  }

  get_gene <- function(
      txdb,
      genome_fa,
      bs_genome) {
    gene_cds <- GenomicFeatures::cdsBy(txdb, by = "gene")
    tx_info <- GenomicFeatures::transcripts(txdb)
    tx_exon <- GenomicFeatures::exonsBy(
      txdb,
      by = "tx",
      use.names = TRUE
    )
    tx_cds_exon <- GenomicFeatures::cdsBy(
      txdb,
      by = "tx",
      use.names = TRUE
    )
    df_cds_id <- data.frame(
      gene_name = rep(
        names(gene_cds),
        IRanges::width(IRanges::PartitioningByEnd(gene_cds))
      ),
      cds_id = gene_cds@unlistData@elementMetadata@listData$cds_id
    )
    tx_strand <- as.factor(as.vector(BiocGenerics::strand(tx_info)))
    names(tx_strand) <- tx_info$tx_name
    genome_fa <- Biostrings::readDNAStringSet(genome_fa)
    names(genome_fa) <- stringr::str_split(
      names(genome_fa), " ",
      simplify = TRUE
    )[, 1]
    return(
      list(
        tx_exon_width = IRanges::width(tx_exon),
        tx_cds_width = IRanges::width(tx_cds_exon),
        tx_seqs = GenomicFeatures::extractTranscriptSeqs(bs_genome, tx_exon),
        tx_cds = GenomicFeatures::transcriptLengths(
          txdb,
          with.utr5_len = TRUE,
          with.cds_len = TRUE,
          with.utr3_len = TRUE
        ),
        gene_cds = gene_cds,
        tx_info = tx_info,
        tx_exon = tx_exon,
        tx_cds_exon = tx_cds_exon,
        df_cds_id = df_cds_id,
        tx_strand = tx_strand,
        genome_fa = genome_fa
      )
    )
  }

  get_cds_position <- function(gene_info) {
    cds_df <- gene_info$tx_cds
    tx_cds_exon <- gene_info$tx_cds_exon
    tx_ids_all <- names(tx_cds_exon)
    if (is.null(tx_ids_all) || any(!nzchar(tx_ids_all))) {
      tx_ids_all <- as.character(seq_along(tx_cds_exon))
      names(tx_cds_exon) <- tx_ids_all
    }

    cds_not_3multiple <- cds_df$tx_name[cds_df$cds_len %% 3L != 0L]
    idx_not3 <- tx_ids_all %in% cds_not_3multiple
    cds_exon_not3 <- tx_cds_exon[idx_not3]
    exon_num <- S4Vectors::elementNROWS(cds_exon_not3)
    gr_not3 <- unlist(cds_exon_not3, use.names = FALSE)

    df_not3 <- list(
      cds_exon = data.frame(
        tx_id = names(cds_exon_not3),
        exon_num = exon_num
      ),
      exon_position = data.frame(
        exon_left = BiocGenerics::start(gr_not3),
        exon_right = BiocGenerics::end(gr_not3)
      )
    )

    cds_exon_is3 <- tx_cds_exon[!idx_not3]
    if (length(gene_info$tx_cds_width) == length(tx_cds_exon)) {
      cds_exon_width <- gene_info$tx_cds_width[!idx_not3]
    } else {
      cds_exon_width <- gene_info$tx_cds_width[names(cds_exon_is3)]
    }
    exon_frame <- lapply(
      cds_exon_width,
      function(x) {
        cumsum(x) %% 3L
      }
    )

    exon_num <- S4Vectors::elementNROWS(cds_exon_is3)
    gr_is3 <- unlist(cds_exon_is3, use.names = FALSE)

    df_is3 <- list(
      cds_exon = data.frame(
        tx_id = names(cds_exon_is3),
        exon_num = exon_num
      ),
      exon_position = data.frame(
        exon_left = BiocGenerics::start(gr_is3),
        exon_right = BiocGenerics::end(gr_is3),
        exon_frame = unlist(
          exon_frame,
          use.names = FALSE
        )
      )
    )

    return(
      list(
        df_not3 = df_not3,
        df_is3 = df_is3
      )
    )
  }

  find_candid_tr <- function(
      min_aa,
      tx_info,
      is_ncrna = FALSE) {
    # Find all stop codons
    cord_stop <- do.call(
      c,
      lapply(
        c("TAG", "TGA", "TAA"),
        function(x) {
          Biostrings::matchPattern(x, tx_info$mix_tx)@ranges
        }
      )
    )
    if (is_ncrna) {
      connect_tx <- IRanges::IRanges(
        start = tx_info$mix_tx_pos$tx_p3,
        width = 2L
      )
    } else {
      connect_tx <- IRanges::IRanges(
        start = tx_info$mix_tx_pos$utr3_p3,
        width = 2L
      )
    }
    # Filter out coordinates that overlap with tx connections
    cord_stop <- cord_stop[
      -IRanges::findOverlaps(
        connect_tx,
        cord_stop,
        minoverlap = 2L
      )@to
    ]
    # Group all stop codons by frame
    if (is_ncrna) {
      utr5_3cord <- c(
        tx_info$mix_tx_pos$tx_p5,
        tx_info$mix_tx_pos$tx_p5 + 1L,
        tx_info$mix_tx_pos$tx_p5 + 2L
      )
    } else {
      utr5_3cord <- c(
        tx_info$mix_tx_pos$utr5_p5,
        tx_info$mix_tx_pos$utr5_p5 + 1L,
        tx_info$mix_tx_pos$utr5_p5 + 2L
      )
    }

    cord_stop_utr5 <- sort(c(utr5_3cord, cord_stop@start))
    idx_stop_utr5 <- cord_stop_utr5 %% 3L

    # Generate all candidate translation regions
    candi_translat <- lapply(
      0:2,
      function(frame_n) {
        # candidate translation regions
        cord_tmp <- cord_stop_utr5[idx_stop_utr5 == frame_n]
        candi_trans <- IRanges::IRanges(
          start = cord_tmp[-length(cord_tmp)],
          end = cord_tmp[-1]
        )
        candi_trans <- candi_trans[
          !(BiocGenerics::end(candi_trans) %in% utr5_3cord)
        ]
        # min aa length
        candi_trans <- candi_trans[
          IRanges::width(candi_trans) > (min_aa - 1L) * 3L
        ]
        return(candi_trans)
      }
    )

    return(candi_translat)
  }

  # Find all candidate ORFs
  find_candid_orf <- function(
      start_codon,
      min_aa,
      tx_info) {
    # Find all start and stop codons
    cord_start <- do.call(
      c,
      lapply(
        start_codon,
        function(x) {
          Biostrings::matchPattern(
            x,
            tx_info$mix_tx
          )@ranges
        }
      )
    )
    cord_stop <- do.call(
      c,
      lapply(
        c("TAG", "TGA", "TAA"),
        function(x) {
          Biostrings::matchPattern(
            x,
            tx_info$mix_tx
          )@ranges
        }
      )
    )
    connect_tx <- IRanges::IRanges(
      start = tx_info$mix_tx_pos$utr3_p3,
      width = 2L
    )
    # Filter out coordinates that overlap with tx connections
    cord_start <- cord_start[
      -IRanges::findOverlaps(
        connect_tx,
        cord_start,
        minoverlap = 2L
      )@to
    ]
    idx_start <- cord_start@start %% 3L
    cord_stop <- cord_stop[
      -IRanges::findOverlaps(
        connect_tx,
        cord_stop,
        minoverlap = 2L
      )@to
    ]
    # Group all stop codons by frame
    utr5_3cord <- c(
      tx_info$mix_tx_pos$utr5_p5,
      tx_info$mix_tx_pos$utr5_p5 + 1L,
      tx_info$mix_tx_pos$utr5_p5 + 2L
    )
    cord_stop_utr5 <- sort(c(utr5_3cord, cord_stop@start))
    idx_stop_utr5 <- cord_stop_utr5 %% 3L
    # Generate all candidate translation regions
    candi_translat <- lapply(
      0:2,
      function(frame_n) {
        # candidate translation regions
        cord_tmp <- cord_stop_utr5[idx_stop_utr5 == frame_n]
        candi_trans <- IRanges::IRanges(
          start = cord_tmp[-length(cord_tmp)],
          end = cord_tmp[-1]
        )
        candi_trans <- candi_trans[
          !(BiocGenerics::end(candi_trans) %in% utr5_3cord)
        ]
        # min aa length
        candi_trans <- candi_trans[
          IRanges::width(candi_trans) > (min_aa - 1L) * 3L
        ]
        # candidate start codon
        cord_tmp_s <- cord_start@start[idx_start == frame_n]
        hit_start <- IRanges::findOverlaps(
          candi_trans,
          IRanges::IRanges(start = cord_tmp_s)
        )
        candi_orf <- IRanges::IRanges(
          start = cord_tmp_s[hit_start@to],
          end = BiocGenerics::end(candi_trans)[hit_start@from]
        )
        # min aa length
        candi_orf <- candi_orf[
          IRanges::width(candi_orf) > (min_aa - 1L) * 3L
        ]
        # Candidate translation regions without initiation codons are removed
        hit_orf <- IRanges::findOverlaps(
          candi_trans,
          candi_orf,
          minoverlap = 3L
        )
        candi_trans <- candi_trans[unique(hit_orf@from)]
        # Find longest ORF that shares stop codon
        candi_longest <- IRanges::reduce(candi_orf)
        return(
          list(
            trans = candi_trans,
            all_orf = candi_orf,
            longset = candi_longest
          )
        )
      }
    )
    # Generate the relationship between the annotated ORF and the candidate ORF
    candi_longest <- do.call(
      c,
      lapply(
        candi_translat,
        function(x) {
          x[[3]]
        }
      )
    )
    # upstream ORF, upstream overlapping orf, n-terminal extension
    idx_utr5 <- tx_info$tx_lens$utr5_len > 10
    rg_utr5 <- IRanges::IRanges(
      start = tx_info$mix_tx_pos$utr5_p5[idx_utr5],
      end = tx_info$mix_tx_pos$utr5_p3[idx_utr5] - 2L
    )
    up_orf <- candi_longest[
      IRanges::findOverlaps(
        candi_longest,
        rg_utr5,
        type = "within"
      )@from
    ]
    # downstream ORF
    idx_utr3 <- tx_info$tx_lens$utr3_len > 10
    rg_utr3 <- IRanges::IRanges(
      start = tx_info$mix_tx_pos$utr3_p5[idx_utr3],
      end = tx_info$mix_tx_pos$utr3_p3[idx_utr3] - 2L
    )
    down_orf <- candi_longest[
      IRanges::findOverlaps(
        candi_longest,
        rg_utr3,
        type = "within"
      )@from
    ]
    # out-of-frame ORF
    rg_cds <- IRanges::IRanges(
      start = tx_info$mix_tx_pos$utr5_p3 + 1L,
      end = tx_info$mix_tx_pos$utr3_p5 - 3L
    )
    outframe_orf <- candi_longest[
      IRanges::findOverlaps(
        candi_longest,
        rg_cds,
        type = "within"
      )@from
    ]
    idx_equal_cds <- unique(
      IRanges::findOverlaps(outframe_orf, rg_cds, type = "equal")@from
    )
    if (length(idx_equal_cds) > 0) {
      outframe_orf <- outframe_orf[-idx_equal_cds]
    }
    # upstream and downstream overlapping ORF, n-terminal extension
    idx_overlap <- IRanges::findOverlaps(
      candi_longest,
      rg_cds
    )@from
    idx_overlap <- idx_overlap[
      !(idx_overlap %in%
        IRanges::findOverlaps(
          candi_longest,
          rg_cds,
          type = "within"
        )@from)
    ]
    rg_stop <- IRanges::IRanges(
      start = tx_info$mix_tx_pos$utr3_p5 - 3L,
      end = tx_info$mix_tx_pos$utr3_p5 - 1L
    )
    rg_start <- IRanges::IRanges(
      start = tx_info$mix_tx_pos$utr5_p3 + 1L
    )
    up_overlap_orf <- candi_longest[idx_overlap][
      IRanges::findOverlaps(
        rg_start,
        candi_longest[idx_overlap],
        type = "within"
      )@to
    ]
    end_orf <- candi_longest[
      IRanges::findOverlaps(
        candi_longest,
        rg_cds,
        type = "end"
      )@from
    ]
    tmp_idx <- IRanges::overlapsAny(up_overlap_orf, end_orf, type = "equal")
    n_extend <- up_overlap_orf[tmp_idx]
    up_overlap_orf <- up_overlap_orf[!tmp_idx]
    down_overlap_orf <- candi_longest[idx_overlap][
      IRanges::findOverlaps(
        rg_stop,
        candi_longest[idx_overlap],
        type = "within"
      )@to
    ]
    orf_relation <- list(
      up_orf = up_orf,
      down_orf = down_orf,
      outframe_orf = outframe_orf,
      n_extend = n_extend,
      up_overlap_orf = up_overlap_orf,
      down_overlap_orf = down_overlap_orf
    )
    return(
      list(
        candi_translat = candi_translat,
        orf_relation = orf_relation
      )
    )
  }

  dir.create(output_dir, recursive = TRUE)
  message("Start build_annotation_index()")
  message(
    paste0(
      "Create BSgenome package for ", scientific_name, " ..."
    )
  )
  create_bs(
    out_dir = output_dir,
    seq_file = genome_file,
    scientific_name = scientific_name,
    anno_name = annotation_version
  )

  message(
    paste0(
      "Create TxDb object for ", scientific_name, " ..."
    )
  )
  df_seqinfo <- GenomeInfoDb::seqinfo(
    rtracklayer::TwoBitFile(
      sub(
        x = genome_file,
        pattern = "\\.fa$|\\.fasta", "\\.2bit"
      )
    )
  )
  txdb <- txdbmaker::makeTxDbFromGFF(
    file = gff_file,
    format = "auto",
    chrominfo = df_seqinfo
  )

  message("Get the transcript containing the longest CDS ...")
  bs_name <- paste0(
    "BSgenome.",
    paste0(
      substr(scientific_name, 1, 1),
      rev(stringr::str_split(scientific_name, pattern = "\\.")[[1]])[1]
    ),
    ".",
    annotation_version,
    ".Ribo"
  )
  require(bs_name,
    character.only = TRUE
  )
  tx_info_lst <- get_tx_info(
    bs = get(bs_name),
    txdb = txdb
  )

  message("Find ORFs on the transcript with the annotated CDS ...")
  candidate_translat_reg <- find_candid_tr(
    min_aa = orf_min_aa,
    tx_info = tx_info_lst$tx_info,
    is_ncrna = FALSE
  )

  candidate_orf <- find_candid_orf(
    start_codon = orf_start_code,
    min_aa = orf_min_aa,
    tx_info = tx_info_lst$tx_info
  )

  message("Find ORFs on ncRNA ...")
  candidate_lnc_reg <- find_candid_tr(
    min_aa = orf_min_aa,
    tx_info = tx_info_lst$lncrna_info,
    is_ncrna = TRUE
  )

  message("Prepare exon position ...")
  gene_info <- get_gene(
    txdb = txdb,
    genome_fa = genome_file,
    bs_genome = get(bs_name)
  )

  cds_exon <- get_cds_position(gene_info = gene_info)

  message("Save R Objects ...")
  dir.create(paste0(output_dir, "/transcript"), recursive = TRUE)
  save(tx_info_lst,
    candidate_translat_reg,
    candidate_lnc_reg,
    candidate_orf,
    gene_info,
    cds_exon,
    file = paste0(output_dir, "/transcript/tx_information.RData")
  )

  Biostrings::writeXStringSet(
    x = tx_info_lst$tx_info$flt_tx_seqs,
    filepath = paste0(output_dir, "/transcript/longest_txs.fa"),
    format = "fasta"
  )

  Biostrings::writeXStringSet(
    x = tx_info_lst$lncrna_info$tx_seqs,
    filepath = paste0(output_dir, "/transcript/lncRNA_txs.fa"),
    format = "fasta"
  )
  message("build_annotation_index() done. Outputs: ", output_dir, "/transcript")
}

#' @export
annotate_prepare <- function(...) {
  .Deprecated("build_annotation_index")
  build_annotation_index(...)
}
