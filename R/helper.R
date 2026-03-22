#' Infer ribosome footprint parameters
#'
#' Estimates ligation and cleavage parameters from mapped ribosome profiling
#' reads. The function iteratively optimises cleavage probabilities and ligation
#' biases before returning a cleaned footprint table together with the learned
#' parameters.
#'
#' @param tx_info Transcript annotation list generated during preprocessing.
#' @param lst_bam List of BAM-derived read tables as produced by [read_bam_alignments()].
#' @param rnase_bias Logical flag indicating whether RNase digestion bias should
#'   be modelled.
#' @param ligat_par Logical flag controlling whether ligation parameters are
#'   re-estimated.
#' @param init_num Number of random initial points used in iterative parameter
#'   optimization.
#' @param maxfeval Maximum function evaluations for each `dfoptim::nmkb` call in
#'   iterative parameter optimization.
#' @param plot_lst_rpf_dist Logical. If `TRUE`, save one bar plot of filtered
#'   read-length distribution used for parameter estimation from
#'   `lst_rpf$rpf_info`. Defaults to `FALSE`.
#' @param plot_file Optional output base path for the filtered read-length plot.
#'   If `NULL`, a timestamped PDF is written to `tempdir()`. If provided, output
#'   file uses `_estimation.pdf` suffix.
#'
#' @return A list with two components: `lst_rpf`, the processed footprint
#'   overlaps, and `par_0`, the learnt parameter set.
#'
#' @export
learn_bias_parameters <- function(tx_info,
                      lst_bam,
                      rnase_bias,
                      RNase,
                      ligat_par,
                      init_num = 20,
                      maxfeval = 5000,
                      plot_lst_rpf_dist = FALSE,
                      plot_file = NULL) {
                        #browser()
  correct_kmer <- function(rpf,
                           tx_info,
                           par_0) {
    correct_rpf <- function(kmer5,
                            kmer3,
                            rpf_num,
                            par_0) {
      corrected_rpf <- exp(
        (rpf_num -
          par_0$eff_intercept -
          par_0$eff_f5[kmer5] -
          par_0$eff_f3[kmer3]
        )
      ) - 1
      return(as.vector(corrected_rpf))
    }

    # 3-mer at the terminal of reads
    rpf_seq <- Biostrings::DNAStringSet(
      x = tx_info$mix_tx
    )[rep(1L, length(rpf$pos))]

    kmer_5 <- as.character(XVector::subseq(
      x = rpf_seq,
      start = rpf$pos,
      width = 3L
    ))

    kmer_3 <- as.character(XVector::subseq(
      x = rpf_seq,
      end = rpf$pos + rpf$qwidth - 1L,
      width = 3L
    ))

    corrected_rpf <- correct_rpf(
      kmer5 = kmer_5,
      kmer3 = kmer_3,
      rpf_num = log(rpf$weight + 1),
      par_0 = par_0
    )

    corrected_rpf[which(corrected_rpf < 0)] <- 0

    if (par_0$eff_intercept != 0) {
      corrected_rpf <- corrected_rpf / (sum(corrected_rpf) / sum(rpf$weight))
    }

    return(corrected_rpf)
  }

  glm2par <- function(fit_coefs,
                      par_0) {
    par_0$eff_intercept[1] <- fit_coefs[1]

    par_0$eff_f5[
      stringr::str_split(
        string = names(fit_coefs)[
          grep(
            pattern = "kmer5",
            x = names(fit_coefs)
          )
        ],
        pattern = "5",
        simplify = TRUE
      )[, 2]
    ] <- fit_coefs[2:64]

    par_0$eff_f3[
      stringr::str_split(
        string = names(fit_coefs)[
          grep(
            pattern = "kmer3",
            x = names(fit_coefs)
          )
        ],
        pattern = "3",
        simplify = TRUE
      )[, 2]
    ] <- fit_coefs[65:127]

    return(
      list(
        par_0 = par_0
      )
    )
  }

  iter_rpf_overlap <- function(potential_rpf_lig2,
                               rpf_overlap) {
    rpf_overlap_iter <- potential_rpf_lig2[potential_rpf_lig2$weight > 0, ]
    rpf_overlap_iter$tag <- rpf_overlap_iter$pos * 100 + rpf_overlap_iter$qwidth
    iter_idx <- match(rpf_overlap_iter$tag, rpf_overlap$tag)
    rpf_overlap_iter$weight <- rpf_overlap$weight[iter_idx]
    return(list(
      iter_idx = iter_idx,
      rpf_overlap_iter = rpf_overlap_iter
    ))
  }

  maintain_prob <- function(cut_seq,
                            prod_hd,
                            bias) {
    prob_mc <- stringr::str_split(
      string = as.vector(cut_seq),
      pattern = "",
      simplify = TRUE
    )

    prob_mp <- matrix(
      data = bias[prob_mc],
      nrow = nrow(prob_mc)
    )

    prob_mp <- prob_mp^(rep(1L, nrow(prob_mp)) %o% prod_hd)

    return(prob_mp)
  }

  get_cleavage_prob <- function(seqs,
                                bias,
                                prob_hd5,
                                prob_hd3) {
    # for 5'
    maintain_prob5 <- maintain_prob(
      cut_seq = seqs$up_seq,
      prod_hd = prob_hd5,
      bias = bias
    )

    cle_p5rev <- 1 - maintain_prob5[, length(prob_hd5):1]

    maintain_cumprod5rev <- matrix(
      data = 1,
      nrow = nrow(maintain_prob5),
      ncol = length(prob_hd5)
    )

    maintain_cumprod5rev[, -1] <- matrixStats::rowCumprods(
      maintain_prob5[, length(prob_hd5):2]
    )

    final_p5 <- maintain_cumprod5rev * cle_p5rev

    final_p5 <- final_p5[, ncol(final_p5):1]

    final_p5 <- final_p5 / Matrix::rowSums(final_p5)

    # for 3'
    maintain_prob3 <- maintain_prob(
      cut_seq = seqs$dn_seq,
      prod_hd = prob_hd3,
      bias = bias
    )

    cle_p3 <- 1 - maintain_prob3

    maintain_cumprod3 <- matrix(
      data = 1,
      nrow = nrow(maintain_prob3),
      ncol = length(prob_hd3)
    )

    maintain_cumprod3[, -1] <- matrixStats::rowCumprods(
      maintain_prob3[, -length(prob_hd3)]
    )

    final_p3 <- maintain_cumprod3 * cle_p3

    final_p3 <- final_p3 / Matrix::rowSums(final_p3)

    return(
      list(
        final_p5 = final_p5,
        final_p3 = final_p3
      )
    )
  }

  get_cleavage_seq <- function(seqs,
                               p_site,
                               ribo_size) {
    seqs <- Biostrings::DNAStringSet(x = seqs)[
      rep(1L, length(p_site))
    ]

    up_seq <- Biostrings::subseq(
      x = seqs,
      start = p_site - sum(ribo_size[1:2]),
      width = ribo_size[1]
    )

    dn_seq <- Biostrings::subseq(
      x = seqs,
      start = p_site + ribo_size[3],
      width = ribo_size[4]
    )

    return(
      list(
        up_seq = up_seq,
        dn_seq = dn_seq
      )
    )
  }

  em_ribo_num <- function(df_rpf,
                          candi_psite,
                          candi_cut5,
                          candi_cut3,
                          tx_info,
                          ribo_size,
                          par_0,
                          pos_psite,
                          iter_num) {
    # filter out read weight 0
    idx_read <- df_rpf$weight > 0L

    df_rpf <- df_rpf[idx_read, ]

    candi_p_weight <- candi_psite <- candi_psite[idx_read, ]

    candi_cut5 <- candi_cut5[idx_read, ]

    candi_cut3 <- candi_cut3[idx_read, ]

    # shrink p-sites
    shrink_p <- as.factor(candi_psite@x)

    shrink_pos <- as.integer(shrink_p)

    index_j <- rank(
      x = shrink_pos,
      ties.method = "first"
    ) -
      rank(
        x = shrink_pos,
        ties.method = "min"
      ) +
      1L

    candi_pos <- as.integer(levels(shrink_p))

    vec_pnum <- rep(1, length(candi_pos))

    # Split read weight
    uni_candi_psite <- unique(candi_psite@x)

    cut_seq <- get_cleavage_seq(
      seqs = tx_info$mix_tx,
      p_site = uni_candi_psite,
      ribo_size = ribo_size
    )

    idx_ij <- match(candi_psite@x, uni_candi_psite)

    # Generating cleavage probabilities for each ribosome terminus
    cut_prob <- get_cleavage_prob(
      seqs = cut_seq,
      bias = par_0$cut_bias$s7,
      prob_hd5 = par_0$prob_hd$p5,
      prob_hd3 = par_0$prob_hd$p3
    )

    base_prob <- cut_prob$final_p5[
      nrow(cut_prob$final_p5) *
        (candi_cut5@x - 1L) +
        idx_ij
    ] *
      cut_prob$final_p3[
        nrow(cut_prob$final_p3) *
          (candi_cut3@x - 1L) +
          idx_ij
      ]

    for (i in 1:iter_num) {
      candi_p_weight@x <- base_prob * vec_pnum[shrink_pos]

      rse_iter <- Matrix::sparseMatrix(
        i = shrink_pos,
        j = index_j,
        x = (df_rpf$weight * candi_p_weight / Matrix::rowSums(candi_p_weight))@x
      )

      vec_pnum[] <- Matrix::rowSums(rse_iter)
    }

    idx_psite <- match(pos_psite, candi_pos)

    vec_pnum <- vec_pnum[idx_psite]

    vec_pnum[is.na(vec_pnum)] <- 0

    return(vec_pnum)
  }

  adjust_ribo <- function(ribo_size,
                          prob_hd,
                          limit1 = 0.05,
                          limit2 = 0.9,
                          limit3 = 0.9) {
    ribo_0 <- ribo_size
    prob_hd <- prob_hd
    p5len <- ribo_size[1]
    p3len <- ribo_size[4]
    chg <- 0
    if (chg == 0) {
      if (prob_hd$p3[1] < limit1) {
        if (prob_hd$p3[2] < limit1) {
          ribo_size[4] <- ribo_size[4] - 1L
          ribo_size[3] <- ribo_size[3] + 1L
          chg <- 1
        }
      } else {
        if (min(prob_hd$p3) == prob_hd$p3[1]) {
          ribo_size[4] <- ribo_size[4] + 1L
          ribo_size[3] <- ribo_size[3] - 1L
          chg <- 1
        }
      }
    }
    if (chg == 0) {
      if (prob_hd$p5[p5len] < limit1) {
        if (prob_hd$p5[p5len - 1] < limit1) {
          ribo_size[1] <- ribo_size[1] - 1L
          ribo_size[2] <- ribo_size[2] + 1L
          chg <- 1
        }
      } else {
        if (min(prob_hd$p5) == prob_hd$p5[p5len]) {
          ribo_size[1] <- ribo_size[1] + 1L
          ribo_size[2] <- ribo_size[2] - 1L
          chg <- 1
        }
      }
    }
    if (chg == 0) {
      if (prob_hd$p3[p3len] > limit2) {
        if (prob_hd$p3[p3len - 1] > limit3) {
          ribo_size[4] <- ribo_size[4] - 1L
          chg <- 1
        }
      } else {
        if (max(prob_hd$p3) != prob_hd$p3[p3len]) {
          ribo_size[4] <- ribo_size[4] - 1L
          chg <- 1
        } else {
          ribo_size[4] <- ribo_size[4] + 1L
          chg <- 1
        }
      }
    }

    if (chg == 0) {
      if (prob_hd$p5[1] > limit2) {
        if (prob_hd$p5[2] > limit3) {
          ribo_size[1] <- ribo_size[1] - 1L
          chg <- 1
        }
      } else {
        if (max(prob_hd$p5) != prob_hd$p5[1]) {
          ribo_size[1] <- ribo_size[1] - 1L
          chg <- 1
        } else {
          ribo_size[1] <- ribo_size[1] + 1L
          chg <- 1
        }
      }
    }

    if (identical(ribo_0, ribo_size)) {
      chg_ribo <- FALSE
    } else {
      (
        chg_ribo <- TRUE
      )
    }
    return(list(ribo_size = ribo_size, chg_ribo = chg_ribo))
  }

  assign_par <- function(out_par,
                         par_0,
                         reference_base_value) {
    bias_len <- length(par_0$cut_bias$s7)

    hindrance_len <- length(par_0$prob_hd$p5)

    out_par <- c(
      out_par[1:(bias_len - 1)],
      reference_base_value,
      out_par[-(1:(bias_len - 1))]
    )

    par_0$cut_bias$s7[
      1:bias_len
    ] <- out_par[1:bias_len]

    par_0$prob_hd$p5 <- out_par[
      (bias_len + 1):(bias_len + hindrance_len)
    ]

    par_0$prob_hd$p3 <- out_par[
      (bias_len + hindrance_len + 1):length(out_par)
    ]

    return(par_0)
  }

  maintain_prob_faster <- function(prob_mc,
                                   base_to_pos,
                                   prod_hd,
                                   bias) {
    prob_mp <- matrix(
      data = bias[prob_mc],
      nrow = nrow(prob_mc)
    )

    prob_mp <- prob_mp^(rep(1L, nrow(prob_mp)) %o% prod_hd)

    return(prob_mp)
  }

  get_cleavage_prob_faster <- function(prob_mc_up,
                                       prob_mc_dn,
                                       base_to_pos,
                                       bias,
                                       prob_hd5,
                                       prob_hd3) {
    # for 5'
    maintain_prob5 <- maintain_prob_faster(
      prob_mc = prob_mc_up,
      base_to_pos = base_to_pos,
      prod_hd = prob_hd5,
      bias = bias
    )

    cle_p5rev <- 1 - maintain_prob5[, length(prob_hd5):1]

    maintain_cumprod5rev <- matrix(
      data = 1,
      nrow = nrow(maintain_prob5),
      ncol = length(prob_hd5)
    )

    maintain_cumprod5rev[, -1] <- matrixStats::rowCumprods(
      x = maintain_prob5[, length(prob_hd5):2]
    )

    final_p5 <- maintain_cumprod5rev * cle_p5rev

    final_p5 <- final_p5[, length(prob_hd5):1]

    final_p5 <- final_p5 / matrixStats::rowSums2(final_p5)

    # for 3'
    maintain_prob3 <- maintain_prob_faster(
      prob_mc = prob_mc_dn,
      base_to_pos = base_to_pos,
      prod_hd = prob_hd3,
      bias = bias
    )

    cle_p3 <- 1 - maintain_prob3

    maintain_cumprod3 <- matrix(
      data = 1,
      nrow = nrow(maintain_prob3),
      ncol = length(prob_hd3)
    )

    maintain_cumprod3[, -1] <- matrixStats::rowCumprods(
      x = maintain_prob3[, -length(prob_hd3)]
    )

    final_p3 <- maintain_cumprod3 * cle_p3

    final_p3 <- final_p3 / matrixStats::rowSums2(final_p3)

    return(
      list(
        final_p5 = final_p5,
        final_p3 = final_p3
      )
    )
  }

  ligation_par <- function(p_num,
                           rpf,
                           candi_psite,
                           candi_cut5,
                           candi_cut3,
                           candi_weight,
                           tx_info,
                           kmer_num,
                           par_0) {
    # 3-mer at the terminal of reads
    rpf_seq <- Biostrings::DNAStringSet(
      x = tx_info$mix_tx
    )[rep(1L, length(rpf$pos))]

    kmer_5 <- as.character(XVector::subseq(
      x = rpf_seq,
      start = rpf$pos,
      width = 3L
    ))

    kmer_3 <- as.character(XVector::subseq(
      x = rpf_seq,
      end = rpf$pos + rpf$qwidth - 1L,
      width = 3L
    ))

    # choose RPFs, sort expected RPF numbers first
    vec_pos <- sort((c(
      sapply(
        names(par_0$eff_f5), function(kmer_i) {
          which(kmer_5 == kmer_i)[1:kmer_num]
        }
      ),
      sapply(
        names(par_0$eff_f3), function(kmer_i) {
          which(kmer_3 == kmer_i)[1:kmer_num]
        }
      )
    )))

    # prepare to calculate expected RPF numbers
    candi_psite <- candi_psite[vec_pos, ]@x

    candi_weight <- candi_weight[vec_pos, ]

    candi_cut5 <- candi_cut5[vec_pos, ]@x

    candi_cut3 <- candi_cut3[vec_pos, ]@x

    uni_candi_psite <- unique(candi_psite)

    idx_ij <- match(candi_psite, uni_candi_psite)

    cut_seq <- get_cleavage_seq(
      seqs = tx_info$mix_tx,
      p_site = uni_candi_psite,
      ribo_size = par_0$ribo_size
    )

    prob_mc_up <- stringr::str_split(
      string = as.character(cut_seq$up_seq),
      pattern = "",
      simplify = TRUE
    )

    prob_mc_dn <- stringr::str_split(
      string = as.character(cut_seq$dn_seq),
      pattern = "",
      simplify = TRUE
    )

    cut_prob <- get_cleavage_prob(
      seqs = cut_seq,
      bias = par_0$cut_bias$s7,
      prob_hd5 = par_0$prob_hd$p5,
      prob_hd3 = par_0$prob_hd$p3
    )

    vec_5 <- nrow(cut_prob$final_p5) * (candi_cut5 - 1L) + idx_ij

    vec_3 <- nrow(cut_prob$final_p3) * (candi_cut3 - 1L) + idx_ij

    obj_num <- rpf$weight[vec_pos]

    vec_weight <- p_num[candi_psite]

    base_to_pos <- 1:4

    bias <- par_0$cut_bias$s7

    hd5_len <- 5:(length(par_0$prob_hd$p5) + 4)

    hd3_len <- seq.int(
      from = length(par_0$prob_hd$p5) + 5,
      to = length(par_0$prob_hd$p5) + length(par_0$prob_hd$p3) + 4
    )

    # calculate expected RPF numbers
    rss_rpf <- function(x) {
      cut_prob <- get_cleavage_prob_faster(
        prob_mc_up = prob_mc_up,
        prob_mc_dn = prob_mc_dn,
        base_to_pos = base_to_pos,
        bias = bias,
        prob_hd5 = x[hd5_len],
        prob_hd3 = x[hd3_len]
      )

      candi_weight@x <- cut_prob$final_p5[vec_5] *
        cut_prob$final_p3[vec_3] *
        vec_weight

      return(Matrix::rowSums(candi_weight))
    }

    exp_rpf_num <- rss_rpf(
      x = c(par_0$cut_bias$s7, par_0$prob_hd$p5, par_0$prob_hd$p3)
    )

    # prepare input matrix of model
    mat_train <- data.frame(
      kmer5 = kmer_5[vec_pos],
      kmer3 = kmer_3[vec_pos],
      exp_num = log(exp_rpf_num + 1),
      obj_num = log(obj_num + 1),
      stringsAsFactors = TRUE
    )

    # fit generalized linear model
    glm_fit <- stats::glm(
      formula = obj_num ~ offset(exp_num) + kmer5 + kmer3,
      data = mat_train,
      family = gaussian(link = "identity")
    )

    return(
      list(
        glm_coef = glm_fit,
        mat_train = mat_train,
        kmer_5 = kmer_5,
        kmer_3 = kmer_3
      )
    )
  }

  iter_estimate <- function(par_0,
                            candi_psite,
                            candi_weight,
                            prob_mc_up,
                            prob_mc_dn,
                            vec_5,
                            vec_3,
                            rpf_group_5,
                            rpf_group_3,
                            obj_num,
                            ribo_num,
                            kmer_cof,
                            par0_unknown,
                            init_num,
                            maxfeval,
                            init_tol,
                            sec_tol,
                            sec_diff,
                            iter_times) {
    obj_num <- obj_num + 1
    safe_nmkb <- function(par, fn, lower, upper, control) {
      res <- try(
        dfoptim::nmkb(
          par = par,
          fn = fn,
          lower = lower,
          upper = upper,
          control = control
        ),
        silent = TRUE
      )
      if (inherits(res, "try-error")) {
        return(NULL)
      }
      if (is.null(res$par) || is.null(res$value) || any(!is.finite(res$par)) || !is.finite(res$value)) {
        return(NULL)
      }
      res
    }

    rss_bias <- function(x) {
      bias[bias_len] <- x[bias_len]

      cut_prob <- get_cleavage_prob_faster(
        prob_mc_up = prob_mc_up,
        prob_mc_dn = prob_mc_dn,
        base_to_pos = base_to_pos,
        bias = bias,
        prob_hd5 = x[hd5_len],
        prob_hd3 = x[hd3_len]
      )

      candi_weight@x <- cut_prob$final_p5[vec_5] *
        cut_prob$final_p3[vec_3] *
        vec_weight

      exp_num <- (Matrix::rowSums(candi_weight) + 1) * kmer_cof

      return(
        sum(
          sqrt(
            (exp_num - obj_num)^2
          ) / (exp_num + exp_num^2 / 10)
        )
      )
    }

    rss_hd <- function(x) {
      bias[bias_len] <- par_0$base_prob

      cut_prob <- get_cleavage_prob_faster(
        prob_mc_up = prob_mc_up,
        prob_mc_dn = prob_mc_dn,
        base_to_pos = base_to_pos,
        bias = bias,
        prob_hd5 = x[hd5_len - bias_lenght + 1],
        prob_hd3 = x[hd3_len - bias_lenght + 1]
      )

      candi_weight@x <- cut_prob$final_p5[vec_5] *
        cut_prob$final_p3[vec_3] *
        vec_weight

      exp_num <- (Matrix::rowSums(candi_weight) + 1) * kmer_cof

      return(
        sum(
          sqrt(
            (exp_num - obj_num)^2 / (exp_num + exp_num^2 / 10)
          )
        )
      )
    }

    bias_lenght <- length(par_0$cut_bias$s7)
    vec_weight <- ribo_num[candi_psite@x]

    bias_len <- 1:(bias_lenght - 1L)

    base_to_pos <- 1:bias_lenght

    bias <- par_0$cut_bias$s7

    bias[bias_lenght] <- par_0$base_prob

    names(base_to_pos) <- names(par_0$cut_bias$s7)

    hd5_len <- bias_lenght:(length(par_0$prob_hd$p5) + bias_lenght - 1L)

    hd3_len <- seq.int(
      from = length(par_0$prob_hd$p5) + bias_lenght,
      to = length(par_0$prob_hd$p5) +
        length(par_0$prob_hd$p3) +
        bias_lenght - 1L
    )

    if (length(unique(par_0$cut_bias$s7)) == 1) {
      lst_pars <- list()
      for (round_i in 1:1) {
        if (par0_unknown) {
          lst_par <- lapply(
            1:init_num, function(i) {
              par0 <- c(
                sort(runif(length(hd5_len), 0.01, 0.9), decreasing = TRUE),
                sort(runif(length(hd3_len), 0.01, 0.9))
              )
              res <- safe_nmkb(
                par = par0,
                fn = rss_hd,
                lower = c(rep(1e-8, length(c(hd3_len, hd5_len)))),
                upper = c(rep(1 - 1e-6, length(c(hd3_len, hd5_len)))),
                control = list(
                  tol = init_tol, # decide search range
                  maxfeval = maxfeval
                )
              )
              if (is.null(res)) NULL else res
            }
          )
          lst_par = lst_par[!vapply(lst_par, is.null, logical(1))]
          if (length(lst_par) > 0L) {
            par_mx = sapply(lst_par, function(opt_out){
              c(opt_out$value, opt_out$par)
            })
            par0_raw <- par_mx[-1, which.min(par_mx[1, ])]
          } else {
            par0_raw <- c(par_0$prob_hd$p5, par_0$prob_hd$p3)
          }
        } else {
          par0_raw <- c(par_0$prob_hd$p5, par_0$prob_hd$p3)
        }
        par0 <- par0_raw

        for (i in seq_along(sec_tol)) {
          diff_value <- 1e4
          value0 <- 0
          round0 <- 1
          par0_raw[which(par0_raw >= (1 - 1e-6))] <- runif(
            length(which(par0_raw >= (1 - 1e-6))),
            (1 - 1e-6 - 2e-7),
            (1 - 1e-6 - 1e-7)
          )
          par0_raw[which(par0_raw <= (1e-8))] <- runif(
            length(which(par0_raw <= (1e-8))), 1e-8 + 1e-9, 1e-8 + 2e-9
          )
          opt_out <- safe_nmkb(
            par = par0_raw,
            fn = rss_hd,
            lower = c(rep(1e-8, length(c(hd3_len, hd5_len)))),
            upper = c(rep(1 - 1e-6, length(c(hd3_len, hd5_len)))),
            control = list(
              tol = sec_tol[i], # decide search range
              maxfeval = maxfeval
            )
          )
          if (!is.null(opt_out)) {
            value0 <- opt_out$value
            par0 <- opt_out$par
          } else {
            next
          }

          while ((abs(diff_value) > sec_diff[i]) & (round0 < iter_times)) {
            diff_value <- value0 - opt_out$value
            value0 <- opt_out$value
            par0 <- opt_out$par
            par0[which(par0 > (1 - 1e-6))] <- runif(
              length(which(par0 > (1 - 1e-6))),
              (1 - 1e-6 - 2e-7),
              (1 - 1e-6 - 1e-7)
            )
            par0[which(par0 < (1e-8))] <- runif(
              length(which(par0 < (1e-8))), 1e-8 + 1e-9, 1e-8 + 2e-9
            )
            round0 <- round0 + 1
            lst_opt <- safe_nmkb(
              par = par0,
              fn = rss_hd,
              lower = c(rep(1e-8, length(c(hd3_len, hd5_len)))),
              upper = c(rep(1 - 1e-6, length(c(hd3_len, hd5_len)))),
              control = list(
                tol = sec_tol[i], # decide search range
                maxfeval = maxfeval
              )
            )
            if (!is.null(lst_opt)) {
              opt_out <- lst_opt
              diff_value <- value0 - opt_out$value
              value0 <- opt_out$value
              par0 <- opt_out$par
            }
          }
          par0_raw <- par0
        }
        lst_pars[[round_i]] <- par0
      }
      lst_pars <- do.call(rbind, lst_pars)
      par0 <- apply(lst_pars, 2, median)
      par0 <- c(rep(par_0$base_prob, bias_lenght - 1L), par0)
    } else {
      if (par0_unknown) {
        lst_par <- lapply(
          1:init_num, function(i) {
            par0 <- c(
              (par_0$cut_bias$s7 +
                par_0$cut_bias$s7 *
                  runif(length(par_0$cut_bias$s7), -0.5, 0.5))[
                -length(par_0$cut_bias$s7)
              ],
              sort(runif(length(hd5_len), 0.01, 0.9), decreasing = TRUE),
              sort(runif(length(hd3_len), 0.01, 0.9))
            )
            res <- safe_nmkb(
              par = par0,
              fn = rss_bias,
              lower = c(
                rep(1e-4, length(bias) - 1),
                rep(1e-8, length(c(hd3_len, hd5_len)))
              ),
              upper = c(
                rep(1 - 1e-2, length(bias) - 1),
                rep(1 - 1e-6, length(c(hd3_len, hd5_len)))
              ),
              control = list(
                tol = init_tol, # decide search range
                maxfeval = maxfeval
              )
            )
            if (is.null(res)) NULL else res
          }
        )
        lst_par = lst_par[!vapply(lst_par, is.null, logical(1))]
        if (length(lst_par) > 0L) {
          par_mx = sapply(lst_par, function(opt_out){
            c(opt_out$value, opt_out$par)
          })
          par0_raw <- par_mx[-1, which.min(par_mx[1, ])]
        } else {
          par0_raw <- c(
            par_0$cut_bias$s7[-bias_lenght],
            par_0$prob_hd$p5, par_0$prob_hd$p3
          )
        }
      } else {
        par0_raw <- c(
          par_0$cut_bias$s7[-bias_lenght],
          par_0$prob_hd$p5, par_0$prob_hd$p3
        )
      }
      par0 <- par0_raw
      for (i in seq_along(sec_tol)) {
        diff_value <- 1e4
        value0 <- 0
        round0 <- 1
        opt_out <- safe_nmkb(
          par = par0_raw,
          fn = rss_bias,
          lower = c(
            rep(1e-4, length(bias) - 1),
            rep(1e-8, length(c(hd3_len, hd5_len)))
          ),
          upper = c(
            rep(1 - 1e-2, length(bias) - 1),
            rep(1 - 1e-6, length(c(hd3_len, hd5_len)))
          ),
          control = list(
            tol = sec_tol[i], # decide search range
            maxfeval = maxfeval
          )
        )
        if (!is.null(opt_out)) {
          value0 <- opt_out$value
          par0 <- opt_out$par
        } else {
          next
        }
        while ((abs(diff_value) > sec_diff[i]) & (round0 < iter_times)) {
          diff_value <- value0 - opt_out$value
          value0 <- opt_out$value
          par0 <- opt_out$par

          par0_bias <- par0[1:(length(bias) - 1)]
          par0_hd <- par0[-(1:(length(bias) - 1))]
          par0_bias[which(par0_bias > (1 - 1e-2))] <- runif(
            length(which(par0_bias > (1 - 1e-2))),
            (1 - 1e-2 - 2e-3),
            (1 - 1e-2 - 1e-3)
          )
          par0_bias[1:(length(bias) - 1)][which(par0_bias < (1e-4))] <- runif(
            length(which(par0_bias < (1e-4))), 1e-4 + 1e-5, 1e-4 + 2e-5
          )
          par0_hd[which(par0_hd > (1 - 1e-6))] <- runif(
            length(which(par0_hd > (1 - 1e-6))),
            (1 - 1e-6 - 2e-7),
            (1 - 1e-6 - 1e-7)
          )
          par0_hd[which(par0_hd < (1e-8))] <- runif(
            length(which(par0_hd < (1e-8))), 1e-8 + 1e-9, 1e-8 + 2e-9
          )
          round0 <- round0 + 1

          lst_opt <- safe_nmkb(
            par = par0,
            fn = rss_bias,
            lower = c(
              rep(1e-4, length(bias) - 1),
              rep(1e-8, length(c(hd3_len, hd5_len)))
            ),
            upper = c(
              rep(1 - 1e-2, length(bias) - 1),
              rep(1 - 1e-6, length(c(hd3_len, hd5_len)))
            ),
            control = list(
              tol = sec_tol[i], # decide search range
              maxfeval = maxfeval
            )
          )
          if (!is.null(lst_opt)) {
            opt_out <- lst_opt
            diff_value <- value0 - opt_out$value
            value0 <- opt_out$value
            par0 <- opt_out$par
          }
        }
        par0_raw <- par0
      }
    }

    par0 <- assign_par(
      out_par = par0,
      par_0 = par_0,
      reference_base_value = par_0$base_prob
    )

    return(par0)
  }

  convert_rpf_psite <- function(ribo_size) {
    convert_idx <- matrix(
      data = sapply(
        seq.int(ribo_size[1]), function(cut5_i) {
          sapply(
            seq.int(ribo_size[4]), function(cut3_i) {
              return(c(
                cut5_i,
                cut3_i,
                ribo_size[1L] - cut5_i + ribo_size[2L] + 1L,
                ribo_size[2L] + ribo_size[3L] + ribo_size[1L] - cut5_i + cut3_i
              ))
            }
          )
        }
      ),
      nrow = 4
    )
    return(convert_idx)
  }

  filter_prob_matrix <- function(candi_psite,
                                 candi_cut5,
                                 candi_cut3,
                                 candi_weight,
                                 rpf_expected,
                                 ribo_num,
                                 tx_info,
                                 par_0,
                                 num_per_group,
                                 ribo_size) {
    candi_frame <- matrixStats::rowMaxs(as.matrix(candi_cut5 %% 3L))

    len_rg <- seq.int(
      from = (par_0$ribo_size[2] + par_0$ribo_size[3] + 1L),
      to = (sum(par_0$ribo_size) - 1L)
    )

    lst_pos <- lapply(
      len_rg, function(len_i) {
        pos_i <- lapply(
          0:2, function(frame_i) {
            pos_i <- which(
              (rpf_expected$qwidth == len_i) & (candi_frame == frame_i)
            )[1:num_per_group]
            if (is.na(pos_i[1])) {
              return(NULL)
            } else {
              return(pos_i)
            }
          }
        )
        return(pos_i)
      }
    )

    vec_pos <- sort(na.omit(unlist(lst_pos)))

    obj_num <- rpf_expected$weight[vec_pos]

    rpf_seq <- Biostrings::DNAStringSet(
      x = tx_info$mix_tx
    )[rep(1L, length(rpf_expected$pos))]

    kmer_5 <- as.character(XVector::subseq(
      x = rpf_seq,
      start = rpf_expected$pos,
      width = 3L
    ))

    kmer_3 <- as.character(XVector::subseq(
      x = rpf_seq,
      end = rpf_expected$pos + rpf_expected$qwidth - 1L,
      width = 3L
    ))

    kmer_cof <- exp(
      as.numeric(
        (par_0$eff_intercept +
          par_0$eff_f5[kmer_5] +
          par_0$eff_f3[kmer_3])[vec_pos]
      )
    )

    candi_psite <- candi_psite[vec_pos, ]

    candi_weight <- candi_weight[vec_pos, ]

    candi_cut5 <- candi_cut5[vec_pos, ]

    candi_cut3 <- candi_cut3[vec_pos, ]

    uni_candi_psite <- sort(unique(candi_psite@x))

    cut_seq <- get_cleavage_seq(
      seqs = tx_info$mix_tx,
      p_site = uni_candi_psite,
      ribo_size = par_0$ribo_size
    )

    idx_ij <- match(candi_psite@x, uni_candi_psite)

    prob_mc_up <- stringr::str_split(
      string = as.character(cut_seq$up_seq),
      pattern = "",
      simplify = TRUE
    )

    prob_mc_dn <- stringr::str_split(
      string = as.character(cut_seq$dn_seq),
      pattern = "",
      simplify = TRUE
    )

    vec_5 <- length(uni_candi_psite) * (candi_cut5@x - 1L) + idx_ij

    vec_3 <- length(uni_candi_psite) * (candi_cut3@x - 1L) + idx_ij

    return(
      list(
        candi_psite = candi_psite,
        candi_weight = candi_weight,
        prob_mc_up = prob_mc_up,
        prob_mc_dn = prob_mc_dn,
        vec_5 = vec_5,
        vec_3 = vec_3,
        obj_num = obj_num,
        kmer_cof = kmer_cof
      )
    )
  }

  prepare_prob_matrix <- function(rpf_expected,
                                  convert_idx,
                                  ribo_num,
                                  rg_cds,
                                  mix_tx,
                                  par_0) {
    qwidth <- as.character(rpf_expected$qwidth)

    psite_num_idx <- as.vector(
      sapply(
        split(
          x = convert_idx[1, ], f = convert_idx[4, ]
        ), length
      )[qwidth]
    )
    psite_num_idx <- as.integer(psite_num_idx)
    valid_idx <- !is.na(psite_num_idx) & psite_num_idx > 0L
    if (!all(valid_idx)) {
      rpf_expected <- rpf_expected[valid_idx, , drop = FALSE]
      qwidth <- qwidth[valid_idx]
      psite_num_idx <- psite_num_idx[valid_idx]
    }

    psite_pos <- rep(x = rpf_expected$pos, times = psite_num_idx) +
      unlist(
        split(
          x = convert_idx[3, ], f = convert_idx[4, ]
        )[qwidth],
        use.names = FALSE
      )

    rg_psite <- IRanges::IRanges(start = psite_pos, width = 1L)

    idx_cds <- IRanges::findOverlaps(query = rg_psite, subject = rg_cds)

    idx_inframe <- idx_cds@from[
      (psite_pos[idx_cds@from] - (rg_cds@start + 1L)[idx_cds@to]) %%
        3L == 0L
    ]

    candi_cut3 <- Matrix::sparseMatrix(
      i = rep(
        x = seq.int(nrow(rpf_expected)), times = psite_num_idx
      )[idx_inframe],
      j = sequence(nvec = psite_num_idx, from = 1L)[idx_inframe],
      x = unlist(
        split(
          x = convert_idx[2, ], f = convert_idx[4, ]
        )[qwidth],
        use.names = FALSE
      )[idx_inframe]
    )

    candi_cut5 <- Matrix::sparseMatrix(
      i = rep(
        x = seq.int(nrow(rpf_expected)), times = psite_num_idx
      )[idx_inframe],
      j = sequence(nvec = psite_num_idx, from = 1L)[idx_inframe],
      x = unlist(
        split(
          x = convert_idx[1, ], f = convert_idx[4, ]
        )[qwidth],
        use.names = FALSE
      )[idx_inframe]
    )

    candi_psite <- Matrix::sparseMatrix(
      i = rep(
        x = seq.int(nrow(rpf_expected)), times = psite_num_idx
      )[idx_inframe],
      j = sequence(nvec = psite_num_idx, from = 1L)[idx_inframe],
      x = psite_pos[idx_inframe]
    )

    candi_p_prob <- Matrix::sparseMatrix(
      i = rep(
        x = seq.int(nrow(rpf_expected)), times = psite_num_idx
      )[idx_inframe],
      j = sequence(nvec = psite_num_idx, from = 1L)[idx_inframe],
      x = 1.1
    )

    uni_candi_psite <- as.integer(sort(unique(candi_psite@x)))

    cut_seq <- get_cleavage_seq(
      seqs = mix_tx,
      p_site = uni_candi_psite,
      ribo_size = par_0$ribo_size
    )

    idx_ij <- match(x = candi_psite@x, table = uni_candi_psite)

    cut_prob <- get_cleavage_prob(
      seqs = cut_seq,
      bias = par_0$cut_bias$s7,
      prob_hd5 = par_0$prob_hd$p5,
      prob_hd3 = par_0$prob_hd$p3
    )

    candi_p_prob@x <- cut_prob$final_p5[
      nrow(cut_prob$final_p5) * (candi_cut5@x - 1L) + idx_ij
    ] *
      cut_prob$final_p3[
        nrow(cut_prob$final_p3) * (candi_cut3@x - 1L) + idx_ij
      ] *
      ribo_num[candi_psite@x]

    idx_order <- order(Matrix::rowSums(candi_p_prob), decreasing = TRUE)

    candi_psite <- candi_psite[idx_order, ]

    candi_p_prob <- candi_p_prob[idx_order, ]

    candi_cut5 <- candi_cut5[idx_order, ]

    candi_cut3 <- candi_cut3[idx_order, ]

    rpf_expected <- rpf_expected[idx_order, ]

    return(list(
      candi_psite = candi_psite,
      candi_cut5 = candi_cut5,
      candi_cut3 = candi_cut3,
      candi_weight = candi_p_prob,
      rpf_expected = rpf_expected,
      uni_candi_psite = uni_candi_psite
    ))
  }

  p_site2rpf <- function(rpf_integrated,
                         p_site_choosed,
                         ribo_size) {
    convert_idx <- convert_rpf_psite(ribo_size = ribo_size)

    # Calculate the reads distribution of selected ribosomes
    df_read <- data.frame(
      pos = rep(
        x = p_site_choosed,
        each = ncol(convert_idx)
      ) - as.integer(convert_idx[3, ]),
      qwidth = rep(
        x = as.integer(convert_idx[4, ]),
        times = length(p_site_choosed)
      )
    )
    read_num_tag <- as.numeric(
      names(
        base::table(df_read$pos * 100 + df_read$qwidth)
      )
    )

    read_0_tag <- base::setdiff(read_num_tag, rpf_integrated$tag)

    df_exp <- rbind(
      data.frame(
        pos = as.integer(read_0_tag %/% 100),
        qwidth = as.integer(read_0_tag %% 100),
        weight = 0L
      ),
      rpf_integrated[rpf_integrated$tag %in% read_num_tag, -4]
    )

    return(df_exp)
  }

  prepare_par0 <- function(par_5add = c(0, 0, 0, 0, 1),
                           eff_f5 = rep(0, 64),
                           eff_f3 = rep(0, 64),
                           eff_intercept = 0,
                           base_prob = 0.1,
                           cut_bias,
                           mix_tx_len,
                           p_init) {
    prob_add5 <- par_5add
    names(prob_add5) <- c("A", "C", "G", "T")

    eff_intercept <- eff_intercept
    eff_f5 <- eff_f5
    eff_f3 <- eff_f3
    names(eff_f5) <-
      names(eff_f3) <-
      names(Biostrings::GENETIC_CODE)

    cut_bias <- cut_bias
    names(cut_bias) <- c("A", "C", "G", "T")

    par_0 <- list(
      prob_add5 = par_5add,
      ribo_size = p_init$ribo_size,
      prob_hd = list(
        p5 = (seq(810, 10, -800 / p_init$ribo_size[1]) / 1000)[
          1:p_init$ribo_size[1]
        ],
        p3 = (seq(10, 810, 800 / p_init$ribo_size[4]) / 1000)[
          1:p_init$ribo_size[4]
        ]
      ),
      base_prob = base_prob,
      cut_bias = list(
        s7 = cut_bias
      ),
      eff_intercept = eff_intercept,
      eff_f5 = eff_f5,
      eff_f3 = eff_f3
    )

    ribo_num <- vector(mode = "integer", length = mix_tx_len)
    ribo_num[p_init$p_pos] <- p_init$p_num

    return(
      list(
        par_0 = par_0,
        ribo_num = ribo_num
      )
    )
  }

  integrate_rpf <- function(rpf_info) {
    read_tag <- base::table(
      rpf_info$pos * 100 + rpf_info$qwidth
    )

    read_num_tag <- as.numeric(names(read_tag))

    read_weight <- data.frame(
      pos = as.integer(read_num_tag %/% 100),
      qwidth = as.integer(read_num_tag %% 100),
      weight = as.vector(read_tag),
      tag = read_num_tag
    )

    return(read_weight)
  }

  estimate_offsets <- function(tx_info,
                               rpf_info,
                               number_min,
                               max_frame_min = 0.3,
                               choose_terminal = "start",
                               limited_range = 8:18) {
    decide_offset <- function(df_read,
                              rg_annot,
                              choose_terminal,
                              limited_range = limited_range) {
      rg_read <- IRanges::IRanges(
        start = df_read$pos,
        width = df_read$qwidth
      )

      tmp_idx <- IRanges::findOverlaps(
        query = rg_annot,
        subject = rg_read,
        minoverlap = 3L,
        type = "within"
      )

      offsets <- rg_annot@start[tmp_idx@from] - df_read$pos[tmp_idx@to]

      freq_v <- base::table(offsets)[as.character(limited_range)]

      offset_p <- as.integer(as.numeric(names((which.max(freq_v)))))

      rg_len <- df_read$qwidth[1]

      if (choose_terminal == "start") {
        offset_p <- offset_p + 1L
      } else {
        offset_p <- rg_len - offset_p - 2L
      }

      return(c(rg_len, offset_p))
    }
    # Convert coordinate information into IRanges
    rg_reads <- IRanges::IRanges(
      start = rpf_info$pos,
      width = rpf_info$qwidth
    )

    rg_cds <- IRanges::IRanges(
      start = tx_info$mix_tx_pos$utr5_p3 + 1L,
      end = tx_info$mix_tx_pos$utr3_p5 - 1L
    )

    rg_codon_start <- IRanges::IRanges(
      start = tx_info$mix_tx_pos$utr5_p3 + 1L,
      width = 3L
    )

    rg_codon_stop <- IRanges::IRanges(
      start = tx_info$mix_tx_pos$utr3_p5 - 3L,
      width = 3L
    )

    idx_cds <- IRanges::findOverlaps(
      query = rg_reads,
      subject = rg_cds,
      minoverlap = 3L
    )

    rpf_info <- rpf_info[idx_cds@from, ]

    rpf_info$frame_idx <-
      (rg_reads@start[idx_cds@from] - rg_cds@start[idx_cds@to]) %% 3L

    lst_rpf <- base::split(x = rpf_info, f = rpf_info$qwidth)

    # filter groups by number
    gp_num <- sapply(lst_rpf, function(x) {
      nrow(x)
    })

    # choose groups
    max_gp <- which.max(gp_num)
    number_min <- sum(gp_num) * (1 - number_min) * 0.5
    if (number_min < 1000) {
      number_min <- 1000
    }
    if ((max_gp - 3) >= 1) {
      tmp1 <- sapply(1:(max_gp - 3), function(i) {
        if ((i == 1 | (i == 2))) {
          FALSE
        } else {
          (gp_num[i] < gp_num[i - 1]) & (gp_num[i] < gp_num[i + 1]) &
            (gp_num[i] > number_min) &
            ((gp_num[i] < gp_num[i - 2]) & (gp_num[i] < gp_num[i + 2]))
        }
      })
    } else {
      tmp1 <- logical(0)
    }
    if ((max_gp + 3) <= length(gp_num)) {
      tmp2 <- sapply((max_gp + 3):length(gp_num), function(i) {
        if ((i == length(gp_num)) | (i == (length(gp_num) - 1))) {
          FALSE
        } else {
          (gp_num[i] < gp_num[i - 1]) & (gp_num[i] < gp_num[i + 1]) &
            (gp_num[i] > number_min) &
            ((gp_num[i] < gp_num[i - 2]) & (gp_num[i] < gp_num[i + 2]))
        }
      })
    } else {
      tmp2 <- logical(0)
    }
    if (sum(tmp1) == 0) {
      rpf_len_5 <- min(which(gp_num > number_min))
    } else {
      rpf_len_5 <- max(which(tmp1))
    }
    if (sum(tmp2) == 0) {
      rpf_len_3 <- max(which(gp_num > number_min))
    } else {
      rpf_len_3 <- min(which(tmp2)) + max_gp
    }

    lst_rpf <- lst_rpf[names(gp_num)[rpf_len_5:rpf_len_3]]

    # filter groups by period
    gp_period <- sapply(
      lst_rpf, function(x) {
        tmp_v <- base::as.vector(base::table(c(0:2, x$frame_idx))) - 1L
        return(tmp_v / sum(tmp_v))
      }
    )

    lst_rpf <- lst_rpf[
      colnames(gp_period)[
        which(apply(gp_period, 2, max) > max_frame_min)
      ]
    ]

    m_offset <- lapply(
      lst_rpf, function(x) {
        if (choose_terminal == "start") {
          offsets <- decide_offset(
            df_read = x,
            rg_annot = rg_codon_start,
            choose_terminal = choose_terminal,
            limited_range = limited_range
          )
        } else {
          offsets <- decide_offset(
            df_read = x,
            rg_annot = rg_codon_stop,
            choose_terminal = choose_terminal,
            limited_range = limited_range
          )
        }
        return(offsets)
      }
    )

    colnames(m_offset) <- NULL

    return(m_offset)
  }

  infer_add5 <- function(df_rpf_add,
                         df_rpf_noadd,
                         tx_info,
                         eff_intercept = 0,
                         eff_f5 = rep(0, 64),
                         eff_f3 = rep(0, 64)) {
    if (unique(eff_f5) == 0) {
      names(eff_f5) <-
        names(eff_f3) <-
        names(Biostrings::GENETIC_CODE)
    }

    # Find the corrected single base added relation matrix
    names_gp <- paste(
      rep(c("A", "C", "G", "T"), each = 4),
      rep(c("A", "C", "G", "T"), 4),
      sep = ""
    )[c(-1, -6, -11, -16)]
    kmer_5 <- Biostrings::xscat(
      Biostrings::DNAStringSet(df_rpf_add$seq),
      Biostrings::subseq(
        Biostrings::DNAStringSet(tx_info$mix_tx)[rep(1L, nrow(df_rpf_add))],
        start = df_rpf_add$pos + 1,
        width = 2L
      )
    )
    kmer_3 <- Biostrings::subseq(
      Biostrings::DNAStringSet(
        tx_info$mix_tx
      )[rep(1L, nrow(df_rpf_add))],
      end = df_rpf_add$pos + df_rpf_add$qwidth - 1,
      width = 3
    )
    base_change <- paste(
      as.character(
        Biostrings::subseq(
          Biostrings::DNAStringSet(
            tx_info$mix_tx
          )[rep(1L, nrow(df_rpf_add))],
          start = df_rpf_add$pos,
          width = 1L
        )
      ),
      df_rpf_add$seq,
      sep = ""
    )

    correction_num <- exp(
      as.numeric(
        (eff_intercept +
          eff_f5[as.character(kmer_5)] +
          eff_f3[as.character(kmer_3)])
      )
    )

    freq_mat <- matrix(
      sapply(
        names_gp, function(name_gp) {
          sum(correction_num[name_gp == base_change])
        }
      ),
      nrow = 3
    )

    # Ratio of adding different bases, suppose add T probability is 1
    rss <- function(x) {
      freq_mat_exp <- matrix(
        c(x[-1], 1, x[-2], 1, x[-3], 1, x),
        nrow = 3
      )
      obj_m <- t(freq_mat) / Matrix::rowSums(t(freq_mat))
      exp_m <- t(freq_mat_exp) / Matrix::rowSums(t(freq_mat_exp))
      return(
        sum(
          (log(obj_m) - log(exp_m))^2
        )
      )
    }
    par0 <- 3:1
    opt_out <- dfoptim::nmkb(
      par = par0,
      fn = rss,
      lower = rep(.01, 3),
      upper = rep(100, 3)
    )
    add_prob <- c(opt_out$par, 1)

    # Ratio without adding bases
    kmer_5 <- Biostrings::subseq(
      Biostrings::DNAStringSet(
        tx_info$mix_tx
      )[rep(1L, nrow(df_rpf_noadd))],
      start = df_rpf_noadd$pos,
      width = 3
    )
    kmer_3 <- Biostrings::subseq(
      Biostrings::DNAStringSet(
        tx_info$mix_tx
      )[rep(1L, nrow(df_rpf_noadd))],
      end = df_rpf_noadd$pos + df_rpf_noadd$qwidth - 1,
      width = 3
    )

    correction_num <- exp(
      as.numeric(
        (eff_intercept +
          eff_f5[as.character(kmer_5)] +
          eff_f3[as.character(kmer_3)])
      )
    )

    num_miss_gp <- sum(
      colMeans(
        freq_mat /
          c(add_prob[-1], add_prob[-2], add_prob[-3], add_prob[-4]) *
          rep(add_prob, each = 3)
      )
    )

    add_prob <- c(
      add_prob,
      (sum(correction_num) - num_miss_gp) /
        (sum(freq_mat) + num_miss_gp) *
        sum(add_prob)
    )
    return(add_prob / sum(add_prob))
  }

  prep_rpf <- function(lst_rpf,
                       tx_info,
                       add5,
                       prob_add5,
                       number_min) {
                        #browser()
    names(prob_add5) <- c("A", "C", "G", "T", "")
    if (!add5) {
      rpf_info <- data.frame(
        pos = lst_rpf$no_add$pos,
        qwidth = lst_rpf$no_add$qwidth
      )

      offsets <- estimate_offsets(
        tx_info = tx_info,
        rpf_info = rpf_info,
        number_min = number_min,
        max_frame_min = 0.3,
        choose_terminal = "start",
        limited_range = 8:18
      )

      tmp_idx <- which(sapply(offsets, length) == 2)
      offsets <- do.call(cbind, offsets[tmp_idx])
      rpf_info <- rpf_info[rpf_info$qwidth %in% offsets[1, ], ]

      rpf_overlap <- integrate_rpf(rpf_info = rpf_info)
    } else {
      rpf_info_add <- data.frame(
        pos = lst_rpf$add5$pos,
        qwidth = lst_rpf$add5$qwidth
      )

      rpf_info_noadd <- data.frame(
        pos = lst_rpf$no_add$pos,
        qwidth = lst_rpf$no_add$qwidth
      )

      rpf_info <- rbind(rpf_info_noadd, rpf_info_add)

      offsets <- estimate_offsets(
        tx_info = tx_info,
        rpf_info = rpf_info,
        number_min = number_min,
        max_frame_min = 0.3,
        choose_terminal = "start",
        limited_range = 10:18
      )
      tmp_idx <- which(sapply(offsets, length) == 2)
      offsets <- do.call(cbind, offsets[tmp_idx])

      rpf_info_add <- rpf_info_add[rpf_info_add$qwidth %in% offsets[1, ], ]
      rpf_info_noadd <- rpf_info_noadd[
        rpf_info_noadd$qwidth %in% offsets[1, ],
      ]

      read_tag_add <- base::table(
        rpf_info_add$pos * 100 + rpf_info_add$qwidth + 99
      )

      read_num_tag_add <- as.numeric(names(read_tag_add))

      read_weight_add <- data.frame(
        pos = as.integer(read_num_tag_add %/% 100),
        qwidth = as.integer(read_num_tag_add %% 100),
        weight = as.vector(read_tag_add),
        tag = read_num_tag_add
      )

      ref_base <- as.character(
        Biostrings::subseq(
          Biostrings::DNAStringSet(
            tx_info$mix_tx
          )[rep(1, nrow(read_weight_add))],
          start = read_weight_add$pos - 1,
          width = 1
        )
      )

      tr_base <- sapply(1:4, function(x) {
        1 + prob_add5[x] / sum(prob_add5[1:4][-x])
      })

      read_weight_add$weight <- read_weight_add$weight * tr_base[ref_base]

      read_tag <- base::table(
        rpf_info_noadd$pos * 100 + rpf_info_noadd$qwidth
      )

      read_num_tag <- as.numeric(names(read_tag))

      read_weight <- data.frame(
        pos = as.integer(read_num_tag %/% 100),
        qwidth = as.integer(read_num_tag %% 100),
        weight = as.vector(read_tag),
        tag = read_num_tag
      )

      ref_base1 <- as.character(
        Biostrings::subseq(
          Biostrings::DNAStringSet(
            tx_info$mix_tx
          )[rep(1, length(read_num_tag))],
          start = read_num_tag %/% 100,
          width = 1
        )
      )

      tr_base2 <- sapply(1:4, function(x) {
        1 - prob_add5[x] /
          (prob_add5[5] +
            prob_add5[x])
      })

      read_weight$weight <- tr_base2[ref_base1] * read_weight$weight

      rpf_overlap <- rbind(read_weight, read_weight_add)
      rpf_info <- rpf_info[rpf_info$qwidth %in% offsets[1, ], ]
    }
    return(list(
      rpf_info = rpf_info,
      rpf_overlap = rpf_overlap,
      offsets = offsets
    ))
  }

  resolve_plot_file <- function(plot_file, suffix, default_prefix) {
    if (is.null(plot_file)) {
      ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
      return(file.path(tempdir(), paste0(default_prefix, "_", ts, "_", suffix, ".pdf")))
    }
    ext <- tools::file_ext(plot_file)
    stem <- if (nzchar(ext)) {
      sub(paste0("\\.", ext, "$"), "", plot_file)
    } else {
      plot_file
    }
    paste0(stem, "_", suffix, ".pdf")
  }

  plot_qwidth_distribution <- function(qwidth_all,
                                       plot_file,
                                       main_title) {
    if (length(qwidth_all) == 0L) {
      warning("No qwidth values available; skip read-length plot.")
      return(NULL)
    }
    counts <- table(as.integer(qwidth_all))
    counts <- counts[order(as.integer(names(counts)))]
    grDevices::pdf(file = plot_file, width = 8.5, height = 6)
    graphics::barplot(
      counts,
      xlab = "Read length (nt)",
      ylab = "Count",
      main = main_title,
      col = "#4E79A7",
      border = NA
    )
    graphics::box()
    grDevices::dev.off()
    plot_file
  }

  plot_filtered_read_length_distribution <- function(lst_rpf,
                                                     plot_file = NULL,
                                                     main_title = "Read length distribution for parameter estimation") {
    if (is.null(lst_rpf$rpf_info) || nrow(lst_rpf$rpf_info) == 0L) {
      warning("lst_rpf$rpf_info is empty; skip estimation read-length plot.")
      return(NULL)
    }
    plot_qwidth_distribution(
      qwidth_all = lst_rpf$rpf_info$qwidth,
      plot_file = plot_file,
      main_title = main_title
    )
  }

  initiate_p_site <- function(tx_info,
                              rg_cds,
                              rpf_info,
                              offsets,
                              ribo_size) {
    # initial p-site number
    p_table <- base::table(
      offsets[2, ][match(rpf_info$qwidth, offsets[1, ])] +
        rpf_info$pos
    )

    p_pos <- vector(
      mode = "integer",
      length = length(tx_info$mix_tx)
    )

    p_pos[as.integer(names(p_table))] <- as.vector(p_table)

    idx_p <- sequence(
      nvec = rg_cds@width,
      from = rg_cds@start
    )

    p_cds <- as.integer(matrixStats::colSums2(
      matrix(
        data = p_pos[idx_p],
        nrow = 3
      )
    ))

    idx_p <- matrix(data = idx_p, nrow = 3)[2, ]

    choose_p <- which(p_cds > 0L)

    p_pos <- idx_p[choose_p]

    p_num <- p_cds[choose_p]

    return(
      list(
        rpf_info = rpf_info,
        p_pos = p_pos,
        p_num = p_num,
        ribo_size = ribo_size
      )
    )
  }

  estimate_ribo_size <- function(offsets_mat, rpf_counts, RNase) {
    min_rpf <- as.numeric(names(rpf_counts[1]))
    if (RNase == "rnase-i") {
      len_gr <- (length(rpf_counts) + 1)
      S5 <- floor(min_rpf / 3)
      L5 <- ceiling((len_gr + 1) / 2)
      L3 <- ceiling((len_gr + 1) / 2)
    } else if (RNase == "mnase") {
      len_gr <- (length(rpf_counts) + 1)
      S5 <- floor(min_rpf / 3)
      L5 <- ceiling((len_gr + 1) / 2)
      L3 <- ceiling((len_gr + 1) / 2)
    } else {
      len_gr <- (length(rpf_counts) + 1)
      S5 <- floor(min_rpf / 3)
      L5 <- ceiling((len_gr + 1) / 2)
      L3 <- ceiling((len_gr + 1) / 2)
    }

    S3 <- min_rpf - S5 - 2
    return(as.integer(c(L5, S5, S3, L3)))
  }

  rg_cds <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p3 + 1L,
    end = tx_info$mix_tx_pos$utr3_p5 - 1L
  )
  trunc_cds <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p3 + 16L,
    end = tx_info$mix_tx_pos$utr3_p5 - 16L
  )
  is_add5 <- length(lst_bam) == 2
  if (is_add5) {
    message("Infer 5' addition parameters ... ", appendLF = FALSE)
    prob_add5 <- infer_add5(
      df_rpf_add = lst_bam$add5,
      df_rpf_noadd = lst_bam$no_add,
      tx_info = tx_info
    )
    message("Done")
  } else {
    prob_add5 <- c(0, 0, 0, 0, 1)
  }
  lst_rpf <- prep_rpf(
    lst_rpf = lst_bam,
    tx_info = tx_info,
    add5 = is_add5,
    prob_add5 = prob_add5,
    number_min = 0.99
  )

  ribo_size <- estimate_ribo_size(
    offsets_mat = lst_rpf$offsets,
    rpf_counts = table(lst_rpf$rpf_info$qwidth),
    RNase = RNase
  )

  init_p_site <- initiate_p_site(
    tx_info = tx_info,
    rg_cds = rg_cds,
    rpf_info = lst_rpf$rpf_info,
    offsets = lst_rpf$offsets,
    ribo_size = ribo_size
  )

  if (!(rnase_bias)) {
    par_0 <- prepare_par0(
      base_prob = 0.1,
      cut_bias = c(0.1, 0.1, 0.1, 0.1),
      mix_tx_len = length(tx_info$mix_tx),
      p_init = init_p_site
    )
  } else {
    par_0 <- prepare_par0(
      base_prob = 0.1,
      cut_bias = c(0.003, 0.7, 0.6, 0.1),
      mix_tx_len = length(tx_info$mix_tx),
      p_init = init_p_site
    )
  }
  convert_idx <- convert_rpf_psite(ribo_size = par_0$par_0$ribo_size)

  choose_p_site <- IRanges::IRanges(
    start = order(par_0$ribo_num, decreasing = TRUE),
    width = 1L
  )

  choosed_p <- choose_p_site@start[
    IRanges::findOverlaps(
      query = choose_p_site,
      subject = trunc_cds
    )@from
  ]

  potential_rpf <- p_site2rpf(
    rpf_integrated = lst_rpf$rpf_overlap,
    p_site_choosed = choosed_p[1:1e4],
    ribo_size = par_0$par_0$ribo_size
  )

  prob_matrix <- prepare_prob_matrix(
    rpf_expected = potential_rpf,
    convert_idx = convert_idx,
    ribo_num = par_0$ribo_num,
    rg_cds = rg_cds,
    mix_tx = tx_info$mix_tx,
    par_0 = par_0$par_0
  )

  prob_matrix_filter <- filter_prob_matrix(
    candi_psite = prob_matrix$candi_psite,
    candi_cut5 = prob_matrix$candi_cut5,
    candi_cut3 = prob_matrix$candi_cut3,
    candi_weight = prob_matrix$candi_weight,
    rpf_expected = prob_matrix$rpf_expected,
    ribo_num = par_0$ribo_num,
    tx_info = tx_info,
    par_0 = par_0$par_0,
    num_per_group = 1e3,
    ribo_size = par_0$par_0$ribo_size
  )

  par_round_1 <- iter_estimate(
    iter_times = 10,
    init_num = init_num,
    maxfeval = maxfeval,
    init_tol = 1e3,
    sec_tol = c(1e2, 1e1, 1e0),
    sec_diff = c(1e0, 1e-1, 1e-2),
    par_0 = par_0$par_0,
    par0_unknown = TRUE,
    candi_psite = prob_matrix_filter$candi_psite,
    candi_weight = prob_matrix_filter$candi_weight,
    prob_mc_up = prob_matrix_filter$prob_mc_up,
    prob_mc_dn = prob_matrix_filter$prob_mc_dn,
    vec_5 = prob_matrix_filter$vec_5,
    vec_3 = prob_matrix_filter$vec_3,
    rpf_group_5 = prob_matrix_filter$rpf_group_5,
    rpf_group_3 = prob_matrix_filter$rpf_group_3,
    obj_num = prob_matrix_filter$obj_num,
    ribo_num = par_0$ribo_num,
    kmer_cof = prob_matrix_filter$kmer_cof
  )
  
  # 大改动 估计——ribonum
  potential_rpf_2 <- p_site2rpf(
    rpf_integrated = lst_rpf$rpf_overlap,
    p_site_choosed = prob_matrix$uni_candi_psite,
    ribo_size = par_round_1$ribo_size
  )

  prob_matrix_2 <- prepare_prob_matrix(
    rpf_expected = potential_rpf_2,
    convert_idx = convert_idx,
    ribo_num = par_0$ribo_num,
    rg_cds = rg_cds,
    mix_tx = tx_info$mix_tx,
    par_0 = par_round_1
  )

  tmp_ribo_num <- em_ribo_num(
    df_rpf = prob_matrix_2$rpf_expected,
    candi_psite = prob_matrix_2$candi_psite,
    candi_cut5 = prob_matrix_2$candi_cut5,
    candi_cut3 = prob_matrix_2$candi_cut3,
    tx_info = tx_info,
    ribo_size = par_round_1$ribo_size,
    par_0 = par_round_1,
    pos_psite = prob_matrix$uni_candi_psite,
    iter_num = 8
  )

  # 第一轮更新ribo_num
  par_0$ribo_num[prob_matrix$uni_candi_psite] <- tmp_ribo_num
  par_0$par_0 <- par_round_1

  # 第二轮估计hd 和 cut_bias
  potential_rpf <- p_site2rpf(
    rpf_integrated = lst_rpf$rpf_overlap,
    p_site_choosed = choosed_p[1:1e4],
    ribo_size = par_0$par_0$ribo_size
  )

  prob_matrix <- prepare_prob_matrix(
    rpf_expected = potential_rpf,
    convert_idx = convert_idx,
    ribo_num = par_0$ribo_num,
    rg_cds = rg_cds,
    mix_tx = tx_info$mix_tx,
    par_0 = par_0$par_0
  )

  prob_matrix_filter <- filter_prob_matrix(
    candi_psite = prob_matrix$candi_psite,
    candi_cut5 = prob_matrix$candi_cut5,
    candi_cut3 = prob_matrix$candi_cut3,
    candi_weight = prob_matrix$candi_weight,
    rpf_expected = prob_matrix$rpf_expected,
    ribo_num = par_0$ribo_num,
    tx_info = tx_info,
    par_0 = par_0$par_0,
    num_per_group = 1e3,
    ribo_size = par_0$par_0$ribo_size
  )

  par_round_2 <- iter_estimate(
    iter_times = 10,
    init_num = init_num,
    maxfeval = maxfeval,
    init_tol = 1e3,
    sec_tol = c(1e2, 1e1, 1e0),
    sec_diff = c(1e0, 1e-1, 1e-2),
    par_0 = par_0$par_0,
    par0_unknown = FALSE,
    candi_psite = prob_matrix_filter$candi_psite,
    candi_weight = prob_matrix_filter$candi_weight,
    prob_mc_up = prob_matrix_filter$prob_mc_up,
    prob_mc_dn = prob_matrix_filter$prob_mc_dn,
    vec_5 = prob_matrix_filter$vec_5,
    vec_3 = prob_matrix_filter$vec_3,
    rpf_group_5 = prob_matrix_filter$rpf_group_5,
    rpf_group_3 = prob_matrix_filter$rpf_group_3,
    obj_num = prob_matrix_filter$obj_num,
    ribo_num = par_0$ribo_num,
    kmer_cof = prob_matrix_filter$kmer_cof
  )
  # 第二轮更新ribo_num
  potential_rpf_2 <- p_site2rpf(
    rpf_integrated = lst_rpf$rpf_overlap,
    p_site_choosed = prob_matrix$uni_candi_psite,
    ribo_size = par_round_2$ribo_size
  )

  prob_matrix_2 <- prepare_prob_matrix(
    rpf_expected = potential_rpf_2,
    convert_idx = convert_idx,
    ribo_num = par_0$ribo_num,
    rg_cds = rg_cds,
    mix_tx = tx_info$mix_tx,
    par_0 = par_round_2
  )

  tmp_ribo_num <- em_ribo_num(
    df_rpf = prob_matrix_2$rpf_expected,
    candi_psite = prob_matrix_2$candi_psite,
    candi_cut5 = prob_matrix_2$candi_cut5,
    candi_cut3 = prob_matrix_2$candi_cut3,
    tx_info = tx_info,
    ribo_size = par_round_2$ribo_size,
    par_0 = par_round_2,
    pos_psite = prob_matrix$uni_candi_psite,
    iter_num = 8
  )
  # 第二轮更新ribo_num后ribo_num最终稳定
  par_0$ribo_num[prob_matrix$uni_candi_psite] <- tmp_ribo_num
  par_0$par_0 <- par_round_2

  if (ligat_par) {
    # 估计ligation参数
    potential_rpf_lig_raw <- p_site2rpf(
      rpf_integrated = lst_rpf$rpf_overlap,
      p_site_choosed = choosed_p[1:1e5],
      ribo_size = par_0$par_0$ribo_size
    )

    prob_matrix_lig_raw <- prepare_prob_matrix(
      rpf_expected = potential_rpf_lig_raw,
      convert_idx = convert_idx,
      ribo_num = par_0$ribo_num,
      rg_cds = rg_cds,
      mix_tx = tx_info$mix_tx,
      par_0 = par_0$par_0
    )

    glm_result <- ligation_par(
      p_num = par_0$ribo_num,
      rpf = prob_matrix_lig_raw$rpf_expected,
      candi_psite = prob_matrix_lig_raw$candi_psite,
      candi_cut5 = prob_matrix_lig_raw$candi_cut5,
      candi_cut3 = prob_matrix_lig_raw$candi_cut3,
      candi_weight = prob_matrix_lig_raw$candi_weight,
      tx_info = tx_info,
      kmer_num = 1e4L,
      par_0 = par_0$par_0
    )

    lig_par <- glm2par(
      fit_coefs = glm_result$glm_coef$coefficients,
      par_0 = par_0$par_0
    )

    rpf_overlap_corrected <- lst_rpf$rpf_overlap
  rpf_overlap_corrected$weight <- correct_kmer(
    rpf = lst_rpf$rpf_overlap,
    tx_info = tx_info,
    par_0 = lig_par$par_0
  )

  names(prob_add5) <- c("A", "C", "G", "T", "")
  lig_par$par_0$prob_add5 <- prob_add5

  qwidth_range <- unique(convert_rpf_psite(lig_par$par_0$ribo_size)[4, ])
  lst_rpf$rpf_overlap <- lst_rpf$rpf_overlap[
    lst_rpf$rpf_overlap$qwidth %in% qwidth_range,
  ]
  lst_rpf$rpf_overlap <- lst_rpf$rpf_overlap[, -4]

  if (isTRUE(plot_lst_rpf_dist)) {
    est_plot_file <- resolve_plot_file(
      plot_file = plot_file,
      suffix = "estimation",
      default_prefix = "riboba_qwidth"
    )
    plot_out <- try(
      plot_filtered_read_length_distribution(
        lst_rpf = lst_rpf,
        plot_file = est_plot_file,
        main_title = "Read-length distribution used for parameter estimation"
      ),
      silent = TRUE
    )
    if (!inherits(plot_out, "try-error") && !is.null(plot_out)) {
      message("Saved filtered read-length distribution: ", plot_out)
    }
  }

  return(list(
    par_0 = lig_par$par_0,
    rpf_corrected = rpf_overlap_corrected[, -4],
    lst_rpf = lst_rpf,
    lig_par = lig_par$par_0
  ))
  }else{
    names(prob_add5) <- c("A", "C", "G", "T", "")
  par_0$par_0$prob_add5 <- prob_add5

  qwidth_range <- unique(convert_rpf_psite(par_0$par_0$ribo_size)[4, ])
  lst_rpf$rpf_overlap <- lst_rpf$rpf_overlap[
    lst_rpf$rpf_overlap$qwidth %in% qwidth_range,
  ]
  lst_rpf$rpf_overlap <- lst_rpf$rpf_overlap[, -4]
  if (isTRUE(plot_lst_rpf_dist)) {
    est_plot_file <- resolve_plot_file(
      plot_file = plot_file,
      suffix = "estimation",
      default_prefix = "riboba_qwidth"
    )
    plot_out <- try(
      plot_filtered_read_length_distribution(
        lst_rpf = lst_rpf,
        plot_file = est_plot_file,
        main_title = "Read-length distribution used for parameter estimation"
      ),
      silent = TRUE
    )
    if (!inherits(plot_out, "try-error") && !is.null(plot_out)) {
      message("Saved filtered read-length distribution: ", plot_out)
    }
  }
  return(list(
    par_0 = par_0$par_0,
    rpf_corrected = NULL,
    lst_rpf = lst_rpf,
    lig_par = NULL
  ))
  }
  
}
#' Plot learned parameter summary
#'
#' Create a compact multi-panel figure for key parameters learned by
#' [learn_bias_parameters()].
#'
#' @param par_result Output list from [learn_bias_parameters()] or a named list
#'   of parameter objects (`par_0`-like lists) from multiple samples.
#' @param output_file Output PDF path. If `NULL`, write to `tempdir()`.
#' @param top_kmer Number of top k-mers (by absolute effect size) to show for
#'   both 5' and 3' ligation effects.
#' @param p0 Reference probability used to normalize steric hindrance
#'   (`-log(p) / -log(p0)`).
#' @param eps Small floor value to avoid `-log(0)`.
#' @param base_family Base font family for plotting device.
#' @param base_size Base font size for plotting device.
#'
#' @return Invisibly returns the plot file path.
#' @export
plot_parameter_summary <- function(par_result,
                                   output_file = NULL,
                                   top_kmer = 15L,
                                   p0 = 0.9,
                                   eps = 1e-3,
                                   base_family = "sans",
                                   base_size = 10) {
  as_entry <- function(x) {
    if (is.list(x) && !is.null(x$par_0)) {
      return(list(par0 = x$par_0, lst_rpf = x$lst_rpf))
    }
    list(par0 = x, lst_rpf = NULL)
  }
  is_par0 <- function(x) {
    is.list(x) &&
      !is.null(x$prob_hd) &&
      !is.null(x$prob_hd$p5) &&
      !is.null(x$prob_hd$p3)
  }

  if (is_par0(as_entry(par_result)$par0)) {
    entry_list <- list(sample1 = as_entry(par_result))
  } else if (is.list(par_result) && length(par_result) > 0L) {
    entry_list <- lapply(par_result, as_entry)
    keep <- vapply(entry_list, function(x) is_par0(x$par0), logical(1))
    entry_list <- entry_list[keep]
    if (length(entry_list) == 0L) {
      stop("No valid parameter object found in par_result.")
    }
    if (is.null(names(entry_list)) || any(names(entry_list) == "")) {
      names(entry_list) <- paste0("sample", seq_along(entry_list))
    }
  } else {
    stop("par_result must be learn_bias_parameters() output or a list of parameter objects.")
  }
  par_list <- lapply(entry_list, `[[`, "par0")
  names(par_list) <- names(entry_list)

  lst_rpf_first <- NULL
  for (nm in names(entry_list)) {
    if (!is.null(entry_list[[nm]]$lst_rpf)) {
      lst_rpf_first <- entry_list[[nm]]$lst_rpf
      break
    }
  }

  if (is.null(output_file)) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_file <- file.path(tempdir(), paste0("riboba_parameter_summary_", ts, ".pdf"))
  }
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  draw_empty <- function(title_txt, body_txt = "No data") {
    graphics::plot.new()
    graphics::title(main = title_txt)
    graphics::text(0.5, 0.5, labels = body_txt)
  }

  infer_enzyme <- function(x) {
    if (stringr::str_detect(x, stringr::regex("RNase\\s*I", ignore_case = TRUE)) ||
      stringr::str_detect(x, stringr::regex("Martinez\\s*,\\s*(Epicentre|TruSeq)", ignore_case = TRUE))) {
      return("RNase I")
    }
    if (stringr::str_detect(x, stringr::regex("\\bMNase\\b", ignore_case = TRUE))) {
      return("MNase")
    }
    if (stringr::str_detect(x, stringr::regex("\\bP1\\b", ignore_case = TRUE))) {
      return("P1")
    }
    "Other"
  }

  build_hindrance_df <- function(par_list, p0, eps) {
    one_sample <- function(par0, id) {
      d5 <- data.frame(
        sample = id, end = "5prime",
        k = seq_along(par0$prob_hd$p5),
        m = as.numeric(par0$prob_hd$p5),
        stringsAsFactors = FALSE
      )
      d3 <- data.frame(
        sample = id, end = "3prime",
        k = seq_along(par0$prob_hd$p3),
        m = as.numeric(par0$prob_hd$p3),
        stringsAsFactors = FALSE
      )
      rbind(d5, d3)
    }
    df_h <- do.call(
      rbind,
      lapply(names(par_list), function(nm) one_sample(par_list[[nm]], nm))
    )
    df_h$m <- pmax(df_h$m, eps)
    df_h$h_norm <- -log(df_h$m) / (-log(p0))

    df_h <- do.call(
      rbind,
      lapply(split(df_h, interaction(df_h$sample, df_h$end, drop = TRUE)), function(d) {
        d <- d[order(d$k), ]
        if (unique(d$end) == "5prime") {
          d$dist <- -rev(seq_len(nrow(d)))
        } else {
          d$dist <- seq_len(nrow(d))
        }
        d
      })
    )

    ribo_size_mat <- sapply(par_list, function(x) x$ribo_size)
    if (is.null(dim(ribo_size_mat))) {
      ribo_size_mat <- matrix(ribo_size_mat, nrow = 4)
      colnames(ribo_size_mat) <- names(par_list)
    }
    shift_5 <- -as.numeric(ribo_size_mat[2, ])
    shift_3 <- as.numeric(ribo_size_mat[3, ])
    names(shift_5) <- colnames(ribo_size_mat)
    names(shift_3) <- colnames(ribo_size_mat)

    idx5 <- df_h$end == "5prime"
    df_h$dist[idx5] <- df_h$dist[idx5] + shift_5[df_h$sample[idx5]]
    df_h$dist[!idx5] <- df_h$dist[!idx5] + shift_3[df_h$sample[!idx5]]
    df_h$enzyme <- vapply(df_h$sample, infer_enzyme, character(1))
    df_h
  }

  draw_hindrance <- function(df_h, title_txt) {
    cols_enzyme <- c("RNase I" = "#0072B2", "MNase" = "#009E73", "P1" = "#D55E00", "Other" = "#666666")
    ends <- c("5prime", "3prime")
    xs <- range(df_h$dist, na.rm = TRUE)
    ys <- range(df_h$h_norm, na.rm = TRUE)
    ys[2] <- max(ys[2], 1)
    graphics::plot.new()
    graphics::plot.window(xlim = xs, ylim = ys)
    graphics::axis(1)
    graphics::axis(2)
    graphics::title(main = title_txt, xlab = "Distance to P-site (nt)", ylab = "Relative steric hindrance")
    graphics::abline(h = pretty(ys), col = "grey90", lwd = 0.8)
    for (ed in ends) {
      df_e <- df_h[df_h$end == ed, ]
      for (sp in unique(df_e$sample)) {
        d <- df_e[df_e$sample == sp, ]
        d <- d[order(d$dist), ]
        col_i <- cols_enzyme[unique(d$enzyme)]
        lty_i <- if (ed == "5prime") 1 else 2
        graphics::lines(d$dist, d$h_norm, col = col_i, lwd = 1.2, lty = lty_i)
        graphics::points(d$dist, d$h_norm, col = col_i, pch = 16, cex = 0.5)
      }
    }
    graphics::legend(
      "topright",
      legend = c("5'", "3'"),
      lty = c(1, 2),
      bty = "n",
      cex = 0.85
    )
  }

  draw_eff <- function(vec_eff, title_txt, top_kmer, col_txt) {
    if (is.null(vec_eff) || length(vec_eff) == 0L) {
      draw_empty(title_txt)
      return(invisible(NULL))
    }
    nz <- vec_eff[abs(vec_eff) > 0]
    if (length(nz) == 0L) {
      draw_empty(title_txt, "All k-mer effects are 0")
      return(invisible(NULL))
    }
    pos <- sort(nz[nz > 0], decreasing = TRUE)
    neg <- sort(nz[nz < 0], decreasing = FALSE)
    n_each <- max(1L, floor(top_kmer / 2L))
    pick_pos <- if (length(pos) > 0L) pos[seq_len(min(n_each, length(pos)))] else numeric(0)
    pick_neg <- if (length(neg) > 0L) neg[seq_len(min(n_each, length(neg)))] else numeric(0)
    pick <- c(pick_pos, pick_neg)
    if (length(pick) < min(top_kmer, length(nz))) {
      rest <- nz[setdiff(names(nz), names(pick))]
      rest <- rest[order(abs(rest), decreasing = TRUE)]
      pick <- c(pick, rest[seq_len(min(min(top_kmer, length(nz)) - length(pick), length(rest)))])
    }
    pick <- pick[order(pick, decreasing = TRUE)]
    cols <- ifelse(pick >= 0, "#D55E00", "#0072B2")
    graphics::barplot(
      pick,
      names.arg = names(pick),
      las = 2,
      cex.names = 0.6,
      col = cols,
      border = NA,
      ylab = "Ligation Bias",
      main = title_txt
    )
    graphics::abline(h = 0, lwd = 0.8, col = "grey40")
  }

  draw_prob_add5 <- function(par0, title_txt = "5' Addition Profile") {
    p <- par0$prob_add5
    if (is.null(p) || length(p) == 0L) {
      draw_empty(title_txt)
      return(invisible(NULL))
    }
    p <- as.numeric(p)
    nms <- names(par0$prob_add5)
    if (is.null(nms) || length(nms) != length(p)) {
      nms <- paste0("V", seq_along(p))
    }
    nms[nms == ""] <- "NoAdd"
    cols <- rep("#B07AA1", length(p))
    cols[nms == "NoAdd"] <- "#9C755F"
    graphics::barplot(
      p,
      names.arg = nms,
      col = cols,
      border = NA,
      ylab = "Probability",
      main = title_txt
    )
  }

  draw_mnase_cut_bias <- function(par0,
                                  title_txt = "Base Cleavage Bias",
                                  normalize = "none") {
    cb <- par0$cut_bias$s7
    if (is.null(cb)) {
      draw_empty(title_txt, "No base-bias data")
      return(invisible(NULL))
    }
    b <- cb[c("A", "C", "G", "T")]
    if (any(is.na(b))) {
      draw_empty(title_txt, "Missing A/C/G/T in cut_bias$s7")
      return(invisible(NULL))
    }
    m <- pmin(pmax(as.numeric(b), 1e-8), 1 - 1e-8)
    haz <- -log(m)
    rel <- switch(normalize,
      sum = haz / sum(haz),
      mean = haz / mean(haz),
      none = haz / min(haz)
    )
    cols <- c(A = "#4E79A7", C = "#59A14F", G = "#E15759", T = "#F28E2B")
    graphics::barplot(
      rel,
      names.arg = c("A", "C", "G", "T"),
      col = cols[c("A", "C", "G", "T")],
      border = "black",
      xlab = "",
      ylab = "Relative MNase cut bias (min = 1)",
      main = title_txt
    )
  }

  draw_qwidth <- function(lst_rpf, title_txt = "Filtered Read-length Distribution") {
    if (is.null(lst_rpf)) {
      draw_empty(title_txt, "No filtered read-length data")
      return(invisible(NULL))
    }
    qwidth_count <- function(qdf) {
      if (is.null(qdf) || !is.data.frame(qdf) || nrow(qdf) == 0L) {
        return(NULL)
      }
      len_col <- intersect(c("qwidth", "width", "read_length", "len", "length"), colnames(qdf))
      if (length(len_col) > 0L) {
        qv <- as.integer(qdf[[len_col[1]]])
      } else if ("seq" %in% colnames(qdf)) {
        qv <- nchar(as.character(qdf$seq))
      } else {
        return(NULL)
      }
      wt_col <- intersect(c("weight", "count", "n"), colnames(qdf))
      w <- if (length(wt_col) > 0L) as.numeric(qdf[[wt_col[1]]]) else rep(1, length(qv))
      keep <- !is.na(qv) & !is.na(w) & w > 0
      if (!any(keep)) {
        return(NULL)
      }
      tapply(w[keep], qv[keep], sum, na.rm = TRUE)
    }

    # Prefer lst_rpf$rpf_info as requested, fallback to rpf_overlap.
    tb <- qwidth_count(lst_rpf$rpf_info)
    if (is.null(tb)) {
      tb <- qwidth_count(lst_rpf$rpf_overlap)
    }
    if (is.null(tb)) {
      draw_empty(title_txt, "No qwidth in filtered reads")
      return(invisible(NULL))
    }
    tb <- tb[!is.na(names(tb))]
    if (length(tb) == 0L) {
      draw_empty(title_txt, "No qwidth in filtered reads")
      return(invisible(NULL))
    }
    tb <- tb[order(as.integer(names(tb)))]
    graphics::barplot(
      as.numeric(tb) / sum(tb),
      names.arg = names(tb),
      col = "#4E79A7",
      border = NA,
      xlab = "Read length (nt)",
      ylab = "Relative abundance",
      main = title_txt
    )
  }

  can_draw_cut_bias <- any(vapply(par_list, function(p) {
    cb <- p$cut_bias$s7
    !is.null(cb) && length(unique(round(cb, 8))) > 1L
  }, logical(1)))

  top_kmer <- as.integer(max(1L, top_kmer))
  grDevices::pdf(file = output_file, width = 13, height = 8.5)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::layout(matrix(1:6, nrow = 2, byrow = TRUE))
  graphics::par(mar = c(6, 4, 3, 1), family = base_family, cex = base_size / 10)

  draw_qwidth(lst_rpf = lst_rpf_first)

  df_h <- build_hindrance_df(par_list = par_list, p0 = p0, eps = eps)
  draw_hindrance(df_h, title_txt = "Ribosome Steric Hindrance")

  if (can_draw_cut_bias) {
    # draw first sample with non-flat cut bias
    idx_cb <- which(vapply(par_list, function(p) {
      cb <- p$cut_bias$s7
      !is.null(cb) && length(unique(round(cb, 8))) > 1L
    }, logical(1)))[1]
    draw_mnase_cut_bias(par0 = par_list[[idx_cb]], title_txt = "Base Cleavage Bias")
  } else {
    draw_empty("Base Cleavage Bias", "No base-bias difference")
  }

  draw_prob_add5(
    par0 = par_list[[1]],
    title_txt = "5' Addition Profile"
  )

  p_first <- par_list[[1]]
  draw_eff(
    vec_eff = p_first$eff_f5,
    title_txt = "Top 5' k-mer Ligation Bias",
    top_kmer = top_kmer,
    col_txt = "#76B7B2"
  )
  draw_eff(
    vec_eff = p_first$eff_f3,
    title_txt = "Top 3' k-mer Ligation Bias",
    top_kmer = top_kmer,
    col_txt = "#EDC948"
  )

  message("Saved parameter summary plot: ", output_file)
  invisible(output_file)
}
#' Plot ORF Biotype Composition Pie Chart
#'
#' Draw a pie chart of ORF counts grouped by `orf_biotype` from
#' [run_riboba_pipeline()] results.
#'
#' @param riboba_result A result object returned by [run_riboba_pipeline()], or a
#'   single-sample result entry, or an `.rds` path saved from that result.
#' @param output_file Output PDF path. If `NULL`, write to `tempdir()`.
#' @param sample_name Optional sample name when `riboba_result` contains multiple
#'   samples. Defaults to the first sample.
#' @param label_min_pct Minimum percentage threshold for drawing in-slice labels.
#'   Slices below this threshold are unlabeled on the pie and shown in legend.
#'
#' @return Invisibly returns the plot file path.
#' @export
plot_orf_biotype_pie <- function(riboba_result,
                                 output_file = NULL,
                                 sample_name = NULL,
                                 label_min_pct = 3) {
  if (is.character(riboba_result) && length(riboba_result) == 1L && file.exists(riboba_result)) {
    riboba_result <- readRDS(riboba_result)
  }

  get_single_sample <- function(x, sample_name = NULL) {
    if (is.list(x) && !is.null(x$orf_info) && is.list(x$orf_info) && !is.null(x$orf_info$ORF)) {
      return(list(sample = ifelse(is.null(sample_name), "sample1", sample_name), res = x))
    }
    if (!is.list(x) || length(x) == 0L) {
      stop("riboba_result must be a run_riboba_pipeline() result or a single-sample result entry.")
    }
    if (is.null(names(x)) || any(names(x) == "")) {
      names(x) <- paste0("sample", seq_along(x))
    }
    if (is.null(sample_name)) {
      sample_name <- names(x)[1]
    }
    if (!(sample_name %in% names(x))) {
      stop("sample_name not found in riboba_result: ", sample_name)
    }
    list(sample = sample_name, res = x[[sample_name]])
  }

  obj <- get_single_sample(riboba_result, sample_name = sample_name)
  sample_name <- obj$sample
  sample_res <- obj$res

  if (is.null(sample_res$orf_info) || is.null(sample_res$orf_info$ORF)) {
    stop("Cannot find orf_info$ORF in selected result.")
  }
  df_orf <- sample_res$orf_info$ORF
  if (!is.data.frame(df_orf) || !("orf_biotype" %in% colnames(df_orf))) {
    stop("orf_info$ORF must be a data.frame containing column 'orf_biotype'.")
  }
  if (nrow(df_orf) == 0L) {
    stop("No ORF rows found in orf_info$ORF.")
  }

  bio <- as.character(df_orf$orf_biotype)
  bio[is.na(bio) | bio == ""] <- "Unknown"
  counts <- sort(table(bio), decreasing = TRUE)

  if (is.null(output_file)) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_file <- file.path(tempdir(), paste0("riboba_orf_biotype_pie_", ts, ".pdf"))
  }
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  cols <- grDevices::hcl.colors(length(counts), palette = "Set 3")
  pct <- as.numeric(counts) / sum(counts) * 100
  pie_labels <- ifelse(
    pct >= label_min_pct,
    paste0(sprintf("%.1f", pct), "%"),
    ""
  )
  legend_labels <- paste0(
    names(counts),
    " (n=", as.integer(counts),
    ", ", sprintf("%.1f", pct), "%)"
  )

  grDevices::pdf(file = output_file, width = 13, height = 7.7)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(2.8, 1.1))
  graphics::par(mar = c(2, 2, 4, 1))
  graphics::pie(
    x = as.numeric(counts),
    labels = pie_labels,
    col = cols,
    border = "white",
    clockwise = TRUE,
    cex = 0.9,
    main = paste0("ORF Biotype Composition: ", sample_name, "\nN = ", sum(counts))
  )
  graphics::par(mar = c(2, 0.5, 4, 2))
  graphics::plot.new()
  graphics::legend(
    "left",
    legend = legend_labels,
    fill = cols,
    border = NA,
    bty = "n",
    cex = 0.9
  )
  grDevices::dev.off()

  message("Saved ORF biotype pie chart: ", output_file)
  invisible(output_file)
}

#' Plot ROC and ORF Biotype in One PDF
#'
#' @param riboba_result A result object returned by [run_riboba_pipeline()], or a
#'   single-sample result entry, or an `.rds` path saved from that result.
#' @param output_file Output PDF path. If `NULL`, write to `tempdir()`.
#' @param sample_name Optional sample name when `riboba_result` contains multiple
#'   samples. Defaults to the first sample.
#' @param label_min_pct Minimum percentage threshold for drawing in-slice labels.
#'
#' @return Invisibly returns the plot file path.
#' @export
plot_orf_roc_biotype_pdf <- function(riboba_result,
                                     output_file = NULL,
                                     sample_name = NULL,
                                     label_min_pct = 3) {
  if (is.character(riboba_result) && length(riboba_result) == 1L && file.exists(riboba_result)) {
    riboba_result <- readRDS(riboba_result)
  }

  get_single_sample <- function(x, sample_name = NULL) {
    if (is.list(x) && !is.null(x[["orf_info", exact = TRUE]]) &&
      is.list(x[["orf_info", exact = TRUE]]) &&
      !is.null(x[["orf_info", exact = TRUE]][["ORF", exact = TRUE]])) {
      return(list(sample = ifelse(is.null(sample_name), "sample1", sample_name), res = x))
    }
    if (!is.list(x) || length(x) == 0L) {
      stop("riboba_result must be a run_riboba_pipeline() result or a single-sample result entry.")
    }
    if (is.null(names(x)) || any(names(x) == "")) {
      names(x) <- paste0("sample", seq_along(x))
    }
    if (is.null(sample_name)) {
      sample_name <- names(x)[1]
    }
    if (!(sample_name %in% names(x))) {
      stop("sample_name not found in riboba_result: ", sample_name)
    }
    list(sample = sample_name, res = x[[sample_name]])
  }

  extract_model_reg <- function(sample_res) {
    safe_get <- function(x, keys) {
      cur <- x
      for (k in keys) {
        if (!is.list(cur) || is.null(cur[[k, exact = TRUE]])) {
          return(NULL)
        }
        cur <- cur[[k, exact = TRUE]]
      }
      cur
    }

    model_reg <- safe_get(sample_res, c("predict_mode", "orf_pred", "model_orf", "model_reg"))
    if (!is.null(model_reg)) {
      return(model_reg)
    }
    train_model <- safe_get(sample_res, c("predict_mode", "orf_pred", "train_model"))
    if (!is.null(train_model)) {
      return(train_model)
    }
    NULL
  }

  calc_roc <- function(score, label01) {
    ord <- order(score, decreasing = TRUE)
    label01 <- label01[ord]
    p <- sum(label01 == 1L)
    n <- sum(label01 == 0L)
    tp <- cumsum(label01 == 1L)
    fp <- cumsum(label01 == 0L)
    tpr <- c(0, tp / p, 1)
    fpr <- c(0, fp / n, 1)
    auc <- sum(diff(fpr) * (head(tpr, -1L) + tail(tpr, -1L)) / 2)
    list(fpr = fpr, tpr = tpr, auc = auc)
  }

  obj <- get_single_sample(riboba_result, sample_name = sample_name)
  sample_name <- obj$sample
  sample_res <- obj$res

  df_orf <- sample_res[["orf_info", exact = TRUE]][["ORF", exact = TRUE]]
  if (!is.data.frame(df_orf) || !("orf_biotype" %in% colnames(df_orf))) {
    stop("orf_info$ORF must be a data.frame containing column 'orf_biotype'.")
  }

  model_reg <- extract_model_reg(sample_res)
  if (is.null(model_reg) || is.null(model_reg$model) || is.null(model_reg$model$predictions)) {
    stop("Cannot find OOB predictions for ROC in result object.")
  }
  pred_mat <- as.matrix(model_reg$model$predictions)
  if (nrow(pred_mat) == 0L || ncol(pred_mat) == 0L) {
    stop("Empty OOB prediction matrix; cannot plot ROC.")
  }
  pos_col <- if ("1" %in% colnames(pred_mat)) "1" else colnames(pred_mat)[ncol(pred_mat)]
  score <- as.numeric(pred_mat[, pos_col])
  y_true <- c(rep(1L, model_reg$n_pos), rep(0L, model_reg$n_neg))
  if (length(score) != length(y_true)) {
    stop("OOB prediction length does not match n_pos + n_neg.")
  }
  roc_obj <- calc_roc(score = score, label01 = y_true)

  bio <- as.character(df_orf$orf_biotype)
  bio[is.na(bio) | bio == ""] <- "Unknown"
  counts <- sort(table(bio), decreasing = TRUE)
  pct <- as.numeric(counts) / sum(counts) * 100
  pie_labels <- ifelse(
    pct >= label_min_pct,
    paste0(sprintf("%.1f", pct), "%"),
    ""
  )
  legend_labels <- paste0(
    names(counts),
    " (n=", as.integer(counts),
    ", ", sprintf("%.1f", pct), "%)"
  )
  cols <- grDevices::hcl.colors(length(counts), palette = "Set 3")

  if (is.null(output_file)) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_file <- file.path(tempdir(), paste0("riboba_orf_roc_biotype_", ts, ".pdf"))
  }
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  grDevices::pdf(file = output_file, width = 14.5, height = 6.8)
  graphics::layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1.8, 1.8, 1.2))

  graphics::par(mar = c(4, 4, 3, 1))
  graphics::plot(
    roc_obj$fpr, roc_obj$tpr,
    type = "l", lwd = 2, col = "#1F77B4",
    xlab = "False Positive Rate",
    ylab = "True Positive Rate",
    main = sprintf("OOB ROC (AUC = %.3f)", roc_obj$auc),
    xlim = c(0, 1), ylim = c(0, 1)
  )
  graphics::abline(0, 1, lty = 2, col = "grey50")

  graphics::par(mar = c(2, 2, 3, 1))
  graphics::pie(
    x = as.numeric(counts),
    labels = pie_labels,
    col = cols,
    border = "white",
    clockwise = TRUE,
    cex = 0.9,
    main = paste0("ORF Biotype: ", sample_name, "\nN = ", sum(counts))
  )

  graphics::par(mar = c(2, 0.5, 3, 2))
  graphics::plot.new()
  graphics::legend(
    "left",
    legend = legend_labels,
    fill = cols,
    border = NA,
    bty = "n",
    cex = 0.85
  )
  grDevices::dev.off()

  message("Saved ROC + ORF biotype plot: ", output_file)
  invisible(output_file)
}

#' Generate Two Standard RiboBA Report PDFs
#'
#' @param par_result Output list from [learn_bias_parameters()] or a named list
#'   of parameter objects.
#' @param riboba_result A result object returned by [run_riboba_pipeline()], or
#'   an `.rds` path saved from that result.
#' @param output_dir Output directory for the two report PDFs.
#' @param sample_name Optional sample name when `riboba_result` contains multiple
#'   samples. Defaults to the first sample.
#' @param file_prefix Filename prefix for the two generated PDFs.
#' @param lst_rpf Optional filtered read object (`par_rpf$lst_rpf`). Use this
#'   when `par_result` does not contain `lst_rpf`.
#'
#' @return A named list with `parameter_pdf` and `roc_pie_pdf`.
#' @export
plot_riboba_two_pdfs <- function(par_result,
                                 riboba_result,
                                 output_dir = file.path(getwd(), "results"),
                                 sample_name = NULL,
                                 file_prefix = "riboba_report",
                                 lst_rpf = NULL) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  parameter_pdf <- file.path(output_dir, paste0(file_prefix, "_parameter_summary_", ts, ".pdf"))
  roc_pie_pdf <- file.path(output_dir, paste0(file_prefix, "_roc_pie_", ts, ".pdf"))

  par_input <- par_result
  if (!is.null(lst_rpf) && is.list(par_result) && !is.null(par_result$par_0) && is.null(par_result$lst_rpf)) {
    par_input <- c(par_result, list(lst_rpf = lst_rpf))
  }

  plot_parameter_summary(
    par_result = par_input,
    output_file = parameter_pdf
  )
  plot_orf_roc_biotype_pdf(
    riboba_result = riboba_result,
    output_file = roc_pie_pdf,
    sample_name = sample_name
  )

  invisible(list(
    parameter_pdf = parameter_pdf,
    roc_pie_pdf = roc_pie_pdf
  ))
}

#' Export RiboBA Result for Downstream Analysis
#'
#' Export one sample into three parts:
#' 1) learned sample parameters (`bias_info`, RDS),
#' 2) ORF annotation object (`orf_info$ORF`, RDS),
#' 3) ORF ribosome distribution (`orf_info$Psite`, RDS).
#'
#' Additionally export:
#' - `orf_info$ORF` as CSV table;
#' - ORF biotype composition pie chart (PDF).
#'
#' @param riboba_result A result object returned by [run_riboba_pipeline()], or
#'   a path to an `.rds` file saved from that result.
#' @param output_dir Output directory.
#' @param sample_name Optional sample name when `riboba_result` contains multiple
#'   samples. Defaults to the first sample.
#' @param file_prefix Filename prefix for exported files.
#' @param par_result Optional `learn_bias_parameters()` output (or a list
#'   containing `lst_rpf`) used to estimate filtered read count.
#' @param rrna_log Optional path to Bowtie rRNA log (`rrna.log`).
#' @param cds_log Optional path to CDS mapping log (`cds.log`).
#' @param lncrna_log Optional path to lncRNA mapping log (`lncrna.log`).
#' @param total_reads Optional total input read count. If `NULL`, parse from
#'   `rrna_log` when available.
#' @param unique_mapped_reads Optional unique mapped read count. If `NULL`, parse
#'   from `cds_log` + `lncrna_log` when available.
#' @param filtered_reads Optional filtered read count. If `NULL`, infer from
#'   `par_result$lst_rpf$rpf_info`.
#'
#' @return A named list with output file paths.
#' @export
export_riboba_downstream <- function(riboba_result,
                                     output_dir = file.path(getwd(), "results"),
                                     sample_name = NULL,
                                     file_prefix = "",
                                     par_result = NULL,
                                     rrna_log = NULL,
                                     cds_log = NULL,
                                     lncrna_log = NULL,
                                     total_reads = NULL,
                                     unique_mapped_reads = NULL,
                                     filtered_reads = NULL) {
  if (is.character(riboba_result) && length(riboba_result) == 1L && file.exists(riboba_result)) {
    riboba_result <- readRDS(riboba_result)
  }

  get_single_sample <- function(x, sample_name = NULL) {
    if (is.list(x) && !is.null(x[["orf_info", exact = TRUE]]) &&
      is.list(x[["orf_info", exact = TRUE]]) &&
      !is.null(x[["orf_info", exact = TRUE]][["ORF", exact = TRUE]])) {
      return(list(sample = ifelse(is.null(sample_name), "sample1", sample_name), res = x))
    }
    if (!is.list(x) || length(x) == 0L) {
      stop("riboba_result must be a run_riboba_pipeline() result or a single-sample result entry.")
    }
    if (is.null(names(x)) || any(names(x) == "")) {
      names(x) <- paste0("sample", seq_along(x))
    }
    if (is.null(sample_name)) {
      sample_name <- names(x)[1]
    }
    if (!(sample_name %in% names(x))) {
      stop("sample_name not found in riboba_result: ", sample_name)
    }
    list(sample = sample_name, res = x[[sample_name]])
  }

  extract_model_reg <- function(sample_res) {
    safe_get <- function(x, keys) {
      cur <- x
      for (k in keys) {
        if (!is.list(cur) || is.null(cur[[k, exact = TRUE]])) {
          return(NULL)
        }
        cur <- cur[[k, exact = TRUE]]
      }
      cur
    }
    model_reg <- safe_get(sample_res, c("predict_mode", "orf_pred", "model_orf", "model_reg"))
    if (!is.null(model_reg)) {
      return(model_reg)
    }
    safe_get(sample_res, c("predict_mode", "orf_pred", "train_model"))
  }

  calc_roc_auc <- function(model_reg) {
    if (is.null(model_reg) || is.null(model_reg$model) || is.null(model_reg$model$predictions) ||
      is.null(model_reg$n_pos) || is.null(model_reg$n_neg)) {
      return(NA_real_)
    }
    pred_mat <- as.matrix(model_reg$model$predictions)
    if (nrow(pred_mat) == 0L || ncol(pred_mat) == 0L) {
      return(NA_real_)
    }
    pos_col <- if ("1" %in% colnames(pred_mat)) "1" else colnames(pred_mat)[ncol(pred_mat)]
    score <- as.numeric(pred_mat[, pos_col])
    y_true <- c(rep(1L, model_reg$n_pos), rep(0L, model_reg$n_neg))
    if (length(score) != length(y_true)) {
      return(NA_real_)
    }
    ord <- order(score, decreasing = TRUE)
    y_true <- y_true[ord]
    p <- sum(y_true == 1L)
    n <- sum(y_true == 0L)
    if (p == 0L || n == 0L) {
      return(NA_real_)
    }
    tp <- cumsum(y_true == 1L)
    fp <- cumsum(y_true == 0L)
    tpr <- c(0, tp / p, 1)
    fpr <- c(0, fp / n, 1)
    sum(diff(fpr) * (head(tpr, -1L) + tail(tpr, -1L)) / 2)
  }

  parse_bowtie_log <- function(log_path) {
    parse_first_int <- function(txt) {
      m <- regmatches(txt, regexpr("[0-9]+", txt))
      if (length(m) == 0L || identical(m, character(0))) {
        return(NA_real_)
      }
      as.numeric(m[1])
    }
    out <- list(reads_processed = NA_real_, reads_aligned = NA_real_)
    if (is.null(log_path) || !nzchar(log_path) || !file.exists(log_path)) {
      return(out)
    }
    x <- readLines(log_path, warn = FALSE)
    proc_line <- grep("^# reads processed:", x, value = TRUE)
    aln_line <- grep("^# reads with at least one alignment:", x, value = TRUE)
    if (length(proc_line) > 0L) {
      out$reads_processed <- parse_first_int(proc_line[1])
    }
    if (length(aln_line) > 0L) {
      out$reads_aligned <- parse_first_int(aln_line[1])
    }
    out
  }

  infer_filtered_reads <- function(par_result) {
    if (is.null(par_result) || !is.list(par_result) || is.null(par_result$lst_rpf)) {
      return(NA_real_)
    }
    rpf_info <- par_result$lst_rpf$rpf_info
    if (is.null(rpf_info) || !is.data.frame(rpf_info) || nrow(rpf_info) == 0L) {
      return(NA_real_)
    }
    if ("weight" %in% colnames(rpf_info)) {
      return(sum(as.numeric(rpf_info$weight), na.rm = TRUE))
    }
    as.numeric(nrow(rpf_info))
  }

  obj <- get_single_sample(riboba_result, sample_name = sample_name)
  sample_name <- obj$sample
  sample_res <- obj$res

  bias_info <- sample_res[["bias_info", exact = TRUE]]
  orf_df <- sample_res[["orf_info", exact = TRUE]][["ORF", exact = TRUE]]
  psite_info <- sample_res[["orf_info", exact = TRUE]][["Psite", exact = TRUE]]

  if (is.null(bias_info) || !is.list(bias_info)) {
    stop("Cannot find list object: bias_info")
  }
  if (is.null(psite_info) || !is.list(psite_info)) {
    stop("Cannot find list object: orf_info$Psite")
  }
  if (!is.data.frame(orf_df)) {
    stop("Cannot find table object: orf_info$ORF")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  sample_safe <- gsub("[^A-Za-z0-9._-]+", "_", sample_name)
  base_name <- if (!is.null(file_prefix) && nzchar(file_prefix)) {
    paste0(file_prefix, "_", sample_safe)
  } else {
    sample_safe
  }

  bias_rds <- file.path(output_dir, paste0(base_name, ".params.rds"))
  orf_rds <- file.path(output_dir, paste0(base_name, ".orf.rds"))
  psite_rds <- file.path(output_dir, paste0(base_name, ".psite.rds"))
  orf_csv <- file.path(output_dir, paste0(base_name, ".orf.csv"))
  biotype_pie_pdf <- file.path(output_dir, paste0(base_name, ".biotype.pdf"))

  saveRDS(bias_info, bias_rds)
  saveRDS(orf_df, orf_rds)
  saveRDS(psite_info, psite_rds)
  utils::write.csv(orf_df, file = orf_csv, row.names = FALSE, quote = TRUE)

  rrna_stat <- parse_bowtie_log(rrna_log)
  cds_stat <- parse_bowtie_log(cds_log)
  lncrna_stat <- parse_bowtie_log(lncrna_log)

  if (is.null(total_reads)) {
    total_reads <- rrna_stat$reads_processed
  }
  if (is.null(unique_mapped_reads)) {
    mapped_parts <- c(cds_stat$reads_aligned, lncrna_stat$reads_aligned)
    if (all(is.na(mapped_parts))) {
      unique_mapped_reads <- NA_real_
    } else {
      unique_mapped_reads <- sum(mapped_parts, na.rm = TRUE)
    }
  }
  if (is.null(filtered_reads)) {
    filtered_reads <- infer_filtered_reads(par_result)
  }

  bio <- as.character(orf_df$orf_biotype)
  annotated_orf_n <- sum(bio %in% c("CDS", "Ext", "Trunc"), na.rm = TRUE)
  ncorf_n <- sum(!(bio %in% c("CDS", "Ext", "Trunc")), na.rm = TRUE)
  roc_auc <- calc_roc_auc(extract_model_reg(sample_res))

  plot_orf_biotype_pie(
    riboba_result = sample_res,
    output_file = biotype_pie_pdf,
    sample_name = sample_name
  )

  message("Saved sample parameters (RDS): ", bias_rds)
  message("Saved ORF info (RDS): ", orf_rds)
  message("Saved ORF table (CSV): ", orf_csv)
  message("Saved ORF ribosome distribution (RDS): ", psite_rds)
  message("Saved ORF biotype pie chart (PDF): ", biotype_pie_pdf)

  invisible(list(
    sample_name = sample_name,
    sample_parameters_rds = bias_rds,
    orf_info_rds = orf_rds,
    orf_table_csv = orf_csv,
    orf_psite_distribution_rds = psite_rds,
    orf_biotype_pie_pdf = biotype_pie_pdf
  ))
}
#' Parse transcriptome-aligned BAM files
#'
#' Reads transcript-aligned BAM files into lightweight data frames containing
#' position, width, sequence, and mapping metadata used throughout the RiboBC
#' pipeline.
#'
#' @param bam_path Character vector with BAM file paths.
#' @param bai_path Character vector with index file paths that correspond to
#'   `bam_path`.
#' @param tx_info Transcript annotation list produced by preprocessing helpers.
#' @param only_uni_map Logical flag; when `TRUE`, only uniquely mapped reads are
#'   retained.
#' @param add_ratio Mismatch tolerance ratio when filtering alignments.
#'
#' @return A list of data frames, one per BAM file.
#'
#' @export
read_bam_alignments <- function(
    bam_path,
    bai_path,
    tx_info,
    only_uni_map,
    add_ratio) {
  # Convert transcript coordinates to merged transcript coordinates
  if (length(names(tx_info)) == 4) {
    pos_shift <- tx_info$mix_tx_pos$utr5_p5 - 1L
    tx_ids <- tx_info$tx_lens$tx_name
  } else {
    pos_shift <- tx_info$mix_tx_pos$tx_p5 - 1L
    tx_ids <- names(tx_info$tx_seqs)
  }
  if (!is.null(tx_ids) && length(tx_ids) == length(pos_shift)) {
    names(pos_shift) <- tx_ids
  }
  # multiple-map and mismatch
  par_scan <- Rsamtools::ScanBamParam(
    what = c("rname", "pos", "qwidth", "seq"),
    tag = c("MD", "XM")
  )
  file_bam <- bam_path
  bai_bam <- bai_path
  lst_bam <- lapply(
    seq_along(file_bam),
    function(i) {
      res <- Rsamtools::scanBam(
        file = file_bam[i],
        index = bai_bam[i],
        param = par_scan
      )[[1]]
      res$multi <- res$tag$XM
      res$tag <- res$tag$MD
      res$seq <- Biostrings::subseq(
        res$seq,
        start = 1,
        width = 1
      )
      res <- data.frame(res, stringsAsFactors = FALSE)
      # if use only uni-map reads
      if (only_uni_map) {
        res <- res[res$multi == 1L, -6]
      }
      shift <- unname(pos_shift[as.character(res$rname)])
      keep <- !is.na(shift)
      if (!all(keep)) {
        res <- res[keep, , drop = FALSE]
        shift <- shift[keep]
      }
      # only use no mismatch and 5' one base addition
      map_type <- ifelse(nchar(res$tag) == 2L, "u", "m")
      map_type[grep("^0", res$tag)] <- "a"
      res$pos <- shift + res$pos
      n_u <- sum(map_type == "u")
      n_a <- sum(map_type == "a")
      if (n_u == 0L) {
        # If exact-match reads are absent, skip add5 fitting to avoid
        # unstable add/no_add split driven by zero denominator.
        idx_no_add <- map_type != "a"
        if (!any(idx_no_add)) {
          idx_no_add <- rep(TRUE, nrow(res))
        }
        tmp_res <- list(
          no_add = res[idx_no_add, c(-1, -4, -5)]
        )
      } else if (n_a < (add_ratio * n_u)) {
        tmp_res <- list(
          no_add = res[map_type == "u", c(-1, -4, -5)]
        )
      } else {
        tmp_res <- list(
          no_add = res[map_type == "u", c(-1, -4, -5)],
          add5 = res[map_type == "a", c(-1, -5)]
        )
      }
      return(tmp_res)
    }
  )
  tmp_names <- stringr::str_split(
    bam_path,
    pattern = "_tx|/"
  )
  names(lst_bam) <- sapply(tmp_names, function(x) x[length(x) - 1])
  return(lst_bam)
}
#'
#' Prepare ncRNA ribosome footprints
#'
#' Applies the learnt parameters from CDS libraries to ncRNA-aligned footprints,
#' returning a harmonised object that can be passed to [predict_orf()].
#'
#' @param lst_rpf Output of [prep_rpf()] for ncRNA alignments.
#' @param tx_info Transcript annotation list for ncRNAs.
#' @param par_learn Parameter list returned by [learn_bias_parameters()] based on CDS data.
#'
#' @return A list containing processed footprints and inherited parameters.
#'
#' @export
prepare_ncrna_rpf <- function(
    lst_rpf,
    tx_info,
    par_learn) {
  # browser()
  integrate_rpf <- function(rpf_info) {
    read_tag <- base::table(
      rpf_info$pos * 100 + rpf_info$qwidth
    )

    read_num_tag <- as.numeric(names(read_tag))

    read_weight <- data.frame(
      pos = as.integer(read_num_tag %/% 100),
      qwidth = as.integer(read_num_tag %% 100),
      weight = as.vector(read_tag),
      tag = read_num_tag
    )

    return(read_weight)
  }

  prob_add5 <- par_learn$par_0$prob_add5
  if (prob_add5[5] == 1) {
    rpf_info <- data.frame(
      pos = lst_rpf$no_add$pos,
      qwidth = lst_rpf$no_add$qwidth
    )

    offsets <- par_learn$lst_rpf$offsets

    rpf_info <- rpf_info[rpf_info$qwidth %in% offsets[1, ], ]

    rpf_overlap <- integrate_rpf(rpf_info = rpf_info)
  } else {
    rpf_info_add <- data.frame(
      pos = lst_rpf$add5$pos,
      qwidth = lst_rpf$add5$qwidth
    )

    rpf_info_noadd <- data.frame(
      pos = lst_rpf$no_add$pos,
      qwidth = lst_rpf$no_add$qwidth
    )

    rpf_info <- rbind(rpf_info_noadd, rpf_info_add)

    offsets <- par_learn$lst_rpf$offsets

    rpf_info_add <- rpf_info_add[rpf_info_add$qwidth %in% offsets[1, ], ]

    rpf_info_noadd <- rpf_info_noadd[
      rpf_info_noadd$qwidth %in% offsets[1, ],
    ]

    read_tag_add <- base::table(
      rpf_info_add$pos * 100 + rpf_info_add$qwidth + 99
    )

    read_num_tag_add <- as.numeric(names(read_tag_add))

    read_weight_add <- data.frame(
      pos = as.integer(read_num_tag_add %/% 100),
      qwidth = as.integer(read_num_tag_add %% 100),
      weight = as.vector(read_tag_add),
      tag = read_num_tag_add
    )

    ref_base <- as.character(
      Biostrings::subseq(
        Biostrings::DNAStringSet(
          tx_info$mix_tx
        )[rep(1, nrow(read_weight_add))],
        start = read_weight_add$pos - 1,
        width = 1
      )
    )

    tr_base <- sapply(1:4, function(x) {
      1 + prob_add5[x] / sum(prob_add5[1:4][-x])
    })

    read_weight_add$weight <- read_weight_add$weight * tr_base[ref_base]

    read_tag <- base::table(
      rpf_info_noadd$pos * 100 + rpf_info_noadd$qwidth
    )

    read_num_tag <- as.numeric(names(read_tag))

    read_weight <- data.frame(
      pos = as.integer(read_num_tag %/% 100),
      qwidth = as.integer(read_num_tag %% 100),
      weight = as.vector(read_tag),
      tag = read_num_tag
    )

    ref_base1 <- as.character(
      Biostrings::subseq(
        Biostrings::DNAStringSet(
          tx_info$mix_tx
        )[rep(1, length(read_num_tag))],
        start = read_num_tag %/% 100,
        width = 1
      )
    )

    tr_base2 <- sapply(1:4, function(x) {
      1 - prob_add5[x] /
        (prob_add5[5] +
          prob_add5[x])
    })

    read_weight$weight <- tr_base2[ref_base1] * read_weight$weight

    rpf_overlap <- rbind(read_weight, read_weight_add)
    rpf_info <- rpf_info[rpf_info$qwidth %in% offsets[1, ], ]
  }

  return(
    list(
      lst_rpf = list(
        rpf_info = rpf_info,
        rpf_overlap = rpf_overlap,
        offsets = offsets
      ),
      par_0 = par_learn$par_0
    )
  )
}
#'
#'
#'
#' Predict ribosome occupancy for candidate ORFs
#'
#' @param orf_candi Candidate ORF ranges, typically output from [predict_orf()].
#' @param par_setup Parameter list (`par_0`) learnt by [learn_bias_parameters()].
#' @param df_rpf Footprint overlap table returned by [learn_bias_parameters()] or
#'   [prepare_ncrna_rpf()].
#' @param tx_info Transcript annotation list.
#' @param iter_num Number of refinement iterations.
#' @param min_ribosome_num Minimum ribosome density retained per codon.
#'
#' @return A list describing predicted P-site occupancies for each ORF class.
#'
#' @export
estimate_orf_psite <- function(
    orf_candi,
    par_setup,
    df_rpf,
    tx_info,
    iter_num,
    min_ribosome_num = 0.5) {
  filter_ribosome <- function(candi_psite,
                              min_ribosome_num) {
    idx_flt_pos <- candi_psite$vec_pnum > min_ribosome_num
    return(
      list(
        vec_pnum = round(candi_psite$vec_pnum[idx_flt_pos], 1),
        candi_pos = candi_psite$candi_pos[idx_flt_pos]
      )
    )
  }

  get_cleavage_seq <- function(seqs,
                               p_site,
                               ribo_size) {
    seqs <- Biostrings::DNAStringSet(x = seqs)[
      rep(1L, length(p_site))
    ]

    up_seq <- Biostrings::subseq(
      x = seqs,
      start = p_site - sum(ribo_size[1:2]),
      width = ribo_size[1]
    )

    dn_seq <- Biostrings::subseq(
      x = seqs,
      start = p_site + ribo_size[3],
      width = ribo_size[4]
    )

    return(
      list(
        up_seq = up_seq,
        dn_seq = dn_seq
      )
    )
  }

  maintain_prob <- function(cut_seq,
                            prod_hd,
                            bias) {
    prob_mc <- stringr::str_split(
      string = as.vector(cut_seq),
      pattern = "",
      simplify = TRUE
    )

    prob_mp <- matrix(
      data = bias[prob_mc],
      nrow = nrow(prob_mc)
    )

    prob_mp <- prob_mp^(rep(1L, nrow(prob_mp)) %o% prod_hd)

    return(prob_mp)
  }

  get_cleavage_prob <- function(seqs,
                                bias,
                                prob_hd5,
                                prob_hd3) {
    # for 5'
    maintain_prob5 <- maintain_prob(
      cut_seq = seqs$up_seq,
      prod_hd = prob_hd5,
      bias = bias
    )

    cle_p5rev <- 1 - maintain_prob5[, length(prob_hd5):1]

    maintain_cumprod5rev <- matrix(
      data = 1,
      nrow = nrow(maintain_prob5),
      ncol = length(prob_hd5)
    )

    maintain_cumprod5rev[, -1] <- matrixStats::rowCumprods(
      maintain_prob5[, length(prob_hd5):2]
    )

    final_p5 <- maintain_cumprod5rev * cle_p5rev

    final_p5 <- final_p5[, ncol(final_p5):1]

    final_p5 <- final_p5 / Matrix::rowSums(final_p5)

    # for 3'
    maintain_prob3 <- maintain_prob(
      cut_seq = seqs$dn_seq,
      prod_hd = prob_hd3,
      bias = bias
    )

    cle_p3 <- 1 - maintain_prob3

    maintain_cumprod3 <- matrix(
      data = 1,
      nrow = nrow(maintain_prob3),
      ncol = length(prob_hd3)
    )

    maintain_cumprod3[, -1] <- matrixStats::rowCumprods(
      maintain_prob3[, -length(prob_hd3)]
    )

    final_p3 <- maintain_cumprod3 * cle_p3

    final_p3 <- final_p3 / Matrix::rowSums(final_p3)

    return(
      list(
        final_p5 = final_p5,
        final_p3 = final_p3
      )
    )
  }
  # iterate n times to reach balance
  pred_psite <- function(df_rpf,
                         candi_psite,
                         candi_cut5,
                         candi_cut3,
                         tx_info,
                         ribo_size,
                         par_0,
                         iter_num) {
    candi_p_weight <- candi_psite # <- candi_psite[idx_read, ]
    # shrink p-sites
    shrink_p <- as.factor(candi_psite@x)
    shrink_pos <- as.integer(shrink_p)
    index_j <- rank(
      shrink_pos,
      ties.method = "first"
    ) -
      rank(shrink_pos, ties.method = "min") + 1L
    candi_pos <- as.integer(levels(shrink_p))
    vec_pnum <- rep(1, length(candi_pos))
    # Split read weight
    uni_candi_psite <- unique(candi_psite@x)
    cut_seq <- get_cleavage_seq(
      seqs = tx_info$mix_tx, p_site = uni_candi_psite, ribo_size = ribo_size
    )
    idx_ij <- match(candi_psite@x, uni_candi_psite)
    # Generating cleavage probabilities for each ribosome terminus
    cut_prob <- get_cleavage_prob(
      seqs = cut_seq, bias = par_0$cut_bias$s7,
      prob_hd5 = par_0$prob_hd$p5, prob_hd3 = par_0$prob_hd$p3
    )
    base_prob <- cut_prob$final_p5[
      nrow(cut_prob$final_p5) * (candi_cut5@x - 1L) + idx_ij
    ] *
      cut_prob$final_p3[nrow(cut_prob$final_p3) * (candi_cut3@x - 1L) + idx_ij]
    prob_cum <- vector(mode = "numeric", length = iter_num)

    for (i in 1:iter_num) {
      candi_p_weight@x <- base_prob * vec_pnum[shrink_pos]
      tmp_sum <- Matrix::rowSums(candi_p_weight)
      prob_cum[i] <- sum(df_rpf$weight * log(tmp_sum))
      sparse_iter <- Matrix::sparseMatrix(
        i = shrink_pos,
        j = index_j,
        x = (df_rpf$weight * candi_p_weight / tmp_sum)@x
      )
      vec_pnum[] <- Matrix::rowSums(sparse_iter)
    }
    return(
      list(
        vec_pnum = vec_pnum,
        candi_pos = candi_pos,
        prob_cum = prob_cum
      )
    )
  }

  ribo_size <- par_setup$ribo_size
  cut5len <- seq.int(ribo_size[1])
  cut3len <- seq.int(ribo_size[4])
  convert_idx <- matrix(data = sapply(cut5len, function(cut5_i) {
    sapply(cut3len, function(cut3_i) {
      return(c(
        cut5_i,
        cut3_i,
        ribo_size[1L] - cut5_i + ribo_size[2L] + 1L,
        ribo_size[2L] + ribo_size[3L] + ribo_size[1L] - cut5_i + cut3_i
      ))
    })
  }), nrow = 4)
  df_rpf <- df_rpf[df_rpf$qwidth %in% unique(convert_idx[4, ]), ]
  qwidth <- as.character(df_rpf$qwidth)
  psite_num_idx <- as.vector(
    sapply(
      split(convert_idx[1, ], convert_idx[4, ]), length
    )[qwidth]
  )
  psite_num_idx <- as.integer(psite_num_idx)
  valid_idx <- !is.na(psite_num_idx) & psite_num_idx > 0L
  if (!all(valid_idx)) {
    df_rpf <- df_rpf[valid_idx, , drop = FALSE]
    qwidth <- qwidth[valid_idx]
    psite_num_idx <- psite_num_idx[valid_idx]
  }
  psite_pos <- rep(
    df_rpf$pos, psite_num_idx
  ) +
    unlist(
      split(convert_idx[3, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )

  x_i <- rep(seq.int(nrow(df_rpf)), psite_num_idx)
  y_i <- sequence(nvec = psite_num_idx, from = 1L)
  candi_psite <- Matrix::sparseMatrix(
    i = x_i,
    j = y_i,
    x = psite_pos
  )

  candi_cut5 <- Matrix::sparseMatrix(
    i = x_i,
    j = y_i,
    x = unlist(
      split(convert_idx[1, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  )

  candi_cut3 <- Matrix::sparseMatrix(
    i = x_i,
    j = y_i,
    x = unlist(
      split(convert_idx[2, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  )

  ribo_num_0 <- pred_psite(
    df_rpf = df_rpf,
    candi_psite = candi_psite,
    candi_cut5 = candi_cut5,
    candi_cut3 = candi_cut3,
    tx_info = tx_info,
    ribo_size = ribo_size,
    par_0 = par_setup,
    iter_num = iter_num
  )

  ribo_num <- filter_ribosome(
    candi_psite = ribo_num_0,
    min_ribosome_num = min_ribosome_num
  )

  return(ribo_num)
}
#'
#'
#' Merge ORF predictions with genomic annotations
#'
#' @param orf.cds Result of [predict_orf()] applied to CDS candidates.
#' @param cds.psite P-site predictions for CDS ORFs, typically from
#'   [estimate_orf_psite()].
#' @param orf.ncrna Result of [predict_orf()] applied to ncRNA candidates.
#' @param ncrna.psite P-site predictions for ncRNA ORFs.
#' @param tx.info List containing CDS and ncRNA transcript annotations.
#' @param gene.info Gene-level metadata object.
#' @param exon.position Data frame describing exon coordinates.
#'
#' @return A unified list linking ORF calls to gene annotations and P-site
#'   profiles.
#'
#' @export
merge_orf_annotations <- function(
    orf.cds,
    cds.psite,
    orf.ncrna,
    ncrna.psite,
    tx.info,
    gene.info,
    exon.position) {
  convert_psite <- function(tx_info,
                            cds_psite,
                            ncrna_psite) {
    cds_tx_range <- IRanges::IRanges(
      start = tx_info$tx_info$mix_tx_pos$utr5_p5,
      end = tx_info$tx_info$mix_tx_pos$utr3_p3
    )

    cds_psite_range <- IRanges::IRanges(
      start = cds_psite$candi_pos,
      width = 1
    )

    cds_psite_hit <- IRanges::findOverlaps(
      cds_psite_range,
      cds_tx_range,
      type = "within"
    )

    cds_psite_pos <- split(
      BiocGenerics::start(cds_psite_range) -
        BiocGenerics::start(cds_tx_range)[cds_psite_hit@to] +
        1L,
      tx_info$tx_info$tx_lens$tx_name[cds_psite_hit@to]
    )

    cds_psite_num <- split(
      cds_psite$vec_pnum,
      tx_info$tx_info$tx_lens$tx_name[cds_psite_hit@to]
    )

    ncrna_tx_range <- IRanges::IRanges(
      start = tx_info$lncrna_info$mix_tx_pos$tx_p5,
      end = tx_info$lncrna_info$mix_tx_pos$tx_p3
    )

    ncrna_psite_range <- IRanges::IRanges(
      start = ncrna_psite$candi_pos,
      width = 1
    )

    ncrna_psite_hit <- IRanges::findOverlaps(
      ncrna_psite_range,
      ncrna_tx_range,
      type = "within"
    )

    ncrna_psite_pos <- split(
      BiocGenerics::start(ncrna_psite_range) -
        BiocGenerics::start(ncrna_tx_range)[ncrna_psite_hit@to] +
        1L,
      names(tx_info$lncrna_info$tx_seqs)[ncrna_psite_hit@to]
    )

    ncrna_psite_num <- split(
      ncrna_psite$vec_pnum,
      names(tx_info$lncrna_info$tx_seqs)[ncrna_psite_hit@to]
    )

    return(
      list(
        psite_position = c(cds_psite_pos, ncrna_psite_pos),
        psite_active = c(cds_psite_num, ncrna_psite_num)
      )
    )
  }

  find_trunc <- function(i_id,
                         tx_orf_cord,
                         tx_exon_width) {
    tx_i <- tx_orf_cord$tid[i_id]
    tx_trunc <- cumsum(tx_exon_width[[tx_i]])
    exon_i <- which(
      (
        (tx_orf_cord$tcoord[i_id] > tx_trunc) +
          (tx_orf_cord$tstop[i_id] > tx_trunc)) == 1L
    )
    if (length(exon_i) == 0L) {
      orf_cord <- c(tx_orf_cord$tcoord[i_id], tx_orf_cord$tstop[i_id])
    } else {
      orf_cord <- c(
        tx_orf_cord$tcoord[i_id],
        t(sapply(exon_i, function(e_i) {
          c(tx_trunc[e_i] + 1L, tx_trunc[e_i])
        })),
        tx_orf_cord$tstop[i_id]
      )
    }
    return(orf_cord)
  }

  tr_start_end <- function(rlocs,
                           strand) {
    starts <- sapply(
      seq_along(rlocs),
      function(i) {
        i_len <- length(rlocs[[i]])
        if (i_len == 2L) {
          if (strand[i] != "-") {
            return(rlocs[[i]][1])
          } else {
            return(rlocs[[i]][2])
          }
        } else {
          if (strand[i] != "-") {
            tmp_v <- vector(mode = "integer", length = i_len - 1L)
            tmp_v[] <- ";"
            idx <- seq.int(i_len / 2)
            tmp_v[idx * 2 - 1] <- rlocs[[i]][idx]
            return(paste0(tmp_v, collapse = ""))
          } else {
            tmp_v <- vector(mode = "integer", length = i_len - 1L)
            tmp_v[] <- ";"
            idx <- seq.int(i_len / 2)
            tmp_v[idx * 2 - 1] <- rlocs[[i]][i_len:(1 + i_len / 2)]
            return(paste0(tmp_v, collapse = ""))
          }
        }
      }
    )
    ends <- sapply(
      seq_along(rlocs),
      function(i) {
        i_len <- length(rlocs[[i]])
        if (i_len == 2L) {
          if (strand[i] != "-") {
            return(rlocs[[i]][2])
          } else {
            return(rlocs[[i]][1])
          }
        } else {
          if (strand[i] != "-") {
            tmp_v <- vector(mode = "integer", length = i_len - 1L)
            tmp_v[] <- ";"
            idx <- seq.int(i_len / 2)
            tmp_v[idx * 2 - 1] <- rlocs[[i]][idx + i_len / 2]
            return(paste0(tmp_v, collapse = ""))
          } else {
            tmp_v <- vector(mode = "integer", length = i_len - 1L)
            tmp_v[] <- ";"
            idx <- seq.int(i_len / 2)
            tmp_v[idx * 2 - 1] <- rlocs[[i]][(i_len / 2):1]
            return(paste0(tmp_v, collapse = ""))
          }
        }
      }
    )
    return(
      list(
        starts = starts,
        ends = ends
      )
    )
  }

  cord_convert <- function(orf_pred,
                           tx_info,
                           gene_info) {
    empty_orf_info <- function() {
      data.frame(
        chrm = character(0),
        strand = character(0),
        starts = character(0),
        ends = character(0),
        gene_id = character(0),
        transcript = character(0),
        tx_start = integer(0),
        tx_end = integer(0),
        orf_length = integer(0),
        start_codon = character(0),
        score = numeric(0),
        orf_sequence = character(0),
        stringsAsFactors = FALSE
      )
    }

    if (length(names(tx_info)) == 4) {
      tx_p5_map <- setNames(tx_info$mix_tx_pos$utr5_p5, tx_info$tx_lens$tx_name)
      tx_p3_map <- setNames(tx_info$mix_tx_pos$utr3_p3, tx_info$tx_lens$tx_name)
      tx_range <- IRanges::IRanges(
        start = tx_info$mix_tx_pos$utr5_p5,
        end = tx_info$mix_tx_pos$utr3_p3
      )

      tx_hit <- IRanges::findOverlaps(
        orf_pred$pred_orf,
        tx_range,
        type = "within"
      )

      tx_orf_cord <- data.frame(
        tid = tx_info$tx_lens$tx_name[tx_hit@to],
        tcoord = BiocGenerics::start(orf_pred$pred_orf) -
          BiocGenerics::start(tx_range)[tx_hit@to] + 1L,
        tstop = BiocGenerics::end(orf_pred$pred_orf) -
          BiocGenerics::start(tx_range)[tx_hit@to] + 1L
      )
    } else {
      tx_p5_map <- setNames(tx_info$mix_tx_pos$tx_p5, names(tx_info$tx_seqs))
      tx_p3_map <- setNames(tx_info$mix_tx_pos$tx_p3, names(tx_info$tx_seqs))
      tx_range <- IRanges::IRanges(
        start = tx_info$mix_tx_pos$tx_p5,
        end = tx_info$mix_tx_pos$tx_p3
      )

      tx_hit <- IRanges::findOverlaps(
        orf_pred$pred_orf,
        tx_range,
        type = "within"
      )

      tx_orf_cord <- data.frame(
        tid = names(tx_info$tx_seqs)[tx_hit@to],
        tcoord = BiocGenerics::start(orf_pred$pred_orf) -
          BiocGenerics::start(tx_range)[tx_hit@to] + 1L,
        tstop = BiocGenerics::end(orf_pred$pred_orf) -
          BiocGenerics::start(tx_range)[tx_hit@to] + 1L
      )
    }
    if (nrow(tx_orf_cord) == 0L) {
      return(empty_orf_info())
    }

    score_vec <- orf_pred$predict_score[tx_hit@from]

    tx_map_len <- vapply(
      tx_orf_cord$tid,
      function(tid) {
        w <- gene_info$tx_exon_width[[tid]]
        if (is.null(w)) return(NA_real_)
        sum(w)
      },
      numeric(1)
    )
    keep_mappable <- !is.na(tx_map_len) &
      tx_orf_cord$tcoord >= 1L &
      tx_orf_cord$tstop >= 1L &
      tx_orf_cord$tcoord <= tx_orf_cord$tstop &
      tx_orf_cord$tstop <= tx_map_len
    if (!all(keep_mappable)) {
      warning(
        sprintf(
          "cord_convert(): %d ORF(s) are outside exon-mappable transcript range; keep for biotype classification with NA genomic coordinates.",
          sum(!keep_mappable)
        )
      )
    }

    orf_info <- data.frame(
      chrm = GenomeInfoDb::seqnames(gene_info$tx_info)[
        match(
          tx_orf_cord$tid,
          gene_info$tx_info@elementMetadata@listData$tx_name
        )
      ],
      strand = BiocGenerics::strand(gene_info$tx_info)[
        match(
          tx_orf_cord$tid,
          gene_info$tx_info@elementMetadata@listData$tx_name
        )
      ]
    )
    orf_info$starts <- NA_character_
    orf_info$ends <- NA_character_

    if (any(keep_mappable)) {
      tx_orf_mappable <- tx_orf_cord[keep_mappable, , drop = FALSE]
      tlocs <- lapply(
        seq_along(tx_orf_mappable$tid), function(i) {
          return(
            find_trunc(
              i_id = i,
              tx_orf_cord = tx_orf_mappable,
              tx_exon_width = gene_info$tx_exon_width
            )
          )
        }
      )

      rlocs <- GenomicFeatures::transcriptLocs2refLocs(
        tlocs,
        BiocGenerics::start(gene_info$tx_exon)[tx_orf_mappable$tid],
        BiocGenerics::end(gene_info$tx_exon)[tx_orf_mappable$tid],
        gene_info$tx_strand[tx_orf_mappable$tid],
        decreasing.rank.on.minus.strand = FALSE
      )

      tx_start_end <- tr_start_end(
        rlocs = rlocs,
        strand = orf_info$strand[keep_mappable]
      )
      orf_info$starts[keep_mappable] <- tx_start_end$start
      orf_info$ends[keep_mappable] <- tx_start_end$ends
    }

    orf_info$transcript <- tx_orf_cord$tid
    orf_info$tx_start <- tx_orf_cord$tcoord
    orf_info$tx_end <- tx_orf_cord$tstop
    orf_info$gene_id <- gene_info$tx_cds$gene_id[
      match(orf_info$transcript, gene_info$tx_cds$tx_name)
    ]
    orf_info$orf_length <- as.integer(
      (orf_info$tx_end - orf_info$tx_start + 1) / 3 - 1
    )

    tx_seq_width <- Biostrings::width(gene_info$tx_seqs[orf_info$transcript])
    tx_len_mix <- as.integer(tx_p3_map[orf_info$transcript] - tx_p5_map[orf_info$transcript] + 1L)
    valid_start_codon <- !is.na(tx_seq_width) &
      orf_info$tx_start >= 1L &
      (orf_info$tx_start + 2L) <= tx_seq_width
    start_codon <- rep(NA_character_, nrow(orf_info))
    if (any(valid_start_codon)) {
      start_codon[valid_start_codon] <- as.character(Biostrings::subseq(
        gene_info$tx_seqs[orf_info$transcript[valid_start_codon]],
        orf_info$tx_start[valid_start_codon],
        width = 3
      ))
    }
    valid_start_codon_mix <- !valid_start_codon &
      !is.na(tx_len_mix) &
      orf_info$tx_start >= 1L &
      (orf_info$tx_start + 2L) <= tx_len_mix
    if (any(valid_start_codon_mix)) {
      abs_start <- tx_p5_map[orf_info$transcript[valid_start_codon_mix]] +
        orf_info$tx_start[valid_start_codon_mix] - 1L
      start_codon[valid_start_codon_mix] <- vapply(
        abs_start,
        function(st) as.character(Biostrings::subseq(tx_info$mix_tx, start = st, width = 3L)),
        character(1)
      )
    }
    orf_info$start_codon <- start_codon
    orf_info$score <- score_vec
    valid_orf_seq <- !is.na(tx_seq_width) &
      tx_orf_cord$tcoord >= 1L &
      tx_orf_cord$tstop >= tx_orf_cord$tcoord &
      tx_orf_cord$tstop <= tx_seq_width
    orf_sequence <- rep(NA_character_, nrow(orf_info))
    if (any(valid_orf_seq)) {
      orf_sequence[valid_orf_seq] <- as.character(
        Biostrings::translate(
          Biostrings::subseq(
            gene_info$tx_seqs[tx_orf_cord$tid[valid_orf_seq]],
            tx_orf_cord$tcoord[valid_orf_seq],
            end = tx_orf_cord$tstop[valid_orf_seq]
          )
        )
      )
    }
    valid_orf_seq_mix <- !valid_orf_seq &
      !is.na(tx_len_mix) &
      tx_orf_cord$tcoord >= 1L &
      tx_orf_cord$tstop >= tx_orf_cord$tcoord &
      tx_orf_cord$tstop <= tx_len_mix
    if (any(valid_orf_seq_mix)) {
      abs_start <- tx_p5_map[tx_orf_cord$tid[valid_orf_seq_mix]] +
        tx_orf_cord$tcoord[valid_orf_seq_mix] - 1L
      abs_end <- tx_p5_map[tx_orf_cord$tid[valid_orf_seq_mix]] +
        tx_orf_cord$tstop[valid_orf_seq_mix] - 1L
      mix_nt <- Biostrings::DNAStringSet(vapply(
        seq_along(abs_start),
        function(i) as.character(Biostrings::subseq(tx_info$mix_tx, start = abs_start[i], end = abs_end[i])),
        character(1)
      ))
      orf_sequence[valid_orf_seq_mix] <- as.character(Biostrings::translate(mix_nt))
    }
    orf_info$orf_sequence <- orf_sequence
    orf_info <- orf_info[
      ,
      c(
        "chrm", "strand", "starts", "ends",
        "gene_id", "transcript", "tx_start", "tx_end",
        "orf_length", "start_codon", "score", "orf_sequence"
      )
    ]
    return(orf_info)
  }

  orf_pos_tx <- function(tx_cds, orf_info) {
    # no cds
    idx <- match(orf_info$transcript, tx_cds$tx_name)
    cds_idx <- which(tx_cds$cds_len[idx] != 0L)
    tmp_info <- orf_info[cds_idx, ]
    tmp_cds <- tx_cds[match(tmp_info$transcript, tx_cds$tx_name), ]
    idx11 <- ifelse(tmp_info$tx_start < (1L + tmp_cds$utr5_len),
      -1,
      ifelse(tmp_info$tx_start == (1L + tmp_cds$utr5_len),
        0, 1
      )
    )
    idx12 <- ifelse(tmp_info$tx_start < (tmp_cds$cds_len + tmp_cds$utr5_len),
      -1, 1
    )
    idx21 <- ifelse(tmp_info$tx_end < (1L + tmp_cds$utr5_len),
      -1, 1
    )
    idx22 <- ifelse(tmp_info$tx_end < (tmp_cds$cds_len + tmp_cds$utr5_len),
      -1,
      ifelse(tmp_info$tx_end == (tmp_cds$cds_len + tmp_cds$utr5_len),
        0, 1
      )
    )
    cl_idx <- paste0(idx11, idx12, idx21, idx22)
    orftype <- vector(length = nrow(tmp_info), mode = "character")
    orftype[] <- "Variant"
    orftype[cl_idx == "-1-1-1-1"] <- "uORF"
    orftype[cl_idx == "-1-11-1"] <- "uoORF"
    orftype[cl_idx == "-1-110"] <- "Ext"
    orftype[cl_idx == "0-110"] <- "CDS"
    orftype[cl_idx == "1-11-1"] <- "intORF"
    orftype[cl_idx == "1-110"] <- "Trunc"
    orftype[cl_idx == "1-111"] <- "doORF"
    orftype[cl_idx == "1111"] <- "dORF"
    orftype2 <- vector(length = nrow(orf_info), mode = "character")
    orftype2[] <- "ncRNA"
    orftype2[cds_idx] <- orftype
    return(orftype2)
  }

  get_orf_exon <- function(orf_info,
                           gene_info,
                           cds_exon) {
    orf_cord_group <- lapply(
      c("+", "-"),
      function(strand_orf) {
        tmp_orf_info <- orf_info[orf_info$strand == strand_orf, ]
        orf_info_not3 <- tmp_orf_info[
          tmp_orf_info$gene_id %in%
            gene_info$tx_cds$gene_id[
              match(
                cds_exon$df_not3$cds_exon$tx_id,
                gene_info$tx_cds$tx_name
              )
            ],
        ]

        orf_info_is3 <- tmp_orf_info[
          tmp_orf_info$gene_id %in%
            gene_info$tx_cds$gene_id[
              match(
                cds_exon$df_is3$cds_exon$tx_id,
                gene_info$tx_cds$tx_name
              )
            ],
        ]
        orf_info_sub <- list(
          not3 = orf_info_not3,
          is3 = orf_info_is3
        )
        orf_cord <- lapply(
          orf_info_sub,
          function(orf_info_i) {
            orf_id <- orf_info_i$orf_id

            exon_left_text <- stringr::str_split(
              orf_info_i$starts, ";"
            )
            exon_right_text <- stringr::str_split(
              orf_info_i$ends, ";"
            )
            exon_left <- as.numeric(unlist(exon_left_text))
            exon_right <- as.numeric(unlist(exon_right_text))

            exon_num <- sapply(exon_left_text, length)
            exon_length <- split(
              exon_right - exon_left + 1L,
              rep(orf_id, exon_num)
            )

            exon_frame <- lapply(
              exon_length,
              function(x) {
                cumsum(x) %% 3L
              }
            )
            exon_frame[orf_info_i$strand == "-"] <- lapply(
              exon_length[orf_info_i$strand == "-"],
              function(x) {
                rev(cumsum(rev(x)) %% 3L)
              }
            )

            orf_exon <- list(
              orf_exon = data.frame(
                orf_id = orf_id,
                exon_num = exon_num
              ),
              exon_position = data.frame(
                exon_left = exon_left,
                exon_right = exon_right,
                exon_frame = as.numeric(unlist(exon_frame, use.names = FALSE))
              )
            )
            return(orf_exon)
          }
        )

        return(orf_cord)
      }
    )
    names(orf_cord_group) <- c("plus", "minus")
    return(orf_cord_group)
  }

  overlap_known_cds <- function(orf_info,
                                gene_info,
                                orf_exon,
                                cds_exon) {
    num_cds <- lapply(
      cds_exon,
      function(cds_exon_i) {
        sapply(
          split(
            cds_exon_i$cds_exon$exon_num,
            gene_info$tx_cds$gene_id[
              match(
                cds_exon_i$cds_exon$tx_id,
                gene_info$tx_cds$tx_name
              )
            ]
          ), sum
        )
      }
    )
    names(num_cds) <- c("not3", "is3")

    orf_ternimal <- lapply(
      orf_exon,
      function(orf_exon_i) {
        orf_cord <- lapply(
          1:2,
          function(terminal_i) {
            tmp_cord <- lapply(
              1:2,
              function(group_i) {
                rep(
                  orf_exon_i[[group_i]]$exon_position[, terminal_i],
                  rep(
                    num_cds[[group_i]][
                      orf_info$gene_id[
                        match(
                          orf_exon_i[[group_i]]$orf_exon$orf_id,
                          orf_info$orf_id
                        )
                      ]
                    ],
                    orf_exon_i[[group_i]]$orf_exon$exon_num
                  )
                )
              }
            )
            names(tmp_cord) <- c("not3", "is3")
            return(tmp_cord)
          }
        )
        names(orf_cord) <- c("left", "right")
        return(orf_cord)
      }
    )

    gene_exon_terminal <- lapply(
      1:2,
      function(not3_i) {
        tmp_i <- lapply(
          1:2,
          function(terminal_i) {
            split(
              cds_exon[[not3_i]]$exon_position[, terminal_i],
              rep(
                gene_info$tx_cds$gene_id[
                  match(
                    cds_exon[[not3_i]]$cds_exon$tx_id,
                    gene_info$tx_cds$tx_name
                  )
                ],
                cds_exon[[not3_i]]$cds_exon$exon_num
              )
            )
          }
        )
        names(tmp_i) <- c("left", "right")
        return(tmp_i)
      }
    )
    names(gene_exon_terminal) <- c("not3", "is3")

    # frame
    gene_exon_frame <- split(
      cds_exon$df_is3$exon_position$exon_frame,
      rep(
        gene_info$tx_cds$gene_id[
          match(
            cds_exon$df_is3$cds_exon$tx_id,
            gene_info$tx_cds$tx_name
          )
        ],
        cds_exon$df_is3$cds_exon$exon_num
      )
    )

    cds_orf_cord <- function(ref_cds,
                             plus_minus,
                             is_not3,
                             cord1,
                             cord2) {
      tmp_i <- orf_exon[[plus_minus]][[is_not3]]$orf_exon
      if (ref_cds) {
        return(
          unlist(
            gene_exon_terminal[[is_not3]][[cord1]][
              orf_info$gene_id[
                rep(
                  match(tmp_i$orf_id, orf_info$orf_id),
                  tmp_i$exon_num
                )
              ]
            ],
            use.names = FALSE
          ) -
            orf_ternimal[[plus_minus]][[cord2]][[is_not3]]
        )
      } else {
        return(
          orf_ternimal[[plus_minus]][[cord1]][[is_not3]] -
            unlist(
              gene_exon_terminal[[is_not3]][[cord2]][
                orf_info$gene_id[
                  rep(
                    match(tmp_i$orf_id, orf_info$orf_id),
                    tmp_i$exon_num
                  )
                ]
              ],
              use.names = FALSE
            )
        )
      }
    }
    cds_num <- lapply(
      1:2,
      function(plus_minus) {
        tmp_num <- lapply(
          1:2,
          function(is_not3) {
            tmp_i <- orf_exon[[plus_minus]][[is_not3]]$orf_exon
            return(sapply(
              gene_exon_terminal[[is_not3]]$left[
                orf_info$gene_id[
                  rep(
                    match(tmp_i$orf_id, orf_info$orf_id),
                    tmp_i$exon_num
                  )
                ]
              ], length
            ))
          }
        )
        names(tmp_num) <- c("not3", "is3")
        return(tmp_num)
      }
    )
    names(cds_num) <- c("plus", "minus")

    plus_not3 <- (cds_orf_cord(TRUE, "plus", "not3", "right", "left") > 0) &
      (cds_orf_cord(FALSE, "plus", "not3", "right", "left") > 0)
    minus_not3 <- (cds_orf_cord(TRUE, "minus", "not3", "right", "left") > 0) &
      (cds_orf_cord(FALSE, "minus", "not3", "right", "left") > 0)
    plus_is3 <- (cds_orf_cord(TRUE, "plus", "is3", "right", "left") > 0) &
      (cds_orf_cord(FALSE, "plus", "is3", "right", "left") > 0)
    minus_is3 <- (cds_orf_cord(TRUE, "minus", "is3", "right", "left") > 0) &
      (cds_orf_cord(FALSE, "minus", "is3", "right", "left") > 0)

    plus_minus <- c("plus", "minus")
    right_left <- c("right", "left")

    is_frame <- lapply(
      1:2,
      function(i) {
        tmp_i <- orf_exon[[plus_minus[i]]]$is3$orf_exon

        frame_distance <- (orf_ternimal[[plus_minus[i]]][[right_left[i]]]$is3 -
          unlist(
            gene_exon_terminal$is3[[right_left[i]]][
              orf_info$gene_id[
                rep(
                  match(tmp_i$orf_id, orf_info$orf_id),
                  tmp_i$exon_num
                )
              ]
            ],
            use.names = FALSE
          )) %% 3L

        frame_orf <- rep(
          orf_exon[[plus_minus[i]]]$is3$exon_position[, 3],
          rep(
            num_cds$is3[
              orf_info$gene_id[
                match(
                  orf_exon[[plus_minus[i]]]$is3$orf_exon$orf_id,
                  orf_info$orf_id
                )
              ]
            ],
            orf_exon[[plus_minus[i]]]$is3$orf_exon$exon_num
          )
        )

        frame_cds <- unlist(
          gene_exon_frame[
            orf_info$gene_id[
              rep(
                match(tmp_i$orf_id, orf_info$orf_id),
                tmp_i$exon_num
              )
            ]
          ],
          use.names = FALSE
        )

        return((frame_orf - frame_distance) %% 3L == frame_cds)
      }
    )
    names(is_frame) <- c("plus", "minus")

    plus_not3_overlap <- data.frame(
      orf_id = orf_exon$plus$not3$orf_exon$orf_id,
      plus_not3 = sapply(
        split(
          sapply(
            split(
              plus_not3,
              rep(
                seq_along(cds_num$plus$not3),
                cds_num$plus$not3
              )
            ), sum
          ), rep(
            seq_along(orf_exon$plus$not3$orf_exon$exon_num),
            orf_exon$plus$not3$orf_exon$exon_num
          )
        ), sum
      ) > 0L
    )

    minus_not3_overlap <- data.frame(
      orf_id = orf_exon$minus$not3$orf_exon$orf_id,
      minus_not3 = sapply(
        split(
          sapply(
            split(
              minus_not3,
              rep(
                seq_along(cds_num$minus$not3),
                cds_num$minus$not3
              )
            ), sum
          ), rep(
            seq_along(orf_exon$minus$not3$orf_exon$exon_num),
            orf_exon$minus$not3$orf_exon$exon_num
          )
        ), sum
      ) > 0L
    )

    plus_is3_overlap <- data.frame(
      orf_id = orf_exon$plus$is3$orf_exon$orf_id,
      plus_is3 = sapply(
        split(
          sapply(
            split(
              plus_is3,
              rep(
                seq_along(cds_num$plus$is3),
                cds_num$plus$is3
              )
            ), sum
          ), rep(
            seq_along(orf_exon$plus$is3$orf_exon$exon_num),
            orf_exon$plus$is3$orf_exon$exon_num
          )
        ), sum
      ) > 0L
    )

    minus_is3_overlap <- data.frame(
      orf_id = orf_exon$minus$is3$orf_exon$orf_id,
      minus_is3 = sapply(
        split(
          sapply(
            split(
              minus_is3,
              rep(
                seq_along(cds_num$minus$is3),
                cds_num$minus$is3
              )
            ), sum
          ), rep(
            seq_along(orf_exon$minus$is3$orf_exon$exon_num),
            orf_exon$minus$is3$orf_exon$exon_num
          )
        ), sum
      ) > 0L
    )

    annotate_overlap <- vector(length = nrow(orf_info))
    annotate_overlap[
      plus_not3_overlap$orf_id[plus_not3_overlap$plus_not3]
    ] <- TRUE
    annotate_overlap[
      plus_is3_overlap$orf_id[plus_is3_overlap$plus_is3]
    ] <- TRUE
    annotate_overlap[
      minus_not3_overlap$orf_id[minus_not3_overlap$minus_not3]
    ] <- TRUE
    annotate_overlap[
      minus_is3_overlap$orf_id[minus_is3_overlap$minus_is3]
    ] <- TRUE

    plus_is3_inframe <- data.frame(
      orf_id = orf_exon$plus$is3$orf_exon$orf_id,
      plus_is3 = sapply(
        split(
          sapply(
            split(
              plus_is3 * 10L + is_frame$plus,
              rep(
                seq_along(cds_num$plus$is3),
                cds_num$plus$is3
              )
            ), max
          ), rep(
            seq_along(orf_exon$plus$is3$orf_exon$exon_num),
            orf_exon$plus$is3$orf_exon$exon_num
          )
        ), max
      ) > 10L
    )

    minus_is3_overlap <- data.frame(
      orf_id = orf_exon$minus$is3$orf_exon$orf_id,
      minus_is3 = sapply(
        split(
          sapply(
            split(
              minus_is3 * 10L + is_frame$minus,
              rep(
                seq_along(cds_num$minus$is3),
                cds_num$minus$is3
              )
            ), max
          ), rep(
            seq_along(orf_exon$minus$is3$orf_exon$exon_num),
            orf_exon$minus$is3$orf_exon$exon_num
          )
        ), max
      ) > 10L
    )

    annotate_inframe <- vector(length = nrow(orf_info))
    annotate_inframe[
      plus_is3_inframe$orf_id[plus_is3_inframe$plus_is3]
    ] <- TRUE
    annotate_inframe[
      minus_is3_overlap$orf_id[minus_is3_overlap$minus_is3]
    ] <- TRUE

    return(
      list(
        annotate_overlap = annotate_overlap,
        annotate_inframe = annotate_inframe
      )
    )
  }

  ribobase_orf <- cord_convert(
    orf_pred = orf.cds,
    tx_info = tx.info$tx_info,
    gene_info = gene.info
  )

  ribobase_orf$orf_biotype <- orf_pos_tx(
    tx_cds = gene.info$tx_cds,
    orf_info = ribobase_orf
  )

  ribobase_orf$orf_id <- seq.int(nrow(ribobase_orf))
  ribobase_orf$is_overlap <- FALSE
  ribobase_orf$is_inframe <- FALSE

  mappable_idx <- which(!is.na(ribobase_orf$starts) & !is.na(ribobase_orf$ends))
  if (length(mappable_idx) > 0L) {
    ribobase_orf_mappable <- ribobase_orf[mappable_idx, , drop = FALSE]
    ribobase_orf_mappable$orf_id <- seq_len(nrow(ribobase_orf_mappable))

    overlap_known <- tryCatch(
      {
        orf_exon <- get_orf_exon(
          orf_info = ribobase_orf_mappable,
          gene_info = gene.info,
          cds_exon = exon.position
        )
        overlap_known_cds(
          orf_info = ribobase_orf_mappable,
          gene_info = gene.info,
          orf_exon = orf_exon,
          cds_exon = exon.position
        )
      },
      error = function(e) {
        warning(
          "Failed to compute overlap/inframe for mapped CDS ORFs; keeping defaults FALSE. Reason: ",
          conditionMessage(e)
        )
        NULL
      }
    )
    if (!is.null(overlap_known)) {
      ribobase_orf$is_overlap[mappable_idx] <- overlap_known$annotate_overlap
      ribobase_orf$is_inframe[mappable_idx] <- overlap_known$annotate_inframe
    }
  }

  if (length(orf.ncrna$pred_orf) == 0L) {
    ncrna_orf <- ribobase_orf[0, , drop = FALSE]
  } else {
    ncrna_orf <- cord_convert(
      orf_pred = orf.ncrna,
      tx_info = tx.info$lncrna_info,
      gene_info = gene.info
    )

    ncrna_orf$orf_biotype <- orf_pos_tx(
      tx_cds = gene.info$tx_cds,
      orf_info = ncrna_orf
    )

    ncrna_orf$orf_id <- seq.int(nrow(ncrna_orf))

    ncrna_orf$is_overlap <- FALSE
    ncrna_orf$is_inframe <- FALSE
  }

  predict_orf <- rbind(ribobase_orf, ncrna_orf)
  predict_orf$orf_id <- seq.int(nrow(predict_orf))

  psite_info <- convert_psite(
    tx_info = tx.info,
    cds_psite = cds.psite,
    ncrna_psite = ncrna.psite
  )
  return(
    list(
      ORF = predict_orf,
      Psite = psite_info
    )
  )
}
#'
#'
#' Run the end-to-end RiboBC analysis pipeline
#'
#' @param transcript.file Path to an `.RData` file produced by the preprocessing
#'   helpers that contains transcript annotations and candidate ORFs.
#' @param bam_dir Directory produced by [prepare_bam()], containing
#'   `*_cds_txsorted.bam(.bai)` and `*_lncrna_txsorted.bam(.bai)`.
#' @param RNase RNase digestion paradigm. Default is `"rnase-i"`;
#'   alternatives are `"mnase"` and `"p1"`.
#' @param min_lnc_reads_for_prediction Minimum ncRNA-aligned read count required
#'   to run lncORF prediction for each sample. If lower than this threshold,
#'   lncORF prediction is skipped with a warning while CDS prediction continues.
#' @param export_downstream Logical; if `TRUE` (default), automatically run
#'   `export_riboba_downstream()` for each sample at the end of pipeline.
#' @param downstream_output_dir Directory used by automatic downstream export.
#'   Defaults to `file.path(getwd(), "results")`.
#' @param downstream_file_prefix File prefix used by automatic downstream export.
#'
#' @return A list of per-sample results including bias parameters, trained
#'   models, and ORF annotations. Automatic export paths are attached as
#'   attribute `downstream_exports` when `export_downstream = TRUE`.
#'
#' @export
run_riboba_pipeline <- function(
    transcript.file,
    bam_dir,
    RNase = "rnase-i",
    min_lnc_reads_for_prediction = 100L,
    export_downstream = TRUE,
    downstream_output_dir = file.path(getwd(), "results"),
    downstream_file_prefix = "") {
  get_prepare_bam_inputs <- function(bam_dir) {
    bam_cds <- sort(list.files(bam_dir, pattern = "_cds_txsorted\\.bam$", full.names = TRUE))
    bam_ncrna <- sort(list.files(bam_dir, pattern = "_lncrna_txsorted\\.bam$", full.names = TRUE))
    bai_cds <- paste0(bam_cds, ".bai")
    bai_ncrna <- paste0(bam_ncrna, ".bai")

    if (length(bam_cds) == 0 || length(bam_ncrna) == 0) {
      stop("No CDS/lncRNA BAM files found in bam_dir: ", bam_dir)
    }
    if (length(bam_cds) != length(bam_ncrna)) {
      stop("CDS and lncRNA BAM counts do not match in bam_dir: ", bam_dir)
    }
    if (!all(file.exists(c(bai_cds, bai_ncrna)))) {
      stop("Missing BAM index (.bai) files in bam_dir: ", bam_dir)
    }

    cds_sample <- sub("_cds_txsorted\\.bam$", "", basename(bam_cds))
    ncrna_sample <- sub("_lncrna_txsorted\\.bam$", "", basename(bam_ncrna))
    if (!identical(cds_sample, ncrna_sample)) {
      stop("Sample names do not match between CDS and lncRNA BAM files in bam_dir: ", bam_dir)
    }

    list(
      bam_cds = bam_cds,
      bai_cds = bai_cds,
      bam_ncrna = bam_ncrna,
      bai_ncrna = bai_ncrna,
      sample_names = cds_sample
    )
  }

  validate_rnase <- function(rnase, n_sample) {
    allowed <- c("rnase-i", "mnase", "p1")
    rnase <- tolower(rnase)
    if (!all(rnase %in% allowed)) {
      stop("RNase must be one of: ", paste(allowed, collapse = ", "))
    }
    if (length(rnase) == 1L) {
      rnase <- rep(rnase, n_sample)
    } else if (length(rnase) != n_sample) {
      stop("RNase length must be 1 or equal to number of samples: ", n_sample)
    }
    rnase
  }

  bam_inputs <- get_prepare_bam_inputs(bam_dir = bam_dir)
  RNase <- validate_rnase(rnase = RNase, n_sample = length(bam_inputs$bam_cds))
  rrna_log_path <- file.path(dirname(bam_dir), "no_rrna", "rrna.log")
  cds_log_path <- file.path(dirname(bam_dir), "map_cds", "cds.log")
  lncrna_log_path <- file.path(dirname(bam_dir), "map_lncrna", "lncrna.log")
  message(
    sprintf(
      "Start RiboBA pipeline for %d sample(s): %s",
      length(bam_inputs$sample_names),
      paste(bam_inputs$sample_names, collapse = ", ")
    )
  )

  message("Loading transcript information ... ", appendLF = FALSE)
  load(transcript.file)
  message("Done")

  message("Loading mapping information ... ", appendLF = FALSE)
  message("CDS ... ", appendLF = FALSE)
  reads_bam <- read_bam_alignments(
    bam_path = bam_inputs$bam_cds,
    bai_path = bam_inputs$bai_cds,
    tx_info = tx_info_lst$tx_info,
    only_uni_map = TRUE,
    add_ratio = 0.1
  )

  message("lncRNA ... ", appendLF = FALSE)
  reads_bam_ncrna <- read_bam_alignments(
    bam_path = bam_inputs$bam_ncrna,
    bai_path = bam_inputs$bai_ncrna,
    tx_info = tx_info_lst$lncrna_info,
    only_uni_map = TRUE,
    add_ratio = 0.1
  )
  message("Done")

  orf_result <- list()
  downstream_exports <- list()
  for (i in seq_along(reads_bam)) {
    sample_name <- bam_inputs$sample_names[i]
    RNase_i = RNase[i]
    message(sprintf("[%d/%d] Processing sample: %s (RNase=%s)", i, length(reads_bam), sample_name, RNase_i))
    if (RNase_i == "mnase") {
      rnase.bias <- TRUE
    } else {
      rnase.bias <- FALSE
    }
    message(
      paste0("Infer bias parameters for sample ", sample_name, " ... "),
      appendLF = FALSE
    )
    par_rpf <- learn_bias_parameters(
      tx_info = tx_info_lst$tx_info,
      lst_bam = reads_bam[[i]],
      rnase_bias = rnase.bias,
      RNase = RNase_i,
      ligat_par = TRUE
    )

    ncrna_rpf <- prepare_ncrna_rpf(
      lst_rpf = reads_bam_ncrna[[i]],
      tx_info = tx_info_lst$lncrna_info,
      par_learn = par_rpf
    )
    n_lnc_reads <- if (!is.null(reads_bam_ncrna[[i]]$no_add)) {
      nrow(reads_bam_ncrna[[i]]$no_add)
    } else {
      0L
    }
    skip_lnc_prediction <- is.na(n_lnc_reads) || n_lnc_reads < as.integer(min_lnc_reads_for_prediction)
    if (skip_lnc_prediction) {
      warning(
        sprintf(
          "Sample %s: lncRNA reads (%s) < min_lnc_reads_for_prediction (%s); skip lncORF prediction.",
          sample_name, n_lnc_reads, as.integer(min_lnc_reads_for_prediction)
        )
      )
    }
    message("Done")

    message(paste0("Predict ORFs for sample ", sample_name, " ... "),
      appendLF = FALSE
    )
    message("CDS ... ", appendLF = FALSE)
    res <- riboba_all(
      candidate_orf = candidate_orf,
      tx_info_lst = tx_info_lst,
      par_rpf = par_rpf,
      candidate_translat_reg = candidate_translat_reg,
      ncrna_rpf = ncrna_rpf,
      candidate_lnc_reg = candidate_lnc_reg,
      rnase.bias = rnase.bias,
      sample_name = sample_name,
      skip_lnc_prediction = skip_lnc_prediction
    )
    message("Done")

    message("Predict ribosome distribution ... ", appendLF = FALSE)
    orf_pred_p <- estimate_orf_psite(
      orf_candi = res$orf_pred$pred_orf,
      par_setup = par_rpf$par_0,
      df_rpf = par_rpf$lst_rpf$rpf_overlap,
      tx_info = tx_info_lst$tx_info,
      iter_num = 5
    )

    if (skip_lnc_prediction || length(res$lnc_orf_pred$pred_orf) == 0L || nrow(ncrna_rpf$lst_rpf$rpf_overlap) == 0L) {
      lncorf_pred_p <- list(candi_pos = integer(0), vec_pnum = numeric(0))
    } else {
      lncorf_pred_p <- estimate_orf_psite(
        orf_candi = res$lnc_orf_pred$pred_orf,
        par_setup = ncrna_rpf$par_0,
        df_rpf = ncrna_rpf$lst_rpf$rpf_overlap,
        tx_info = tx_info_lst$lncrna_info,
        iter_num = 5
      )
    }
    message("Done")

    message("Positional relationship ncORFs and CDS ... ", appendLF = FALSE)
    orf_info <- merge_orf_annotations(
      orf.cds = res$orf_pred,
      cds.psite = orf_pred_p,
      ncrna.psite = lncorf_pred_p,
      orf.ncrna = res$lnc_orf_pred,
      tx.info = tx_info_lst,
      gene.info = gene_info,
      exon.position = cds_exon
    )
    message("Done")

    orf_result[[i]] <- list(
      bias_info = par_rpf$lig_par,
      predict_mode = res,
      orf_info = orf_info
    )
    if (isTRUE(export_downstream)) {
      message("Export downstream files ... ", appendLF = FALSE)
      export_i <- tryCatch(
        export_riboba_downstream(
          riboba_result = orf_result[[i]],
          output_dir = downstream_output_dir,
          sample_name = sample_name,
          file_prefix = downstream_file_prefix,
          par_result = par_rpf,
          rrna_log = rrna_log_path,
          cds_log = cds_log_path,
          lncrna_log = lncrna_log_path
        ),
        error = function(e) {
          warning(
            sprintf(
              "Sample %s: automatic export_riboba_downstream failed: %s",
              sample_name, conditionMessage(e)
            )
          )
          NULL
        }
      )
      downstream_exports[[sample_name]] <- export_i
      message("Done")
    }
    message(
      sprintf(
        "Sample %s finished. Predicted ORFs: %d",
        sample_name,
        nrow(orf_info$ORF)
      )
    )

  }

  names(orf_result) <- bam_inputs$sample_names
  if (isTRUE(export_downstream)) {
    attr(orf_result, "downstream_exports") <- downstream_exports
  }
  message("RiboBA pipeline finished.")
  return(orf_result)
}

#' @export
infer_par <- function(...) {
  .Deprecated("learn_bias_parameters")
  learn_bias_parameters(...)
}

#' @export
input_bam <- function(...) {
  .Deprecated("read_bam_alignments")
  read_bam_alignments(...)
}

#' @export
prep_rpf_ncrna <- function(...) {
  .Deprecated("prepare_ncrna_rpf")
  prepare_ncrna_rpf(...)
}

#' @export
pred_p_candi <- function(...) {
  .Deprecated("estimate_orf_psite")
  estimate_orf_psite(...)
}

#' @export
cord_unify <- function(...) {
  .Deprecated("merge_orf_annotations")
  merge_orf_annotations(...)
}

#' @export
ribo_base <- function(...) {
  .Deprecated("run_riboba_pipeline")
  run_riboba_pipeline(...)
}

riboba_all <- function(
    candidate_orf,
    tx_info_lst,
    par_rpf,
    candidate_translat_reg,
    ncrna_rpf,
    candidate_lnc_reg,
    rnase.bias,
    sample_name = NULL,
    skip_lnc_prediction = FALSE) {
      get_cleavage_seq <- function(seqs,
                             p_site,
                             ribo_size) {
  seqs <- Biostrings::DNAStringSet(x = seqs)[
    rep(1L, length(p_site))
  ]

  up_seq <- Biostrings::subseq(
    x = seqs,
    start = p_site - sum(ribo_size[1:2]),
    width = ribo_size[1]
  )

  dn_seq <- Biostrings::subseq(
    x = seqs,
    start = p_site + ribo_size[3],
    width = ribo_size[4]
  )

  return(
    list(
      up_seq = up_seq,
      dn_seq = dn_seq
    )
  )
}

maintain_prob <- function(cut_seq,
                          prod_hd,
                          bias) {
  prob_mc <- stringr::str_split(
    string = as.vector(cut_seq),
    pattern = "",
    simplify = TRUE
  )

  prob_mp <- matrix(
    data = bias[prob_mc],
    nrow = nrow(prob_mc)
  )

  prob_mp <- prob_mp^(rep(1L, nrow(prob_mp)) %o% prod_hd)

  return(prob_mp)
}

get_cleavage_prob <- function(seqs,
                              bias,
                              prob_hd5,
                              prob_hd3) {
  # for 5'
  maintain_prob5 <- maintain_prob(
    cut_seq = seqs$up_seq,
    prod_hd = prob_hd5,
    bias = bias
  )

  cle_p5rev <- 1 - maintain_prob5[, length(prob_hd5):1]

  maintain_cumprod5rev <- matrix(
    data = 1,
    nrow = nrow(maintain_prob5),
    ncol = length(prob_hd5)
  )

  maintain_cumprod5rev[, -1] <- matrixStats::rowCumprods(
    maintain_prob5[, length(prob_hd5):2]
  )

  final_p5 <- maintain_cumprod5rev * cle_p5rev

  final_p5 <- final_p5[, ncol(final_p5):1]

  final_p5 <- final_p5 / Matrix::rowSums(final_p5)

  # for 3'
  maintain_prob3 <- maintain_prob(
    cut_seq = seqs$dn_seq,
    prod_hd = prob_hd3,
    bias = bias
  )

  cle_p3 <- 1 - maintain_prob3

  maintain_cumprod3 <- matrix(
    data = 1,
    nrow = nrow(maintain_prob3),
    ncol = length(prob_hd3)
  )

  maintain_cumprod3[, -1] <- matrixStats::rowCumprods(
    maintain_prob3[, -length(prob_hd3)]
  )

  final_p3 <- maintain_cumprod3 * cle_p3

  final_p3 <- final_p3 / Matrix::rowSums(final_p3)

  return(
    list(
      final_p5 = final_p5,
      final_p3 = final_p3
    )
  )
}

infer_ribo_num_train <- function(
    par_lst,
    tx_info,
    is_mnase = FALSE) {
  par_rpf <- par_lst
  ribo_size <- par_rpf$par_0$ribo_size
  cut5len <- seq.int(ribo_size[1])
  cut3len <- seq.int(ribo_size[4])
  convert_idx <- matrix(data = sapply(cut5len, function(cut5_i) {
    sapply(cut3len, function(cut3_i) {
      return(c(
        cut5_i,
        cut3_i,
        ribo_size[1L] - cut5_i + ribo_size[2L] + 1L,
        ribo_size[2L] + ribo_size[3L] + ribo_size[1L] - cut5_i + cut3_i
      ))
    })
  }), nrow = 4)
  par_rpf$lst_rpf$rpf_overlap <- par_rpf$lst_rpf$rpf_overlap[
    par_rpf$lst_rpf$rpf_overlap$qwidth %in% unique(convert_idx[4, ]),
  ]

  qwidth <- as.character(par_rpf$lst_rpf$rpf_overlap$qwidth)
  psite_num_idx <- as.vector(
    sapply(
      split(convert_idx[1, ], convert_idx[4, ]), length
    )[qwidth]
  )
  psite_num_idx <- as.integer(psite_num_idx)
  valid_idx <- !is.na(psite_num_idx) & psite_num_idx > 0L
  if (!all(valid_idx)) {
    par_rpf$lst_rpf$rpf_overlap <- par_rpf$lst_rpf$rpf_overlap[valid_idx, , drop = FALSE]
    qwidth <- qwidth[valid_idx]
    psite_num_idx <- psite_num_idx[valid_idx]
  }
  psite_pos <- rep(
    par_rpf$lst_rpf$rpf_overlap$pos, psite_num_idx
  ) +
    unlist(
      split(convert_idx[3, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  candi_psite <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = psite_pos
  )
  candi_cut5 <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = unlist(
      split(convert_idx[1, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  )
  candi_cut3 <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = unlist(
      split(convert_idx[2, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  )
  df_rpf <- par_rpf$lst_rpf$rpf_overlap
  par_0 <- par_rpf$par_0
  rg_cds <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p3 + 1L,
    end = tx_info$mix_tx_pos$utr3_p5 - 1L
  )

  frame_sum <- lapply(par_rpf$lst_rpf$offsets[1, ], function(len_i) {
    idx_i <- df_rpf$qwidth == len_i
    candi_psite <- candi_psite[idx_i, ]
    candi_cut5 <- candi_cut5[idx_i, ]
    candi_cut3 <- candi_cut3[idx_i, ]
    candi_p_weight <- candi_psite
    # shrink p-sites
    shrink_p <- as.factor(candi_psite@x)
    shrink_pos <- as.integer(shrink_p)
    index_j <- rank(
      shrink_pos,
      ties.method = "first"
    ) -
      rank(shrink_pos, ties.method = "min") + 1L
    candi_pos <- as.integer(levels(shrink_p))
    vec_pnum <- rep(1, length(candi_pos))
    # Split read weight
    uni_candi_psite <- unique(candi_psite@x)
    if (is_mnase == TRUE) {
      cut_seq <- get_cleavage_seq(
        seqs = tx_info$mix_tx, p_site = uni_candi_psite, ribo_size = ribo_size
      )
      idx_ij <- match(candi_psite@x, uni_candi_psite)
      # Generating cleavage probabilities for each ribosome terminus
      cut_prob <- get_cleavage_prob(
        seqs = cut_seq, bias = par_0$cut_bias$s7,
        prob_hd5 = par_0$prob_hd$p5, prob_hd3 = par_0$prob_hd$p3
      )
      base_prob <- cut_prob$final_p5[
        nrow(cut_prob$final_p5) * (candi_cut5@x - 1L) + idx_ij
      ] *
        cut_prob$final_p3[nrow(cut_prob$final_p3) * (candi_cut3@x - 1L) + idx_ij]
      candi_p_weight@x <- base_prob * vec_pnum[shrink_pos]
    } else {
      cut_seq <- get_cleavage_seq(
        seqs = tx_info$mix_tx, p_site = uni_candi_psite[1:2], ribo_size = ribo_size
      )
      idx_ij <- match(candi_psite@x, uni_candi_psite)
      # Generating cleavage probabilities for each ribosome terminus
      cut_prob <- get_cleavage_prob(
        seqs = cut_seq, bias = par_0$cut_bias$s7,
        prob_hd5 = par_0$prob_hd$p5, prob_hd3 = par_0$prob_hd$p3
      )
      cut_prob$final_p5 <- matrix(
        rep(cut_prob$final_p5[1, ], length(uni_candi_psite)),
        ncol = ncol(cut_prob$final_p5), byrow = T
      )
      cut_prob$final_p3 <- matrix(
        rep(cut_prob$final_p3[1, ], length(uni_candi_psite)),
        ncol = ncol(cut_prob$final_p3), byrow = T
      )
      base_prob <- cut_prob$final_p5[
        nrow(cut_prob$final_p5) * (candi_cut5@x - 1L) + idx_ij
      ] *
        cut_prob$final_p3[nrow(cut_prob$final_p3) * (candi_cut3@x - 1L) + idx_ij]
      candi_p_weight@x <- base_prob * vec_pnum[shrink_pos]
    }
    M <- candi_p_weight
    sx <- Matrix::summary(M)
    ord <- order(sx$i, -sx$x)
    sx2 <- sx[ord, ]
    sx2 <- sx2[!duplicated(sx2$i), ]
    sx$x <- 0L
    sx$x[as.numeric(rownames(sx2))] <- 1L
    M_argmax <- Matrix::sparseMatrix(
      i    = sx$i,
      j    = sx$j,
      x    = sx$x,
      dims = dim(M)
    )
    sparse_iter <- Matrix::sparseMatrix(
      i = shrink_pos,
      j = index_j,
      x = (df_rpf[idx_i, ]$weight * M_argmax)@x
    )
    vec_pnum[] <- Matrix::rowSums(sparse_iter)
    ribo_num <- list(
      candi_pos = candi_pos,
      vec_pnum = vec_pnum
    )
    tmp_psite_gr <- IRanges::IRanges(start = ribo_num$candi_pos, width = 1L)
    idx_cds_psite <- IRanges::findOverlaps(
      query = tmp_psite_gr,
      subject = rg_cds,
      minoverlap = 1L
    )
    psite_cds <- tmp_psite_gr@start[idx_cds_psite@from] - rg_cds@start[idx_cds_psite@to]
    cds_frame_pnum <- split(ribo_num$vec_pnum[idx_cds_psite@from], psite_cds %% 3L)
    frame_sum <- sapply(cds_frame_pnum, sum)
    frame_prob = sort(frame_sum / sum(frame_sum), decreasing = T)
    len_weight = (log(frame_prob[1] / (1 / 3)) - log(frame_prob[2] / (1 / 3)))
    frame_max <- which.max(frame_sum)
    if (frame_max == 1) {
      ribo_num$candi_pos <- ribo_num$candi_pos + 1L
    } else if (frame_max == 3) {
      ribo_num$candi_pos <- ribo_num$candi_pos - 1L
    }
    if (len_weight < 0.1) {
      ribo_num$vec_pnum <- ribo_num$vec_pnum * 0
    }
    return(list(
      ribo_num = ribo_num,
      frame_max = frame_max,
      len_weight = len_weight
    ))
  })
  ribo_num_pos <- colSums(do.call(rbind, lapply(frame_sum, function(x) {
    ribo_num_0 <- vector("numeric", length = length(tx_info_lst$tx_info$mix_tx))
    ribo_num_0[x$ribo_num$candi_pos] <- x$ribo_num$vec_pnum
    return(ribo_num_0)
  })))
  frame_max = sapply(frame_sum, function(x) x$frame_max)
  len_weight = sapply(frame_sum, function(x) x$len_weight)
  return(list(
    ribo_num_pos = ribo_num_pos,
    frame_max = frame_max,
    len_weight = len_weight
  ))
}

shift_psite_new <- function(tx_info,
                            rpf_info,
                            offsets,
                            shift_pos,
                            ribo_num = NULL) {
  # create negative RPF
  rpf_info_neg <- rpf_info
  rpf_info_neg$pos <- rpf_info$pos + sample(shift_pos, nrow(rpf_info), replace = TRUE)

  # initial p-site number (POS)
  p_table <- base::table(
    offsets[2, ][match(rpf_info$qwidth, offsets[1, ])] +
      rpf_info$pos
  )
  # NEG
  p_table_neg <- base::table(
    offsets[2, ][match(rpf_info_neg$qwidth, offsets[1, ])] +
      rpf_info_neg$pos
  )

  p_pos <- p_pos_neg <- integer(length(tx_info$mix_tx))
  p_pos[as.integer(names(p_table))] <- as.vector(p_table)
  p_pos_neg[as.integer(names(p_table_neg))] <- as.vector(p_table_neg)

  choose_p <- which(p_pos > 0L)
  choose_p_neg <- which(p_pos_neg > 0L)

  ribo_num_pos <- list(candi_pos = choose_p, vec_pnum = p_pos[choose_p])
  ribo_num_neg <- list(candi_pos = choose_p_neg, vec_pnum = p_pos_neg[choose_p_neg])
  if (is.null(ribo_num)) {
    return(list(
      ribo_num = ribo_num_pos,
      ribo_num_neg = ribo_num_neg
    ))
  } else {
    p_pos <- ribo_num
    choose_p <- which(p_pos > 0L)
    ribo_num_pos <- list(candi_pos = choose_p, vec_pnum = p_pos[choose_p])
    return(list(
      ribo_num = ribo_num_pos,
      ribo_num_neg = ribo_num_neg
    ))
  }
}

train_features <- function(ribo_num_train,
                           candi_translat,
                           tx_info,
                           min_pos_num,
                           eff_ribo_num,
                           bin_num) {
  library(IRanges)
  library(S4Vectors)
  library(Matrix)
  library(Biostrings)
  library(randomForest)
  ## ===================== Helper1: 基础 spa 矩阵 =====================

  build_spa_base <- function(rg_p, candi_rg, p_num,
                             eff_ribo_num, min_pos_num) {
    hits <- IRanges::findOverlaps(rg_p, candi_rg)
    if (length(hits) == 0L) {
      return(NULL)
    }

    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)

    # 每个 ORF 的 P-site 命中数；过滤 read 数太少的 ORF
    orf_counts <- tabulate(sh, nbins = length(candi_rg))
    keep_orf <- which(orf_counts >= min_pos_num)
    hits <- hits[sh %in% keep_orf]
    if (length(hits) == 0L) {
      return(NULL)
    }

    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)

    # ORF index → 行号 1..n_orf
    orf_levels <- sort(unique(sh))
    orf_factor <- factor(sh, levels = orf_levels)
    idx_orf <- as.integer(orf_factor)
    n_orf <- length(orf_levels)

    # 与 ORF 起点的距离（nt）
    dist_p <- IRanges::start(rg_p)[qh] - IRanges::start(candi_rg)[sh]
    idx_pos <- (dist_p %/% 3L) + 1L
    idx_frame <- dist_p %% 3L
    max_pos <- max(idx_pos)

    w <- p_num[qh] # 每个 P-site 的权重（核糖体数）

    make_spa_frame <- function(which_frame) {
      sel <- idx_frame == which_frame
      Matrix::sparseMatrix(
        i    = c(idx_orf[sel], n_orf + 1L),
        j    = c(idx_pos[sel], max_pos),
        x    = c(w[sel], 1),
        dims = c(n_orf + 1L, max_pos)
      )[1:n_orf, , drop = FALSE]
    }

    spa_frame0 <- make_spa_frame(0L)
    spa_frame1 <- make_spa_frame(1L)
    spa_frame2 <- make_spa_frame(2L)

    spa_all <- spa_frame0 + spa_frame1 + spa_frame2
    spa_eff <- spa_all > eff_ribo_num # 有效 codon

    # 三帧补零对齐到 spa_all 的非零位置
    spa0 <- spa_all
    spa0@x <- rep(0, length(spa0@x))
    spa_frame0 <- spa_frame0 + spa0
    spa_frame1 <- spa_frame1 + spa0
    spa_frame2 <- spa_frame2 + spa0

    eps <- 1e-12
    frac_from <- function(num_mat, den_mat) {
      out <- den_mat
      out@x <- (num_mat@x + eps) / (den_mat@x + 3 * eps)
      out
    }
    spa_frame0_fq <- frac_from(spa_frame0, spa_all)
    spa_frame1_fq <- frac_from(spa_frame1, spa_all)
    spa_frame2_fq <- frac_from(spa_frame2, spa_all)

    orf_aa <- Biostrings::width(candi_rg[orf_levels]) / 3L

    list(
      orf_id        = orf_levels,
      orf_aa        = orf_aa,
      spa_all       = spa_all,
      spa_eff       = spa_eff,
      spa_frame0    = spa_frame0,
      spa_frame1    = spa_frame1,
      spa_frame2    = spa_frame2,
      spa_frame0_fq = spa_frame0_fq,
      spa_frame1_fq = spa_frame1_fq,
      spa_frame2_fq = spa_frame2_fq
    )
  }

  ## ===================== Helper2: 加帧相关特征 =====================

  add_frame_features <- function(spa_lst, cds_frame = NULL) {
    eps <- 1e-12

    # log-odds
    logit_from <- function(fq_mat) {
      out <- fq_mat
      p <- fq_mat@x
      out@x <- log((p + eps) / (1 - p + eps))
      out
    }
    bias_logodds_0 <- logit_from(spa_lst$spa_frame0_fq)
    bias_logodds_1 <- logit_from(spa_lst$spa_frame1_fq)
    bias_logodds_2 <- logit_from(spa_lst$spa_frame2_fq)

    # 熵
    entropy_from <- function(fq_mat) {
      out <- fq_mat
      p <- fq_mat@x
      out@x <- -(p + eps) * log(p + eps)
      out
    }
    entropy_0 <- entropy_from(spa_lst$spa_frame0_fq)
    entropy_1 <- entropy_from(spa_lst$spa_frame1_fq)
    entropy_2 <- entropy_from(spa_lst$spa_frame2_fq)

    # 均匀 (1/3,1/3,1/3) 的 L2 距离
    dist_uniform_1 <- spa_lst$spa_all
    dist_uniform_1@x <- sqrt(
      (spa_lst$spa_frame0_fq@x - 1 / 3)^2 +
        (spa_lst$spa_frame1_fq@x - 1 / 3)^2 +
        (spa_lst$spa_frame2_fq@x - 1 / 3)^2
    )

    # 如给了 cds_frame，则生成更多“目标帧分布”的距离
    if (!is.null(cds_frame)) {
      cds_frame <- cds_frame / sum(cds_frame)
      cds_frame_ahead <- cds_frame[c(2, 3, 1)]
      cds_frame_behind <- cds_frame[c(3, 1, 2)]

      cf1 <- cds_frame + cds_frame_ahead
      cds_frame1 <- cf1 / sum(cf1)

      cf2 <- cds_frame + cds_frame_behind
      cds_frame2 <- cf2 / sum(cf2)

      cf3 <- cds_frame * (1 / 3) + cds_frame_ahead * (2 / 3)
      cds_frame3 <- cf3 / sum(cf3)

      cf4 <- cds_frame * (1 / 3) + cds_frame_behind * (2 / 3)
      cds_frame4 <- cf4 / sum(cf4)

      mk_dist <- function(target_vec) {
        m <- spa_lst$spa_all
        m@x <- sqrt(
          (spa_lst$spa_frame0_fq@x - target_vec[1])^2 +
            (spa_lst$spa_frame1_fq@x - target_vec[2])^2 +
            (spa_lst$spa_frame2_fq@x - target_vec[3])^2
        )
        m
      }

      dist_uniform_2 <- mk_dist(cds_frame_ahead)
      dist_uniform_3 <- mk_dist(cds_frame_behind)
      dist_uniform_4 <- mk_dist(cds_frame)
      dist_uniform_5 <- mk_dist(cds_frame1)
      dist_uniform_6 <- mk_dist(cds_frame2)
      dist_uniform_7 <- mk_dist(cds_frame3)
      dist_uniform_8 <- mk_dist(cds_frame4)
    } else {
      dist_uniform_2 <- dist_uniform_3 <- dist_uniform_4 <- NULL
      dist_uniform_5 <- dist_uniform_6 <- dist_uniform_7 <- dist_uniform_8 <- NULL
    }

    spa_lst$bias_logodds_0 <- bias_logodds_0
    spa_lst$bias_logodds_1 <- bias_logodds_1
    spa_lst$bias_logodds_2 <- bias_logodds_2
    spa_lst$entropy_0 <- entropy_0
    spa_lst$entropy_1 <- entropy_1
    spa_lst$entropy_2 <- entropy_2
    spa_lst$dist_uniform_1 <- dist_uniform_1
    spa_lst$dist_uniform_2 <- dist_uniform_2
    spa_lst$dist_uniform_3 <- dist_uniform_3
    spa_lst$dist_uniform_4 <- dist_uniform_4
    spa_lst$dist_uniform_5 <- dist_uniform_5
    spa_lst$dist_uniform_6 <- dist_uniform_6
    spa_lst$dist_uniform_7 <- dist_uniform_7
    spa_lst$dist_uniform_8 <- dist_uniform_8

    spa_lst
  }

  ## ===================== Helper3: ORF 级特征矩阵 =====================

  get_feature <- function(spa_lst) {
    bin_eff <- function(spa_eff) {
      Matrix::rowSums(spa_eff)
    }
    bin_prop_sum <- function(spa_prop, spa_eff) {
      Matrix::rowSums(spa_prop * spa_eff)
    }
    bin_weight_sum <- function(spa_ribo, spa_prop, spa_eff) {
      Matrix::rowSums(spa_eff * (spa_ribo * spa_prop))
    }

    # 按 codon 位置线性权重（递增 / 递减）
    row_linear_fill <- function(M, from = 0, to = 1) {
      stopifnot(inherits(M, "dgCMatrix"))
      S <- Matrix::summary(M)
      if (nrow(S) == 0L) {
        return(Matrix::sparseMatrix(
          i = integer(0),
          j = integer(0),
          x = numeric(0),
          dims = dim(M),
          dimnames = dimnames(M)
        ))
      }
      idx_by_row <- split(seq_len(nrow(S)), S$i)
      newx <- numeric(nrow(S))
      for (ind in idx_by_row) {
        k <- length(ind)
        if (k == 0L) next
        ord <- order(S$j[ind])
        idx <- ind[ord]
        newx[idx] <- seq(from, to, length.out = k)
      }
      Matrix::sparseMatrix(
        i = S$i,
        j = S$j,
        x = newx,
        dims = dim(M),
        dimnames = dimnames(M)
      )
    }

    max_frame <- function(f1, f2, f3) {
      data.frame(
        max1 = Matrix::rowSums((f1 < f2) & (f2 < f3)),
        max2 = Matrix::rowSums((f2 < f1) & (f1 < f3)),
        max3 = Matrix::rowSums((f2 < f3) & (f3 < f1)),
        max4 = Matrix::rowSums((f3 < f2) & (f2 < f1)),
        max5 = Matrix::rowSums((f1 < f3) & (f3 < f2)),
        max6 = Matrix::rowSums((f3 < f1) & (f1 < f2))
      )
    }

    spa_all_w <- spa_lst$spa_all
    spa_all_inc <- row_linear_fill(spa_all_w, from = 0, to = 1)
    spa_all_dec <- row_linear_fill(spa_all_w, from = 1, to = 0)

    feature_mat <- cbind(
      orf_aa       = spa_lst$orf_aa,
      rpf_all      = Matrix::rowSums(spa_lst$spa_all),
      rpf_f0       = Matrix::rowSums(spa_lst$spa_frame0),
      rpf_f1       = Matrix::rowSums(spa_lst$spa_frame1),
      rpf_f2       = Matrix::rowSums(spa_lst$spa_frame2),
      eff_codons   = bin_eff(spa_lst$spa_eff),
      f0_prop_sum  = bin_prop_sum(spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      f1_prop_sum  = bin_prop_sum(spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      f2_prop_sum  = bin_prop_sum(spa_lst$spa_frame2_fq, spa_lst$spa_eff),
      f0_weighted  = bin_weight_sum(spa_all_w, spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      f1_weighted  = bin_weight_sum(spa_all_w, spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      f2_weighted  = bin_weight_sum(spa_all_w, spa_lst$spa_frame2_fq, spa_lst$spa_eff),
      f0_inc_w     = bin_weight_sum(spa_all_inc, spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      f1_inc_w     = bin_weight_sum(spa_all_inc, spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      f2_inc_w     = bin_weight_sum(spa_all_inc, spa_lst$spa_frame2_fq, spa_lst$spa_eff),
      f0_dec_w     = bin_weight_sum(spa_all_dec, spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      f1_dec_w     = bin_weight_sum(spa_all_dec, spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      f2_dec_w     = bin_weight_sum(spa_all_dec, spa_lst$spa_frame2_fq, spa_lst$spa_eff),
      logodds0_sum = bin_prop_sum(spa_lst$bias_logodds_0, spa_lst$spa_eff),
      logodds1_sum = bin_prop_sum(spa_lst$bias_logodds_1, spa_lst$spa_eff),
      logodds2_sum = bin_prop_sum(spa_lst$bias_logodds_2, spa_lst$spa_eff),
      entropy0_sum = bin_prop_sum(spa_lst$entropy_0, spa_lst$spa_eff),
      entropy1_sum = bin_prop_sum(spa_lst$entropy_1, spa_lst$spa_eff),
      entropy2_sum = bin_prop_sum(spa_lst$entropy_2, spa_lst$spa_eff),
      dist_u1_sum  = bin_prop_sum(spa_lst$dist_uniform_1, spa_lst$spa_eff),
      dist_u2_sum  = if (!is.null(spa_lst$dist_uniform_2)) bin_prop_sum(spa_lst$dist_uniform_2, spa_lst$spa_eff) else 0,
      dist_u3_sum  = if (!is.null(spa_lst$dist_uniform_3)) bin_prop_sum(spa_lst$dist_uniform_3, spa_lst$spa_eff) else 0,
      dist_u4_sum  = if (!is.null(spa_lst$dist_uniform_4)) bin_prop_sum(spa_lst$dist_uniform_4, spa_lst$spa_eff) else 0,
      dist_u5_sum  = if (!is.null(spa_lst$dist_uniform_5)) bin_prop_sum(spa_lst$dist_uniform_5, spa_lst$spa_eff) else 0,
      dist_u6_sum  = if (!is.null(spa_lst$dist_uniform_6)) bin_prop_sum(spa_lst$dist_uniform_6, spa_lst$spa_eff) else 0,
      dist_u7_sum  = if (!is.null(spa_lst$dist_uniform_7)) bin_prop_sum(spa_lst$dist_uniform_7, spa_lst$spa_eff) else 0,
      dist_u8_sum  = if (!is.null(spa_lst$dist_uniform_8)) bin_prop_sum(spa_lst$dist_uniform_8, spa_lst$spa_eff) else 0
    )

    feature_mat <- cbind(
      feature_mat,
      max_frame(spa_lst$spa_frame0, spa_lst$spa_frame1, spa_lst$spa_frame2)
    )

    colnames(feature_mat) <- paste0("feat_", seq_len(ncol(feature_mat)))
    feature_mat
  }

  ## ===================== Helper4: 短 ORF 分 bin 特征 =====================

  build_norf_segments <- function(spa_lst,
                                  cds_frame,
                                  bin_num = 10,
                                  aa_min = 6,
                                  rpf_min = 1,
                                  eff_pos_min = 3) {
    stopifnot(inherits(spa_lst$spa_all, "dgCMatrix"))

    spa_all <- spa_lst$spa_all
    dims <- dim(spa_all)

    # 完全没信号直接空矩阵
    if (length(spa_all@x) == 0L) {
      return(matrix(
        numeric(0),
        nrow = 0,
        ncol = ncol(get_feature(add_frame_features(spa_lst, cds_frame)))
      ))
    }

    of <- Matrix::summary(spa_all) # columns: i(row), j(col), x
    aa_vec <- spa_lst$orf_aa

    # codon → bin 1..bin_num
    bin_idx <- floor((of$j - 0.5) / aa_vec[of$i] * bin_num) + 1L
    bin_idx[bin_idx < 1L] <- 1L
    bin_idx[bin_idx > bin_num] <- bin_num

    res_list <- vector("list", bin_num)

    for (b in seq_len(bin_num)) {
      sel <- (bin_idx == b)
      if (!any(sel)) next

      # 当前 bin 的 codon mask（0/1）
      mask <- Matrix::sparseMatrix(
        i    = of$i[sel],
        j    = of$j[sel],
        x    = 1,
        dims = dims
      )

      seg <- spa_lst

      seg$spa_all <- spa_lst$spa_all * mask
      seg$spa_frame0 <- spa_lst$spa_frame0 * mask
      seg$spa_frame1 <- spa_lst$spa_frame1 * mask
      seg$spa_frame2 <- spa_lst$spa_frame2 * mask

      seg$spa_eff <- seg$spa_all > 0

      eps <- 1e-12
      frac_from <- function(num_mat, den_mat) {
        out <- den_mat
        out@x <- (num_mat@x + eps) / (den_mat@x + 3 * eps)
        out
      }
      seg$spa_frame0_fq <- frac_from(seg$spa_frame0, seg$spa_all)
      seg$spa_frame1_fq <- frac_from(seg$spa_frame1, seg$spa_all)
      seg$spa_frame2_fq <- frac_from(seg$spa_frame2, seg$spa_all)

      # 这个 bin 中每条 ORF 实际有多少 codon
      seg$orf_aa <- Matrix::rowSums(mask > 0)

      keep_orf <- seg$orf_aa > 0
      if (!any(keep_orf)) next

      # 截掉全 0 行
      seg$orf_aa <- seg$orf_aa[keep_orf]
      seg$spa_all <- seg$spa_all[keep_orf, , drop = FALSE]
      seg$spa_eff <- seg$spa_eff[keep_orf, , drop = FALSE]
      seg$spa_frame0 <- seg$spa_frame0[keep_orf, , drop = FALSE]
      seg$spa_frame1 <- seg$spa_frame1[keep_orf, , drop = FALSE]
      seg$spa_frame2 <- seg$spa_frame2[keep_orf, , drop = FALSE]
      seg$spa_frame0_fq <- seg$spa_frame0_fq[keep_orf, , drop = FALSE]
      seg$spa_frame1_fq <- seg$spa_frame1_fq[keep_orf, , drop = FALSE]
      seg$spa_frame2_fq <- seg$spa_frame2_fq[keep_orf, , drop = FALSE]

      # 加 log-odds / entropy / dist_uniform_x
      seg <- add_frame_features(seg, cds_frame = cds_frame)

      # 和主模型一样抽特征
      feat <- get_feature(seg)

      # 用统一语义过滤：feat_1=orf_aa, feat_2=rpf_all, feat_6=eff_codons
      orf_aa <- feat[, "feat_1"]
      rpf_all <- feat[, "feat_2"]
      eff_cod <- feat[, "feat_6"]

      keep_seg <- (orf_aa >= aa_min &
        rpf_all >= rpf_min &
        eff_cod >= eff_pos_min)

      res_list[[b]] <- feat[keep_seg, , drop = FALSE]
    }

    feature_mat <- do.call(
      rbind,
      res_list[!vapply(res_list, is.null, logical(1))]
    )

    if (is.null(feature_mat)) {
      feature_mat <- matrix(
        numeric(0),
        nrow = 0,
        ncol = ncol(get_feature(add_frame_features(spa_lst, cds_frame)))
      )
    }

    feature_mat
  }


  ## ===================== Helper5: overlap ORF 特征 ==========================

  build_overlap_features <- function(back_frame,
                                     head_frame,
                                     cds_frame,
                                     bin_num = 10,
                                     aa_min = 6,
                                     rpf_min = 40,
                                     eff_pos_min = 13) {
    # back_frame / head_frame 都应该来自 build_spa_base() 或
    # build_spa_base() + add_frame_features()
    stopifnot(
      inherits(back_frame$spa_all, "dgCMatrix"),
      inherits(head_frame$spa_all, "dgCMatrix"),
      nrow(back_frame$spa_all) == nrow(head_frame$spa_all),
      length(back_frame$orf_aa) == length(head_frame$orf_aa)
    )

    ## ---------- 1) 小工具：按列数补零对齐 ----------

    pad_ncol <- function(M, target_ncol) {
      if (is.null(M)) {
        return(NULL)
      }
      if (ncol(M) == target_ncol) {
        return(M)
      }
      if (ncol(M) > target_ncol) {
        stop("pad_ncol: matrix has more columns than target_ncol.")
      }
      extra <- target_ncol - ncol(M)
      if (extra <= 0L) {
        return(M)
      }

      zero_block <- Matrix::Matrix(
        0,
        nrow   = nrow(M),
        ncol   = extra,
        sparse = TRUE
      )
      cbind(M, zero_block)
    }

    # 只对真正用来叠加的几个矩阵做列对齐
    max_cols <- max(ncol(back_frame$spa_all), ncol(head_frame$spa_all))
    for (nm in c("spa_all", "spa_frame0", "spa_frame1", "spa_frame2")) {
      back_frame[[nm]] <- pad_ncol(back_frame[[nm]], max_cols)
      head_frame[[nm]] <- pad_ncol(head_frame[[nm]], max_cols)
    }

    ## ---------- 2) 构造 overlap 的基础 spa（只是叠加） ----------

    eps <- 1e-12
    frac_from <- function(num_mat, den_mat) {
      out <- den_mat
      out@x <- (num_mat@x + eps) / (den_mat@x + 3 * eps)
      out
    }

    ov <- list()
    ov$orf_aa <- back_frame$orf_aa

    ov$spa_all <- back_frame$spa_all + head_frame$spa_all
    ov$spa_frame0 <- back_frame$spa_frame0 + head_frame$spa_frame0
    ov$spa_frame1 <- back_frame$spa_frame1 + head_frame$spa_frame1
    ov$spa_frame2 <- back_frame$spa_frame2 + head_frame$spa_frame2

    ov$spa_eff <- ov$spa_all > 0
    ov$spa_frame0_fq <- frac_from(ov$spa_frame0, ov$spa_all)
    ov$spa_frame1_fq <- frac_from(ov$spa_frame1, ov$spa_all)
    ov$spa_frame2_fq <- frac_from(ov$spa_frame2, ov$spa_all)

    ## ---------- 3) 直接当成新的 spa_lst，走“短 ORF 分段”逻辑 ----------
    ##     build_norf_segments 内部会：
    ##       - 调用 add_frame_features(ov, cds_frame)
    ##       - 调用 get_feature(…) 得到和主模型一致的特征列
    ##       - 再按 bin_num 做分段，筛选 aa_min / rpf_min / eff_pos_min

    build_norf_segments(
      spa_lst     = ov,
      cds_frame   = cds_frame,
      bin_num     = bin_num,
      aa_min      = aa_min,
      rpf_min     = rpf_min,
      eff_pos_min = eff_pos_min
    )
  }

  ## ===================== Helper6: RF模型 ==========================
  train_model_cds_rf <- function(m_pos, m_neg,
                                 sample_name = NULL,
                                 ntree = 500,
                                 mtry = NULL,
                                 classwt = NULL) {
    # -------- 组装训练数据 --------
    X <- rbind(m_pos, m_neg)
    y <- factor(
      c(
        rep(1L, nrow(m_pos)), # 正类：翻译 CDS
        rep(0L, nrow(m_neg))
      ), # 负类：非翻译
      levels = c(0L, 1L)
    )

    train_df <- data.frame(y = y, X)

    # -------- 默认 mtry：sqrt(p) --------
    if (is.null(mtry)) {
      mtry <- max(1L, floor(sqrt(ncol(X))))
    }

    # -------- 训练 randomForest（OOB 预测可直接用于验证 ROC） --------
    rf <- ranger::ranger(
      formula = y ~ .,
      data = train_df,
      num.trees = ntree,
      mtry = mtry,
      min.node.size = 10,
      sample.fraction = 1,
      max.depth = NULL,
      splitrule = "gini",
      probability = TRUE,  # 输出类别概率
      respect.unordered.factors = "partition",
      importance = "impurity",
      num.threads = max(1L, parallel::detectCores() - 1L),
      oob.error = TRUE
    )

    calc_roc <- function(score, label01) {
      ord <- order(score, decreasing = TRUE)
      label01 <- label01[ord]
      p <- sum(label01 == 1L)
      n <- sum(label01 == 0L)
      tp <- cumsum(label01 == 1L)
      fp <- cumsum(label01 == 0L)
      tpr <- c(0, tp / p, 1)
      fpr <- c(0, fp / n, 1)
      auc <- sum(diff(fpr) * (head(tpr, -1L) + tail(tpr, -1L)) / 2)
      list(fpr = fpr, tpr = tpr, auc = auc)
    }

    roc_val <- list(auc = NA_real_)
    roc_plot_file <- NA_character_
    if (!is.null(rf$predictions) && nrow(as.matrix(rf$predictions)) == nrow(train_df)) {
      pred_oob <- as.matrix(rf$predictions)[, "1"]
      y_true <- as.integer(as.character(train_df$y))
      roc_val <- calc_roc(pred_oob, y_true)
      roc_out_dir <- file.path(getwd(), "results")
      if (!dir.exists(roc_out_dir)) {
        dir.create(roc_out_dir, recursive = TRUE, showWarnings = FALSE)
      }
      roc_plot_file <- file.path(
        roc_out_dir,
        if (!is.null(sample_name) && nzchar(sample_name)) {
          paste0(gsub("[^A-Za-z0-9._-]+", "_", sample_name), "_oob_roc.pdf")
        } else {
          paste0("oob_roc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
        }
      )
      grDevices::pdf(file = roc_plot_file, width = 6.5, height = 6)
      graphics::plot(
        roc_val$fpr, roc_val$tpr,
        type = "l", lwd = 2, col = "#1F77B4",
        xlab = "False Positive Rate",
        ylab = "True Positive Rate",
        main = sprintf("OOB ROC (AUC = %.3f)", roc_val$auc),
        xlim = c(0, 1), ylim = c(0, 1)
      )
      graphics::abline(0, 1, lty = 2, col = "grey50")
      grDevices::dev.off()
      message("Saved OOB ROC curve: ", roc_plot_file)
    } else {
      warning("No usable OOB predictions in ranger model; skip ROC plot.")
    }

    list(
      model = rf,
      n_pos = nrow(m_pos),
      n_neg = nrow(m_neg),
      val_auc = roc_val$auc,
      val_roc_plot = roc_plot_file
    )
  }

  classify_orf_relative_to_cds <- function(candi_ribo_trans, tx_info) {
    # 1) Build transcript intervals
    rg_tx <- IRanges::IRanges(
      start = tx_info$mix_tx_pos$utr5_p5,
      end   = tx_info$mix_tx_pos$utr3_p3
    )

    # 2) Match candidate ORFs to their transcripts
    hits <- IRanges::findOverlaps(candi_ribo_trans, rg_tx, select = "first")
    n <- length(candi_ribo_trans)

    # Default everything to Variant
    orftype <- rep("Variant", n)

    ok <- !is.na(hits)
    if (any(ok)) {
      # 3) Only compute coordinates for the matched entries
      s <- IRanges::start(candi_ribo_trans)[ok]
      e <- IRanges::end(candi_ribo_trans)[ok]

      # Define UTR/CDS boundaries
      bp5 <- tx_info$mix_tx_pos$utr5_p3[hits[ok]] + 1L # CDS start
      bp3 <- tx_info$mix_tx_pos$utr3_p5[hits[ok]] - 1L # CDS end

      # 4) Four region indicators (-1/0/1)
      idx11 <- ifelse(s < bp5, -1L, ifelse(s == bp5, 0L, 1L))
      idx12 <- ifelse(s < bp3, -1L, 1L) # no 0 branch
      idx21 <- ifelse(e < bp5, -1L, 1L) # no 0 branch
      idx22 <- ifelse(e < bp3, -1L, ifelse(e == bp3, 0L, 1L))

      key <- paste0(idx11, idx12, idx21, idx22)

      # 5) Lookup table keeps the rules centralized and easy to read
      lut <- c(
        "-1-1-1-1" = "uORF",
        "-1-11-1"  = "uoORF",
        "-1-110"   = "Ext",
        "0-110"    = "CDS",
        "1-11-1"   = "intORF",
        "1-110"    = "Trunc",
        "1-111"    = "doORF",
        "1111"     = "dORF"
      )

      mapped <- lut[match(key, names(lut))]
      orftype[ok] <- ifelse(is.na(mapped), "Variant", mapped)
    }

    # Fix the factor levels for downstream consistency
    factor(
      orftype,
      levels = c("CDS", "uORF", "uoORF", "intORF", "Ext", "Trunc", "doORF", "dORF", "Variant")
    )
  }

  ## ===================== 1. P-site & ORF 定义 ==============================

  rg_p_pos <- IRanges::IRanges(
    start = ribo_num_train$ribo_num$candi_pos,
    width = 1L
  )
  rg_p_neg <- IRanges::IRanges(
    start = ribo_num_train$ribo_num_neg$candi_pos,
    width = 1L
  )

  # CDS 区域
  rg_cds <- IRanges::IRanges(
    start = tx_info$mix_tx_pos$utr5_p3 + 1L,
    end   = tx_info$mix_tx_pos$utr3_p5 - 1L
  )

  # 所有候选翻译区域
  candi_all_orf <- do.call(c, candi_translat)
  candi_all_orf@start <- candi_all_orf@start + 3L
  candi_all_orf@width <- candi_all_orf@width - 1L
  candi_all_orf_back <- candi_all_orf

  # CDS stop-stop 对应 ORF
  cds_idx <- IRanges::end(candi_all_orf) %in% IRanges::end(rg_cds)
  candi_orf <- candi_all_orf[cds_idx]

  ## ===================== 2. 全 ORF + 全局 cds_frame =======================

  spa_all_pos <- build_spa_base(
    rg_p         = rg_p_pos,
    candi_rg     = candi_all_orf,
    p_num        = ribo_num_train$ribo_num$vec_pnum,
    eff_ribo_num = eff_ribo_num,
    min_pos_num  = min_pos_num
  )
  if (is.null(spa_all_pos)) {
    stop("No ORF passed min_pos_num in positive set (all ORFs).")
  }

  cds_rows <- sort(unique(
    IRanges::findOverlaps(
      query = candi_all_orf_back[spa_all_pos$orf_id],
      subject = candi_orf,
      type = "equal"
    )@from
  ))

  if (length(cds_rows) > 0L) {
    cds_frame_all <- c(
      sum(spa_all_pos$spa_frame0[cds_rows, ]),
      sum(spa_all_pos$spa_frame1[cds_rows, ]),
      sum(spa_all_pos$spa_frame2[cds_rows, ])
    )
    cds_frame_all <- cds_frame_all / sum(spa_all_pos$spa_all[cds_rows, ])
  } else {
    cds_frame_all <- c(
      sum(spa_all_pos$spa_frame0),
      sum(spa_all_pos$spa_frame1),
      sum(spa_all_pos$spa_frame2)
    )
    cds_frame_all <- cds_frame_all / sum(spa_all_pos$spa_all)
  }
  spa_all_pos <- add_frame_features(spa_all_pos, cds_frame_all)
  feature_mat_all <- get_feature(spa_all_pos)
  candi_all_orf <- candi_all_orf_back[spa_all_pos$orf_id]

  ## ===================== 3. CDS 正/负集 + ahead/behind =====================

  spa_frame_pos <- build_spa_base(
    rg_p         = rg_p_pos,
    candi_rg     = candi_orf,
    p_num        = ribo_num_train$ribo_num$vec_pnum,
    eff_ribo_num = eff_ribo_num,
    min_pos_num  = min_pos_num
  )
  spa_frame_neg <- build_spa_base(
    rg_p         = rg_p_neg,
    candi_rg     = candi_orf,
    p_num        = ribo_num_train$ribo_num_neg$vec_pnum,
    eff_ribo_num = eff_ribo_num,
    min_pos_num  = min_pos_num
  )
  if (is.null(spa_frame_pos) || is.null(spa_frame_neg)) {
    stop("Not enough positive/negative CDS ORFs for training.")
  }

  cds_frame_cds <- c(
    sum(spa_frame_pos$spa_frame0),
    sum(spa_frame_pos$spa_frame1),
    sum(spa_frame_pos$spa_frame2)
  )
  cds_frame_cds <- cds_frame_cds / sum(spa_frame_pos$spa_all)

  spa_frame_pos <- add_frame_features(spa_frame_pos, cds_frame_cds)
  spa_frame_neg <- add_frame_features(spa_frame_neg, cds_frame_cds)
# 增加去掉不符合frame周期性的CDS
  feature_mat_p1 <- get_feature(spa_lst = spa_frame_pos)
  feature_mat_neg <- get_feature(spa_lst = spa_frame_neg)

  change_pos_lable <- function(feature_mat_p1, feature_mat_neg) {
    X <- rbind(feature_mat_p1, feature_mat_neg)
    y <- factor(c(rep(1L, nrow(feature_mat_p1)), rep(0L, nrow(feature_mat_neg))), levels = c(0L, 1L))
    train_df <- data.frame(y = y, X)
    # ---- Train Random Forest (ranger) ----
    rf <- ranger::ranger(
      formula = y ~ .,
      data = train_df,
      num.trees = 500,
      mtry = floor(sqrt(ncol(X))),
      min.node.size = 10,
      sample.fraction = 1,
      max.depth = NULL,
      splitrule = "gini",
      probability = TRUE,
      respect.unordered.factors = "partition",
      importance = "impurity",
      num.threads = max(1L, parallel::detectCores() - 1L)
    )
    p_val <- predict(rf, data = data.frame(feature_mat_p1))$predictions[, "0"]
    feature_mat_p2 <- feature_mat_p1[p_val < 0.5, ]
    X <- rbind(feature_mat_p2, feature_mat_neg)
    y <- factor(c(rep(1L, nrow(feature_mat_p2)), rep(0L, nrow(feature_mat_neg))), levels = c(0L, 1L))
    train_df <- data.frame(y = y, X)
    # ---- Train Random Forest (ranger) ----
    rf_2 <- ranger::ranger(
      formula = y ~ .,
      data = train_df,
      num.trees = 500,
      mtry = floor(sqrt(ncol(X))),
      min.node.size = 10,
      sample.fraction = 1,
      max.depth = NULL,
      splitrule = "gini",
      probability = TRUE,
      respect.unordered.factors = "partition",
      importance = "impurity",
      num.threads = max(1L, parallel::detectCores() - 1L)
    )
    p_val <- predict(rf_2, data = data.frame(feature_mat_p1))$predictions[, "0"]
    return(p_val)
  }

  prob_cds_false <- change_pos_lable(feature_mat_p1 = feature_mat_p1, feature_mat_neg = feature_mat_neg)

  # 得到清洗过后的CDS frame阳性和阴性集
  spa_frame_pos <- c(
    lapply(spa_frame_pos[c('orf_id', 'orf_aa')], function(x) x[prob_cds_false < 0.5]),
    lapply(spa_frame_pos[-which(names(spa_frame_pos) %in% c('orf_id', 'orf_aa'))], function(x) x[prob_cds_false < 0.5, ])
  )
  # 错框（+1 / -1）
  spa_frame_ahead <- spa_frame_pos
  spa_frame_behind <- spa_frame_pos

  spa_frame_ahead$spa_frame0 <- spa_frame_pos$spa_frame1
  spa_frame_ahead$spa_frame1 <- spa_frame_pos$spa_frame2
  spa_frame_ahead$spa_frame2 <- spa_frame_pos$spa_frame0
  spa_frame_ahead$spa_frame0_fq <- spa_frame_pos$spa_frame1_fq
  spa_frame_ahead$spa_frame1_fq <- spa_frame_pos$spa_frame2_fq
  spa_frame_ahead$spa_frame2_fq <- spa_frame_pos$spa_frame0_fq

  spa_frame_behind$spa_frame2 <- spa_frame_pos$spa_frame1
  spa_frame_behind$spa_frame1 <- spa_frame_pos$spa_frame0
  spa_frame_behind$spa_frame0 <- spa_frame_pos$spa_frame2
  spa_frame_behind$spa_frame2_fq <- spa_frame_pos$spa_frame1_fq
  spa_frame_behind$spa_frame1_fq <- spa_frame_pos$spa_frame0_fq
  spa_frame_behind$spa_frame0_fq <- spa_frame_pos$spa_frame2_fq

  spa_frame_ahead <- add_frame_features(spa_frame_ahead, cds_frame_cds)
  spa_frame_behind <- add_frame_features(spa_frame_behind, cds_frame_cds)

  feat_m_pos <- get_feature(spa_frame_pos)
  feat_m_ahead <- get_feature(spa_frame_ahead)
  feat_m_behind <- get_feature(spa_frame_behind)
  feature_mat_neg <- get_feature(spa_frame_neg)

  ## ===================== 4. 主 CDS 模型（randomForest） ===================


  ## ===================== 5. 短 ORF 模型（ncORF，u/d/lnc） ====================

  feat_norf_pos <- build_norf_segments(spa_frame_pos, cds_frame = cds_frame_cds, bin_num = bin_num)
  feat_norf_neg <- build_norf_segments(spa_frame_neg, cds_frame = cds_frame_cds, bin_num = bin_num)
  feat_norf_ahead <- build_norf_segments(spa_frame_ahead, cds_frame = cds_frame_cds, bin_num = bin_num)
  feat_norf_behind <- build_norf_segments(spa_frame_behind, cds_frame = cds_frame_cds, bin_num = bin_num)


  ## ===================== 6. overlap ORF 模型 ===============================

  # overlap 正例：背景 = 错框 CDS，头 = 真 CDS
  feat_oorf_p1 <- build_overlap_features(
    back_frame = spa_frame_ahead,
    head_frame = spa_frame_pos,
    cds_frame  = cds_frame_cds,
    bin_num    = 10
  )
  feat_oorf_p2 <- build_overlap_features(
    back_frame = spa_frame_behind,
    head_frame = spa_frame_pos,
    cds_frame  = cds_frame_cds,
    bin_num    = 10
  )

  overlap_orf_pos <- rbind(feat_oorf_p1, feat_oorf_p2)
  uorf_pos <- feat_norf_pos
  uorf_neg <- rbind(feat_norf_neg, feat_norf_ahead, feat_norf_behind)
  cds_pos <- feat_m_pos
  cds_neg <- rbind(feature_mat_neg, feat_m_ahead, feat_m_behind)
  modle_orf <- train_model_cds_rf(
    m_pos = rbind(cds_pos, uorf_pos, overlap_orf_pos),
    m_neg = rbind(cds_neg, uorf_neg),
    sample_name = sample_name
  )
  ## ===================== 7. 返回结果 ========================================
  orf_type <- classify_orf_relative_to_cds(candi_all_orf, tx_info)
  score <- predict(modle_orf$model, data = data.frame(feature_mat_all))$predictions[, "1"]
  idx_score = ifelse(orf_type %in% c('CDS','Ext','Trunc','Variant'), 0.5, 
  ifelse(orf_type %in% c('uoORF','uORF','ncRNA'), 0.3, ifelse(orf_type %in% c('dORF'), 0.7, ifelse(orf_type %in% c('doORF'), 0.6, 0.75))))

  idx_orf = score >= idx_score 
  list(
    train_model = modle_orf,
    candidate_rg = candi_all_orf[idx_orf],
    orf_type = orf_type[idx_orf],
    score = score[idx_orf],
    cds_frame_cds = cds_frame_cds
  )
}

ribo_num_pos <- infer_ribo_num_train(
  par_lst = par_rpf,
  tx_info = tx_info_lst$tx_info,
  is_mnase = rnase.bias
)

psite_pos <- shift_psite_new(
    tx_info = tx_info_lst$tx_info,
    rpf_info = par_rpf$lst_rpf$rpf_info,
    offsets = par_rpf$lst_rpf$offsets,
    shift_pos = c(-6:-1, 1:6),
    ribo_num = ribo_num_pos$ribo_num_pos
  )

model_pos <- train_features(
    ribo_num_train = psite_pos,
    candi_translat = candidate_translat_reg,
    tx_info        = tx_info_lst$tx_info,
    min_pos_num    = 10,
    eff_ribo_num   = 1,
    bin_num = 10
  )

  predict_orfs <- function(
      tx_info,
      model_info,
      ribo_num) {
        orf_type = model_info$orf_type
        candidate_rg = model_info$candidate_rg
        y_pred_prob = model_info$score
        train_model = model_info$train_model
    idx_cds <- (orf_type %in% c("Ext", "CDS", "Trunc")) & (y_pred_prob > 0.5)
    # 起始位点训练模型（用真实数据的CDS起始位点为阳性，其余内部atg为阴性）
    train_start_features <- function(start_codon,
                                     real_ribo_num,
                                     candid_translat_orf,
                                     tx_info) {
      if (length(start_codon) == 1) {
        atg_start <- Biostrings::matchPattern(start_codon, tx_info$mix_tx)@ranges
      } else {
        atg_start <- do.call(
          c, lapply(
            start_codon, function(x) {
              Biostrings::matchPattern(x, tx_info$mix_tx)@ranges
            }
          )
        )
        atg_start <- BiocGenerics::sort(atg_start)
      }

      rg_cds <- IRanges::IRanges(
        start = tx_info$mix_tx_pos$utr5_p3 + 1L,
        end   = tx_info$mix_tx_pos$utr3_p5 - 1L
      )
      hits <- IRanges::findOverlaps(candid_translat_orf, rg_cds)
      rg_cds_pos <- rg_cds[hits@to]

      hit_start <- IRanges::findOverlaps(
        atg_start,
        rg_cds_pos,
        type = "within"
      )
      hit_start <- hit_start[
        (atg_start@start[hit_start@from] - rg_cds_pos@start[hit_start@to]) %%
          3L == 0L
      ]
      middle_rg <- IRanges::IRanges(
        start = atg_start@start[hit_start@from] - 15L,
        width = 30L
      )

      rg_p <- IRanges::IRanges(start = real_ribo_num$candi_pos, width = 1L)
      hit_p <- IRanges::findOverlaps(rg_p, middle_rg)
      # Use start-codon anchored coordinates: keep biological window [-15, +14]
      anchor_start <- atg_start@start[hit_start@from][hit_p@to]
      rel_to_start <- rg_p@start[hit_p@from] - anchor_start
      idx_orf <- as.integer(as.factor(hit_p@to))
      orf_id <- as.integer(levels(as.factor(hit_p@to)))
      valid_mid <- !is.na(rel_to_start) & rel_to_start >= -15L & rel_to_start <= 14L
      dist_p <- rel_to_start[valid_mid] + 16L
      idx_orf_mid <- idx_orf[valid_mid]
      x_mid <- real_ribo_num$vec_pnum[hit_p@from][valid_mid]
        spa_mid <- Matrix::sparseMatrix(
          i = dist_p,
          j = idx_orf_mid,
          x = x_mid,
          dims = c(30L, length(orf_id))
        )

        start_atg_info <- matrix(0, ncol = length(middle_rg), nrow = 30)
        spa_mid_mat <- as.matrix(spa_mid)
        n_fill <- min(length(orf_id), ncol(spa_mid_mat))
        if (n_fill > 0L) {
          orf_id2 <- orf_id[seq_len(n_fill)]
          for (k in seq_len(n_fill)) {
            col_id <- orf_id2[k]
            if (!is.na(col_id) && col_id >= 1L && col_id <= ncol(start_atg_info)) {
              start_atg_info[, col_id] <- spa_mid_mat[, k]
            }
          }
        }
        start_atg_info <- as.data.frame(t(start_atg_info))
      near_start <- stringr::str_split(
        string = as.character(Biostrings::subseq(
          x = Biostrings::DNAStringSet(
            tx_info$mix_tx
          )[rep(1L, length(middle_rg))],
          start = middle_rg@start + 10L,
          width = 10L
        )),
        pattern = "",
        simplify = TRUE
      )

      keep <- setdiff(seq_len(ncol(near_start)), 6:8)
      df_seq <- as.data.frame(near_start[, keep, drop = FALSE],
        stringsAsFactors = FALSE
      )

      # clean to A/C/G/T/N and set unified levels per position
      lvl <- c("A", "C", "G", "T", "N")
      for (j in seq_along(df_seq)) {
        x <- df_seq[[j]]
        x[!(x %in% lvl)] <- "N"
        df_seq[[j]] <- factor(x, levels = lvl)
      }
      colnames(df_seq) <- paste0("pos", keep)

      # sparse one-hot: one block of 5 columns per position
      X_seq <- Matrix::sparse.model.matrix(~ . - 1, data = df_seq)
      X <- cbind(as.matrix(start_atg_info), X_seq)
      y <- rep(0, nrow(X))
      y[BiocGenerics::"%in%"(atg_start[hit_start@from]@start, rg_cds@start)] <- 1L

      idx_pos <- which(y == 1)
      idx_neg <- which(y == 0)

      tr_pos <- sample(idx_pos, floor(0.8 * length(idx_pos))) # 80% 训练
      va_pos <- setdiff(idx_pos, tr_pos) # 20% 验证

      tr_neg <- sample(idx_neg, floor(0.8 * length(idx_neg)))
      va_neg <- setdiff(idx_neg, tr_neg)

      idx_train <- c(tr_pos, tr_neg)
      idx_val <- c(va_pos, va_neg)

      ## ---- Class weight（全局类不平衡）----
      pos_w <- sum(y[idx_train] == 0) / sum(y[idx_train] == 1)
      w_train <- ifelse(y[idx_train] == 1, pos_w, 1)
      w_val <- rep(1, length(idx_val)) # 验证集一般不用权重

      ## ---- DMatrix（把权重放进来）----
      dtrain <- xgboost::xgb.DMatrix(
        data   = X[idx_train, , drop = FALSE],
        label  = y[idx_train],
        weight = w_train
      )
      dval <- xgboost::xgb.DMatrix(
        data   = X[idx_val, , drop = FALSE],
        label  = y[idx_val],
        weight = w_val
      )

      ## ---- Train ----
      params <- list(
        objective        = "binary:logistic",
        eval_metric      = "aucpr",
        eta              = 0.2,
        max_depth        = 6,
        subsample        = 0.8,
        colsample_bytree = 0.8
      )

      bst <- xgboost::xgb.train(
        params = params,
        data = dtrain,
        nrounds = 1000,
        evals = list(train = dtrain, val = dval),
        early_stopping_rounds = 100,
        verbose = 0
      )
      eval_log <- bst$evaluation_log
      best_iter <- bst$best_iteration
      if (!is.null(eval_log) && nrow(eval_log) >= best_iter && "val_aucpr" %in% colnames(eval_log)) {
        message(
          sprintf(
            "Initiation model trained. best_iteration=%d, val_aucpr=%.6f",
            best_iter,
            as.numeric(eval_log$val_aucpr[best_iter])
          )
        )
      } else {
        message(sprintf("Initiation model trained. best_iteration=%d", best_iter))
      }

      ## ---- Predict ----
      pred_prob_val <- predict(bst, dval)


      pred <- predict(bst, dval)
      ths <- seq(0, 1, by = 0.001)
      f1 <- function(t, p, y) {
        prd <- as.integer(p >= t)
        TP <- sum(prd == 1 & y == 1)
        FP <- sum(prd == 1 & y == 0)
        FN <- sum(prd == 0 & y == 1)
        prec <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
        rec <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
        ifelse(prec + rec == 0, 0, 2 * prec * rec / (prec + rec))
      }
      t_opt <- ths[which.max(sapply(ths, f1, p = pred, y = y[idx_val]))]

      return(list(
        bst = bst,
        t_opt = t_opt
      ))
    }
    message("Start train initiation site model...")
    model_start <- train_start_features(
      start_codon = "ATG",
      tx_info = tx_info,
      real_ribo_num = ribo_num,
      candid_translat_orf = candidate_rg[idx_cds]
    )

    start_features <- function(start_codon,
                               real_ribo_num,
                               candid_translat_orf,
                               tx_info) {
      if (length(start_codon) == 1) {
        rg_start <- Biostrings::matchPattern(start_codon, tx_info$mix_tx)@ranges
      } else {
        rg_start <- do.call(
          c, lapply(
            start_codon, function(x) {
              Biostrings::matchPattern(x, tx_info$mix_tx)@ranges
            }
          )
        )
        rg_start <- BiocGenerics::sort(rg_start)
      }

      hit_start <- IRanges::findOverlaps(
        rg_start,
        candid_translat_orf,
        type = "within"
      )
      hit_start <- hit_start[
        (rg_start@start[hit_start@from] - candid_translat_orf@start[hit_start@to]) %%
          3L == 0L
      ]

      if (length(hit_start) == 0L) {
        return(list(
          start_atg_info = matrix(numeric(0), nrow = 0L, ncol = 0L),
          atg_start = rg_start,
          hit_start = hit_start
        ))
      }

      middle_rg <- IRanges::IRanges(
        start = rg_start@start[hit_start@from] - 15L,
        width = 30L
      )

      rg_p <- IRanges::IRanges(start = real_ribo_num$candi_pos, width = 1L)
      hit_p <- IRanges::findOverlaps(rg_p, middle_rg)
      # Use start-codon anchored coordinates: keep biological window [-15, +14]
      anchor_start <- rg_start@start[hit_start@from][hit_p@to]
      rel_to_start <- rg_p@start[hit_p@from] - anchor_start
      idx_orf <- as.integer(as.factor(hit_p@to))
      orf_id <- as.integer(levels(as.factor(hit_p@to)))
      valid_mid <- !is.na(rel_to_start) & rel_to_start >= -15L & rel_to_start <= 14L
      dist_p <- rel_to_start[valid_mid] + 16L
      idx_orf_mid <- idx_orf[valid_mid]
      x_mid <- real_ribo_num$vec_pnum[hit_p@from][valid_mid]
        spa_mid <- Matrix::sparseMatrix(
          i = dist_p,
          j = idx_orf_mid,
          x = x_mid,
          dims = c(30L, length(orf_id))
        )

        start_atg_info <- matrix(0, ncol = length(middle_rg), nrow = 30)
        spa_mid_mat <- as.matrix(spa_mid)
        n_fill <- min(length(orf_id), ncol(spa_mid_mat))
        if (n_fill > 0L) {
          orf_id2 <- orf_id[seq_len(n_fill)]
          for (k in seq_len(n_fill)) {
            col_id <- orf_id2[k]
            if (!is.na(col_id) && col_id >= 1L && col_id <= ncol(start_atg_info)) {
              start_atg_info[, col_id] <- spa_mid_mat[, k]
            }
          }
        }
        start_atg_info <- as.data.frame(t(start_atg_info))
      near_start <- stringr::str_split(
        string = as.character(Biostrings::subseq(
          x = Biostrings::DNAStringSet(
            tx_info$mix_tx
          )[rep(1L, length(middle_rg))],
          start = middle_rg@start + 10L,
          width = 10L
        )),
        pattern = "",
        simplify = TRUE
      )

      keep <- setdiff(seq_len(ncol(near_start)), 6:8)
      df_seq <- as.data.frame(near_start[, keep, drop = FALSE],
        stringsAsFactors = FALSE
      )

      # clean to A/C/G/T/N and set unified levels per position
      lvl <- c("A", "C", "G", "T", "N")
      for (j in seq_along(df_seq)) {
        x <- df_seq[[j]]
        x[!(x %in% lvl)] <- "N"
        df_seq[[j]] <- factor(x, levels = lvl)
      }
      colnames(df_seq) <- paste0("pos", keep)

      # sparse one-hot: one block of 5 columns per position
      X_seq <- Matrix::sparse.model.matrix(~ . - 1, data = df_seq)
      X <- cbind(as.matrix(start_atg_info), X_seq)

      return(list(
        start_atg_info = X,
        atg_start = rg_start,
        hit_start = hit_start
      ))
    }
    # pred near atg position

    pred_near_atg <- function(model_start,
                              lst_nonstart) {
      if (length(lst_nonstart$hit_start) == 0L ||
          nrow(lst_nonstart$start_atg_info) == 0L) {
        return(NULL)
      }
      near_atg <- predict(
        model_start$bst,
        lst_nonstart$start_atg_info
      )
      hit_start <- lst_nonstart$hit_start[near_atg >= model_start$t_opt]
      near_atg <- sapply(
        split(
          hit_start@from,
          factor(hit_start@to, levels = unique(hit_start@to))
        ),
        function(x) {
          x[1]
        }
      )
      if (length(near_atg) > 0) {
        return(
          list(
            idx_posi = as.integer(names(near_atg)),
            candi_atg = lst_nonstart$atg_start[near_atg]
          )
        )
      } else {
        return(NULL)
      }
    }
    message("Start predict orfs...")
    # 翻译区域起始位点预测 ####
    stop_pos_rg <- candidate_rg
    predict_score <- y_pred_prob

    lst_start <- start_features(
      start_codon = "ATG",
      tx_info = tx_info,
      real_ribo_num = ribo_num,
      candid_translat_orf = stop_pos_rg
    )

    hit_start <- lst_start$hit_start
    atg_start <- lst_start$atg_start

    if (length(hit_start) == 0) {
      rg_orf_pred <- NULL
    } else {
      start_prob <- predict(
        model_start$bst,
        lst_start$start_atg_info
      )
      grp <- factor(hit_start@to, levels = unique(hit_start@to))
      row_id <- seq_along(hit_start@to)
      pick_first <- tapply(
        row_id, grp,
        function(ix) {
          lab <- ifelse(start_prob[ix] >= model_start$t_opt, "y", "n")
          if (all(lab == "n")) ix[1L] else ix[which(lab == "y")[1L]]
        }
      )
      candi_atg_rows <- as.integer(pick_first)
      rg_candi_atg <- atg_start[hit_start@from[candi_atg_rows]]

      # 对于有ATG的翻译区域，预测ATG之前的non-ATG起始位点概率
      rg_nonatg <- IRanges::IRanges(
        start = stop_pos_rg[hit_start@to[candi_atg_rows]]@start,
        end = rg_candi_atg@start - 1L
      )

      rg_no_atg <- stop_pos_rg[-hit_start@to]

      lst_nonstart <- start_features(
        start_codon = c(
          "CTG", "GTG", "TTG",
          "AAG", "ACG", "AGG",
          "ATA", "ATC", "ATT"
        ),
        tx_info = tx_info,
        real_ribo_num = ribo_num,
        candid_translat_orf = rg_nonatg
      )

      lst_nonstart2 <- start_features(
        start_codon = c(
          "CTG", "GTG", "TTG",
          "AAG", "ACG", "AGG",
          "ATA", "ATC", "ATT"
        ),
        tx_info = tx_info,
        real_ribo_num = ribo_num,
        candid_translat_orf = rg_no_atg
      )

      candi_near <- pred_near_atg(
        model_start = model_start,
        lst_nonstart = lst_nonstart
      )

      if (is.null(candi_near)) {
        rg_orf_pred <- IRanges::IRanges(
          start = rg_candi_atg@start,
          end = BiocGenerics::end(
            stop_pos_rg[
              !BiocGenerics::"%in%"(stop_pos_rg, rg_no_atg)
            ]
          )
        )
      } else {
        rg_candi_atg[candi_near$idx_posi] <- candi_near$candi_atg

        rg_orf_pred <- IRanges::IRanges(
          start = rg_candi_atg@start,
          end = BiocGenerics::end(stop_pos_rg[hit_start@to[candi_atg_rows]])
        )
      }

      if (length(lst_nonstart2$hit_start) == 0) {
        rg_orf_pred2 <- NULL
      } else {
        near_atg_prob2 <- predict(
          model_start$bst,
          lst_nonstart2$start_atg_info
        )

        probs <- as.numeric(near_atg_prob2)
        t_opt <- model_start$t_opt
        to_vec <- lst_nonstart2$hit_start@to # group id per row
        from_vec <- lst_nonstart2$hit_start@from # atg_start index per row

        # Lock group order to first appearance
        grp <- factor(to_vec, levels = unique(to_vec))
        rows <- seq_along(to_vec)

        # For each group: pick first ≥ t_opt; if none, pick argmax
        pick_row <- tapply(rows, grp, function(ix) {
          ix <- ix[!is.na(probs[ix])]
          if (length(ix) == 0L) {
            return(NA_integer_)
          }
          ok <- which(probs[ix] >= t_opt)
          if (length(ok) > 0L) {
            ix[ok[1]] # first that meets threshold
          } else {
            ix[which.max(probs[ix])] # fallback: max within group
          }
        })

        # Drop empty groups (defensive)
        pick_row <- unlist(pick_row, use.names = TRUE)
        keep <- !is.na(pick_row)
        pick_row <- pick_row[keep]

        # Map back to ranges
        picked_to_id <- as.integer(names(pick_row)) # group ids (must index rg_no_atg)
        picked_from <- from_vec[pick_row] # index into atg_start

        rg_candi_atg2 <- lst_nonstart2$atg_start[picked_from]

        rg_orf_pred2 <- IRanges::IRanges(
          start = BiocGenerics::start(rg_candi_atg2),
          end   = BiocGenerics::end(rg_no_atg[picked_to_id])
        )
      }
      rg_orf_pred <- c(rg_orf_pred, rg_orf_pred2)

      predict_score <- predict_score[
        match(
          BiocGenerics::end(rg_orf_pred),
          BiocGenerics::end(stop_pos_rg)
        )
      ]
    }
    return(list(
      model_orf = list(
        model_start = model_start,
        model_reg = train_model
      ),
      pred_orf = rg_orf_pred,
      predict_score = predict_score
    ))
  }
  message("Start predict ORFs on real data...")
  orf_pred <- predict_orfs(
    tx_info = tx_info_lst$tx_info,
    model_info = model_pos,
    ribo_num = psite_pos$ribo_num
  )
predict_orfs_lnc <- function(par_lnc,
  tx_info,
                               ribo_num_pos,
                               is_mnase,
                               candidate_translat_reg,
                               model_learned,
                               model_pos,
                               min_pos_num,
                               eff_ribo_num) {
  infer_ribo_lncrna <- function(
    par_lst,
    ribo_num_pos,
    lncrna_info,
    is_mnase) {
  par_rpf <- par_lst
  ribo_size <- par_rpf$par_0$ribo_size
  cut5len <- seq.int(ribo_size[1])
  cut3len <- seq.int(ribo_size[4])
  convert_idx <- matrix(data = sapply(cut5len, function(cut5_i) {
    sapply(cut3len, function(cut3_i) {
      return(c(
        cut5_i,
        cut3_i,
        ribo_size[1L] - cut5_i + ribo_size[2L] + 1L,
        ribo_size[2L] + ribo_size[3L] + ribo_size[1L] - cut5_i + cut3_i
      ))
    })
  }), nrow = 4)
  par_rpf$lst_rpf$rpf_overlap <- par_rpf$lst_rpf$rpf_overlap[
    par_rpf$lst_rpf$rpf_overlap$qwidth %in% unique(convert_idx[4, ]),
  ]

  qwidth <- as.character(par_rpf$lst_rpf$rpf_overlap$qwidth)
  psite_num_idx <- as.vector(
    sapply(
      split(convert_idx[1, ], convert_idx[4, ]), length
    )[qwidth]
  )
  psite_num_idx <- as.integer(psite_num_idx)
  valid_idx <- !is.na(psite_num_idx) & psite_num_idx > 0L
  if (!all(valid_idx)) {
    par_rpf$lst_rpf$rpf_overlap <- par_rpf$lst_rpf$rpf_overlap[valid_idx, , drop = FALSE]
    qwidth <- qwidth[valid_idx]
    psite_num_idx <- psite_num_idx[valid_idx]
  }
  psite_pos <- rep(
    par_rpf$lst_rpf$rpf_overlap$pos, psite_num_idx
  ) +
    unlist(
      split(convert_idx[3, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  candi_psite <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = psite_pos
  )
  candi_cut5 <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = unlist(
      split(convert_idx[1, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  )
  candi_cut3 <- Matrix::sparseMatrix(
    i = rep(seq.int(nrow(par_rpf$lst_rpf$rpf_overlap)), psite_num_idx),
    j = sequence(nvec = psite_num_idx, from = 1L),
    x = unlist(
      split(convert_idx[2, ], convert_idx[4, ])[qwidth],
      use.names = FALSE
    )
  )
  df_rpf <- par_rpf$lst_rpf$rpf_overlap
  par_0 <- par_rpf$par_0
  len_levels <- par_rpf$lst_rpf$offsets[1, ]
  len_counts <- vapply(
    len_levels,
    function(x) sum(df_rpf$qwidth == x),
    integer(1)
  )
  min_reads_per_len <- max(20L, floor(0.05 * nrow(df_rpf)))
  keep_len <- len_counts >= min_reads_per_len
  if (!any(keep_len)) {
    keep_len <- len_counts > 0L
    warning(
      "No read-length group reaches min_reads_per_len in infer_ribo_lncrna; fallback to all non-empty lengths."
    )
  }

  frame_sum = lapply(seq_along(len_levels), function(len_i){
    idx_i <- df_rpf$qwidth == len_levels[len_i]
    if(!keep_len[len_i] || sum(idx_i) == 0L){
      return(NULL)
    }else{
      candi_psite <- candi_psite[idx_i, ]
    candi_cut5 <- candi_cut5[idx_i, ]
    candi_cut3 <- candi_cut3[idx_i, ]
    candi_p_weight <- candi_psite
    # shrink p-sites
    shrink_p <- as.factor(candi_psite@x)
    shrink_pos <- as.integer(shrink_p)
    index_j <- rank(
      shrink_pos,
      ties.method = "first"
    ) -
      rank(shrink_pos, ties.method = "min") + 1L
    candi_pos <- as.integer(levels(shrink_p))
    vec_pnum <- rep(1, length(candi_pos))
    # Split read weight
    uni_candi_psite <- unique(candi_psite@x)
    if (is_mnase == TRUE) {
      cut_seq <- get_cleavage_seq(
        seqs = lncrna_info$mix_tx, p_site = uni_candi_psite, ribo_size = ribo_size
      )
      idx_ij <- match(candi_psite@x, uni_candi_psite)
      # Generating cleavage probabilities for each ribosome terminus
      cut_prob <- get_cleavage_prob(
        seqs = cut_seq, bias = par_0$cut_bias$s7,
        prob_hd5 = par_0$prob_hd$p5, prob_hd3 = par_0$prob_hd$p3
      )
      base_prob <- cut_prob$final_p5[
        nrow(cut_prob$final_p5) * (candi_cut5@x - 1L) + idx_ij
      ] *
        cut_prob$final_p3[nrow(cut_prob$final_p3) * (candi_cut3@x - 1L) + idx_ij]
      candi_p_weight@x <- base_prob * vec_pnum[shrink_pos]
    } else {
      cut_seq <- get_cleavage_seq(
        seqs = lncrna_info$mix_tx, p_site = uni_candi_psite[1:2], ribo_size = ribo_size
      )
      idx_ij <- match(candi_psite@x, uni_candi_psite)
      # Generating cleavage probabilities for each ribosome terminus
      cut_prob <- get_cleavage_prob(
        seqs = cut_seq, bias = par_0$cut_bias$s7,
        prob_hd5 = par_0$prob_hd$p5, prob_hd3 = par_0$prob_hd$p3
      )
      cut_prob$final_p5 <- matrix(
        rep(cut_prob$final_p5[1, ], length(uni_candi_psite)),
        ncol = ncol(cut_prob$final_p5), byrow = T
      )
      cut_prob$final_p3 <- matrix(
        rep(cut_prob$final_p3[1, ], length(uni_candi_psite)),
        ncol = ncol(cut_prob$final_p3), byrow = T
      )
      base_prob <- cut_prob$final_p5[
        nrow(cut_prob$final_p5) * (candi_cut5@x - 1L) + idx_ij
      ] *
        cut_prob$final_p3[nrow(cut_prob$final_p3) * (candi_cut3@x - 1L) + idx_ij]
      candi_p_weight@x <- base_prob * vec_pnum[shrink_pos]
    }
    M <- candi_p_weight
    sx <- Matrix::summary(M)
    ord <- order(sx$i, -sx$x)
    sx2 <- sx[ord, ]
    sx2 <- sx2[!duplicated(sx2$i), ]
    sx$x <- 0L
    sx$x[as.numeric(rownames(sx2))] <- 1L
    M_argmax <- Matrix::sparseMatrix(
      i    = sx$i,
      j    = sx$j,
      x    = sx$x,
      dims = dim(M)
    )
    sparse_iter <- Matrix::sparseMatrix(
      i = shrink_pos,
      j = index_j,
      x = (df_rpf[idx_i, ]$weight * M_argmax)@x
    )
    vec_pnum[] <- Matrix::rowSums(sparse_iter)
    ribo_num <- list(
      candi_pos = candi_pos,
      vec_pnum = vec_pnum
    )

    frame_max <- ribo_num_pos$frame_max[len_i]
    len_weight = ribo_num_pos$len_weight[len_i]
    if (frame_max == 1) {
      ribo_num$candi_pos <- ribo_num$candi_pos + 1L
    } else if (frame_max == 3) {
      ribo_num$candi_pos <- ribo_num$candi_pos - 1L
    }
    if (len_weight < 0.1) {
      ribo_num$vec_pnum <- ribo_num$vec_pnum * 0
    }
    return(ribo_num)
    }
  })
  idx <- vapply(
    frame_sum,
    function(x) {
      is.list(x) &&
        length(x$candi_pos) > 0L &&
        length(x$vec_pnum) > 0L
    },
    logical(1)
  )
  if (!any(idx)) {
    warning("No valid length group remained in infer_ribo_lncrna; returning zero ribo_num for lncRNA.")
    return(rep(0, length(lncrna_info$mix_tx)))
  }
  ribo_mat <- do.call(rbind, lapply(frame_sum[idx], function(x) {
    ribo_num_0 <- vector("numeric", length = length(lncrna_info$mix_tx))
    ribo_num_0[x$candi_pos] <- x$vec_pnum
    return(ribo_num_0)
  }))
  if (is.null(dim(ribo_mat))) {
    ribo_mat <- matrix(ribo_mat, nrow = 1L)
  }
  ribo_num_lnc <- colSums(ribo_mat)
  return(ribo_num_lnc)
}

get_features_lncrna <- function(ribo_num_lnc,
                           candi_translat,
                           tx_info,
                           cds_frame_all,
                           min_pos_num,
                           eff_ribo_num) {
  library(IRanges)
  library(S4Vectors)
  library(Matrix)
  library(Biostrings)
  library(randomForest)
  ## ===================== Helper1: 基础 spa 矩阵 =====================

  build_spa_base <- function(rg_p, candi_rg, p_num,
                             eff_ribo_num, min_pos_num) {
    hits <- IRanges::findOverlaps(rg_p, candi_rg)
    if (length(hits) == 0L) {
      return(NULL)
    }

    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)

    # 每个 ORF 的 P-site 命中数；过滤 read 数太少的 ORF
    orf_counts <- tabulate(sh, nbins = length(candi_rg))
    keep_orf <- which(orf_counts >= min_pos_num)
    hits <- hits[sh %in% keep_orf]
    if (length(hits) == 0L) {
      return(NULL)
    }

    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)

    # ORF index → 行号 1..n_orf
    orf_levels <- sort(unique(sh))
    orf_factor <- factor(sh, levels = orf_levels)
    idx_orf <- as.integer(orf_factor)
    n_orf <- length(orf_levels)

    # 与 ORF 起点的距离（nt）
    dist_p <- IRanges::start(rg_p)[qh] - IRanges::start(candi_rg)[sh]
    idx_pos <- (dist_p %/% 3L) + 1L
    idx_frame <- dist_p %% 3L
    max_pos <- max(idx_pos)

    w <- p_num[qh] # 每个 P-site 的权重（核糖体数）

    make_spa_frame <- function(which_frame) {
      sel <- idx_frame == which_frame
      Matrix::sparseMatrix(
        i    = c(idx_orf[sel], n_orf + 1L),
        j    = c(idx_pos[sel], max_pos),
        x    = c(w[sel], 1),
        dims = c(n_orf + 1L, max_pos)
      )[1:n_orf, , drop = FALSE]
    }

    spa_frame0 <- make_spa_frame(0L)
    spa_frame1 <- make_spa_frame(1L)
    spa_frame2 <- make_spa_frame(2L)

    spa_all <- spa_frame0 + spa_frame1 + spa_frame2
    spa_eff <- spa_all > eff_ribo_num # 有效 codon

    # 三帧补零对齐到 spa_all 的非零位置
    spa0 <- spa_all
    spa0@x <- rep(0, length(spa0@x))
    spa_frame0 <- spa_frame0 + spa0
    spa_frame1 <- spa_frame1 + spa0
    spa_frame2 <- spa_frame2 + spa0

    eps <- 1e-12
    frac_from <- function(num_mat, den_mat) {
      out <- den_mat
      out@x <- (num_mat@x + eps) / (den_mat@x + 3 * eps)
      out
    }
    spa_frame0_fq <- frac_from(spa_frame0, spa_all)
    spa_frame1_fq <- frac_from(spa_frame1, spa_all)
    spa_frame2_fq <- frac_from(spa_frame2, spa_all)

    orf_aa <- Biostrings::width(candi_rg[orf_levels]) / 3L

    list(
      orf_id        = orf_levels,
      orf_aa        = orf_aa,
      spa_all       = spa_all,
      spa_eff       = spa_eff,
      spa_frame0    = spa_frame0,
      spa_frame1    = spa_frame1,
      spa_frame2    = spa_frame2,
      spa_frame0_fq = spa_frame0_fq,
      spa_frame1_fq = spa_frame1_fq,
      spa_frame2_fq = spa_frame2_fq
    )
  }

  ## ===================== Helper2: 加帧相关特征 =====================

  add_frame_features <- function(spa_lst, cds_frame = NULL) {
    eps <- 1e-12

    # log-odds
    logit_from <- function(fq_mat) {
      out <- fq_mat
      p <- fq_mat@x
      out@x <- log((p + eps) / (1 - p + eps))
      out
    }
    bias_logodds_0 <- logit_from(spa_lst$spa_frame0_fq)
    bias_logodds_1 <- logit_from(spa_lst$spa_frame1_fq)
    bias_logodds_2 <- logit_from(spa_lst$spa_frame2_fq)

    # 熵
    entropy_from <- function(fq_mat) {
      out <- fq_mat
      p <- fq_mat@x
      out@x <- -(p + eps) * log(p + eps)
      out
    }
    entropy_0 <- entropy_from(spa_lst$spa_frame0_fq)
    entropy_1 <- entropy_from(spa_lst$spa_frame1_fq)
    entropy_2 <- entropy_from(spa_lst$spa_frame2_fq)

    # 均匀 (1/3,1/3,1/3) 的 L2 距离
    dist_uniform_1 <- spa_lst$spa_all
    dist_uniform_1@x <- sqrt(
      (spa_lst$spa_frame0_fq@x - 1 / 3)^2 +
        (spa_lst$spa_frame1_fq@x - 1 / 3)^2 +
        (spa_lst$spa_frame2_fq@x - 1 / 3)^2
    )

    # 如给了 cds_frame，则生成更多“目标帧分布”的距离
    if (!is.null(cds_frame)) {
      cds_frame <- cds_frame / sum(cds_frame)
      cds_frame_ahead <- cds_frame[c(2, 3, 1)]
      cds_frame_behind <- cds_frame[c(3, 1, 2)]

      cf1 <- cds_frame + cds_frame_ahead
      cds_frame1 <- cf1 / sum(cf1)

      cf2 <- cds_frame + cds_frame_behind
      cds_frame2 <- cf2 / sum(cf2)

      cf3 <- cds_frame * (1 / 3) + cds_frame_ahead * (2 / 3)
      cds_frame3 <- cf3 / sum(cf3)

      cf4 <- cds_frame * (1 / 3) + cds_frame_behind * (2 / 3)
      cds_frame4 <- cf4 / sum(cf4)

      mk_dist <- function(target_vec) {
        m <- spa_lst$spa_all
        m@x <- sqrt(
          (spa_lst$spa_frame0_fq@x - target_vec[1])^2 +
            (spa_lst$spa_frame1_fq@x - target_vec[2])^2 +
            (spa_lst$spa_frame2_fq@x - target_vec[3])^2
        )
        m
      }

      dist_uniform_2 <- mk_dist(cds_frame_ahead)
      dist_uniform_3 <- mk_dist(cds_frame_behind)
      dist_uniform_4 <- mk_dist(cds_frame)
      dist_uniform_5 <- mk_dist(cds_frame1)
      dist_uniform_6 <- mk_dist(cds_frame2)
      dist_uniform_7 <- mk_dist(cds_frame3)
      dist_uniform_8 <- mk_dist(cds_frame4)
    } else {
      dist_uniform_2 <- dist_uniform_3 <- dist_uniform_4 <- NULL
      dist_uniform_5 <- dist_uniform_6 <- dist_uniform_7 <- dist_uniform_8 <- NULL
    }

    spa_lst$bias_logodds_0 <- bias_logodds_0
    spa_lst$bias_logodds_1 <- bias_logodds_1
    spa_lst$bias_logodds_2 <- bias_logodds_2
    spa_lst$entropy_0 <- entropy_0
    spa_lst$entropy_1 <- entropy_1
    spa_lst$entropy_2 <- entropy_2
    spa_lst$dist_uniform_1 <- dist_uniform_1
    spa_lst$dist_uniform_2 <- dist_uniform_2
    spa_lst$dist_uniform_3 <- dist_uniform_3
    spa_lst$dist_uniform_4 <- dist_uniform_4
    spa_lst$dist_uniform_5 <- dist_uniform_5
    spa_lst$dist_uniform_6 <- dist_uniform_6
    spa_lst$dist_uniform_7 <- dist_uniform_7
    spa_lst$dist_uniform_8 <- dist_uniform_8

    spa_lst
  }

  ## ===================== Helper3: ORF 级特征矩阵 =====================

  get_feature <- function(spa_lst) {
    bin_eff <- function(spa_eff) {
      Matrix::rowSums(spa_eff)
    }
    bin_prop_sum <- function(spa_prop, spa_eff) {
      Matrix::rowSums(spa_prop * spa_eff)
    }
    bin_weight_sum <- function(spa_ribo, spa_prop, spa_eff) {
      Matrix::rowSums(spa_eff * (spa_ribo * spa_prop))
    }

    # 按 codon 位置线性权重（递增 / 递减）
    row_linear_fill <- function(M, from = 0, to = 1) {
      stopifnot(inherits(M, "dgCMatrix"))
      S <- Matrix::summary(M)
      if (nrow(S) == 0L) {
        return(Matrix::sparseMatrix(
          i = integer(0),
          j = integer(0),
          x = numeric(0),
          dims = dim(M),
          dimnames = dimnames(M)
        ))
      }
      idx_by_row <- split(seq_len(nrow(S)), S$i)
      newx <- numeric(nrow(S))
      for (ind in idx_by_row) {
        k <- length(ind)
        if (k == 0L) next
        ord <- order(S$j[ind])
        idx <- ind[ord]
        newx[idx] <- seq(from, to, length.out = k)
      }
      Matrix::sparseMatrix(
        i = S$i,
        j = S$j,
        x = newx,
        dims = dim(M),
        dimnames = dimnames(M)
      )
    }

    max_frame <- function(f1, f2, f3) {
      data.frame(
        max1 = Matrix::rowSums((f1 < f2) & (f2 < f3)),
        max2 = Matrix::rowSums((f2 < f1) & (f1 < f3)),
        max3 = Matrix::rowSums((f2 < f3) & (f3 < f1)),
        max4 = Matrix::rowSums((f3 < f2) & (f2 < f1)),
        max5 = Matrix::rowSums((f1 < f3) & (f3 < f2)),
        max6 = Matrix::rowSums((f3 < f1) & (f1 < f2))
      )
    }

    spa_all_w <- spa_lst$spa_all
    spa_all_inc <- row_linear_fill(spa_all_w, from = 0, to = 1)
    spa_all_dec <- row_linear_fill(spa_all_w, from = 1, to = 0)

    feature_mat <- cbind(
      orf_aa       = spa_lst$orf_aa,
      rpf_all      = Matrix::rowSums(spa_lst$spa_all),
      rpf_f0       = Matrix::rowSums(spa_lst$spa_frame0),
      rpf_f1       = Matrix::rowSums(spa_lst$spa_frame1),
      rpf_f2       = Matrix::rowSums(spa_lst$spa_frame2),
      eff_codons   = bin_eff(spa_lst$spa_eff),
      f0_prop_sum  = bin_prop_sum(spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      f1_prop_sum  = bin_prop_sum(spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      f2_prop_sum  = bin_prop_sum(spa_lst$spa_frame2_fq, spa_lst$spa_eff),
      f0_weighted  = bin_weight_sum(spa_all_w, spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      f1_weighted  = bin_weight_sum(spa_all_w, spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      f2_weighted  = bin_weight_sum(spa_all_w, spa_lst$spa_frame2_fq, spa_lst$spa_eff),
      f0_inc_w     = bin_weight_sum(spa_all_inc, spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      f1_inc_w     = bin_weight_sum(spa_all_inc, spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      f2_inc_w     = bin_weight_sum(spa_all_inc, spa_lst$spa_frame2_fq, spa_lst$spa_eff),
      f0_dec_w     = bin_weight_sum(spa_all_dec, spa_lst$spa_frame0_fq, spa_lst$spa_eff),
      f1_dec_w     = bin_weight_sum(spa_all_dec, spa_lst$spa_frame1_fq, spa_lst$spa_eff),
      f2_dec_w     = bin_weight_sum(spa_all_dec, spa_lst$spa_frame2_fq, spa_lst$spa_eff),
      logodds0_sum = bin_prop_sum(spa_lst$bias_logodds_0, spa_lst$spa_eff),
      logodds1_sum = bin_prop_sum(spa_lst$bias_logodds_1, spa_lst$spa_eff),
      logodds2_sum = bin_prop_sum(spa_lst$bias_logodds_2, spa_lst$spa_eff),
      entropy0_sum = bin_prop_sum(spa_lst$entropy_0, spa_lst$spa_eff),
      entropy1_sum = bin_prop_sum(spa_lst$entropy_1, spa_lst$spa_eff),
      entropy2_sum = bin_prop_sum(spa_lst$entropy_2, spa_lst$spa_eff),
      dist_u1_sum  = bin_prop_sum(spa_lst$dist_uniform_1, spa_lst$spa_eff),
      dist_u2_sum  = if (!is.null(spa_lst$dist_uniform_2)) bin_prop_sum(spa_lst$dist_uniform_2, spa_lst$spa_eff) else 0,
      dist_u3_sum  = if (!is.null(spa_lst$dist_uniform_3)) bin_prop_sum(spa_lst$dist_uniform_3, spa_lst$spa_eff) else 0,
      dist_u4_sum  = if (!is.null(spa_lst$dist_uniform_4)) bin_prop_sum(spa_lst$dist_uniform_4, spa_lst$spa_eff) else 0,
      dist_u5_sum  = if (!is.null(spa_lst$dist_uniform_5)) bin_prop_sum(spa_lst$dist_uniform_5, spa_lst$spa_eff) else 0,
      dist_u6_sum  = if (!is.null(spa_lst$dist_uniform_6)) bin_prop_sum(spa_lst$dist_uniform_6, spa_lst$spa_eff) else 0,
      dist_u7_sum  = if (!is.null(spa_lst$dist_uniform_7)) bin_prop_sum(spa_lst$dist_uniform_7, spa_lst$spa_eff) else 0,
      dist_u8_sum  = if (!is.null(spa_lst$dist_uniform_8)) bin_prop_sum(spa_lst$dist_uniform_8, spa_lst$spa_eff) else 0
    )

    feature_mat <- cbind(
      feature_mat,
      max_frame(spa_lst$spa_frame0, spa_lst$spa_frame1, spa_lst$spa_frame2)
    )

    colnames(feature_mat) <- paste0("feat_", seq_len(ncol(feature_mat)))
    feature_mat
  }

  ## ===================== 1. P-site & ORF 定义 ==============================

  rg_p_pos <- IRanges::IRanges(
    start = ribo_num_lnc$candi_pos,
    width = 1L
  )

  # 所有候选翻译区域
  candi_all_orf <- do.call(c, candi_translat)
  candi_all_orf@start <- candi_all_orf@start + 3L
  candi_all_orf@width <- candi_all_orf@width - 1L

  ## ===================== 2. 全 ORF + 全局 cds_frame =======================

  spa_all_pos <- build_spa_base(
    rg_p         = rg_p_pos,
    candi_rg     = candi_all_orf,
    p_num        = ribo_num_lnc$vec_pnum,
    eff_ribo_num = eff_ribo_num,
    min_pos_num  = min_pos_num
  )
  if (is.null(spa_all_pos)) {
    stop("No ORF passed min_pos_num in positive set (all ORFs).")
  }

  spa_all_pos <- add_frame_features(spa_all_pos, cds_frame_all)
  feature_mat_all <- get_feature(spa_all_pos)
    candi_all_orf <- candi_all_orf[spa_all_pos$orf_id]

  return(list(
    feature_mat_all = feature_mat_all,
    candi_orf = candi_all_orf
  ) )
}

    start_features <- function(start_codon,
                               real_ribo_num,
                               candid_translat_orf,
                               tx_info) {
      if (length(start_codon) == 1) {
        rg_start <- Biostrings::matchPattern(start_codon, tx_info$mix_tx)@ranges
      } else {
        rg_start <- do.call(
          c, lapply(
            start_codon, function(x) {
              Biostrings::matchPattern(x, tx_info$mix_tx)@ranges
            }
          )
        )
        rg_start <- BiocGenerics::sort(rg_start)
      }

      hit_start <- IRanges::findOverlaps(
        rg_start,
        candid_translat_orf,
        type = "within"
      )
      hit_start <- hit_start[
        (rg_start@start[hit_start@from] - candid_translat_orf@start[hit_start@to]) %%
          3L == 0L
      ]

      if (length(hit_start) == 0L) {
        return(list(
          start_atg_info = matrix(numeric(0), nrow = 0L, ncol = 0L),
          atg_start = rg_start,
          hit_start = hit_start
        ))
      }

      middle_rg <- IRanges::IRanges(
        start = rg_start@start[hit_start@from] - 15L,
        width = 30L
      )

      rg_p <- IRanges::IRanges(start = real_ribo_num$candi_pos, width = 1L)
      hit_p <- IRanges::findOverlaps(rg_p, middle_rg)
      # Use start-codon anchored coordinates: keep biological window [-15, +14]
      anchor_start <- rg_start@start[hit_start@from][hit_p@to]
      rel_to_start <- rg_p@start[hit_p@from] - anchor_start
      idx_orf <- as.integer(as.factor(hit_p@to))
      orf_id <- as.integer(levels(as.factor(hit_p@to)))
      valid_mid <- !is.na(rel_to_start) & rel_to_start >= -15L & rel_to_start <= 14L
      dist_p <- rel_to_start[valid_mid] + 16L
      idx_orf_mid <- idx_orf[valid_mid]
      x_mid <- real_ribo_num$vec_pnum[hit_p@from][valid_mid]
        spa_mid <- Matrix::sparseMatrix(
          i = dist_p,
          j = idx_orf_mid,
          x = x_mid,
          dims = c(30L, length(orf_id))
        )

        start_atg_info <- matrix(0, ncol = length(middle_rg), nrow = 30)
        spa_mid_mat <- as.matrix(spa_mid)
        n_fill <- min(length(orf_id), ncol(spa_mid_mat))
        if (n_fill > 0L) {
          orf_id2 <- orf_id[seq_len(n_fill)]
          for (k in seq_len(n_fill)) {
            col_id <- orf_id2[k]
            if (!is.na(col_id) && col_id >= 1L && col_id <= ncol(start_atg_info)) {
              start_atg_info[, col_id] <- spa_mid_mat[, k]
            }
          }
        }
        start_atg_info <- as.data.frame(t(start_atg_info))
 
      near_start <- stringr::str_split(
        string = as.character(Biostrings::subseq(
          x = Biostrings::DNAStringSet(
            tx_info$mix_tx
          )[rep(1L, length(middle_rg))],
          start = middle_rg@start + 10L,
          width = 10L
        )),
        pattern = "",
        simplify = TRUE
      )

      keep <- setdiff(seq_len(ncol(near_start)), 6:8)
      df_seq <- as.data.frame(near_start[, keep, drop = FALSE],
        stringsAsFactors = FALSE
      )

      # clean to A/C/G/T/N and set unified levels per position
      lvl <- c("A", "C", "G", "T", "N")
      for (j in seq_along(df_seq)) {
        x <- df_seq[[j]]
        x[!(x %in% lvl)] <- "N"
        df_seq[[j]] <- factor(x, levels = lvl)
      }
      colnames(df_seq) <- paste0("pos", keep)

      # sparse one-hot: one block of 5 columns per position
      X_seq <- Matrix::sparse.model.matrix(~ . - 1, data = df_seq)
      X <- cbind(as.matrix(start_atg_info), X_seq)

      return(list(
        start_atg_info = X,
        atg_start = rg_start,
        hit_start = hit_start
      ))
    }
    # pred near atg position

    pred_near_atg <- function(model_start,
                              lst_nonstart) {
      if (length(lst_nonstart$hit_start) == 0L ||
          nrow(lst_nonstart$start_atg_info) == 0L) {
        return(NULL)
      }
      near_atg <- predict(
        model_start$bst,
        lst_nonstart$start_atg_info
      )
      hit_start <- lst_nonstart$hit_start[near_atg >= model_start$t_opt]
      near_atg <- sapply(
        split(
          hit_start@from,
          factor(hit_start@to, levels = unique(hit_start@to))
        ),
        function(x) {
          x[1]
        }
      )
      if (length(near_atg) > 0) {
        return(
          list(
            idx_posi = as.integer(names(near_atg)),
            candi_atg = lst_nonstart$atg_start[near_atg]
          )
        )
      } else {
        return(NULL)
      }
    }

ribo_num = infer_ribo_lncrna(
  par_lst = par_lnc,
    ribo_num_pos = ribo_num_pos,
    lncrna_info = tx_info,
    is_mnase = is_mnase
)
    p_pos <- ribo_num
    choose_p <- which(p_pos > 0L)
    ribo_num_lnc <- list(candi_pos = choose_p, vec_pnum = p_pos[choose_p])

    message("Start extract translation features for real ribo-seq data...")
    translation_features <- get_features_lncrna(
      ribo_num_lnc = ribo_num_lnc,
    candi_translat = candidate_translat_reg,
    tx_info        = tx_info,
    cds_frame_all = model_pos$cds_frame_cds,
    min_pos_num    = 10,
    eff_ribo_num   = 1
    )

      predict_score <- predict(
      model_learned$model_orf$model_reg$model, 
      data = data.frame(translation_features$feature_mat_all))$predictions[, "1"]
      stop_pos_rg <- translation_features$candi_orf
      idx_orf = predict_score >= 0.5
      stop_pos_rg = stop_pos_rg[idx_orf]
      predict_score = predict_score[idx_orf]

    message("Start predict orfs...")
    # 翻译区域起始位点预测 ####
    lst_start <- start_features(
      start_codon = "ATG",
      tx_info = tx_info,
      real_ribo_num = ribo_num_lnc,
      candid_translat_orf = stop_pos_rg
    )

    hit_start <- lst_start$hit_start
    atg_start <- lst_start$atg_start

    if (length(hit_start) == 0) {
      rg_orf_pred <- NULL
    } else {
      start_prob <- predict(
        model_learned$model_orf$model_start$bst,
        lst_start$start_atg_info
      )
      grp <- factor(hit_start@to, levels = unique(hit_start@to))
      row_id <- seq_along(hit_start@to)
      pick_first <- tapply(
        row_id, grp,
        function(ix) {
          lab <- ifelse(start_prob[ix] >= model_learned$model_start$t_opt, "y", "n")
          if (all(lab == "n")) ix[1L] else ix[which(lab == "y")[1L]]
        }
      )
      candi_atg_rows <- as.integer(pick_first)
      rg_candi_atg <- atg_start[hit_start@from[candi_atg_rows]]

      # 对于有ATG的翻译区域，预测ATG之前的non-ATG起始位点概率
      rg_nonatg <- IRanges::IRanges(
        start = stop_pos_rg[hit_start@to[candi_atg_rows]]@start,
        end = rg_candi_atg@start - 1L
      )

      rg_no_atg <- stop_pos_rg[-hit_start@to]

      lst_nonstart <- start_features(
        start_codon = c(
          "CTG", "GTG", "TTG",
          "AAG", "ACG", "AGG",
          "ATA", "ATC", "ATT"
        ),
        tx_info = tx_info,
        real_ribo_num = ribo_num_lnc,
        candid_translat_orf = rg_nonatg
      )

      lst_nonstart2 <- start_features(
        start_codon = c(
          "CTG", "GTG", "TTG",
          "AAG", "ACG", "AGG",
          "ATA", "ATC", "ATT"
        ),
        tx_info = tx_info,
        real_ribo_num = ribo_num_lnc,
        candid_translat_orf = rg_no_atg
      )

      candi_near <- pred_near_atg(
        model_start = model_learned$model_orf$model_start,
        lst_nonstart = lst_nonstart
      )

      if (is.null(candi_near)) {
        rg_orf_pred <- IRanges::IRanges(
          start = rg_candi_atg@start,
          end = BiocGenerics::end(
            stop_pos_rg[
              !BiocGenerics::"%in%"(stop_pos_rg, rg_no_atg)
            ]
          )
        )
      } else {
        rg_candi_atg[candi_near$idx_posi] <- candi_near$candi_atg

        rg_orf_pred <- IRanges::IRanges(
          start = rg_candi_atg@start,
          end = BiocGenerics::end(stop_pos_rg[hit_start@to[candi_atg_rows]])
        )
      }

      if (length(lst_nonstart2$hit_start) == 0) {
        rg_orf_pred2 <- NULL
      } else {
        near_atg_prob2 <- predict(
          model_learned$model_orf$model_start$bst,
          lst_nonstart2$start_atg_info
        )

        probs <- as.numeric(near_atg_prob2)
        t_opt <- model_learned$model_orf$model_start$t_opt
        to_vec <- lst_nonstart2$hit_start@to # group id per row
        from_vec <- lst_nonstart2$hit_start@from # atg_start index per row

        # Lock group order to first appearance
        grp <- factor(to_vec, levels = unique(to_vec))
        rows <- seq_along(to_vec)

        # For each group: pick first ≥ t_opt; if none, pick argmax
        pick_row <- tapply(rows, grp, function(ix) {
          ix <- ix[!is.na(probs[ix])]
          if (length(ix) == 0L) {
            return(NA_integer_)
          }
          ok <- which(probs[ix] >= t_opt)
          if (length(ok) > 0L) {
            ix[ok[1]] # first that meets threshold
          } else {
            ix[which.max(probs[ix])] # fallback: max within group
          }
        })

        # Drop empty groups (defensive)
        pick_row <- unlist(pick_row, use.names = TRUE)
        keep <- !is.na(pick_row)
        pick_row <- pick_row[keep]

        # Map back to ranges
        picked_to_id <- as.integer(names(pick_row)) # group ids (must index rg_no_atg)
        picked_from <- from_vec[pick_row] # index into atg_start

        rg_candi_atg2 <- lst_nonstart2$atg_start[picked_from]

        rg_orf_pred2 <- IRanges::IRanges(
          start = BiocGenerics::start(rg_candi_atg2),
          end   = BiocGenerics::end(rg_no_atg[picked_to_id])
        )
      }
      rg_orf_pred <- c(rg_orf_pred, rg_orf_pred2)

      predict_score <- predict_score[
        match(
          BiocGenerics::end(rg_orf_pred),
          BiocGenerics::end(stop_pos_rg)
        )
      ]
    }
    return(list(
      pred_orf = rg_orf_pred,
      predict_score = predict_score
    ))
  }

if (skip_lnc_prediction) {
  lnc_orf_pred <- list(
    pred_orf = IRanges::IRanges(),
    predict_score = numeric(0)
  )
} else {
  lnc_orf_pred = predict_orfs_lnc(
    par_lnc = ncrna_rpf,
    tx_info= tx_info_lst$lncrna_info,
                                 ribo_num_pos=ribo_num_pos,
                                 is_mnase= rnase.bias,
                                 candidate_translat_reg = candidate_lnc_reg,
                                 model_learned = orf_pred,
                                 model_pos = model_pos,
                                 min_pos_num = 10,
                                 eff_ribo_num = 1
  )
}
  #  Return everything requested
  return(list(
    orf_pred      = orf_pred,
    lnc_orf_pred  = lnc_orf_pred
  ))
}
