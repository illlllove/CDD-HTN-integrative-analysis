library(data.table)
library(coloc)

harm_to_eqtl <-
function(
        eq,
        gw,
        gwas_name = "GWAS",
        freq_tol = 0.20
) {
    
    eq <- data.table::copy(eq)
    gw <- data.table::copy(gw)
    
    eq[, `:=`(
        A1 = toupper(A1),
        A2 = toupper(A2)
    )]
    
    gw[, `:=`(
        A1 = toupper(A1),
        A2 = toupper(A2)
    )]
    
    data.table::setnames(
        gw,
        old = c(
            "A1", "A2", "freq",
            "b", "se", "p"
        ),
        new = c(
            "GW_A1", "GW_A2", "GW_freq",
            "GW_b", "GW_se", "GW_p"
        )
    )
    
    m <- merge(
        eq,
        gw,
        by = "SNP",
        all = FALSE
    )
    
    m[, palindromic :=
          (A1 == "A" & A2 == "T") |
          (A1 == "T" & A2 == "A") |
          (A1 == "C" & A2 == "G") |
          (A1 == "G" & A2 == "C")
    ]
    
    m[, same :=
          A1 == GW_A1 &
          A2 == GW_A2
    ]
    
    m[, swapped :=
          A1 == GW_A2 &
          A2 == GW_A1
    ]
    
    m[, mismatch :=
          !palindromic &
          !same &
          !swapped
    ]
    
    # 先初始化，避免子集赋值时出现不稳定行为
    m[, `:=`(
        GW_b_aligned = NA_real_,
        GW_freq_aligned = NA_real_
    )]
    
    # allele方向完全一致
    m[
        same == TRUE,
        `:=`(
            GW_b_aligned = GW_b,
            GW_freq_aligned = GW_freq
        )
    ]
    
    # allele方向相反
    m[
        swapped == TRUE,
        `:=`(
            GW_b_aligned = -GW_b,
            GW_freq_aligned = 1 - GW_freq
        )
    ]
    
    m[, freq_diff :=
          abs(Freq - GW_freq_aligned)
    ]
    
    m[, freq_fail :=
          !palindromic &
          !mismatch &
          (
              !is.finite(Freq) |
                  !is.finite(GW_freq_aligned) |
                  freq_diff > freq_tol
          )
    ]
    
    keep <- m[
        palindromic == FALSE &
            mismatch == FALSE &
            freq_fail == FALSE &
            is.finite(b) &
            is.finite(SE) &
            SE > 0 &
            is.finite(GW_b_aligned) &
            is.finite(GW_se) &
            GW_se > 0
    ]
    
    qc <- data.table(
        GWAS = gwas_name,
        overlap = nrow(m),
        palindromic = sum(m$palindromic, na.rm = TRUE),
        mismatch = sum(m$mismatch, na.rm = TRUE),
        freq_fail = sum(m$freq_fail, na.rm = TRUE),
        kept = nrow(keep)
    )
    
    list(
        data = keep,
        qc = qc
    )
}
quiet_coloc_pair <-
function(
        dat,
        p12,
        pair_name
) {
    
    d_eqtl <- list(
        beta = dat$b,
        varbeta = dat$SE^2,
        snp = dat$SNP,
        type = "quant",
        sdY = 1
    )
    
    d_gwas <- list(
        beta = dat$GW_b_aligned,
        varbeta = dat$GW_se^2,
        snp = dat$SNP,
        type = "cc"
    )
    
    invisible(
        capture.output(
            res <- coloc::coloc.abf(
                dataset1 = d_eqtl,
                dataset2 = d_gwas,
                p1 = 1e-4,
                p2 = 1e-4,
                p12 = p12
            )
        )
    )
    
    s <- res$summary
    
    list(
        summary = data.table(
            pair = pair_name,
            p12 = p12,
            nsnps = as.integer(s["nsnps"]),
            H0 = as.numeric(s["PP.H0.abf"]),
            H1 = as.numeric(s["PP.H1.abf"]),
            H2 = as.numeric(s["PP.H2.abf"]),
            H3 = as.numeric(s["PP.H3.abf"]),
            H4 = as.numeric(s["PP.H4.abf"])
        ),
        result = res
    )
}
run_pair_priors <-
function(
        dat,
        pair_name,
        priors = c(1e-6, 1e-5, 5e-5)
) {
    
    fits <- lapply(
        priors,
        function(x) {
            quiet_coloc_pair(
                dat = dat,
                p12 = x,
                pair_name = pair_name
            )
        }
    )
    
    list(
        summary = rbindlist(
            lapply(
                fits,
                function(x) x$summary
            )
        ),
        fits = fits
    )
}
harm_disease_pair <-
function(
        htn,
        cdd,
        freq_tol = 0.20
) {
    
    h <- copy(htn)
    c <- copy(cdd)
    
    h[, `:=`(
        A1 = toupper(A1),
        A2 = toupper(A2)
    )]
    
    c[, `:=`(
        A1 = toupper(A1),
        A2 = toupper(A2)
    )]
    
    setnames(
        c,
        old = c(
            "A1", "A2", "freq",
            "b", "se", "p"
        ),
        new = c(
            "CDD_A1", "CDD_A2", "CDD_freq",
            "CDD_b", "CDD_se", "CDD_p"
        )
    )
    
    m <- merge(
        h,
        c,
        by = "SNP",
        all = FALSE
    )
    
    m[, palindromic :=
          (A1 == "A" & A2 == "T") |
          (A1 == "T" & A2 == "A") |
          (A1 == "C" & A2 == "G") |
          (A1 == "G" & A2 == "C")
    ]
    
    m[, same :=
          A1 == CDD_A1 &
          A2 == CDD_A2
    ]
    
    m[, swapped :=
          A1 == CDD_A2 &
          A2 == CDD_A1
    ]
    
    m[, mismatch :=
          palindromic == FALSE &
          same == FALSE &
          swapped == FALSE
    ]
    
    m[, `:=`(
        CDD_b_aligned = NA_real_,
        CDD_freq_aligned = NA_real_
    )]
    
    m[
        same == TRUE,
        `:=`(
            CDD_b_aligned = CDD_b,
            CDD_freq_aligned = CDD_freq
        )
    ]
    
    m[
        swapped == TRUE,
        `:=`(
            CDD_b_aligned = -CDD_b,
            CDD_freq_aligned = 1 - CDD_freq
        )
    ]
    
    m[, freq_diff :=
          abs(freq - CDD_freq_aligned)
    ]
    
    m[, freq_fail :=
          palindromic == FALSE &
          mismatch == FALSE &
          (
              !is.finite(freq) |
                  !is.finite(CDD_freq_aligned) |
                  freq_diff > freq_tol
          )
    ]
    
    keep <- m[
        palindromic == FALSE &
            mismatch == FALSE &
            freq_fail == FALSE &
            is.finite(b) &
            is.finite(se) &
            se > 0 &
            is.finite(CDD_b_aligned) &
            is.finite(CDD_se) &
            CDD_se > 0
    ]
    
    qc <- data.table(
        HTN_region_raw = nrow(h),
        CDD_overlap = nrow(m),
        palindromic = sum(m$palindromic, na.rm = TRUE),
        mismatch = sum(m$mismatch, na.rm = TRUE),
        freq_fail = sum(m$freq_fail, na.rm = TRUE),
        kept = nrow(keep)
    )
    
    list(
        data = keep,
        qc = qc
    )
}
quiet_coloc_dd <-
function(
        dat,
        p12
) {
    
    d_htn <- list(
        beta = dat$b,
        varbeta = dat$se^2,
        snp = dat$SNP,
        type = "cc"
    )
    
    d_cdd <- list(
        beta = dat$CDD_b_aligned,
        varbeta = dat$CDD_se^2,
        snp = dat$SNP,
        type = "cc"
    )
    
    invisible(
        capture.output(
            res <- coloc::coloc.abf(
                dataset1 = d_htn,
                dataset2 = d_cdd,
                p1 = 1e-4,
                p2 = 1e-4,
                p12 = p12
            )
        )
    )
    
    s <- res$summary
    
    list(
        summary = data.table(
            pair = "HTN_CDD",
            p12 = p12,
            nsnps = as.integer(s["nsnps"]),
            H0 = as.numeric(s["PP.H0.abf"]),
            H1 = as.numeric(s["PP.H1.abf"]),
            H2 = as.numeric(s["PP.H2.abf"]),
            H3 = as.numeric(s["PP.H3.abf"]),
            H4 = as.numeric(s["PP.H4.abf"])
        ),
        result = res
    )
}
run_dd_priors <-
function(
        dat,
        priors = c(1e-6, 1e-5, 5e-5)
) {
    
    fits <- lapply(
        priors,
        function(x) {
            quiet_coloc_dd(
                dat = dat,
                p12 = x
            )
        }
    )
    
    list(
        summary = rbindlist(
            lapply(
                fits,
                function(x) x$summary
            )
        ),
        fits = fits
    )
}
formal_probe_analysis <-
function(
        probe_id,
        freeze_table,
        eqtl_data,
        htn_pair_data,
        cdd_pair_data,
        htn_full,
        cdd_full,
        priors = c(1e-6, 1e-5, 5e-5)
) {
    
    info <- freeze_table[
        probeID == probe_id
    ]
    
    stopifnot(nrow(info) == 1)
    
    chr_now <- as.integer(info$ProbeChr)
    center_now <- as.integer(info$Probe_bp)
    
    region_start <- center_now - 1000000L
    region_end   <- center_now + 1000000L
    
    # -------------------------
    # 1. eQTL
    # -------------------------
    
    eq <- eqtl_data[
        Probe == probe_id
    ]
    
    if (nrow(eq) == 0) {
        stop(
            paste(
                "No eQTL records:",
                probe_id
            )
        )
    }
    
    # -------------------------
    # 2. eQTL–HTN
    # -------------------------
    
    harm_eh <- harm_to_eqtl(
        eq = eq,
        gw = htn_pair_data,
        gwas_name = "HTN"
    )
    
    fit_eh <- run_pair_priors(
        dat = harm_eh$data,
        pair_name = "EQTL_HTN",
        priors = priors
    )
    
    # -------------------------
    # 3. eQTL–CDD
    # -------------------------
    
    harm_ec <- harm_to_eqtl(
        eq = eq,
        gw = cdd_pair_data,
        gwas_name = "CDD"
    )
    
    fit_ec <- run_pair_priors(
        dat = harm_ec$data,
        pair_name = "EQTL_CDD",
        priors = priors
    )
    
    # -------------------------
    # 4. HTN ±1 Mb region
    # -------------------------
    
    htn_region_raw <- htn_full[
        chromosome == chr_now &
            base_pair_location >= region_start &
            base_pair_location <= region_end
    ]
    
    htn_region <- htn_region_raw[
        ,
        .(
            SNP  = variant_id,
            A1   = toupper(effect_allele),
            A2   = toupper(other_allele),
            freq = effect_allele_frequency,
            b    = beta,
            se   = standard_error,
            p    = p_value,
            chr  = chromosome,
            pos  = base_pair_location
        )
    ]
    
    # -------------------------
    # 5. HTN–CDD
    # -------------------------
    
    harm_dd <- harm_disease_pair(
        htn = htn_region,
        cdd = cdd_full
    )
    
    fit_dd <- run_dd_priors(
        dat = harm_dd$data,
        priors = priors
    )
    
    # -------------------------
    # 6. summary
    # -------------------------
    
    summary_all <- rbindlist(
        list(
            fit_eh$summary,
            fit_ec$summary,
            fit_dd$summary
        )
    )
    
    summary_all[
        ,
        `:=`(
            probeID = probe_id,
            ProbeChr = chr_now,
            Probe_bp = center_now,
            region_start = region_start,
            region_end = region_end
        )
    ]
    
    setcolorder(
        summary_all,
        c(
            "probeID",
            "ProbeChr",
            "Probe_bp",
            "region_start",
            "region_end",
            "pair",
            "p12",
            "nsnps",
            "H0",
            "H1",
            "H2",
            "H3",
            "H4"
        )
    )
    
    # -------------------------
    # 7. QC table
    # -------------------------
    
    qc_all <- rbindlist(
        list(
            
            data.table(
                probeID = probe_id,
                pair = "EQTL_HTN",
                input_or_region = nrow(eq),
                overlap = harm_eh$qc$overlap,
                palindromic = harm_eh$qc$palindromic,
                mismatch = harm_eh$qc$mismatch,
                freq_fail = harm_eh$qc$freq_fail,
                kept = harm_eh$qc$kept
            ),
            
            data.table(
                probeID = probe_id,
                pair = "EQTL_CDD",
                input_or_region = nrow(eq),
                overlap = harm_ec$qc$overlap,
                palindromic = harm_ec$qc$palindromic,
                mismatch = harm_ec$qc$mismatch,
                freq_fail = harm_ec$qc$freq_fail,
                kept = harm_ec$qc$kept
            ),
            
            data.table(
                probeID = probe_id,
                pair = "HTN_CDD",
                input_or_region =
                    harm_dd$qc$HTN_region_raw,
                overlap =
                    harm_dd$qc$CDD_overlap,
                palindromic =
                    harm_dd$qc$palindromic,
                mismatch =
                    harm_dd$qc$mismatch,
                freq_fail =
                    harm_dd$qc$freq_fail,
                kept =
                    harm_dd$qc$kept
            )
        )
    )
    
    list(
        info = info,
        summary = summary_all,
        qc = qc_all,
        
        harmonized = list(
            EQTL_HTN = harm_eh$data,
            EQTL_CDD = harm_ec$data,
            HTN_CDD = harm_dd$data
        ),
        
        fits = list(
            EQTL_HTN = fit_eh,
            EQTL_CDD = fit_ec,
            HTN_CDD = fit_dd
        )
    )
}
