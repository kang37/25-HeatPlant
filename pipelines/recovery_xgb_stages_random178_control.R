# 随机178站对照: 从全部站里随机抽178站, 比较CCM的178站是否更好/更差。参数同 recovery_xgb_shap_phi_stages_all.R (KB YEAR0 YEAR1)
src <- readLines("pipelines/recovery_xgb_shap_phi_stages_all.R", encoding = "UTF-8")
idx <- grep("^res178 <- run_set", src)[1]
invisible(capture.output(eval(parse(text = paste(src[1:(idx - 1)], collapse = "\n")))))
NR <- 20L
ctrl <- rbindlist(lapply(seq_len(NS), function(s) {
  dd <- build_stage(s, dt18$stat_id); feats <- attr(dd, "feats")
  out <- rbindlist(mclapply(seq_len(NR), function(i) {
    set.seed(7000 + 100 * s + i); sub <- dd[sample(.N, 178)]
    a <- one_rep(sub, feats, "random", r = i, baselines = FALSE); b <- one_rep(sub, feats, "spatial", r = i, baselines = FALSE)
    data.table(i = i, r2_random = a$r2_w, r2_spatial = b$r2_w) }, mc.cores = NCORE))
  out[, stage := s]; out }))
cc <- fread(sprintf("data_proc/output_loose_classify/recovery_xgb_stages_ccm178%s_cv_summary.csv", SUF))
res <- ctrl[, .(rand178_random_mean = mean(r2_random), rand178_random_q10 = quantile(r2_random, .1), rand178_random_q90 = quantile(r2_random, .9),
                rand178_spatial_mean = mean(r2_spatial), rand178_spatial_q10 = quantile(r2_spatial, .1), rand178_spatial_q90 = quantile(r2_spatial, .9)), by = stage]
res <- merge(res, cc[, .(stage, years, ccm178_random = r2w_random, ccm178_spatial = r2w_spatial)], by = "stage")
res[, pct_random := sapply(stage, function(s) mean(ctrl[stage == s]$r2_random <= cc[stage == s]$r2w_random))]
res[, pct_spatial := sapply(stage, function(s) mean(ctrl[stage == s]$r2_spatial <= cc[stage == s]$r2w_spatial))]
print(res[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(res, sprintf("data_proc/output_loose_classify/recovery_xgb_stages_random178_control%s.csv", SUF))
cat("完成\n")
