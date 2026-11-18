#!/usr/bin/env Rscript
# =============================================================================
# recovery_xgb_shap_phi_stages_all.R
#   不管因果关系: 把dt18里全部协变量齐全的904站纳入, 2001-2021年(21年), 每3年一个阶段
#   (2001-03,...,2019-21共7段), 每阶段一个XGBoost+SHAP模型; 同一套流程、同年份在178站(CCM判据(1))
#   子集上重跑一遍作对照, 以便把"是否做CCM筛选"的影响和"年份/阶段"的影响分开。
#   φ按阶段重新估计: 每站每阶段3年窗口合并AR1(年内lag), 去季节用2001-2021气候态。
#   DV加权 1/(SE²+τ²); X=阶段内均值动态6个+静态8个, invest缺失保留(阶段1只有2003年有invest)。
#   评估: 随机5折 + 空间分块5折(kmeans地理聚类整簇留出), ridge/RF基线, 置换零分布, 去invest消融,
#         全样本TreeSHAP, 站点bootstrap名次稳定性。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(xgboost); library(ggplot2)
                                library(glmnet); library(ranger); library(parallel)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")
NCORE <- max(1L, min(6L, detectCores() - 1L))
ARGS <- commandArgs(trailingOnly = TRUE)   # 可选: KB YEAR0 YEAR1, 默认 3 2001 2021
KB <- if (length(ARGS) >= 1) as.integer(ARGS[1]) else 3L
YEAR0 <- if (length(ARGS) >= 2) as.integer(ARGS[2]) else 2001L
YEAR1 <- if (length(ARGS) >= 3) as.integer(ARGS[3]) else 2021L
SUF <- if (KB == 3L && YEAR1 == 2021L) "" else sprintf("_k%d_%d_%d", KB, YEAR0, YEAR1)
REPS <- 5L; K <- 5L; NPERM <- 10L; NBOOT <- 60L
NS <- (YEAR1 - YEAR0 + 1L) %/% KB

env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]; cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18 <- env$dt
STA <- c("imperv", "grass", "water", "ever_needle", "ever_broad", "deci_broad", "elev", "rsds_sd")
DYN <- c("tavg", "rh", "precip", "rsds", "ntl", "invest")
FE_ALL <- c(DYN, STA)
stat_tab <- dt18[, c("stat_id", STA), with = FALSE]
co <- fread("data_raw/WS土地利用比例_1000m.csv", select = 1:3); setnames(co, c("stat_id", "lon", "lat"))
pc <- fread(file.path(OUT, "recovery_panel_covariates_all.csv"))
nanmean <- function(x) { x <- x[!is.na(x)]; if (length(x)) mean(x) else NA_real_ }

# ---- 178站清单(与之前一致) ----------------------------------------------------
neg <- fread("data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv")
neg[, drho := rho - rho_min]
LCB <- "data_proc/output_hcsif_buf1000/landcover_break_stations.csv"
n_break <- NA_integer_
if (file.exists(LCB)) { lcb <- fread(LCB); n_break <- sum(dt18$stat_id %in% lcb[has_break == TRUE]$stat_id); neg <- neg[!meteo_stat %in% lcb[has_break == TRUE]$stat_id] }
sr <- neg[tp >= 0 & tp <= 8]; sr[, pass1 := (drho > 0) & (p_surr < 0.05)]
sel178 <- unique(sr[pass1 == TRUE]$meteo_stat)
cat("全部站:", nrow(dt18), "| 178站:", length(sel178), "| 全部站中的地表覆盖断点站:", n_break, "\n")

# ---- 全部站: 阶段φ(3年窗口合并AR1) -------------------------------------------
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
d   <- rbindlist(lapply(fs, fread), fill = TRUE)[meteo_stat %in% dt18$stat_id & year >= YEAR0 & year <= YEAR1, .(meteo_stat, year, doy, SIF_buf1000)]
setorder(d, meteo_stat, year, doy)
deseason_z <- function(x, doy) { clim <- tapply(x, doy, mean, na.rm = TRUE); a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x))); (a - mean(a, na.rm = TRUE)) / s }
d[, sif_anom := deseason_z(SIF_buf1000, doy), by = meteo_stat]
d[, lag1 := shift(sif_anom, 1), by = .(meteo_stat, year)]
d <- d[!is.na(lag1) & !is.na(sif_anom)]
d[, stage := (year - YEAR0) %/% KB + 1L]
fit_blk <- function(sub) { if (nrow(sub) < 15) return(NULL)
  co_ <- summary(lm(sif_anom ~ lag1, data = sub))$coefficients
  list(phi = co_["lag1", "Estimate"], se = co_["lag1", "Std. Error"], n = nrow(sub)) }
blocks <- d[, fit_blk(.SD), by = .(meteo_stat, stage)]
cat("阶段φ: ", nrow(blocks), "个站-阶段; 站数", uniqueN(blocks$meteo_stat), "; 每阶段中位站数", median(blocks[, .N, by = stage]$N), "\n")

# ---- 工具 ----------------------------------------------------------------------
GRID <- expand.grid(max_depth = c(2, 4), min_child_weight = c(5, 10))
base_p <- list(objective = "reg:squarederror", eval_metric = "rmse", eta = 0.05, subsample = 0.8, colsample_bytree = 0.8, lambda = 5, nthread = 1)
rand_folds <- function(n, nf, seed) { set.seed(seed); sample(rep(seq_len(nf), length.out = n)) }
tune_xgb <- function(X, y, w, seed) {
  fid <- rand_folds(nrow(X), 5, seed); folds <- lapply(1:5, function(f) which(fid == f))
  dtr <- xgb.DMatrix(X, label = y, weight = w); best <- list(score = Inf)
  for (i in seq_len(nrow(GRID))) {
    p <- c(base_p, as.list(GRID[i, ]), list(seed = seed))
    cv <- xgb.cv(params = p, data = dtr, nrounds = 600, folds = folds, early_stopping_rounds = 30, verbose = 0)
    el <- as.data.frame(cv$evaluation_log); j <- which.min(el$test_rmse_mean)
    if (el$test_rmse_mean[j] < best$score) best <- list(score = el$test_rmse_mean[j], p = p, n = j)
  }
  best$fit <- xgb.train(params = best$p, data = dtr, nrounds = best$n); best
}
r2u <- function(y, p) 1 - sum((y - p)^2) / sum((y - mean(y))^2)
r2w <- function(y, p, w) { yb <- sum(w * y) / sum(w); 1 - sum(w * (y - p)^2) / sum(w * (y - yb)^2) }
impute <- function(Xtr, Xte) { med <- apply(Xtr, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(Xtr))) { Xtr[is.na(Xtr[, j]), j] <- med[j]; Xte[is.na(Xte[, j]), j] <- med[j] }
  list(tr = Xtr, te = Xte) }
make_folds <- function(dd, scheme, r) {
  if (scheme == "random") return(rand_folds(nrow(dd), K, 1000 + r))
  set.seed(2000 + r); kc <- if (nrow(dd) > 400) 40L else 20L
  km <- kmeans(scale(as.matrix(dd[, .(lon, lat)])), centers = kc, nstart = 5)
  cf <- sample(rep(seq_len(K), length.out = kc)); cf[km$cluster] }
one_rep <- function(dd, feats, scheme, r, baselines, y = dd$phi, w = dd$w) {
  X <- as.matrix(dd[, feats, with = FALSE]); n <- nrow(X); fid <- make_folds(dd, scheme, r)
  px <- pr <- pf <- rep(NA_real_, n)
  for (f in seq_len(K)) {
    tr <- fid != f; te <- !tr
    m <- tune_xgb(X[tr, , drop = FALSE], y[tr], w[tr], seed = 100 * r + f)
    px[te] <- predict(m$fit, X[te, , drop = FALSE])
    if (baselines) {
      im <- impute(X[tr, , drop = FALSE], X[te, , drop = FALSE])
      mu <- colMeans(im$tr); sdv <- apply(im$tr, 2, sd); sdv[sdv == 0] <- 1
      cvr <- cv.glmnet(scale(im$tr, mu, sdv), y[tr], alpha = 0, nfolds = 5, weights = w[tr], standardize = FALSE)
      pr[te] <- as.numeric(predict(cvr, scale(im$te, mu, sdv), s = "lambda.1se"))
      rf <- ranger(x = as.data.frame(im$tr), y = y[tr], num.trees = 300, case.weights = w[tr], seed = 1)
      pf[te] <- predict(rf, as.data.frame(im$te))$predictions
    }
  }
  data.table(rep = r, r2_u = r2u(y, px), r2_w = r2w(y, px, w), rmse = sqrt(mean((y - px)^2)),
             ridge_r2w = if (baselines) r2w(y, pr, w) else NA_real_, rf_r2w = if (baselines) r2w(y, pf, w) else NA_real_)
}
cv_eval <- function(dd, feats, scheme, baselines = TRUE)
  rbindlist(mclapply(seq_len(REPS), function(r) one_rep(dd, feats, scheme, r, baselines), mc.cores = NCORE))

build_stage <- function(s, stn) {
  yrs <- YEAR0 + (s - 1L) * KB + 0:(KB - 1L)
  b <- blocks[stage == s & meteo_stat %in% stn]
  dyn <- pc[year %in% yrs & stat_id %in% stn, lapply(.SD, nanmean), by = stat_id, .SDcols = DYN]
  dd <- merge(b, dyn, by.x = "meteo_stat", by.y = "stat_id")
  dd <- merge(dd, stat_tab, by.x = "meteo_stat", by.y = "stat_id")
  dd <- merge(dd, co, by.x = "meteo_stat", by.y = "stat_id")
  dd <- dd[is.finite(phi) & is.finite(se)]
  feats <- FE_ALL[sapply(FE_ALL, function(v) any(!is.na(dd[[v]])))]
  dd <- dd[complete.cases(dd[, setdiff(feats, "invest"), with = FALSE])]
  tau2 <- max(var(dd$phi) - mean(dd$se^2), 0.25 * var(dd$phi))
  dd[, w := 1 / (se^2 + tau2)]; dd[, w := w / mean(w)]
  attr(dd, "feats") <- feats; attr(dd, "years") <- paste(range(yrs), collapse = "-"); attr(dd, "noise_share") <- mean(dd$se^2) / var(dd$phi)
  dd
}

run_set <- function(stn, tag) {
  cat(sprintf("\n################ 站点集: %s (候选%d站) ################\n", tag, length(stn)))
  cv_tab <- list(); imp_all <- list(); shap_all <- list()
  for (s in seq_len(NS)) {
    dd <- build_stage(s, stn); feats <- attr(dd, "feats"); yrs <- attr(dd, "years"); hasinv <- "invest" %in% feats
    cat(sprintf("\n=== [%s] 阶段%d (%s): n=%d, 特征=%d, invest缺失=%s, SE噪声占比=%.2f ===\n", tag, s, yrs, nrow(dd), length(feats),
                if (hasinv) sum(is.na(dd$invest)) else "全缺", attr(dd, "noise_share")))
    rr <- cv_eval(dd, feats, "random"); rs <- cv_eval(dd, feats, "spatial")
    noinv <- if (hasinv) cv_eval(dd, setdiff(feats, "invest"), "random", FALSE) else NULL
    nullr <- unlist(mclapply(seq_len(NPERM), function(b) { set.seed(5000 + b); ix <- sample(nrow(dd))
      one_rep(dd, feats, "random", r = b, baselines = FALSE, y = dd$phi[ix], w = dd$w[ix])$r2_w }, mc.cores = NCORE))
    p_null <- (1 + sum(nullr >= mean(rr$r2_w))) / (NPERM + 1)
    cat(sprintf("随机5折 R2w=%.3f±%.3f | ridge=%.3f RF=%.3f\n空间分块 R2w=%.3f±%.3f | ridge=%.3f RF=%.3f\n置换max=%.3f p=%.3f%s\n",
                mean(rr$r2_w), sd(rr$r2_w), mean(rr$ridge_r2w), mean(rr$rf_r2w), mean(rs$r2_w), sd(rs$r2_w), mean(rs$ridge_r2w), mean(rs$rf_r2w),
                max(nullr), p_null, if (!is.null(noinv)) sprintf(" | 去invest R2w=%.3f", mean(noinv$r2_w)) else ""))
    cv_tab[[s]] <- data.table(set = tag, stage = s, years = yrs, n = nrow(dd), has_invest = hasinv, noise_share = attr(dd, "noise_share"),
      r2w_random = mean(rr$r2_w), r2w_random_sd = sd(rr$r2_w), r2w_spatial = mean(rs$r2_w), r2w_spatial_sd = sd(rs$r2_w),
      ridge_random = mean(rr$ridge_r2w), rf_random = mean(rr$rf_r2w), ridge_spatial = mean(rs$ridge_r2w), rf_spatial = mean(rs$rf_r2w),
      null_max = max(nullr), perm_p = p_null, noinvest_r2w_random = if (!is.null(noinv)) mean(noinv$r2_w) else NA_real_)

    X <- as.matrix(dd[, feats, with = FALSE]); y <- dd$phi
    m <- tune_xgb(X, y, dd$w, seed = 7)
    sh <- predict(m$fit, X, predcontrib = TRUE)[, feats, drop = FALSE]
    imp <- data.table(set = tag, stage = s, years = yrs, var = feats, mean_abs_shap = colMeans(abs(sh)))
    imp[, share := mean_abs_shap / sum(mean_abs_shap)]; imp[, rank := rank(-mean_abs_shap)]
    imp[, shap_vs_value_spearman := sapply(var, function(v) { ok <- !is.na(X[, v]); suppressWarnings(cor(X[ok, v], sh[ok, v], method = "spearman")) })]
    set.seed(99 + s)
    bt <- mclapply(seq_len(NBOOT), function(b) { idx <- sample(nrow(X), replace = TRUE)
      mb <- xgb.train(params = c(m$p, list(seed = b)), data = xgb.DMatrix(X[idx, , drop = FALSE], label = y[idx], weight = dd$w[idx]), nrounds = m$n)
      colMeans(abs(predict(mb, X, predcontrib = TRUE)[, feats, drop = FALSE])) }, mc.cores = NCORE)
    bt <- do.call(rbind, bt); rk <- t(apply(-bt, 1, rank))
    stab <- data.table(var = feats, rank_q10 = apply(rk, 2, quantile, .1), rank_q90 = apply(rk, 2, quantile, .9), top5_freq = colMeans(rk <= 5))
    imp <- merge(imp, stab, by = "var")[order(rank)]
    print(imp[, .(var, share = round(share, 3), rank, spearman = round(shap_vs_value_spearman, 2), rank_q10, rank_q90, top5_freq)])
    imp_all[[s]] <- imp
    shap_all[[s]] <- data.table(set = tag, stage = s, years = yrs, meteo_stat = rep(dd$meteo_stat, length(feats)), var = rep(feats, each = nrow(dd)),
                                shap = as.vector(sh), value = as.vector(X))
  }
  cvt <- rbindlist(cv_tab); imp <- rbindlist(imp_all); shp <- rbindlist(shap_all)
  fwrite(cvt, file.path(OUT, sprintf("recovery_xgb_stages_%s_cv_summary.csv", tag)))
  fwrite(imp, file.path(OUT, sprintf("recovery_xgb_stages_%s_importance.csv", tag)))
  fwrite(shp, file.path(OUT, sprintf("recovery_xgb_stages_%s_shap_long.csv", tag)))
  summ <- imp[, .(mean_share = mean(share), min_share = min(share), max_share = max(share), n_top5 = sum(rank <= 5), n_stage = .N,
                  n_pos = sum(shap_vs_value_spearman > 0.3, na.rm = TRUE), n_neg = sum(shap_vs_value_spearman < -0.3, na.rm = TRUE),
                  mean_spearman = mean(shap_vs_value_spearman, na.rm = TRUE)), by = var][order(-mean_share)]
  fwrite(summ, file.path(OUT, sprintf("recovery_xgb_stages_%s_cross_stage_summary.csv", tag)))
  cat(sprintf("\n=== [%s] 各阶段CV汇总 ===\n", tag))
  print(cvt[, .(stage, years, n, noise = round(noise_share, 2), r2w_random = round(r2w_random, 3), r2w_spatial = round(r2w_spatial, 3),
                rf_random = round(rf_random, 3), ridge_spatial = round(ridge_spatial, 3), perm_p = round(perm_p, 3), noinvest = round(noinvest_r2w_random, 3))])
  cat(sprintf("\n=== [%s] 跨阶段稳定性 ===\n", tag)); print(summ[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

  lab_s <- function(s) { y <- YEAR0 + (s - 1L) * KB; sprintf("%d-%02d", y, (y + KB - 1L) %% 100) }
  lv <- sapply(seq_len(NS), lab_s)
  imp[, stage_lab := factor(lab_s(stage), levels = lv)]; ord <- summ$var; imp[, var := factor(var, levels = rev(ord))]
  p1 <- ggplot(imp, aes(stage_lab, var, fill = share)) + geom_tile(color = "white") + geom_text(aes(label = rank), size = 3) +
    scale_fill_gradient(low = "#F7FBFF", high = "#C1121F", name = "Share of\nmean |SHAP|") +
    labs(title = sprintf("phi drivers by 3-year stage [%s]: SHAP share, number = rank", tag), x = "Stage", y = NULL) + theme_bw(base_size = 11)
  ggsave(file.path(OUT, sprintf("recovery_xgb_stages_%s_share_heatmap.png", tag)), p1, width = 9, height = 6, dpi = 300)
  p2 <- ggplot(imp, aes(stage_lab, var, fill = shap_vs_value_spearman)) + geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", shap_vs_value_spearman)), size = 2.8) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#C1121F", midpoint = 0, limits = c(-1, 1), name = "Spearman\n(value, SHAP)") +
    labs(title = sprintf("Direction by stage [%s] (red: higher value -> higher phi = slower recovery)", tag), x = "Stage", y = NULL) + theme_bw(base_size = 11)
  ggsave(file.path(OUT, sprintf("recovery_xgb_stages_%s_direction_heatmap.png", tag)), p2, width = 9, height = 6, dpi = 300)
  shp[, stage_lab := factor(lab_s(stage), levels = lv)]
  iv <- shp[var == "invest" & !is.na(value)]
  if (nrow(iv)) {
    p4 <- ggplot(iv, aes(value, shap)) + geom_hline(yintercept = 0, color = "grey70") + geom_point(alpha = 0.3, size = 0.7, color = "#2166AC") +
      geom_smooth(method = "loess", se = FALSE, color = "#C1121F", linewidth = 0.8, formula = y ~ x) + facet_wrap(~stage_lab, ncol = 3) +
      labs(title = sprintf("invest SHAP dependence by stage [%s] (shared y-axis; 1 dot = 1 station)", tag), x = "invest (log1p, block mean)", y = "SHAP value") + theme_bw(base_size = 10)
    ggsave(file.path(OUT, sprintf("recovery_xgb_stages_%s_invest_dependence.png", tag)), p4, width = 10, height = 6.5, dpi = 300)
  }
  list(cv = cvt, imp = imp, summ = summ)
}

res178 <- run_set(sel178, paste0("ccm178", SUF))
resall <- run_set(dt18$stat_id, paste0("all904", SUF))

cat("\n################ 对照: 178站 vs 全部904站 (同年份、同流程) ################\n")
cmp <- merge(res178$cv[, .(stage, years, n178 = n, r2w_random_178 = r2w_random, r2w_spatial_178 = r2w_spatial, noise178 = noise_share)],
             resall$cv[, .(stage, n_all = n, r2w_random_all = r2w_random, r2w_spatial_all = r2w_spatial, noise_all = noise_share,
                           rf_random_all = rf_random, ridge_spatial_all = ridge_spatial, perm_p_all = perm_p, noinvest_all = noinvest_r2w_random)], by = "stage")
print(cmp[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(cmp, file.path(OUT, sprintf("recovery_xgb_stages_all_vs_178_cv%s.csv", SUF)))
sc <- merge(res178$summ[, .(var, share178 = mean_share, top5_178 = n_top5, spear178 = mean_spearman)],
            resall$summ[, .(var, share_all = mean_share, top5_all = n_top5, pos_all = n_pos, neg_all = n_neg, spear_all = mean_spearman)], by = "var")[order(-share_all)]
print(sc[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(sc, file.path(OUT, sprintf("recovery_xgb_stages_all_vs_178_shap%s.csv", SUF)))
cat("\n完成。\n")
