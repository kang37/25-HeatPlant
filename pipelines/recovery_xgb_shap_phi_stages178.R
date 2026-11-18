#!/usr/bin/env Rscript
# =============================================================================
# recovery_xgb_shap_phi_stages178.R
#   参考 Yuan et al. 2024 (J. Cleaner Prod. 476:143640, XGBoost-SHAP, 分时期分别建模再比较SHAP):
#   178站(CCM判据(1)), 每3年一个阶段(2000-02,...,2018-20; 2021-22不足3年丢弃), 每阶段一个模型:
#     DV = 该阶段该站的φ(3年窗口合并AR1, 见recovery_phi_block_pooling.R), 按1/(SE²+τ²)加权
#     X  = 阶段内平均的动态协变量(tavg/rh/precip/rsds/ntl/invest) + 静态8个(含invest, 缺失保留)
#   阶段1(2000-02)没有任何invest数据(源数据始于2003), 该阶段只能不含invest。
#   每阶段: 嵌套重复CV(随机5折+空间分块5折, 内层调参), ridge/RF基线, 置换零分布,
#   全样本TreeSHAP, 站点bootstrap名次稳定性; 最后跨阶段比较SHAP份额/名次/方向。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(xgboost); library(ggplot2)
                                library(glmnet); library(ranger); library(parallel)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")
NCORE <- max(1L, min(6L, detectCores() - 1L))
REPS <- 10L; K <- 5L; NPERM <- 20L; NBOOT <- 100L; KB <- 3L; YEAR0 <- 2000L; NYEAR <- 23L
NS <- NYEAR %/% KB

blocks <- fread(file.path(OUT, "recovery_phi_block_pooling_blocks.csv"))[k == KB]
pc <- fread(file.path(OUT, "recovery_panel_covariates_178.csv"))
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]; cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
rsds_sd <- env$dt[, .(stat_id, rsds_sd)]
co <- fread("data_raw/WS土地利用比例_1000m.csv", select = 1:3); setnames(co, c("stat_id", "lon", "lat"))
DYN <- c("tavg", "rh", "precip", "rsds", "ntl", "invest")
STA <- c("imperv", "grass", "water", "ever_needle", "ever_broad", "deci_broad", "elev", "rsds_sd")
FE_ALL <- c(DYN, STA)
nanmean <- function(x) { x <- x[!is.na(x)]; if (length(x)) mean(x) else NA_real_ }

build_stage <- function(s) {
  yrs <- YEAR0 + (s - 1L) * KB + 0:(KB - 1L)
  b <- blocks[blk == s]
  dyn <- pc[year %in% yrs, lapply(.SD, nanmean), by = stat_id, .SDcols = DYN]
  sta <- pc[, lapply(.SD, function(x) x[1]), by = stat_id, .SDcols = setdiff(STA, "rsds_sd")]
  d <- merge(b, dyn, by.x = "meteo_stat", by.y = "stat_id")
  d <- merge(d, sta, by.x = "meteo_stat", by.y = "stat_id")
  d <- merge(d, rsds_sd, by.x = "meteo_stat", by.y = "stat_id")
  d <- merge(d, co, by.x = "meteo_stat", by.y = "stat_id")
  d <- d[is.finite(phi) & is.finite(se)]
  feats <- FE_ALL[sapply(FE_ALL, function(v) any(!is.na(d[[v]])))]
  d <- d[complete.cases(d[, setdiff(feats, "invest"), with = FALSE])]
  tau2 <- max(var(d$phi) - mean(d$se^2), 0.25 * var(d$phi))
  d[, w := 1 / (se^2 + tau2)]; d[, w := w / mean(w)]
  attr(d, "feats") <- feats; attr(d, "years") <- paste(range(yrs), collapse = "-")
  attr(d, "noise_share") <- mean(d$se^2) / var(d$phi)
  d
}

GRID <- expand.grid(max_depth = c(2, 3), min_child_weight = c(3, 6))
base_p <- list(objective = "reg:squarederror", eval_metric = "rmse", eta = 0.03,
               subsample = 0.8, colsample_bytree = 0.8, lambda = 5, nthread = 1)
rand_folds <- function(n, nf, seed) { set.seed(seed); sample(rep(seq_len(nf), length.out = n)) }
tune_xgb <- function(X, y, w, seed) {
  fid <- rand_folds(nrow(X), 5, seed); folds <- lapply(1:5, function(f) which(fid == f))
  dtr <- xgb.DMatrix(X, label = y, weight = w); best <- list(score = Inf)
  for (i in seq_len(nrow(GRID))) {
    p <- c(base_p, as.list(GRID[i, ]), list(seed = seed))
    cv <- xgb.cv(params = p, data = dtr, nrounds = 800, folds = folds, early_stopping_rounds = 30, verbose = 0)
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
  set.seed(2000 + r); kc <- 20L
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
      rf <- ranger(x = as.data.frame(im$tr), y = y[tr], num.trees = 500, case.weights = w[tr], seed = 1)
      pf[te] <- predict(rf, as.data.frame(im$te))$predictions
    }
  }
  data.table(rep = r, r2_u = r2u(y, px), r2_w = r2w(y, px, w), rmse = sqrt(mean((y - px)^2)),
             ridge_r2w = if (baselines) r2w(y, pr, w) else NA_real_, rf_r2w = if (baselines) r2w(y, pf, w) else NA_real_)
}
cv_eval <- function(dd, feats, scheme, baselines = TRUE)
  rbindlist(mclapply(seq_len(REPS), function(r) one_rep(dd, feats, scheme, r, baselines), mc.cores = NCORE))

cv_tab <- list(); imp_all <- list(); stab_all <- list(); shap_all <- list()
for (s in seq_len(NS)) {
  dd <- build_stage(s); feats <- attr(dd, "feats"); yrs <- attr(dd, "years")
  cat(sprintf("\n=========== 阶段%d (%s): n=%d, 特征=%d%s, invest缺失=%d, SE噪声占比=%.2f ===========\n", s, yrs, nrow(dd), length(feats),
              if (!"invest" %in% feats) " (无invest数据)" else "", if ("invest" %in% feats) sum(is.na(dd$invest)) else NA_integer_, attr(dd, "noise_share")))
  rr <- cv_eval(dd, feats, "random"); rs <- cv_eval(dd, feats, "spatial")
  noinv <- if ("invest" %in% feats) cv_eval(dd, setdiff(feats, "invest"), "random", FALSE) else NULL
  nullr <- unlist(mclapply(seq_len(NPERM), function(b) { set.seed(5000 + b); ix <- sample(nrow(dd))
    one_rep(dd, feats, "random", r = b, baselines = FALSE, y = dd$phi[ix], w = dd$w[ix])$r2_w }, mc.cores = NCORE))
  p_null <- (1 + sum(nullr >= mean(rr$r2_w))) / (NPERM + 1)
  cat(sprintf("随机5折  R2w=%.3f±%.3f R2u=%.3f | ridge=%.3f RF=%.3f\n空间分块 R2w=%.3f±%.3f | ridge=%.3f RF=%.3f\n置换零分布 max=%.3f p=%.3f%s\n",
              mean(rr$r2_w), sd(rr$r2_w), mean(rr$r2_u), mean(rr$ridge_r2w), mean(rr$rf_r2w),
              mean(rs$r2_w), sd(rs$r2_w), mean(rs$ridge_r2w), mean(rs$rf_r2w), max(nullr), p_null,
              if (!is.null(noinv)) sprintf(" | 去invest R2w=%.3f", mean(noinv$r2_w)) else ""))
  cv_tab[[s]] <- data.table(stage = s, years = yrs, n = nrow(dd), n_feats = length(feats), has_invest = "invest" %in% feats,
    noise_share = attr(dd, "noise_share"), r2w_random = mean(rr$r2_w), r2w_random_sd = sd(rr$r2_w), r2u_random = mean(rr$r2_u),
    r2w_spatial = mean(rs$r2_w), r2w_spatial_sd = sd(rs$r2_w), ridge_random = mean(rr$ridge_r2w), rf_random = mean(rr$rf_r2w),
    ridge_spatial = mean(rs$ridge_r2w), rf_spatial = mean(rs$rf_r2w), null_max = max(nullr), perm_p = p_null,
    noinvest_r2w_random = if (!is.null(noinv)) mean(noinv$r2_w) else NA_real_)

  X <- as.matrix(dd[, feats, with = FALSE]); y <- dd$phi
  m <- tune_xgb(X, y, dd$w, seed = 7)
  sh <- predict(m$fit, X, predcontrib = TRUE)[, feats, drop = FALSE]
  imp <- data.table(stage = s, years = yrs, var = feats, mean_abs_shap = colMeans(abs(sh)))
  imp[, share := mean_abs_shap / sum(mean_abs_shap)]; imp[, rank := rank(-mean_abs_shap)]
  imp[, shap_vs_value_spearman := sapply(var, function(v) { ok <- !is.na(X[, v]); suppressWarnings(cor(X[ok, v], sh[ok, v], method = "spearman")) })]
  set.seed(99 + s)
  bt <- mclapply(seq_len(NBOOT), function(b) { idx <- sample(nrow(X), replace = TRUE)
    mb <- xgb.train(params = c(m$p, list(seed = b)), data = xgb.DMatrix(X[idx, , drop = FALSE], label = y[idx], weight = dd$w[idx]), nrounds = m$n)
    colMeans(abs(predict(mb, X, predcontrib = TRUE)[, feats, drop = FALSE])) }, mc.cores = NCORE)
  bt <- do.call(rbind, bt); rk <- t(apply(-bt, 1, rank))
  stab <- data.table(stage = s, years = yrs, var = feats, rank_median = apply(rk, 2, median),
                     rank_q10 = apply(rk, 2, quantile, .1), rank_q90 = apply(rk, 2, quantile, .9), top5_freq = colMeans(rk <= 5))
  imp <- merge(imp, stab[, .(var, rank_q10, rank_q90, top5_freq)], by = "var")[order(rank)]
  print(imp[, .(var, share = round(share, 3), rank, spearman = round(shap_vs_value_spearman, 2), rank_q10, rank_q90, top5_freq)])
  imp_all[[s]] <- imp
  shap_all[[s]] <- data.table(stage = s, years = yrs, meteo_stat = rep(dd$meteo_stat, length(feats)), var = rep(feats, each = nrow(dd)),
                              shap = as.vector(sh), value = as.vector(X))
}
cvt <- rbindlist(cv_tab); imp <- rbindlist(imp_all); shp <- rbindlist(shap_all)
fwrite(cvt, file.path(OUT, "recovery_xgb_stages178_cv_summary.csv"))
fwrite(imp, file.path(OUT, "recovery_xgb_stages178_importance.csv"))
fwrite(shp, file.path(OUT, "recovery_xgb_stages178_shap_long.csv"))
cat("\n=== 各阶段CV汇总 ===\n"); print(cvt[, .(stage, years, n, has_invest, noise_share = round(noise_share, 2), r2w_random = round(r2w_random, 3),
  r2w_spatial = round(r2w_spatial, 3), rf_random = round(rf_random, 3), ridge_spatial = round(ridge_spatial, 3), perm_p = round(perm_p, 3),
  noinvest_r2w_random = round(noinvest_r2w_random, 3))])

# ---- 跨阶段比较 ----------------------------------------------------------------
summ <- imp[, .(mean_share = mean(share), min_share = min(share), max_share = max(share), n_top5 = sum(rank <= 5), n_stage = .N,
                n_pos = sum(shap_vs_value_spearman > 0.3), n_neg = sum(shap_vs_value_spearman < -0.3),
                mean_spearman = mean(shap_vs_value_spearman)), by = var][order(-mean_share)]
cat("\n=== 跨阶段稳定性(阶段间) ===\n"); print(summ[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(summ, file.path(OUT, "recovery_xgb_stages178_cross_stage_summary.csv"))

lab_s <- function(s) { y <- YEAR0 + (s - 1L) * KB; sprintf("%d-%02d", y, (y + KB - 1L) %% 100) }
imp[, stage_lab := factor(lab_s(stage), levels = sapply(seq_len(NS), lab_s))]
ord <- summ$var; imp[, var := factor(var, levels = rev(ord))]
p1 <- ggplot(imp, aes(stage_lab, var, fill = share)) + geom_tile(color = "white") +
  geom_text(aes(label = rank), size = 3, color = "black") +
  scale_fill_gradient(low = "#F7FBFF", high = "#C1121F", name = "Share of\nmean |SHAP|", na.value = "grey90") +
  labs(title = "phi drivers by 3-year stage (178 stations): SHAP share, number = rank", x = "Stage", y = NULL,
       subtitle = "blank = variable not available in that stage (no invest data before 2003)") + theme_bw(base_size = 11)
ggsave(file.path(OUT, "recovery_xgb_stages178_share_heatmap.png"), p1, width = 9, height = 6, dpi = 300)
p2 <- ggplot(imp, aes(stage_lab, var, fill = shap_vs_value_spearman)) + geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", shap_vs_value_spearman)), size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#C1121F", midpoint = 0, limits = c(-1, 1), name = "Spearman\n(value, SHAP)") +
  labs(title = "Direction of effect by stage (red: higher value -> higher phi = slower recovery)", x = "Stage", y = NULL) + theme_bw(base_size = 11)
ggsave(file.path(OUT, "recovery_xgb_stages178_direction_heatmap.png"), p2, width = 9, height = 6, dpi = 300)

shp[, stage_lab := factor(lab_s(stage), levels = sapply(seq_len(NS), lab_s))]
top_vars <- as.character(head(ord, 6))
sb <- shp[var %in% top_vars]; sb[, var := factor(var, levels = rev(top_vars))]
sb[, vs := { rg <- range(value, na.rm = TRUE); if (diff(rg) == 0) rep(0.5, .N) else (value - rg[1]) / diff(rg) }, by = .(stage, var)]
p3 <- ggplot(sb, aes(shap, var, color = vs)) + geom_vline(xintercept = 0, color = "grey70") + geom_jitter(height = 0.2, size = 0.7, alpha = 0.7) +
  scale_color_gradient(low = "#2166AC", high = "#C1121F", na.value = "grey60", name = "Value") + facet_wrap(~stage_lab, ncol = 4) +
  labs(title = "SHAP beeswarm by stage (top-6 variables overall)", x = "SHAP value", y = NULL) + theme_bw(base_size = 9)
ggsave(file.path(OUT, "recovery_xgb_stages178_beeswarm_top6.png"), p3, width = 13, height = 6.5, dpi = 300)
iv <- shp[var == "invest" & !is.na(value)]
if (nrow(iv)) {
  p4 <- ggplot(iv, aes(value, shap)) + geom_hline(yintercept = 0, color = "grey70") + geom_point(alpha = 0.6, size = 1, color = "#2166AC") +
    geom_smooth(method = "loess", se = FALSE, color = "#C1121F", linewidth = 0.8, formula = y ~ x) + facet_wrap(~stage_lab, ncol = 3) +
    labs(title = "invest: SHAP dependence by stage (shared y-axis; 1 dot = 1 station)", x = "invest (log1p, block mean)", y = "SHAP value") + theme_bw(base_size = 10)
  ggsave(file.path(OUT, "recovery_xgb_stages178_invest_dependence.png"), p4, width = 10, height = 6.5, dpi = 300)
}
cat("\n完成。\n")
