#!/usr/bin/env Rscript
# =============================================================================
# recovery_xgb_shap_phi_panel.R
#   φ的站-块面板版: 把相邻k年(k=3主分析, k=6敏感性)合并估一个φ(见
#   recovery_phi_block_pooling.R), 配块内平均的动态协变量(tavg/rh/precip/rsds/ntl/invest),
#   用XGBoost+TreeSHAP建模。动态协变量面板文件只覆盖117站, 所以这里限定117站。
#   评估: 嵌套重复CV, 外层与内层都按"站点分组"切(同一站的不同块绝不跨训练/测试);
#         另做空间分块CV(kmeans地理聚类, 整簇留出)。按 1/(SE^2+tau^2) 加权。
#   对照: ridge/RF基线; 置换y零分布; 去invest消融; 加block中点年份(检验ntl/invest是否只是时间趋势代理)。
#   解释: SHAP(含invest) + 按站点整体重抽样的bootstrap名次稳定性。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(xgboost); library(ggplot2)
                                library(glmnet); library(ranger); library(parallel)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")
NCORE <- max(1L, min(6L, detectCores() - 1L))
REPS <- 10L; K <- 5L; NPERM <- 20L; NBOOT <- 100L

# ---- 数据 --------------------------------------------------------------------
blocks <- fread(file.path(OUT, "recovery_phi_block_pooling_blocks.csv"))
pc <- fread(file.path(OUT, "recovery_liu2021_panel_covariates.csv"))
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]; cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
rsds_sd <- env$dt[, .(stat_id, rsds_sd)]
co <- fread("data_raw/WS土地利用比例_1000m.csv", select = 1:3); setnames(co, c("stat_id", "lon", "lat"))
DYN <- c("tavg", "rh", "precip", "rsds", "ntl", "invest")
STA <- c("imperv", "grass", "water", "ever_needle", "ever_broad", "deci_broad", "elev", "rsds_sd")
YEAR0 <- 2000L; NYEAR <- 23L
nanmean <- function(x) { x <- x[!is.na(x)]; if (length(x)) mean(x) else NA_real_ }

build_panel <- function(kk) {
  b <- blocks[blocks$k == kk & meteo_stat %in% unique(pc$stat_id)]
  nb <- floor(NYEAR / kk); stopifnot(max(b$blk) == nb)
  p <- copy(pc); p[, blk := (year - YEAR0) %/% kk + 1L]; p <- p[blk <= nb]
  dyn <- p[, lapply(.SD, nanmean), by = .(stat_id, blk), .SDcols = DYN]
  sta <- p[, lapply(.SD, function(x) x[1]), by = stat_id, .SDcols = setdiff(STA, "rsds_sd")]
  yr  <- p[, .(year_mid = mean(year)), by = .(stat_id, blk)]
  d <- merge(b, dyn, by.x = c("meteo_stat", "blk"), by.y = c("stat_id", "blk"))
  d <- merge(d, sta, by.x = "meteo_stat", by.y = "stat_id")
  d <- merge(d, yr, by.x = c("meteo_stat", "blk"), by.y = c("stat_id", "blk"))
  d <- merge(d, rsds_sd, by.x = "meteo_stat", by.y = "stat_id")
  d <- merge(d, co, by.x = "meteo_stat", by.y = "stat_id")
  d <- d[is.finite(phi) & is.finite(se)]
  d <- d[complete.cases(d[, c(setdiff(DYN, "invest"), STA), with = FALSE])]
  tau2 <- max(var(d$phi) - mean(d$se^2), 0.25 * var(d$phi))
  d[, w := 1 / (se^2 + tau2)]; d[, w := w / mean(w)]
  attr(d, "noise_share") <- mean(d$se^2) / var(d$phi)
  d
}

# ---- 评估工具 -----------------------------------------------------------------
GRID <- expand.grid(max_depth = c(2, 3, 4), min_child_weight = c(5, 10))
base_p <- list(objective = "reg:squarederror", eval_metric = "rmse", eta = 0.03,
               subsample = 0.8, colsample_bytree = 0.8, lambda = 5, nthread = 1)
group_folds <- function(g, nf, seed) {
  set.seed(seed); ug <- unique(g); f <- sample(rep(seq_len(nf), length.out = length(ug))); names(f) <- ug
  f[as.character(g)] }
tune_xgb <- function(X, y, w, g, seed) {
  fid <- group_folds(g, 5, seed); folds <- lapply(1:5, function(f) which(fid == f))
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
station_folds <- function(dd, scheme, r) {
  st <- unique(dd[, .(meteo_stat, lon, lat)])
  if (scheme == "group") { f <- group_folds(st$meteo_stat, K, 1000 + r)
  } else { set.seed(2000 + r); kc <- 20L
    km <- kmeans(scale(as.matrix(st[, .(lon, lat)])), centers = kc, nstart = 5)
    cf <- sample(rep(seq_len(K), length.out = kc)); f <- cf[km$cluster]; names(f) <- st$meteo_stat }
  f[as.character(dd$meteo_stat)] }

one_rep <- function(dd, feats, scheme, r, baselines, y = dd$phi, w = dd$w) {
  X <- as.matrix(dd[, feats, with = FALSE]); g <- dd$meteo_stat; n <- nrow(X)
  fid <- station_folds(dd, scheme, r)
  px <- pr <- pf <- rep(NA_real_, n)
  for (f in seq_len(K)) {
    tr <- fid != f; te <- !tr
    m <- tune_xgb(X[tr, , drop = FALSE], y[tr], w[tr], g[tr], seed = 100 * r + f)
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
  agg <- data.table(g = g, y = y, p = px)[, .(y = mean(y), p = mean(p)), by = g]
  data.table(rep = r, r2_u = r2u(y, px), r2_w = r2w(y, px, w), rmse = sqrt(mean((y - px)^2)),
             r2_station = r2u(agg$y, agg$p),
             ridge_r2w = if (baselines) r2w(y, pr, w) else NA_real_,
             rf_r2w    = if (baselines) r2w(y, pf, w) else NA_real_)
}
cv_eval <- function(dd, feats, scheme, baselines = TRUE, reps = REPS)
  rbindlist(mclapply(seq_len(reps), function(r) one_rep(dd, feats, scheme, r, baselines), mc.cores = NCORE))
sm <- function(r) sprintf("R2w=%.3f±%.3f R2u=%.3f R2_stn=%.3f RMSE=%.3f | ridge R2w=%.3f RF R2w=%.3f",
  mean(r$r2_w), sd(r$r2_w), mean(r$r2_u), mean(r$r2_station), mean(r$rmse), mean(r$ridge_r2w), mean(r$rf_r2w))

# ---- 主流程 -------------------------------------------------------------------
FE <- c(DYN, STA); FE_NOINV <- setdiff(FE, "invest"); FE_TIME <- c(FE, "year_mid")
lab <- "phi block (AR1 persistence; higher = slower recovery)"
all_tab <- list(); all_imp <- list(); all_stab <- list()
for (k in c(3L, 6L)) {
  dd <- build_panel(k)
  ns <- attr(dd, "noise_share")
  cat(sprintf("\n================ k=%d: rows=%d stations=%d invest缺失行=%d | SE噪声占比=%.2f | 站间ICC见summary ================\n",
              k, nrow(dd), uniqueN(dd$meteo_stat), sum(is.na(dd$invest)), ns))
  res <- list(base_group = cv_eval(dd, FE, "group"), base_spatial = cv_eval(dd, FE, "spatial"),
              noinv_group = cv_eval(dd, FE_NOINV, "group", FALSE), noinv_spatial = cv_eval(dd, FE_NOINV, "spatial", FALSE),
              time_group = cv_eval(dd, FE_TIME, "group", FALSE))
  for (nm in names(res)) cat(sprintf("[%-13s] %s\n", nm, sm(res[[nm]])))
  nullr <- unlist(mclapply(seq_len(NPERM), function(b) { set.seed(5000 + b); ix <- sample(nrow(dd))
    one_rep(dd, FE, "group", r = b, baselines = FALSE, y = dd$phi[ix], w = dd$w[ix])$r2_w }, mc.cores = NCORE))
  p_null <- (1 + sum(nullr >= mean(res$base_group$r2_w))) / (NPERM + 1)
  cat(sprintf("置换零分布 R2w: mean=%.3f max=%.3f -> p=%.3f\n", mean(nullr), max(nullr), p_null))
  all_tab[[as.character(k)]] <- data.table(k = k, rows = nrow(dd), stations = uniqueN(dd$meteo_stat), noise_share = ns,
    scheme = c("group", "spatial", "group_noinvest", "spatial_noinvest", "group_withtime"),
    r2_w = sapply(res, function(r) mean(r$r2_w)), r2_w_sd = sapply(res, function(r) sd(r$r2_w)),
    r2_u = sapply(res, function(r) mean(r$r2_u)), r2_station = sapply(res, function(r) mean(r$r2_station)),
    ridge_r2w = sapply(res, function(r) mean(r$ridge_r2w)), rf_r2w = sapply(res, function(r) mean(r$rf_r2w)),
    null_r2w_mean = mean(nullr), null_r2w_max = max(nullr), perm_p = p_null)

  # ---- 全样本模型 + SHAP ----
  X <- as.matrix(dd[, FE, with = FALSE]); y <- dd$phi
  m <- tune_xgb(X, y, dd$w, dd$meteo_stat, seed = 7)
  cat("全样本超参: depth=", m$p$max_depth, " min_child_weight=", m$p$min_child_weight, " nrounds=", m$n, "\n")
  sh <- predict(m$fit, X, predcontrib = TRUE)[, FE, drop = FALSE]
  imp <- data.table(k = k, var = FE, mean_abs_shap = colMeans(abs(sh)))
  imp[, share := mean_abs_shap / sum(mean_abs_shap)]; imp[, rank := rank(-mean_abs_shap)]
  imp[, shap_vs_value_spearman := sapply(var, function(v) { ok <- !is.na(X[, v]); suppressWarnings(cor(X[ok, v], sh[ok, v], method = "spearman")) })]
  imp[, type := ifelse(var %in% DYN, "dynamic", "static")]
  imp <- imp[order(rank)]; cat("\nSHAP全局重要性:\n"); print(imp)
  inv_na <- is.na(X[, "invest"])
  if (any(inv_na)) cat(sprintf("invest缺失行的invest-SHAP均值=%.4f (非缺失行 mean|SHAP|=%.4f)\n", mean(sh[inv_na, "invest"]), mean(abs(sh[!inv_na, "invest"]))))

  stns <- unique(dd$meteo_stat); set.seed(99)
  bt <- mclapply(seq_len(NBOOT), function(b) {
    sb <- sample(stns, replace = TRUE); idx <- unlist(lapply(sb, function(s) which(dd$meteo_stat == s)))
    mb <- xgb.train(params = c(m$p, list(seed = b)), data = xgb.DMatrix(X[idx, , drop = FALSE], label = y[idx], weight = dd$w[idx]), nrounds = m$n)
    colMeans(abs(predict(mb, X, predcontrib = TRUE)[, FE, drop = FALSE])) }, mc.cores = NCORE)
  bt <- do.call(rbind, bt); rk <- t(apply(-bt, 1, rank))
  stab <- data.table(k = k, var = FE, boot_mean = colMeans(bt), rank_median = apply(rk, 2, median),
                     rank_q10 = apply(rk, 2, quantile, .1), rank_q90 = apply(rk, 2, quantile, .9), top5_freq = colMeans(rk <= 5))[order(rank_median)]
  cat("\n站点整体bootstrap稳定性(", NBOOT, "次):\n"); print(stab)
  all_imp[[as.character(k)]] <- imp; all_stab[[as.character(k)]] <- stab
  fwrite(cbind(dd[, .(meteo_stat, blk, phi, se, w)], as.data.table(sh)), file.path(OUT, sprintf("recovery_xgb_panel_phi_k%d_shap_values.csv", k)))

  # ---- 图 ----
  sl <- data.table(var = rep(FE, each = nrow(sh)), shap = as.vector(sh), value = as.vector(X))
  sl[, var := factor(var, levels = rev(imp$var))]
  sl[, vs := { rg <- range(value, na.rm = TRUE); if (diff(rg) == 0) rep(0.5, .N) else (value - rg[1]) / diff(rg) }, by = var]
  p1 <- ggplot(sl, aes(shap, var, color = vs)) + geom_vline(xintercept = 0, color = "grey70") +
    geom_jitter(height = 0.2, size = 1.1, alpha = 0.75) +
    scale_color_gradient(low = "#2166AC", high = "#C1121F", na.value = "grey60", name = "Feature value\n(low -> high)") +
    labs(title = sprintf("SHAP summary: phi panel, k=%d-year blocks", k),
         subtitle = sprintf("%s | rows=%d, stations=%d | station-grouped CV R2w=%.3f | grey = missing", lab, nrow(dd), uniqueN(dd$meteo_stat), mean(res$base_group$r2_w)),
         x = "SHAP value", y = NULL) + theme_bw(base_size = 10)
  ggsave(file.path(OUT, sprintf("recovery_xgb_panel_phi_k%d_beeswarm.png", k)), p1, width = 9, height = 6.5, dpi = 300)
  imp2 <- copy(imp); imp2[, var := factor(var, levels = rev(var))]
  p2 <- ggplot(imp2, aes(var, mean_abs_shap, fill = type)) + geom_col() + coord_flip() +
    scale_fill_manual(values = c(dynamic = "#C1121F", static = "#2166AC")) +
    labs(title = sprintf("Mean |SHAP|: phi panel, k=%d", k), subtitle = "red = dynamic (block-mean) covariates, blue = static", x = NULL, y = "mean |SHAP|", fill = NULL) + theme_bw(base_size = 11)
  ggsave(file.path(OUT, sprintf("recovery_xgb_panel_phi_k%d_importance.png", k)), p2, width = 7.5, height = 5.5, dpi = 300)
  top6 <- unique(c("invest", imp$var[1:5]))[1:6]
  dep <- sl[var %in% top6 & !is.na(value)]; dep[, var := factor(var, levels = top6)]
  p3 <- ggplot(dep, aes(value, shap)) + geom_hline(yintercept = 0, color = "grey70") + geom_point(alpha = 0.5, size = 1, color = "#2166AC") +
    geom_smooth(method = "loess", se = FALSE, color = "#C1121F", linewidth = 0.8, formula = y ~ x) + facet_wrap(~var, scales = "free", ncol = 3) +
    labs(title = sprintf("SHAP dependence (note: free y-axes): phi panel k=%d", k), x = "feature value", y = "SHAP value") + theme_bw(base_size = 10)
  ggsave(file.path(OUT, sprintf("recovery_xgb_panel_phi_k%d_dependence.png", k)), p3, width = 10, height = 6, dpi = 300)
}
fwrite(rbindlist(all_tab),  file.path(OUT, "recovery_xgb_panel_phi_cv_summary.csv"))
fwrite(rbindlist(all_imp),  file.path(OUT, "recovery_xgb_panel_phi_importance.csv"))
fwrite(rbindlist(all_stab), file.path(OUT, "recovery_xgb_panel_phi_boot_stability.csv"))
cat("\n完成。\n")
