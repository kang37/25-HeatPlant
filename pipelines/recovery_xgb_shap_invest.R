#!/usr/bin/env Rscript
# =============================================================================
# recovery_xgb_shap_invest.R
#   此前恢复力driver分析为避免覆盖率问题/共线性,把invest(城市园林绿化投资)剔除了。
#   这里强制保留invest, 用XGBoost+TreeSHAP分别建φ(恢复慢度,越大越差)与
#   logOmega(Isbell抵抗力Ω=Yn/|Ye-Yn|取对数,越大=偏离越小=抵抗力越强)两个模型。xgboost原生处理缺失值,invest
#   缺失的站(117站里20个,57站里10个)不丢弃。
#
#   评估: 嵌套重复K折(内层xgb.cv选max_depth/min_child_weight/nrounds, 外层只评估),
#         对照 ridge / RF 基线、置换y的零分布、去掉invest的消融。
#   解释: TreeSHAP(predict(predcontrib=TRUE)); bootstrap稳定性(invest名次)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(xgboost); library(ggplot2)
                                library(glmnet); library(ranger); library(parallel)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")
NCORE <- max(1L, min(6L, detectCores() - 1L))
REPS_SCREEN <- 5L; REPS_FINAL <- 10L; K <- 5L; NPERM <- 20L; NBOOT <- 100L

# ---- 协变量 -----------------------------------------------------------------
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18  <- env$dt
VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))
VAR10 <- setdiff(VAR15, c("mixedleaf", "deci_needle", "tavg", "rh", "rsds_mean"))
VPD2  <- c("vpd_mean_stn", "vpd_cv_stn")
vc    <- fread(file.path(OUT, "vpd_station_covariate.csv"))[, .(meteo_stat, vpd_mean_stn, vpd_cv_stn)]

# ---- 站点清单 ----------------------------------------------------------------
sel117 <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")$meteo_stat
neg <- fread("data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv")
neg[, drho := rho - rho_min]
LCB <- "data_proc/output_hcsif_buf1000/landcover_break_stations.csv"
if (file.exists(LCB)) { lcb <- fread(LCB); neg <- neg[!meteo_stat %in% lcb[has_break == TRUE]$stat_id] }
sr <- neg[tp >= 0 & tp <= 8]; sr[, pass1 := (drho > 0) & (p_surr < 0.05)]
sel178 <- unique(sr[pass1 == TRUE]$meteo_stat)

# ---- DV1: phi(纯AR1) ---------------------------------------------------------
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)[meteo_stat %in% sel178, .(meteo_stat, year, doy, SIF_buf1000)]
d   <- copy(sif); setorder(d, meteo_stat, year, doy)
deseason_z <- function(x, doy) { clim <- tapply(x, doy, mean, na.rm = TRUE); a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x))); (a - mean(a, na.rm = TRUE)) / s }
d[, sif_anom := deseason_z(SIF_buf1000, doy), by = meteo_stat]
d[, lag1 := shift(sif_anom, 1), by = .(meteo_stat, year)]
fit_one <- function(sub) { sub <- sub[!is.na(lag1) & !is.na(sif_anom)]; if (nrow(sub) < 30) return(NULL)
  list(phi = unname(coef(lm(sif_anom ~ lag1, data = sub))["lag1"])) }
phi <- d[, fit_one(.SD), by = meteo_stat]
dphi <- merge(phi, dt18, by.x = "meteo_stat", by.y = "stat_id")
dphi <- merge(dphi, vc, by = "meteo_stat", all.x = TRUE)

# ---- DV2: logOmega(方向修正, 57站) -------------------------------------------
ev  <- fread(file.path(OUT, "recovery_liu2021_metrics_dirfix_117_events.csv"))
ev  <- ev[is.finite(Omega) & Omega > 0]
agg <- ev[, .(logOmega = as.numeric(mean(log(Omega), na.rm = TRUE)), direction_2v = direction_2v[1]), by = meteo_stat]
dom <- merge(agg, dt18, by.x = "meteo_stat", by.y = "stat_id")
dom <- merge(dom, vc, by = "meteo_stat")
dom[, is_promote := as.integer(direction_2v == "Promote")]

# ---- 配置: 全部强制含invest ---------------------------------------------------
CFG <- list(
  list(id = "phi_117_V15inv",    dv = "phi", data = dphi[meteo_stat %in% sel117], x = c(VAR15, "invest")),
  list(id = "phi_178_V15inv",    dv = "phi", data = dphi,                          x = c(VAR15, "invest")),
  list(id = "phi_117_V10vpdinv", dv = "phi", data = dphi[meteo_stat %in% sel117], x = c(VAR10, VPD2, "invest")),
  list(id = "phi_178_V10vpdinv", dv = "phi", data = dphi,                          x = c(VAR10, VPD2, "invest")),
  list(id = "omg_57_V10vpdinvP", dv = "logOmega", data = dom, x = c(VAR10, VPD2, "is_promote", "invest")),
  list(id = "omg_57_V15invP",    dv = "logOmega", data = dom, x = c(VAR15, "is_promote", "invest")),
  list(id = "omg_57_V10invP",    dv = "logOmega", data = dom, x = c(VAR10, "is_promote", "invest"))
)
for (i in seq_along(CFG)) {
  cc <- CFG[[i]]; dd <- cc$data
  non_inv <- setdiff(cc$x, "invest")
  dd <- dd[complete.cases(dd[, c(cc$dv, non_inv), with = FALSE])]
  CFG[[i]]$data <- dd
  cat(sprintf("%-20s n=%3d  invest缺失=%2d  p=%d\n", cc$id, nrow(dd), sum(is.na(dd$invest)), length(cc$x)))
}

# ---- xgboost: 内层调参 + 外层评估 ---------------------------------------------
GRID <- expand.grid(max_depth = c(2, 3), min_child_weight = c(1, 3))
base_p <- list(objective = "reg:squarederror", eval_metric = "rmse", eta = 0.03,
               subsample = 0.8, colsample_bytree = 0.8, lambda = 5, nthread = 1)
tune_xgb <- function(X, y, seed) {
  dtr <- xgb.DMatrix(X, label = y); best <- list(score = Inf)
  for (i in seq_len(nrow(GRID))) {
    p <- c(base_p, as.list(GRID[i, ]), list(seed = seed))
    cv <- xgb.cv(params = p, data = dtr, nrounds = 800, nfold = 5, early_stopping_rounds = 30, verbose = 0)
    el <- as.data.frame(cv$evaluation_log); j <- which.min(el$test_rmse_mean)
    if (el$test_rmse_mean[j] < best$score) best <- list(score = el$test_rmse_mean[j], p = p, n = j)
  }
  best$fit <- xgb.train(params = best$p, data = dtr, nrounds = best$n)
  best
}
r2f <- function(y, p) 1 - sum((y - p)^2) / sum((y - mean(y))^2)
impute <- function(Xtr, Xte) { med <- apply(Xtr, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(Xtr))) { Xtr[is.na(Xtr[, j]), j] <- med[j]; Xte[is.na(Xte[, j]), j] <- med[j] }
  list(tr = Xtr, te = Xte) }

one_rep <- function(X, y, r, baselines) {
  set.seed(1000 + r); n <- nrow(X); folds <- sample(rep(seq_len(K), length.out = n))
  o <- list(xgb = rep(NA_real_, n), ridge = rep(NA_real_, n), rf = rep(NA_real_, n))
  for (k in seq_len(K)) {
    tr <- folds != k; te <- !tr
    m <- tune_xgb(X[tr, , drop = FALSE], y[tr], seed = 100 * r + k)
    o$xgb[te] <- predict(m$fit, X[te, , drop = FALSE])
    if (baselines) {
      im <- impute(X[tr, , drop = FALSE], X[te, , drop = FALSE])
      mu <- colMeans(im$tr); sdv <- apply(im$tr, 2, sd); sdv[sdv == 0] <- 1
      Ztr <- scale(im$tr, mu, sdv); Zte <- scale(im$te, mu, sdv)
      cvr <- cv.glmnet(Ztr, y[tr], alpha = 0, nfolds = 5, standardize = FALSE)
      o$ridge[te] <- as.numeric(predict(cvr, Zte, s = "lambda.1se"))
      rf <- ranger(x = as.data.frame(im$tr), y = y[tr], num.trees = 500, seed = 1)
      o$rf[te] <- predict(rf, as.data.frame(im$te))$predictions
    }
  }
  data.table(rep = r, r2_xgb = r2f(y, o$xgb), rmse_xgb = sqrt(mean((y - o$xgb)^2)),
             rho_xgb = suppressWarnings(cor(y, o$xgb, method = "spearman")),
             r2_ridge = if (baselines) r2f(y, o$ridge) else NA_real_,
             r2_rf    = if (baselines) r2f(y, o$rf)    else NA_real_)
}
cv_eval <- function(X, y, reps, baselines = TRUE)
  rbindlist(mclapply(seq_len(reps), function(r) one_rep(X, y, r, baselines), mc.cores = NCORE))

# ---- 阶段1: 全部配置筛选 ------------------------------------------------------
screen <- rbindlist(lapply(CFG, function(cc) {
  X <- as.matrix(cc$data[, cc$x, with = FALSE]); y <- cc$data[[cc$dv]]
  t0 <- Sys.time(); r <- cv_eval(X, y, REPS_SCREEN)
  cat(sprintf("[screen] %-20s xgb R2=%.3f | ridge=%.3f | rf=%.3f  (%.0fs)\n", cc$id,
              mean(r$r2_xgb), mean(r$r2_ridge), mean(r$r2_rf), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  data.table(id = cc$id, dv = cc$dv, n = nrow(X), p = ncol(X), n_invest_na = sum(is.na(X[, "invest"])),
             xgb_r2 = mean(r$r2_xgb), xgb_r2_sd = sd(r$r2_xgb), xgb_rmse = mean(r$rmse_xgb),
             xgb_spearman = mean(r$rho_xgb), ridge_r2 = mean(r$r2_ridge), rf_r2 = mean(r$r2_rf))
}))
fwrite(screen, file.path(OUT, "recovery_xgb_shap_invest_screen.csv"))
print(screen)

# ---- 阶段2: 每个DV取最优配置 --------------------------------------------------
win <- screen[, .SD[which.max(xgb_r2)], by = dv]
cat("\n入选:\n"); print(win[, .(dv, id, n, xgb_r2)])

lab <- c(phi = "phi (AR1 persistence; higher = slower recovery)",
         logOmega = "logOmega (higher = smaller deviation = HIGHER resistance)")
final_tab <- list(); shap_long_all <- list(); imp_all <- list(); stab_all <- list()
for (w in seq_len(nrow(win))) {
  cc <- CFG[[which(sapply(CFG, `[[`, "id") == win$id[w])]]
  X <- as.matrix(cc$data[, cc$x, with = FALSE]); y <- cc$data[[cc$dv]]; dv <- cc$dv
  cat("\n==============", cc$id, "==============\n")

  r <- cv_eval(X, y, REPS_FINAL)
  Xn <- X[, setdiff(colnames(X), "invest"), drop = FALSE]
  r_noinv <- cv_eval(Xn, y, REPS_FINAL, baselines = FALSE)
  nullr2 <- unlist(mclapply(seq_len(NPERM), function(b) { set.seed(5000 + b)
    one_rep(X, sample(y), r = b, baselines = FALSE)$r2_xgb }, mc.cores = NCORE))
  p_null <- (1 + sum(nullr2 >= mean(r$r2_xgb))) / (NPERM + 1)
  cat(sprintf("嵌套CV(%d重x%d折): xgb R2=%.3f±%.3f RMSE=%.4f Spearman=%.3f | ridge=%.3f | RF=%.3f\n",
              REPS_FINAL, K, mean(r$r2_xgb), sd(r$r2_xgb), mean(r$rmse_xgb), mean(r$rho_xgb), mean(r$r2_ridge), mean(r$r2_rf)))
  cat(sprintf("置换零分布 R2: mean=%.3f max=%.3f -> p=%.3f | 去掉invest R2=%.3f (ΔR2=%+.3f)\n",
              mean(nullr2), max(nullr2), p_null, mean(r_noinv$r2_xgb), mean(r$r2_xgb) - mean(r_noinv$r2_xgb)))
  final_tab[[w]] <- data.table(dv = dv, id = cc$id, n = nrow(X), p = ncol(X), n_invest_na = sum(is.na(X[, "invest"])),
    xgb_r2 = mean(r$r2_xgb), xgb_r2_sd = sd(r$r2_xgb), xgb_rmse = mean(r$rmse_xgb), xgb_spearman = mean(r$rho_xgb),
    y_sd = sd(y), ridge_r2 = mean(r$r2_ridge), rf_r2 = mean(r$r2_rf),
    null_r2_mean = mean(nullr2), null_r2_max = max(nullr2), perm_p = p_null,
    noinvest_r2 = mean(r_noinv$r2_xgb), delta_r2_invest = mean(r$r2_xgb) - mean(r_noinv$r2_xgb))

  m <- tune_xgb(X, y, seed = 7)
  cat("全样本选定超参: depth=", m$p$max_depth, " min_child_weight=", m$p$min_child_weight, " nrounds=", m$n, "\n")
  sh <- predict(m$fit, X, predcontrib = TRUE); sh <- sh[, colnames(X), drop = FALSE]
  imp <- data.table(dv = dv, var = colnames(sh), mean_abs_shap = colMeans(abs(sh)))
  imp[, share := mean_abs_shap / sum(mean_abs_shap)]
  imp[, rank := rank(-mean_abs_shap)]
  imp[, shap_vs_value_spearman := sapply(var, function(v) { ok <- !is.na(X[, v])
    suppressWarnings(cor(X[ok, v], sh[ok, v], method = "spearman")) })]
  imp <- imp[order(rank)]
  cat("\nSHAP全局重要性:\n"); print(imp)
  inv_na <- is.na(X[, "invest"])
  if (any(inv_na)) cat(sprintf("invest缺失站的invest-SHAP: mean=%.4f (非缺失站 mean|SHAP|=%.4f)\n",
                                mean(sh[inv_na, "invest"]), mean(abs(sh[!inv_na, "invest"]))))

  set.seed(99)
  bt <- mclapply(seq_len(NBOOT), function(b) { idx <- sample(nrow(X), replace = TRUE)
    mb <- xgb.train(params = c(m$p, list(seed = b)), data = xgb.DMatrix(X[idx, , drop = FALSE], label = y[idx]), nrounds = m$n)
    s <- predict(mb, X, predcontrib = TRUE); colMeans(abs(s[, colnames(X), drop = FALSE])) }, mc.cores = NCORE)
  bt <- do.call(rbind, bt); rk <- t(apply(-bt, 1, rank))
  stab <- data.table(dv = dv, var = colnames(X), boot_mean = colMeans(bt),
                     rank_median = apply(rk, 2, median), rank_q10 = apply(rk, 2, quantile, .1), rank_q90 = apply(rk, 2, quantile, .9),
                     top5_freq = colMeans(rk <= 5))[order(rank_median)]
  cat("\nbootstrap稳定性(", NBOOT, "次):\n"); print(stab)
  stab_all[[w]] <- stab; imp_all[[w]] <- imp

  sl <- data.table(dv = dv, station = rep(cc$data$meteo_stat, ncol(sh)), var = rep(colnames(sh), each = nrow(sh)),
                   shap = as.vector(sh), value = as.vector(X))
  shap_long_all[[w]] <- sl
  fwrite(cbind(meteo_stat = cc$data$meteo_stat, y = y, pred = predict(m$fit, X), as.data.table(sh)),
         file.path(OUT, sprintf("recovery_xgb_shap_invest_%s_shap_values.csv", dv)))

  # ---- 图 ----
  ord <- imp$var; sl[, var := factor(var, levels = rev(ord))]
  sl[, vs := { v <- value; rg <- range(v, na.rm = TRUE); if (diff(rg) == 0) rep(0.5, .N) else (v - rg[1]) / diff(rg) }, by = var]
  sl[, is_inv := var == "invest"]
  p1 <- ggplot(sl, aes(shap, var, color = vs)) +
    geom_vline(xintercept = 0, color = "grey70") +
    geom_jitter(height = 0.2, size = 1.4, alpha = 0.8) +
    scale_color_gradient(low = "#2166AC", high = "#C1121F", na.value = "grey60", name = "Feature value\n(low -> high)") +
    labs(title = sprintf("SHAP summary: %s", dv), subtitle = sprintf("%s | n=%d | nested-CV R2=%.3f | grey = missing", lab[[dv]], nrow(X), mean(r$r2_xgb)),
         x = "SHAP value (impact on prediction)", y = NULL) + theme_bw(base_size = 11)
  ggsave(file.path(OUT, sprintf("recovery_xgb_shap_invest_%s_beeswarm.png", dv)), p1, width = 8.5, height = 6.5, dpi = 300)

  imp[, var := factor(var, levels = rev(var))]
  p2 <- ggplot(imp, aes(var, mean_abs_shap, fill = var == "invest")) + geom_col() + coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#C1121F", `FALSE` = "#2166AC"), guide = "none") +
    labs(title = sprintf("Mean |SHAP|: %s", dv), subtitle = "red = invest", x = NULL, y = "mean |SHAP|") + theme_bw(base_size = 11)
  ggsave(file.path(OUT, sprintf("recovery_xgb_shap_invest_%s_importance.png", dv)), p2, width = 7, height = 5.5, dpi = 300)

  top6 <- unique(c("invest", as.character(imp$var[imp$rank <= 6])))[1:6]
  dep <- sl[var %in% top6 & !is.na(value)]; dep[, var := factor(var, levels = top6)]
  p3 <- ggplot(dep, aes(value, shap)) + geom_hline(yintercept = 0, color = "grey70") + geom_point(alpha = 0.6, size = 1.3, color = "#2166AC") +
    geom_smooth(method = "loess", se = FALSE, color = "#C1121F", linewidth = 0.8, formula = y ~ x) +
    facet_wrap(~var, scales = "free", ncol = 3) +
    labs(title = sprintf("SHAP dependence: %s (invest + top features)", dv), x = "feature value", y = "SHAP value") + theme_bw(base_size = 18) + theme(strip.text = element_text(size = 19), plot.title = element_text(size = 20))
  ggsave(file.path(OUT, sprintf("recovery_xgb_shap_invest_%s_dependence.png", dv)), p3, width = 13, height = 8, dpi = 300)

}
fwrite(rbindlist(final_tab), file.path(OUT, "recovery_xgb_shap_invest_final_models.csv"))
fwrite(rbindlist(imp_all),   file.path(OUT, "recovery_xgb_shap_invest_importance.csv"))
fwrite(rbindlist(stab_all),  file.path(OUT, "recovery_xgb_shap_invest_boot_stability.csv"))
cat("\n完成。\n")
