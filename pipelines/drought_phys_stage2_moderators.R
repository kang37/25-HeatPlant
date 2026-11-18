#!/usr/bin/env Rscript
# =============================================================================
# drought_phys_stage2_moderators.R
#   两阶段法·第二阶段: 以第一阶段的站点敏感性 S = -beta (S越大=VPD抑制SIF越强; 控制了气温/
#   辐射/降水距平)为因变量, 回归调节变量, 其中invest为必含项。
#   权重 w = 1/(se_beta^2 + tau^2)(再归一化), 抵消第一阶段估计误差(噪声占比~0.36)。
#   方法: (1) WLS + VIF + 城市聚类稳健标准误  (2) 岭回归 + 城市整簇bootstrap
#         (3) XGBoost+TreeSHAP, 嵌套空间分块CV, 置换零分布, 去invest消融, 城市整簇bootstrap名次稳定性
#   注意: invest是城市层面变量, 同城站共用同一个值, 聚类/bootstrap都以城市为单位。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(car); library(glmnet)
                                library(ranger); library(xgboost); library(parallel)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")
NCORE <- max(1L, min(6L, detectCores() - 1L)); REPS <- 10L; K <- 5L; NPERM <- 20L; NBOOT_R <- 500L; NBOOT_X <- 100L

env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]; cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18 <- env$dt; smap <- env$smap
VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))
VAR10 <- setdiff(VAR15, c("mixedleaf", "deci_needle", "tavg", "rh", "rsds_mean"))
PRED  <- c(VAR10, "vpd_mean_stn", "vpd_cv_stn", "invest")

b <- fread(file.path(OUT, "drought_phys_arx_station_beta.csv"))
vp <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))
vc <- vp[, .(vpd_mean_stn = mean(vpd_mean, na.rm = TRUE), vpd_cv_stn = sd(vpd_mean, na.rm = TRUE) / mean(vpd_mean, na.rm = TRUE)), by = meteo_stat]
d <- merge(b, dt18[, c("stat_id", VAR10, "invest"), with = FALSE], by.x = "meteo_stat", by.y = "stat_id")
d <- merge(d, vc, by = "meteo_stat"); d <- merge(d, smap[, .(stat_id, city)], by.x = "meteo_stat", by.y = "stat_id", all.x = TRUE)
d[is.na(city), city := paste0("stn", meteo_stat)]
d <- d[complete.cases(d[, c(VAR10, "vpd_mean_stn", "vpd_cv_stn", "lon", "lat"), with = FALSE])]
d[, S := -beta]
tau2 <- max(var(d$S) - mean(d$beta_se^2), 0.25 * var(d$S)); d[, w := 1 / (beta_se^2 + tau2)]; d[, w := w / mean(w)]
cat(sprintf("站点数=%d | invest缺失=%d | 城市数=%d | S均值=%.3f SD=%.3f | 噪声占比=%.2f\n", nrow(d), sum(is.na(d$invest)),
            uniqueN(d[!is.na(invest)]$city), mean(d$S), sd(d$S), mean(d$beta_se^2) / var(d$S)))

r2w <- function(y, p, w) { yb <- sum(w * y) / sum(w); 1 - sum(w * (y - p)^2) / sum(w * (y - yb)^2) }
r2u <- function(y, p) 1 - sum((y - p)^2) / sum((y - mean(y))^2)

# ---- (0) invest的朴素关系 ------------------------------------------------------
dc <- d[!is.na(invest)]
dc[, tert := cut(invest, quantile(invest, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("low", "mid", "high"))]
cat("\n== invest与S的朴素关系 ==\nSpearman(S, invest) =", round(cor(dc$S, dc$invest, method = "spearman"), 3), "\n")
print(dc[, .(n = .N, S_mean = round(weighted.mean(S, w), 3), se = round(sqrt(sum(w^2 * (S - weighted.mean(S, w))^2)) / sum(w), 3)), by = tert][order(tert)])
cat("invest与其他调节变量Spearman:\n"); print(round(sapply(c(VAR10, "vpd_mean_stn", "vpd_cv_stn"), function(v) cor(dc$invest, dc[[v]], method = "spearman")), 2))

# ---- (1) WLS + VIF + 城市聚类稳健SE -------------------------------------------
Xz <- scale(as.matrix(dc[, PRED, with = FALSE])); y <- dc$S; w <- dc$w
fitw <- lm(S ~ ., data = data.frame(S = y, Xz), weights = w)
vif_v <- sort(vif(fitw), decreasing = TRUE); cat("\nVIF(降序):\n"); print(round(vif_v, 2))
X1 <- cbind(1, Xz); colnames(X1)[1] <- "(Intercept)"
cl_rob <- function(X1, y, w, g) {
  fit <- lm.wfit(X1, y, w); e <- fit$residuals; bread <- solve(crossprod(X1 * sqrt(w)))
  S_g <- rowsum(X1 * (w * e), g); meat <- crossprod(S_g); G <- nrow(S_g); n <- nrow(X1); kk <- ncol(X1)
  V <- bread %*% meat %*% bread * (G / (G - 1)) * ((n - 1) / (n - kk))
  se <- sqrt(diag(V)); tt <- fit$coefficients / se; list(coef = fit$coefficients, se = se, p = 2 * pt(-abs(tt), df = G - 1), G = G) }
rb <- cl_rob(X1, y, w, dc$city)
ols <- data.table(var = names(rb$coef), coef = rb$coef, se_cluster = rb$se, p = rb$p)[var != "(Intercept)"]
ols[, `:=`(ci_lo = coef - qt(.975, rb$G - 1) * se_cluster, ci_hi = coef + qt(.975, rb$G - 1) * se_cluster,
           sig = ifelse(p < .01, "**", ifelse(p < .05, "*", ifelse(p < .1, ".", ""))))]
cat(sprintf("\n== WLS(标准化X; 因变量S单位=SD; 聚类=%d个城市) R2w=%.3f ==\n", rb$G, r2w(y, as.numeric(X1 %*% rb$coef), w)))
print(ols[order(p)][, .(var, coef = round(coef, 4), se = round(se_cluster, 4), p = round(p, 4), sig)])
fwrite(ols, file.path(OUT, "drought_phys_stage2_wls_cluster.csv"))

# ---- (2) 岭回归 + 城市整簇bootstrap --------------------------------------------
set.seed(42); cvr <- cv.glmnet(Xz, y, alpha = 0, nfolds = 10, weights = w, standardize = FALSE)
rc <- as.numeric(coef(cvr, s = "lambda.1se"))[-1]; names(rc) <- PRED
ridge_r2 <- r2w(y, as.numeric(predict(cvr, Xz, s = "lambda.1se")), w)
cities <- unique(dc$city); idx_by_city <- split(seq_len(nrow(dc)), dc$city)
set.seed(1)
bm <- do.call(cbind, mclapply(seq_len(NBOOT_R), function(bi) { cs <- sample(cities, replace = TRUE); ix <- unlist(idx_by_city[cs])
  tryCatch({ m <- cv.glmnet(Xz[ix, ], y[ix], alpha = 0, nfolds = 5, weights = w[ix], standardize = FALSE); as.numeric(coef(m, s = "lambda.1se"))[-1] },
           error = function(e) rep(NA_real_, length(PRED))) }, mc.cores = NCORE))
rownames(bm) <- PRED
ridge <- data.table(var = PRED, ridge_coef = rc, ci_lo = apply(bm, 1, quantile, .025, na.rm = TRUE), ci_hi = apply(bm, 1, quantile, .975, na.rm = TRUE))
ridge[, excl0 := (ci_lo > 0) | (ci_hi < 0)]
cat(sprintf("\n== 岭回归(1se) R2w=%.3f; 城市整簇bootstrap 95%%CI ==\n", ridge_r2)); print(ridge[order(-abs(ridge_coef))][, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
fwrite(ridge, file.path(OUT, "drought_phys_stage2_ridge_boot.csv"))
ridge[, var := factor(var, levels = rev(ridge[order(-abs(ridge_coef))]$var))]
ggsave(file.path(OUT, "drought_phys_stage2_ridge_ci.png"),
  ggplot(ridge, aes(var, ridge_coef, color = excl0)) + geom_hline(yintercept = 0, color = "grey70") + geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi)) +
    coord_flip() + scale_color_manual(values = c(`TRUE` = "#C1121F", `FALSE` = "grey50"), guide = "none") +
    labs(title = "Moderators of VPD sensitivity S (=-beta): ridge + city-cluster bootstrap 95% CI",
         subtitle = sprintf("n=%d stations, %d cities | ridge R2w=%.3f | red = CI excludes 0 | S higher = stronger VPD suppression of SIF", nrow(dc), length(cities), ridge_r2),
         x = NULL, y = "Standardized ridge coefficient") + theme_bw(base_size = 11), width = 9, height = 6, dpi = 300)

# ---- (3) XGBoost + SHAP(空间分块CV; invest缺失保留) ----------------------------
GRID <- expand.grid(max_depth = c(2, 3, 4), min_child_weight = c(5, 10))
base_p <- list(objective = "reg:squarederror", eval_metric = "rmse", eta = 0.03, subsample = 0.8, colsample_bytree = 0.8, lambda = 5, nthread = 1)
spatial_folds <- function(lon, lat, nf, ncl, seed) { set.seed(seed); km <- kmeans(scale(cbind(lon, lat)), centers = ncl, nstart = 3)
  cf <- sample(rep(seq_len(nf), length.out = ncl)); cf[km$cluster] }
tune_xgb <- function(X, yv, wv, lon, lat, seed) {
  fid <- spatial_folds(lon, lat, 5, 12, seed); folds <- lapply(1:5, function(f) which(fid == f))
  dtr <- xgb.DMatrix(X, label = yv, weight = wv); best <- list(score = Inf)
  for (i in seq_len(nrow(GRID))) { p <- c(base_p, as.list(GRID[i, ]), list(seed = seed))
    cv <- xgb.cv(params = p, data = dtr, nrounds = 800, folds = folds, early_stopping_rounds = 30, verbose = 0)
    el <- as.data.frame(cv$evaluation_log); j <- which.min(el$test_rmse_mean)
    if (el$test_rmse_mean[j] < best$score) best <- list(score = el$test_rmse_mean[j], p = p, n = j) }
  best$fit <- xgb.train(params = best$p, data = dtr, nrounds = best$n); best }
impute <- function(Xtr, Xte) { med <- apply(Xtr, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(Xtr))) { Xtr[is.na(Xtr[, j]), j] <- med[j]; Xte[is.na(Xte[, j]), j] <- med[j] }; list(tr = Xtr, te = Xte) }
one_rep <- function(dd, feats, r, baselines, yv = dd$S, wv = dd$w) {
  X <- as.matrix(dd[, feats, with = FALSE]); fid <- spatial_folds(dd$lon, dd$lat, K, 30, 2000 + r)
  px <- pr <- pf <- rep(NA_real_, nrow(X))
  for (f in seq_len(K)) { tr <- fid != f; te <- !tr
    m <- tune_xgb(X[tr, , drop = FALSE], yv[tr], wv[tr], dd$lon[tr], dd$lat[tr], seed = 100 * r + f)
    px[te] <- predict(m$fit, X[te, , drop = FALSE])
    if (baselines) { im <- impute(X[tr, , drop = FALSE], X[te, , drop = FALSE]); mu <- colMeans(im$tr); sdv <- apply(im$tr, 2, sd); sdv[sdv == 0] <- 1
      cvm <- cv.glmnet(scale(im$tr, mu, sdv), yv[tr], alpha = 0, nfolds = 5, weights = wv[tr], standardize = FALSE)
      pr[te] <- as.numeric(predict(cvm, scale(im$te, mu, sdv), s = "lambda.1se"))
      rf <- ranger(x = as.data.frame(im$tr), y = yv[tr], num.trees = 500, case.weights = wv[tr], seed = 1)
      pf[te] <- predict(rf, as.data.frame(im$te))$predictions } }
  data.table(rep = r, r2_w = r2w(yv, px, wv), r2_u = r2u(yv, px), rmse = sqrt(mean((yv - px)^2)),
             ridge_r2w = if (baselines) r2w(yv, pr, wv) else NA_real_, rf_r2w = if (baselines) r2w(yv, pf, wv) else NA_real_) }
cv_eval <- function(dd, feats, baselines = TRUE) rbindlist(mclapply(seq_len(REPS), function(r) one_rep(dd, feats, r, baselines), mc.cores = NCORE))

rA <- cv_eval(d, PRED); rN <- cv_eval(d, setdiff(PRED, "invest"), FALSE)
nullr <- unlist(mclapply(seq_len(NPERM), function(bi) { set.seed(5000 + bi); ix <- sample(nrow(d)); one_rep(d, PRED, bi, FALSE, yv = d$S[ix], wv = d$w[ix])$r2_w }, mc.cores = NCORE))
p_null <- (1 + sum(nullr >= mean(rA$r2_w))) / (NPERM + 1)
cat(sprintf("\n== XGBoost 空间分块CV(%d重x%d折, n=%d, invest缺失保留) ==\nR2w=%.3f±%.3f R2u=%.3f RMSE=%.3f | ridge R2w=%.3f RF R2w=%.3f\n置换零分布 max=%.3f p=%.3f | 去invest R2w=%.3f (invest增益 %+.3f)\n",
            REPS, K, nrow(d), mean(rA$r2_w), sd(rA$r2_w), mean(rA$r2_u), mean(rA$rmse), mean(rA$ridge_r2w), mean(rA$rf_r2w),
            max(nullr), p_null, mean(rN$r2_w), mean(rA$r2_w) - mean(rN$r2_w)))
fwrite(data.table(n = nrow(d), r2w = mean(rA$r2_w), r2w_sd = sd(rA$r2_w), ridge_r2w = mean(rA$ridge_r2w), rf_r2w = mean(rA$rf_r2w),
                  null_max = max(nullr), perm_p = p_null, noinvest_r2w = mean(rN$r2_w)), file.path(OUT, "drought_phys_stage2_xgb_cv.csv"))

X <- as.matrix(d[, PRED, with = FALSE]); m <- tune_xgb(X, d$S, d$w, d$lon, d$lat, seed = 7)
sh <- predict(m$fit, X, predcontrib = TRUE)[, PRED, drop = FALSE]
imp <- data.table(var = PRED, mean_abs_shap = colMeans(abs(sh)))
imp[, `:=`(share = mean_abs_shap / sum(mean_abs_shap), rank = rank(-mean_abs_shap))]
imp[, shap_vs_value_spearman := sapply(var, function(v) { ok <- !is.na(X[, v]); suppressWarnings(cor(X[ok, v], sh[ok, v], method = "spearman")) })]
city_idx <- split(seq_len(nrow(d)), d$city); cs_all <- names(city_idx); set.seed(99)
bt <- do.call(rbind, mclapply(seq_len(NBOOT_X), function(bi) { ix <- unlist(city_idx[sample(cs_all, replace = TRUE)])
  mb <- xgb.train(params = c(m$p, list(seed = bi)), data = xgb.DMatrix(X[ix, , drop = FALSE], label = d$S[ix], weight = d$w[ix]), nrounds = m$n)
  colMeans(abs(predict(mb, X, predcontrib = TRUE)[, PRED, drop = FALSE])) }, mc.cores = NCORE))
rk <- t(apply(-bt, 1, rank))
imp[, `:=`(rank_q10 = apply(rk, 2, quantile, .1)[var], rank_q90 = apply(rk, 2, quantile, .9)[var], top5_freq = colMeans(rk <= 5)[var])]
imp <- imp[order(rank)]; cat("\nSHAP全局重要性(城市整簇bootstrap名次):\n"); print(imp[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(imp, file.path(OUT, "drought_phys_stage2_xgb_importance.csv"))

sl <- data.table(var = rep(PRED, each = nrow(sh)), shap = as.vector(sh), value = as.vector(X)); sl[, var := factor(var, levels = rev(imp$var))]
sl[, vs := { rg <- range(value, na.rm = TRUE); if (diff(rg) == 0) rep(0.5, .N) else (value - rg[1]) / diff(rg) }, by = var]
ggsave(file.path(OUT, "drought_phys_stage2_shap_beeswarm.png"),
  ggplot(sl, aes(shap, var, color = vs)) + geom_vline(xintercept = 0, color = "grey70") + geom_jitter(height = 0.2, size = 0.8, alpha = 0.7) +
    scale_color_gradient(low = "#2166AC", high = "#C1121F", na.value = "grey60", name = "Feature value\n(low -> high)") +
    labs(title = "SHAP: moderators of VPD sensitivity S (=-beta; higher = stronger suppression)",
         subtitle = sprintf("n=%d | spatial-block CV R2w=%.3f | grey = missing (invest)", nrow(d), mean(rA$r2_w)), x = "SHAP value", y = NULL) + theme_bw(base_size = 10),
  width = 9, height = 6.5, dpi = 300)
top6 <- unique(c("invest", as.character(imp$var[1:5])))[1:6]; dep <- sl[var %in% top6 & !is.na(value)]; dep[, var := factor(var, levels = top6)]
ggsave(file.path(OUT, "drought_phys_stage2_shap_dependence.png"),
  ggplot(dep, aes(value, shap)) + geom_hline(yintercept = 0, color = "grey70") + geom_point(alpha = 0.4, size = 0.8, color = "#2166AC") +
    geom_smooth(method = "loess", se = FALSE, color = "#C1121F", linewidth = 0.8, formula = y ~ x) + facet_wrap(~var, scales = "free", ncol = 3) +
    labs(title = "SHAP dependence (invest + top moderators; free axes)", x = "feature value", y = "SHAP value") + theme_bw(base_size = 10),
  width = 10, height = 6, dpi = 300)
cat("\n完成。\n")
