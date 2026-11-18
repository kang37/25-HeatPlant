#!/usr/bin/env Rscript
# =============================================================================
# recovery_phi_best_model.R
#   8个φ版本里按驱动模型R²选出最好的: adjA+ED(SIF) = 气候调整(气温/降水/辐射, Guo 2026式)
#   + 遮掉SIF扰动期后的φ, 178站(气候齐全配对>=100者, n=175), 15变量。
#   方法: OLS(+VIF) / RF(置换重要性, OOB) / 岭回归+bootstrap(500)  y与X均标准化;
#   性能: RF OOB R², 岭回归10折CV R²(20次重复), 空间分块CV R²(kmeans 20簇, 5折, 20次);
#   敏感性: 加入站点平均SIF(信号强度, 检验"干旱站SIF信号弱、噪声压低自相关"这一解释)。
#   输出: 总结图(系数森林图 + RF重要性 + 实测vs OOB预测)与CSV。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(car); library(glmnet); library(ranger); library(patchwork)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- Sys.getenv("PHI_OUT", file.path(PROJ, "data_proc/output_loose_classify")); NBOOT <- 500L
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]; cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18 <- env$dt; VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))
CLIM <- c("tavg", "rh", "precip", "rsds_mean", "rsds_sd", "elev")
ph <- fread(file.path(OUT, "recovery_phi_variants_178.csv"))[, .(meteo_stat, phi = phi_adjA_ED_sif)]
co <- fread("data_raw/WS土地利用比例_1000m.csv", select = 1:3); setnames(co, c("meteo_stat", "lon", "lat"))
fs <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sifm <- rbindlist(lapply(fs, fread, select = c("meteo_stat", "doy", "SIF_buf1000")), fill = TRUE)[doy >= 124 & doy <= 268, .(sif_mean = mean(SIF_buf1000, na.rm = TRUE)), by = meteo_stat]
dd <- Reduce(function(a, b) merge(a, b, by = "meteo_stat"), list(ph, co, sifm, dt18[, c("stat_id", VAR15), with = FALSE][, meteo_stat := stat_id][, stat_id := NULL]))
dd <- dd[is.finite(phi) & complete.cases(dd[, VAR15, with = FALSE])]
cat("建模站数:", nrow(dd), "\n"); n <- nrow(dd)
dd[, yy := as.numeric(scale(phi))]
f <- as.formula(paste("yy ~", paste(VAR15, collapse = "+")))

# ---- OLS + VIF ---------------------------------------------------------------
Xz <- dd[, lapply(.SD, function(z) as.numeric(scale(z))), .SDcols = VAR15]; Xz[, yy := dd$yy]
ols <- lm(f, data = Xz); vv <- vif(ols); sm <- summary(ols)$coefficients
cat(sprintf("OLS R2=%.3f | 最大VIF=%.1f(%s)\n", summary(ols)$r.squared, max(vv), names(which.max(vv))))
# ---- RF ----------------------------------------------------------------------
set.seed(42); rf <- ranger(f, data = dd, num.trees = 2000, importance = "permutation", seed = 42)
cat(sprintf("RF OOB R2=%.3f\n", rf$r.squared))
# ---- 岭回归 + bootstrap ---------------------------------------------------------
X <- scale(as.matrix(dd[, VAR15, with = FALSE])); y <- dd$yy
set.seed(42); cv <- cv.glmnet(X, y, alpha = 0, nfolds = 10, standardize = FALSE)
rc <- as.numeric(coef(cv, s = "lambda.1se"))[-1]; names(rc) <- VAR15
set.seed(1); bm <- replicate(NBOOT, { ix <- sample(n, replace = TRUE)
  tryCatch(as.numeric(coef(cv.glmnet(X[ix, , drop = FALSE], y[ix], alpha = 0, nfolds = 5, standardize = FALSE), s = "lambda.1se"))[-1], error = function(e) rep(NA_real_, length(VAR15))) })
lo <- apply(bm, 1, quantile, .025, na.rm = TRUE); hi <- apply(bm, 1, quantile, .975, na.rm = TRUE)
res <- data.table(var = VAR15, ridge = rc, lo = lo, hi = hi, robust = (lo > 0) | (hi < 0),
                  ols = sm[VAR15, 1], ols_p = sm[VAR15, 4], vif = vv[VAR15], rf_imp = rf$variable.importance[VAR15])
res[, type := fifelse(var %in% CLIM, "climate/terrain", fifelse(var %in% c("ntl", "imperv"), "urbanization", "land cover"))]
res <- res[order(-abs(ridge))]; print(res[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])
cat("气候类变量占|岭回归系数|比例:", round(sum(abs(res[type == "climate/terrain"]$ridge)) / sum(abs(res$ridge)), 3), "\n")

# ---- 性能: 随机10折 / 空间分块CV ------------------------------------------------
r2 <- function(obs, pr) 1 - sum((obs - pr)^2) / sum((obs - mean(obs))^2)
cv_run <- function(folds_fun, reps) {
  rr <- lapply(seq_len(reps), function(r) { fid <- folds_fun(r); pr <- pf <- rep(NA_real_, n)
    for (k in unique(fid)) { tr <- fid != k; te <- !tr
      mu <- colMeans(X[tr, , drop = FALSE]); sdv <- apply(X[tr, , drop = FALSE], 2, sd); sdv[sdv == 0] <- 1
      cvr <- cv.glmnet(scale(X[tr, , drop = FALSE], mu, sdv), y[tr], alpha = 0, nfolds = 5, standardize = FALSE)
      pr[te] <- as.numeric(predict(cvr, scale(X[te, , drop = FALSE], mu, sdv), s = "lambda.1se"))
      pf[te] <- predict(ranger(x = as.data.frame(X[tr, , drop = FALSE]), y = y[tr], num.trees = 500, seed = 1), as.data.frame(X[te, , drop = FALSE]))$predictions }
    c(ridge = r2(y, pr), rf = r2(y, pf)) })
  do.call(rbind, rr) }
rand_f <- function(r) { set.seed(100 + r); sample(rep(1:10, length.out = n)) }
spat_f <- function(r) { set.seed(200 + r); km <- kmeans(scale(as.matrix(dd[, .(lon, lat)])), centers = 20, nstart = 5); cf <- sample(rep(1:5, length.out = 20)); cf[km$cluster] }
pr_rand <- cv_run(rand_f, 20); pr_spat <- cv_run(spat_f, 20)
perf <- data.table(scheme = c("random 10-fold CV", "spatial-block CV (20 clusters, 5 folds)"),
  ridge_r2 = c(mean(pr_rand[, "ridge"]), mean(pr_spat[, "ridge"])), ridge_sd = c(sd(pr_rand[, "ridge"]), sd(pr_spat[, "ridge"])),
  rf_r2 = c(mean(pr_rand[, "rf"]), mean(pr_spat[, "rf"])), rf_sd = c(sd(pr_rand[, "rf"]), sd(pr_spat[, "rf"])))
perf <- rbind(data.table(scheme = "in-sample / OOB", ridge_r2 = 1 - sum((y - predict(cv, X, s = "lambda.1se"))^2) / sum((y - mean(y))^2), ridge_sd = NA, rf_r2 = rf$r.squared, rf_sd = NA), perf)
perf[, ols_r2 := c(summary(ols)$r.squared, NA, NA)]
print(perf[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])

# ---- 敏感性: 加入站点平均SIF -----------------------------------------------------
cat(sprintf("\nSpearman(phi, 站点平均SIF)=%.2f | (precip, 平均SIF)=%.2f | (rh, 平均SIF)=%.2f\n", cor(dd$phi, dd$sif_mean, method = "spearman"),
            cor(dd$precip, dd$sif_mean, method = "spearman"), cor(dd$rh, dd$sif_mean, method = "spearman")))
V2 <- c(VAR15, "sif_mean"); X2 <- scale(as.matrix(dd[, V2, with = FALSE]))
set.seed(42); cv2 <- cv.glmnet(X2, y, alpha = 0, nfolds = 10, standardize = FALSE); rc2 <- as.numeric(coef(cv2, s = "lambda.1se"))[-1]; names(rc2) <- V2
set.seed(1); bm2 <- replicate(NBOOT, { ix <- sample(n, replace = TRUE)
  tryCatch(as.numeric(coef(cv.glmnet(X2[ix, , drop = FALSE], y[ix], alpha = 0, nfolds = 5, standardize = FALSE), s = "lambda.1se"))[-1], error = function(e) rep(NA_real_, length(V2))) })
lo2 <- apply(bm2, 1, quantile, .025, na.rm = TRUE); hi2 <- apply(bm2, 1, quantile, .975, na.rm = TRUE)
sens <- data.table(var = V2, ridge_with_sifmean = rc2, lo = lo2, hi = hi2, robust = (lo2 > 0) | (hi2 < 0))
sens <- merge(res[, .(var, ridge_base = ridge, robust_base = robust)], sens, by = "var", all = TRUE)[order(-abs(ridge_with_sifmean))]
cat("加入站点平均SIF后的岭回归系数(对照原模型):\n"); print(sens[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])
fwrite(res, file.path(OUT, "recovery_phi_best_model_effects.csv")); fwrite(perf, file.path(OUT, "recovery_phi_best_model_performance.csv"))
fwrite(sens, file.path(OUT, "recovery_phi_best_model_sifmean_sensitivity.csv"))

# ---- 图 ----------------------------------------------------------------------
pal <- c("climate/terrain" = "#C1121F", "urbanization" = "#E08A00", "land cover" = "#2166AC")
res[, var := factor(var, levels = rev(res$var))]
pA <- ggplot(res, aes(var, ridge, color = type, alpha = robust)) + geom_hline(yintercept = 0, color = "grey70") +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.55) + coord_flip() + scale_color_manual(values = pal, name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  labs(title = "A  Effect on phi (ridge + bootstrap 95% CI)", subtitle = "standardized; solid = CI excludes 0; positive = slower recovery", x = NULL, y = "Standardized ridge coefficient") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 7.5, color = "grey35"))
imp <- copy(res)[, var := factor(var, levels = rev(res[order(-rf_imp)]$var))]
pB <- ggplot(imp, aes(var, rf_imp, fill = type)) + geom_col() + coord_flip() + scale_fill_manual(values = pal, guide = "none") +
  labs(title = "B  Random forest importance", subtitle = sprintf("OOB R2 = %.2f (permutation importance)", rf$r.squared), x = NULL, y = "Importance") +
  theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size = 7.5, color = "grey35"))
set.seed(42); rf0 <- ranger(as.formula(paste("phi ~", paste(VAR15, collapse = "+"))), data = dd, num.trees = 2000, seed = 42)
pc <- data.table(obs = dd$phi, pred = rf0$predictions)
pC <- ggplot(pc, aes(pred, obs)) + geom_abline(color = "grey60") + geom_point(alpha = 0.7, color = "#2166AC", size = 1.4) +
  geom_smooth(method = "lm", se = FALSE, color = "#C1121F", linewidth = 0.7, formula = y ~ x) +
  labs(title = "C  Observed vs out-of-bag predicted phi", subtitle = sprintf("n = %d | OOB R2 = %.2f | spatial-block CV R2 (RF/ridge) = %.2f/%.2f", n, rf$r.squared, perf[3]$rf_r2, perf[3]$ridge_r2),
       x = "RF out-of-bag prediction", y = "Observed phi (climate-adjusted, event-masked)") + theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size = 7.5, color = "grey35"))
p <- (pA | (pB / pC)) + plot_layout(widths = c(1.25, 1)) +
  plot_annotation(title = "Best phi driver model: climate-adjusted, SIF-event-masked phi (178 CCM-screened stations, 15 variables)",
                  theme = theme(plot.title = element_text(size = 11, face = "bold")))
ggsave(file.path(OUT, "recovery_phi_best_model_summary.png"), p, width = 12, height = 7.5, dpi = 300)
cat("\n完成。\n")
