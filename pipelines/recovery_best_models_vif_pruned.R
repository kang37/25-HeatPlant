#!/usr/bin/env Rscript
# =============================================================================
# recovery_best_models_vif_pruned.R
#   对phi(117x15var)和logOmega(57x15var+VPD+direction)两个模型做VIF驱动的变量
#   剪枝: 先固定丢近常数变量(mixedleaf/deci_needle, 91-97%站点近零, 已由
#   tp_lag_var_screening_117.R诊断过, 数值不稳定、不是真实效应来源), 再对剩下
#   的变量做逐步VIF剪枝(每次丢VIF最高的一个, 直到全部<5), 剪枝后用Conley空间
#   HAC重新做显著性检验+RF重要性对照。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ranger); library(car)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

conley_lm <- function(fit, lon, lat, cut = 200) {
  X <- model.matrix(fit); u <- residuals(fit)
  S <- X * u
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  dy <- outer(lat, lat, "-") * 111.32
  K  <- pmax(0, 1 - sqrt(dx^2 + dy^2) / cut); dim(K) <- c(length(lon), length(lon))
  bread <- solve(crossprod(X))
  V <- bread %*% (crossprod(S, K %*% S)) %*% bread
  n <- nrow(X); k <- ncol(X); V * (n / (n - k))
}
run_conley <- function(fit, d) {
  Vc <- conley_lm(fit, d$longitude, d$latitude, cut = 200)
  b <- coef(fit); se_c <- sqrt(diag(Vc))
  tval <- b / se_c; pval <- 2 * pt(-abs(tval), df = nrow(d) - length(b))
  res <- data.table(var = names(b), coef = round(b, 4), se_conley = round(se_c, 4),
                     t = round(tval, 2), p_conley = round(pval, 4))[var != "(Intercept)"]
  res[, sig := fifelse(p_conley < .01, "***", fifelse(p_conley < .05, "**", fifelse(p_conley < .1, "*", "")))]
  res[, direction := fifelse(coef > 0, "positive(+)", "negative(-)")]
  setorder(res, p_conley); res
}
vif_prune <- function(pred, data, thresh = 5) {
  repeat {
    f <- as.formula(paste("y ~", paste(pred, collapse = "+")))
    v <- vif(lm(f, data = data))
    if (max(v) < thresh || length(pred) <= 3) break
    drop_var <- names(which.max(v))
    cat("  VIF剪枝: 丢", drop_var, "(VIF=", round(max(v),1), ")\n")
    pred <- setdiff(pred, drop_var)
  }
  pred
}

env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18 <- env$dt; st <- fread("data_raw/hcsif/stations_924.csv")
vc <- fread(file.path(OUT, "vpd_station_covariate.csv"))[, .(meteo_stat, vpd_mean_stn, vpd_cv_stn)]
VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))
VAR13 <- setdiff(VAR15, c("mixedleaf", "deci_needle"))   # 先固定丢近常数变量

# =============================================================================
# 模型1: phi, 117站
# =============================================================================
sel117 <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")$meteo_stat
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)[meteo_stat %in% sel117, .(meteo_stat, year, doy, date, SIF_buf1000)]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[meteo_stat %in% sel117]
d117 <- merge(sif, vpd, by = c("meteo_stat", "year", "doy")); setorder(d117, meteo_stat, year, doy)
deseason_z <- function(x, doy) {
  clim <- tapply(x, doy, mean, na.rm = TRUE); a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x))); (a - mean(a, na.rm=TRUE)) / s
}
d117[, sif_anom := deseason_z(SIF_buf1000, doy), by = meteo_stat]
d117[, lag1 := shift(sif_anom, 1), by = .(meteo_stat, year)]
fit_one <- function(sub) {
  sub <- sub[!is.na(lag1) & !is.na(sif_anom)]; if (nrow(sub) < 30) return(NULL)
  m <- tryCatch(lm(sif_anom ~ lag1, data = sub), error = function(e) NULL); if (is.null(m)) return(NULL)
  co <- summary(m)$coefficients; if (!"lag1" %in% rownames(co)) return(NULL)
  list(phi = co["lag1", "Estimate"])
}
phi <- d117[, fit_one(.SD), by = meteo_stat]
dm1 <- merge(merge(phi, dt18, by.x = "meteo_stat", by.y = "stat_id"), st, by = "meteo_stat")
dm1 <- dm1[complete.cases(dm1[, VAR13, with = FALSE])]
dm1[, y := phi]
cat("=== 模型1: phi, 剪枝前13变量(已去mixedleaf/deci_needle), n=", nrow(dm1), "===\n")
PRED1 <- vif_prune(VAR13, dm1, thresh = 5)
cat("剪枝后剩", length(PRED1), "个变量:", paste(PRED1, collapse=", "), "\n")

f1 <- as.formula(paste("phi ~", paste(PRED1, collapse = "+")))
ols1 <- lm(f1, data = dm1)
res1 <- run_conley(ols1, dm1)
set.seed(42); rf1 <- ranger(f1, data = dm1, num.trees = 2000, importance = "permutation", seed = 42)
imp1 <- data.table(var = names(rf1$variable.importance), rf_importance = rf1$variable.importance)
sm1 <- merge(res1, imp1, by = "var"); setorder(sm1, p_conley)
cat("\n模型1最终: OLS R2=", round(summary(ols1)$r.squared,3), " RF OOB R2=", round(rf1$r.squared,3), "\n")
print(sm1[, .(var, coef, direction, p_conley, sig, rf_importance = round(rf_importance,4))])
fwrite(sm1, file.path(OUT, "recovery_best_phi_vifpruned_conley.csv"))

# =============================================================================
# 模型2: logOmega, 57站 (加VPD+direction, 先丢vpd_mean_stn保留vpd_cv_stn,
#   因为上次RF重要性里vpd_cv_stn远高于vpd_mean_stn)
# =============================================================================
evO <- fread(file.path(OUT, "recovery_liu2021_metrics_dirfix_117_events.csv"))
aggO <- evO[, .(logOmega = as.numeric(mean(log(Omega), na.rm=TRUE)), direction_2v = direction_2v[1]), by = meteo_stat]
dm2 <- merge(merge(aggO, dt18, by.x = "meteo_stat", by.y = "stat_id"), vc, by = "meteo_stat")
dm2 <- merge(dm2, st, by = "meteo_stat")
dm2[, is_promote := as.integer(direction_2v == "Promote")]
VAR13_vpd <- c(VAR13, "vpd_cv_stn", "is_promote")   # 只留vpd_cv_stn, 丢vpd_mean_stn
dm2 <- dm2[complete.cases(dm2[, VAR13_vpd, with = FALSE])]
dm2[, y := logOmega]
cat("\n=== 模型2: logOmega, 剪枝前", length(VAR13_vpd), "变量(已去mixedleaf/deci_needle/vpd_mean_stn), n=", nrow(dm2), "===\n")
PRED2 <- vif_prune(VAR13_vpd, dm2, thresh = 5)
cat("剪枝后剩", length(PRED2), "个变量:", paste(PRED2, collapse=", "), "\n")

f2 <- as.formula(paste("logOmega ~", paste(PRED2, collapse = "+")))
ols2 <- lm(f2, data = dm2)
res2 <- run_conley(ols2, dm2)
set.seed(42); rf2 <- ranger(f2, data = dm2, num.trees = 2000, importance = "permutation", seed = 42)
imp2 <- data.table(var = names(rf2$variable.importance), rf_importance = rf2$variable.importance)
sm2 <- merge(res2, imp2, by = "var"); setorder(sm2, p_conley)
cat("\n模型2最终: OLS R2=", round(summary(ols2)$r.squared,3), " RF OOB R2=", round(rf2$r.squared,3), "\n")
print(sm2[, .(var, coef, direction, p_conley, sig, rf_importance = round(rf_importance,4))])
fwrite(sm2, file.path(OUT, "recovery_best_logOmega_vifpruned_conley.csv"))

# ---- 图 -----------------------------------------------------------------------
plot_conley <- function(sm, title, fname) {
  sm2p <- copy(sm); sm2p[, var := factor(var, levels = rev(var[order(coef)]))]
  p <- ggplot(sm2p, aes(var, coef)) +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_point(aes(color = sig != ""), size = 2.5) +
    geom_errorbar(aes(ymin = coef - 1.96*se_conley, ymax = coef + 1.96*se_conley, color = sig != ""), width = .3) +
    scale_color_manual(values = c(`TRUE`="#C1121F", `FALSE`="grey50"), guide = "none") +
    coord_flip() + labs(title = title, x = NULL, y = "Coefficient (Conley spatial-HAC 95% CI)") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT, fname), p, width = 8, height = 6, dpi = 300)
}
plot_conley(sm1, "phi, VIF-pruned (117 stations)", "recovery_best_phi_vifpruned_conley.png")
plot_conley(sm2, "logOmega, VIF-pruned (57 stations, +direction)", "recovery_best_logOmega_vifpruned_conley.png")
cat("\n-> recovery_best_{phi,logOmega}_vifpruned_conley.png\n")
