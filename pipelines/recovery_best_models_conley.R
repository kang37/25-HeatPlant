#!/usr/bin/env Rscript
# =============================================================================
# recovery_best_models_conley.R
#   给目前两个最优模型(phi~117站x15变量; logOmega~57站x15变量+VPD+direction)
#   补上效应方向+显著性——用OLS+Conley空间HAC标准误(200km Bartlett核,
#   沿用tp_lag_drivers_117_final.R里已经验证过的conley_lm()), 不是普通OLS
#   (中国气象站空间自相关强, 普通OLS显著性会虚高)。和已经算好的RF重要性
#   拼在一张表里, 三件事(方向/显著性/重要性)一次看全。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ranger); library(car)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

conley_lm <- function(fit, lon, lat, cut = 200) {
  X <- model.matrix(fit); u <- residuals(fit)
  S <- X * u
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  dy <- outer(lat, lat, "-") * 111.32
  K  <- pmax(0, 1 - sqrt(dx^2 + dy^2) / cut)
  dim(K) <- c(length(lon), length(lon))
  bread <- solve(crossprod(X))
  V <- bread %*% (crossprod(S, K %*% S)) %*% bread
  n <- nrow(X); k <- ncol(X)
  V * (n / (n - k))
}
run_conley <- function(fit, d) {
  Vc <- conley_lm(fit, d$longitude, d$latitude, cut = 200)
  b <- coef(fit); se_c <- sqrt(diag(Vc))
  tval <- b / se_c; pval <- 2 * pt(-abs(tval), df = nrow(d) - length(b))
  res <- data.table(var = names(b), coef = round(b, 4), se_conley = round(se_c, 4),
                     t = round(tval, 2), p_conley = round(pval, 4))[var != "(Intercept)"]
  res[, sig := fifelse(p_conley < .01, "***", fifelse(p_conley < .05, "**", fifelse(p_conley < .1, "*", "")))]
  res[, direction := fifelse(coef > 0, "positive(+)", "negative(-)")]
  setorder(res, p_conley)
  res
}

env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18 <- env$dt
st   <- fread("data_raw/hcsif/stations_924.csv")
vc   <- fread(file.path(OUT, "vpd_station_covariate.csv"))[, .(meteo_stat, vpd_mean_stn, vpd_cv_stn)]
VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))

# =============================================================================
# 模型1: phi ~ 117站 x 15变量
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
dm1 <- merge(phi, dt18, by.x = "meteo_stat", by.y = "stat_id")
dm1 <- merge(dm1, st, by = "meteo_stat")
dm1 <- dm1[complete.cases(dm1[, VAR15, with = FALSE])]
cat("模型1 (phi, 117x15var) 建模站数:", nrow(dm1), "\n")

f1 <- as.formula(paste("phi ~", paste(VAR15, collapse = "+")))
ols1 <- lm(f1, data = dm1)
vif1 <- sort(vif(ols1), decreasing = TRUE)
cat("VIF(最大3个):", paste(names(vif1)[1:3], round(vif1[1:3],1), collapse=", "), "\n")
res1 <- run_conley(ols1, dm1)
set.seed(42); rf1 <- ranger(f1, data = dm1, num.trees = 2000, importance = "permutation", seed = 42)
imp1 <- data.table(var = names(rf1$variable.importance), rf_importance = rf1$variable.importance)
sm1 <- merge(res1, imp1, by = "var"); setorder(sm1, p_conley)
cat("\n=== 模型1: phi ~ 117站 x 15变量 (OLS R2=", round(summary(ols1)$r.squared,3),
    ", RF OOB R2=", round(rf1$r.squared,3), ") ===\n")
print(sm1[, .(var, coef, direction, p_conley, sig, rf_importance = round(rf_importance,4))])
fwrite(sm1, file.path(OUT, "recovery_best_phi_117_15var_conley.csv"))

# =============================================================================
# 模型2: logOmega ~ 57站 x 15变量+VPD+direction
# =============================================================================
evO <- fread(file.path(OUT, "recovery_liu2021_metrics_dirfix_117_events.csv"))
aggO <- evO[, .(logOmega = as.numeric(mean(log(Omega), na.rm=TRUE)), direction_2v = direction_2v[1]), by = meteo_stat]
dm2 <- merge(aggO, dt18, by.x = "meteo_stat", by.y = "stat_id")
dm2 <- merge(dm2, vc, by = "meteo_stat")
dm2 <- merge(dm2, st, by = "meteo_stat")
VAR2 <- c(VAR15, "vpd_mean_stn", "vpd_cv_stn")
dm2 <- dm2[complete.cases(dm2[, VAR2, with = FALSE])]
dm2[, is_promote := as.integer(direction_2v == "Promote")]
cat("\n模型2 (logOmega, 57x15var+VPD+direction) 建模站数:", nrow(dm2), "\n")

f2 <- as.formula(paste("logOmega ~", paste(c(VAR2, "is_promote"), collapse = "+")))
ols2 <- lm(f2, data = dm2)
vif2 <- sort(vif(ols2), decreasing = TRUE)
cat("VIF(最大3个):", paste(names(vif2)[1:3], round(vif2[1:3],1), collapse=", "), "\n")
res2 <- run_conley(ols2, dm2)
set.seed(42); rf2 <- ranger(f2, data = dm2, num.trees = 2000, importance = "permutation", seed = 42)
imp2 <- data.table(var = names(rf2$variable.importance), rf_importance = rf2$variable.importance)
sm2 <- merge(res2, imp2, by = "var"); setorder(sm2, p_conley)
cat("\n=== 模型2: logOmega ~ 57站 x 15变量+VPD+direction (OLS R2=", round(summary(ols2)$r.squared,3),
    ", RF OOB R2=", round(rf2$r.squared,3), ") ===\n")
print(sm2[, .(var, coef, direction, p_conley, sig, rf_importance = round(rf_importance,4))])
fwrite(sm2, file.path(OUT, "recovery_best_logOmega_57_conley.csv"))

# ---- 图: 系数+置信区间, 两个模型各一张 ---------------------------------------
plot_conley <- function(sm, title, fname) {
  sm2p <- copy(sm); sm2p[, var := factor(var, levels = rev(var[order(coef)]))]
  p <- ggplot(sm2p, aes(var, coef)) +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_point(aes(color = sig != ""), size = 2.5) +
    geom_errorbar(aes(ymin = coef - 1.96*se_conley, ymax = coef + 1.96*se_conley, color = sig != ""), width = .3) +
    scale_color_manual(values = c(`TRUE`="#C1121F", `FALSE`="grey50"), guide = "none") +
    coord_flip() + labs(title = title, x = NULL, y = "Coefficient (Conley spatial-HAC 95% CI)") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT, fname), p, width = 8, height = 6.5, dpi = 300)
}
plot_conley(sm1, "phi ~ 117 stations x 15 variables (Conley spatial-HAC SE)", "recovery_best_phi_117_15var_conley.png")
plot_conley(sm2, "logOmega ~ 57 stations x 15var+VPD+direction (Conley spatial-HAC SE)", "recovery_best_logOmega_57_conley.png")
cat("\n-> recovery_best_phi_117_15var_conley.png, recovery_best_logOmega_57_conley.png\n")
