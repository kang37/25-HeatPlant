#!/usr/bin/env Rscript
# =============================================================================
# recovery_omega_178_bydirection.R
#   方案: 用median_coef_2v符号(不用75%阈值)给178站分Promote/Inhibit两类(不留
#   Ambiguous), 用|frac_positive_2v-0.5|做置信度权重(离0.5越远=方向越不像
#   抛硬币=权重越大), 事件构造按方向翻转(同recovery_liu2021_metrics_dirfix_117.R
#   逻辑), 10变量精简集+VPD, 促进/抑制分开跑加权RF+加权岭回归+加权bootstrap。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ranger); library(glmnet)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. 178站方向(sign(median))+ 置信度权重 --------------------------------
allsig <- fread("data_proc/smap_bivar_134/smap_allsigtp_178.csv")
opt <- allsig[order(meteo_stat, -rho)][, .SD[1], by = meteo_stat]
opt[, direction_median := ifelse(median_coef_2v > 0, "Promote", "Inhibit")]
opt[, conf_weight := abs(frac_positive_2v - 0.5)]
stn_dir <- opt[, .(meteo_stat, direction_median, conf_weight, tp)]
cat("=== 178站新口径方向分布 ===\n"); print(stn_dir[, .N, by = direction_median])
cat("置信度权重分布:\n"); print(summary(stn_dir$conf_weight))

# ---- 2. HCSIF+VPD距平, 事件构造(按方向翻转, 同dirfix逻辑) -------------------
sel <- stn_dir$meteo_stat
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)[meteo_stat %in% sel, .(meteo_stat, year, doy, date, SIF_buf1000)]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[meteo_stat %in% sel]
d   <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
d   <- merge(d, stn_dir[, .(meteo_stat, direction_median)], by = "meteo_stat")
setorder(d, meteo_stat, year, doy)

clim_by_doy <- function(x, doy) tapply(x, doy, mean, na.rm = TRUE)
deseason_z <- function(x, doy) {
  clim <- clim_by_doy(x, doy); a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (a - mean(a, na.rm = TRUE)) / s
}
d[, sif_anom := deseason_z(SIF_buf1000, doy), by = meteo_stat]
d[, vpd_z    := as.numeric(scale(vpd_mean)), by = meteo_stat]
d[, sif_clim_doy := clim_by_doy(SIF_buf1000, doy)[as.character(doy)], by = meteo_stat]

calc_event <- function(sub, dirn1) {
  yy <- sub[order(doy)]; n <- nrow(yy)
  i0 <- which.max(yy$vpd_z)
  above <- yy$vpd_z > 0
  if (is.na(above[i0]) || !above[i0]) return(NULL)
  lo <- i0; while (lo > 1 && isTRUE(above[lo - 1])) lo <- lo - 1
  hi <- i0; while (hi < n && isTRUE(above[hi + 1])) hi <- hi + 1
  onset <- lo; end <- hi
  seg1 <- onset:n; seg2 <- end:n
  if (all(is.na(yy$sif_anom[seg1])) || all(is.na(yy$sif_anom[seg2]))) return(NULL)
  if (dirn1 == "Inhibit") {
    severe_i <- seg1[which.min(yy$sif_anom[seg1])]; recov_i <- seg2[which.max(yy$sif_anom[seg2])]
  } else {
    severe_i <- seg1[which.max(yy$sif_anom[seg1])]; recov_i <- seg2[which.min(yy$sif_anom[seg2])]
  }
  Yn_severe <- yy$sif_clim_doy[severe_i]; Ye <- yy$SIF_buf1000[severe_i]
  Yn_recov  <- yy$sif_clim_doy[recov_i];  Ye1 <- yy$SIF_buf1000[recov_i]
  dev_e <- Ye - Yn_severe; dev_e1 <- Ye1 - Yn_recov
  list(Omega = as.numeric(Yn_severe / abs(dev_e)))
}
ev <- d[, calc_event(.SD, direction_median[1]), by = .(meteo_stat, year, direction_median)]
ev <- ev[!is.null(Omega) & is.finite(Omega) & Omega > 0]
agg <- ev[, .(logOmega = as.numeric(mean(log(Omega), na.rm = TRUE))), by = .(meteo_stat, direction_median)]
agg <- merge(agg, stn_dir[, .(meteo_stat, conf_weight)], by = "meteo_stat")
cat("\n站-年事件数:", nrow(ev), " | 聚合后站数:", nrow(agg), "\n")

# ---- 3. 协变量: 10变量精简集 + VPD -------------------------------------------
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18 <- env$dt
VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))
VAR10 <- setdiff(VAR15, c("mixedleaf", "deci_needle", "tavg", "rh", "rsds_mean"))
vc <- fread(file.path(OUT, "vpd_station_covariate.csv"))[, .(meteo_stat, vpd_mean_stn, vpd_cv_stn)]
PRED <- c(VAR10, "vpd_mean_stn", "vpd_cv_stn")

dall <- merge(agg, dt18, by.x = "meteo_stat", by.y = "stat_id")
dall <- merge(dall, vc, by = "meteo_stat")
dall <- dall[complete.cases(dall[, PRED, with = FALSE])]
cat("协变量完整可建模站数:", nrow(dall), "\n"); print(dall[, .N, by = direction_median])
fwrite(dall, file.path(OUT, "recovery_omega_178_bydirection_data.csv"))

# ---- 4. 促进/抑制分开: 加权RF + 加权岭回归 + 加权bootstrap ------------------
run_group <- function(grp) {
  sub <- dall[direction_median == grp]
  w <- sub$conf_weight
  f <- as.formula(paste("logOmega ~", paste(PRED, collapse = "+")))
  set.seed(42)
  rf <- ranger(f, data = sub, case.weights = w, num.trees = 2000, importance = "permutation", seed = 42)
  imp <- data.table(group = grp, var = names(rf$variable.importance), importance = rf$variable.importance)

  X <- scale(as.matrix(sub[, PRED, with = FALSE])); y <- sub$logOmega
  set.seed(42); cv <- cv.glmnet(X, y, weights = w, alpha = 0, nfolds = 10, standardize = FALSE)
  ridge_coef <- as.numeric(coef(cv, s = "lambda.1se"))[-1]; names(ridge_coef) <- PRED
  pred <- predict(cv, newx = X, s = "lambda.1se")
  ridge_r2 <- 1 - sum(w * (y - pred)^2) / sum(w * (y - weighted.mean(y, w))^2)

  set.seed(1)
  boot <- replicate(500, {
    idx <- sample(nrow(sub), replace = TRUE)
    Xb <- X[idx, , drop = FALSE]; yb <- y[idx]; wb <- w[idx]
    tryCatch({
      cvb <- cv.glmnet(Xb, yb, weights = wb, alpha = 0, nfolds = 5, standardize = FALSE)
      as.numeric(coef(cvb, s = "lambda.1se"))[-1]
    }, error = function(e) rep(NA_real_, length(PRED)))
  })
  rownames(boot) <- PRED
  ci_lo <- apply(boot, 1, quantile, .025, na.rm = TRUE); ci_hi <- apply(boot, 1, quantile, .975, na.rm = TRUE)
  res <- data.table(group = grp, var = PRED, ridge_coef = round(ridge_coef, 4),
                     boot_ci_lo = round(ci_lo, 4), boot_ci_hi = round(ci_hi, 4),
                     ci_excludes_0 = (ci_lo > 0) | (ci_hi < 0))
  list(rf_r2 = rf$r.squared, ridge_r2 = ridge_r2, imp = imp, res = res[order(-abs(ridge_coef))], n = nrow(sub))
}

out_p <- run_group("Promote")
out_i <- run_group("Inhibit")

cat("\n=== Promote (n=", out_p$n, ") RF OOB R2=", round(out_p$rf_r2,3), " 岭回归R2=", round(out_p$ridge_r2,3), " ===\n")
print(out_p$res)
cat("\n=== Inhibit (n=", out_i$n, ") RF OOB R2=", round(out_i$rf_r2,3), " 岭回归R2=", round(out_i$ridge_r2,3), " ===\n")
print(out_i$res)

fwrite(rbind(out_p$res, out_i$res), file.path(OUT, "recovery_omega_178_bydirection_inference.csv"))
fwrite(rbind(out_p$imp, out_i$imp), file.path(OUT, "recovery_omega_178_bydirection_importance.csv"))

# ---- 5. 图: 两组系数+CI并排 --------------------------------------------------
res_all <- rbind(out_p$res, out_i$res)
res_all[, group := factor(group, levels = c("Promote", "Inhibit"))]
res_all[, var := factor(var, levels = rev(unique(res_all[order(-abs(ridge_coef))]$var))), by = group]
p <- ggplot(res_all, aes(var, ridge_coef, color = ci_excludes_0)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_pointrange(aes(ymin = boot_ci_lo, ymax = boot_ci_hi)) + coord_flip() +
  facet_wrap(~group, scales = "free_x", labeller = as_labeller(c(
    Promote = sprintf("Promote (n=%d, RF R2=%.2f)", out_p$n, out_p$rf_r2),
    Inhibit = sprintf("Inhibit (n=%d, RF R2=%.2f)", out_i$n, out_i$rf_r2)))) +
  scale_color_manual(values = c(`TRUE` = "#C1121F", `FALSE` = "grey50"), guide = "none") +
  labs(title = "logOmega drivers by direction (sign(median) 2-class, confidence-weighted)",
       subtitle = "10 decorrelated variables + VPD | weight = |frac_positive_2v - 0.5| | red = bootstrap 95% CI excludes 0",
       x = NULL, y = "Weighted standardized ridge coefficient") +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size = 8, color = "grey35"))
ggsave(file.path(OUT, "recovery_omega_178_bydirection_coef_ci.png"), p, width = 11, height = 6.5, dpi = 300)
cat("\n-> recovery_omega_178_bydirection_coef_ci.png\n")
