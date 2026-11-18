#!/usr/bin/env Rscript
# =============================================================================
# recovery_omega_drivers_unscreened_allstations.R
#   logOmega最佳(干净)模型是57站(CCM筛选117站中方向明确的Promote/Inhibit)x
#   10变量精简集+VPD+is_promote方向哑变量。is_promote来自二变量S-map, 而
#   二变量S-map的状态空间锚定在每站CCM确认的optimal_tp上——这个量对"未经CCM
#   筛选"的站点根本不存在, 所以57站最佳模型不能直接搬到全量站点上比较。
#
#   能做到"同一套指标构造、同一套协变量、只改筛选与否"的公平对照是:
#   不用方向修正(dirfix)版本, 而用最初"假设VPD在所有站都抑制SIF"的直接版本
#   (recovery_liu2021_metrics_117.R的构造, 不需要S-map/CCM方向标签), 分别在
#   (a) 117 CCM筛选站 (b) 全量~904未筛选站 上都算一遍, 用10变量+VPD(不含
#   is_promote, 因为未筛选站没有这个量)跑同一套三方法, 比较筛选前后差异。
#   注意: 这个对照隔离的是"筛选"这一个变量, 不含方向修正的贡献(已在其他脚本
#   记录: dirfix让logOmega R2从0.33左右提升到0.57左右, 是另一个独立的效应)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(car); library(glmnet); library(ranger)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18  <- env$dt
VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))
VAR10 <- setdiff(VAR15, c("mixedleaf", "deci_needle", "tavg", "rh", "rsds_mean"))
cat("10变量精简集:", paste(VAR10, collapse=", "), "\n")

vc <- fread(file.path(OUT, "vpd_station_covariate.csv"))[, .(meteo_stat, vpd_mean_stn, vpd_cv_stn)]
ALLPRED <- c(VAR10, "vpd_mean_stn", "vpd_cv_stn")

# vpd_station_covariate.csv 只覆盖117/178站并集(178站), 对"未筛选全量站"这里
# 需要用vpd_8day_cache.rds里的全部站点重新算一份vpd_mean_stn/vpd_cv_stn,
# 构造与vpd_station_covariate.R完全一致, 只是不做sel117/sel178子集限制。
vpd_full <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))
vc_all <- vpd_full[, .(vpd_mean_stn = mean(vpd_mean, na.rm = TRUE),
                        vpd_cv_stn   = sd(vpd_mean, na.rm = TRUE) / mean(vpd_mean, na.rm = TRUE)),
                    by = meteo_stat]
cat("全量VPD协变量站数:", nrow(vc_all), "\n")

# ---- (a) 117 CCM筛选站: 复用已有的非方向修正事件表 --------------------------
ev117 <- fread(file.path(OUT, "recovery_liu2021_metrics_117_events.csv"))
agg117 <- ev117[, .(logOmega = as.numeric(mean(log(Omega), na.rm=TRUE))), by = meteo_stat]
do117 <- merge(agg117, dt18, by.x = "meteo_stat", by.y = "stat_id")
do117 <- merge(do117, vc, by = "meteo_stat")
do117 <- do117[complete.cases(do117[, ALLPRED, with = FALSE])]
cat("(a) 117站(CCM筛选, 非方向修正): 可建模站数", nrow(do117), "\n")

# ---- (b) 全量~904未筛选站: 同一套非方向修正构造 ------------------------------
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)[, .(meteo_stat, year, doy, date, SIF_buf1000)]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))
d   <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
setorder(d, meteo_stat, year, doy)
cat("SIF∩VPD 全量站数(未筛选):", uniqueN(d$meteo_stat), "\n")

clim_by_doy <- function(x, doy) tapply(x, doy, mean, na.rm = TRUE)
deseason_z <- function(x, doy) {
  clim <- clim_by_doy(x, doy); a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (a - mean(a, na.rm = TRUE)) / s
}
d[, sif_anom := deseason_z(SIF_buf1000, doy), by = meteo_stat]
d[, vpd_z    := as.numeric(scale(vpd_mean)), by = meteo_stat]
d[, sif_clim_doy := clim_by_doy(SIF_buf1000, doy)[as.character(doy)], by = meteo_stat]

calc_event <- function(sub) {
  yy <- sub[order(doy)]; n <- nrow(yy)
  i0 <- which.max(yy$vpd_z)
  above <- yy$vpd_z > 0
  if (is.na(above[i0]) || !above[i0]) return(NULL)
  lo <- i0; while (lo > 1 && isTRUE(above[lo - 1])) lo <- lo - 1
  hi <- i0; while (hi < n && isTRUE(above[hi + 1])) hi <- hi + 1
  onset <- lo; end <- hi
  seg1 <- onset:n
  if (all(is.na(yy$sif_anom[seg1]))) return(NULL)
  severe_i <- seg1[which.min(yy$sif_anom[seg1])]
  seg2 <- end:n
  if (all(is.na(yy$sif_anom[seg2]))) return(NULL)
  recov_i <- seg2[which.max(yy$sif_anom[seg2])]
  Yn_severe <- yy$sif_clim_doy[severe_i]; Ye <- yy$SIF_buf1000[severe_i]
  dev_e  <- Ye - Yn_severe
  Omega <- as.numeric(Yn_severe / abs(dev_e))
  list(Omega = Omega)
}
ev_all <- d[, calc_event(.SD), by = .(meteo_stat, year)]
cat("全量未筛选: 站-年事件数", nrow(ev_all), "| 站数", uniqueN(ev_all$meteo_stat), "\n")
n_inf <- sum(!is.finite(ev_all$Omega))
if (n_inf > 0) cat("剔除Omega非有限值(分母趋近0, 117筛选站未遇到过这个边界情况):", n_inf, "站-年\n")
ev_all <- ev_all[is.finite(Omega) & Omega > 0]

agg_all <- ev_all[, .(logOmega = as.numeric(mean(log(Omega), na.rm=TRUE))), by = meteo_stat]
do_all <- merge(agg_all, dt18, by.x = "meteo_stat", by.y = "stat_id")
do_all <- merge(do_all, vc_all, by = "meteo_stat")
do_all <- do_all[complete.cases(do_all[, ALLPRED, with = FALSE])]
cat("(b) 全量未筛选(非方向修正): 可建模站数", nrow(do_all), "\n")

run_battery <- function(dd, ALLPRED) {
  f <- as.formula(paste("logOmega ~", paste(ALLPRED, collapse = "+")))
  ols <- lm(f, data = dd)
  vif_v <- sort(vif(ols), decreasing = TRUE)
  co <- as.data.table(summary(ols)$coefficients, keep.rownames = "var")
  setnames(co, c("var", "ols_coef", "se", "t", "p")); co <- co[var != "(Intercept)"]
  co[, sig := ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,".","")))]
  ols_r2 <- summary(ols)$r.squared

  set.seed(42)
  rf <- ranger(f, data = dd, num.trees = 2000, importance = "permutation", seed = 42)
  imp <- data.table(var = names(rf$variable.importance), importance = rf$variable.importance)

  X <- scale(as.matrix(dd[, ALLPRED, with = FALSE])); y <- dd$logOmega
  set.seed(42); cv <- cv.glmnet(X, y, alpha = 0, nfolds = 10, standardize = FALSE)
  ridge_coef <- as.numeric(coef(cv, s = "lambda.1se"))[-1]; names(ridge_coef) <- ALLPRED
  pred <- predict(cv, X, s = "lambda.1se")
  ridge_r2 <- 1 - sum((y-pred)^2)/sum((y-mean(y))^2)

  set.seed(1)
  boot_mat <- replicate(500, {
    idx <- sample(nrow(dd), replace = TRUE)
    Xb <- X[idx,,drop=FALSE]; yb <- y[idx]
    tryCatch({ cvb<-cv.glmnet(Xb,yb,alpha=0,nfolds=5,standardize=FALSE); as.numeric(coef(cvb,s="lambda.1se"))[-1] },
             error=function(e) rep(NA_real_, length(ALLPRED)))
  })
  rownames(boot_mat) <- ALLPRED
  ci_lo <- apply(boot_mat, 1, quantile, .025, na.rm = TRUE)
  ci_hi <- apply(boot_mat, 1, quantile, .975, na.rm = TRUE)

  res <- data.table(var = ALLPRED, ridge_coef = round(ridge_coef, 4),
                     boot_ci_lo = round(ci_lo, 4), boot_ci_hi = round(ci_hi, 4),
                     ci_excludes_0 = (ci_lo>0)|(ci_hi<0))
  res <- merge(res, co[, .(var, ols_coef, p, sig)], by = "var")
  res <- res[order(-abs(ridge_coef))]
  list(n = nrow(dd), ols_r2 = ols_r2, rf_r2 = rf$r.squared, ridge_r2 = ridge_r2,
       vif = vif_v, res = res, imp = imp[order(-importance)])
}

out_117 <- run_battery(do117, ALLPRED)
out_all <- run_battery(do_all, ALLPRED)

cat("\n=== (a) 117站(CCM筛选, 非方向修正) ===\n")
cat("OLS R2 =", round(out_117$ols_r2,3), " | RF OOB R2 =", round(out_117$rf_r2,3),
    " | 岭回归(1se) R2 =", round(out_117$ridge_r2,3), "\n")
cat("VIF(降序):\n"); print(round(out_117$vif, 1))
cat("岭回归+bootstrap+OLS:\n"); print(out_117$res)

cat("\n=== (b) 全量", nrow(do_all), "站(未筛选, 非方向修正) ===\n")
cat("OLS R2 =", round(out_all$ols_r2,3), " | RF OOB R2 =", round(out_all$rf_r2,3),
    " | 岭回归(1se) R2 =", round(out_all$ridge_r2,3), "\n")
cat("VIF(降序):\n"); print(round(out_all$vif, 1))
cat("岭回归+bootstrap+OLS:\n"); print(out_all$res)

r2_tab <- data.table(set = c(sprintf("CCM-screened, 117 stations (n=%d)", out_117$n),
                              sprintf("Unscreened (n=%d)", out_all$n)),
                      ols_r2 = c(out_117$ols_r2, out_all$ols_r2),
                      rf_r2  = c(out_117$rf_r2,  out_all$rf_r2),
                      ridge_r2 = c(out_117$ridge_r2, out_all$ridge_r2))
cat("\n=== R2 对照(logOmega, 非方向修正版, 10var+VPD) ===\n"); print(r2_tab)

cmp <- merge(
  out_117$res[, .(var, ridge_coef_screened = ridge_coef, ci_excl0_screened = ci_excludes_0,
                   ols_p_screened = p, ols_sig_screened = sig)],
  out_all$res[, .(var, ridge_coef_unscreened = ridge_coef, ci_excl0_unscreened = ci_excludes_0,
                   ols_p_unscreened = p, ols_sig_unscreened = sig)],
  by = "var", all = TRUE
)
cmp[, sign_agree := sign(ridge_coef_screened) == sign(ridge_coef_unscreened)]
cmp <- cmp[order(-abs(ridge_coef_screened))]
cat("\n=== 筛选(117) vs 未筛选 岭回归系数对照 ===\n"); print(cmp)

fwrite(out_117$res, file.path(OUT, "recovery_omega_nodirfix_117_10var_vpd_inference.csv"))
fwrite(out_all$res, file.path(OUT, "recovery_omega_nodirfix_unscreened_10var_vpd_inference.csv"))
fwrite(r2_tab, file.path(OUT, "recovery_omega_screened_vs_unscreened_R2.csv"))
fwrite(cmp,    file.path(OUT, "recovery_omega_screened_vs_unscreened_compare.csv"))

# ---- 图 -----------------------------------------------------------------------
plot_dt <- rbindlist(list(
  out_all$res[, .(var, ridge_coef, boot_ci_lo, boot_ci_hi, ci_excludes_0, set = sprintf("Unscreened (n=%d)", out_all$n))],
  out_117$res[, .(var, ridge_coef, boot_ci_lo, boot_ci_hi, ci_excludes_0, set = sprintf("CCM-screened, 117 stations (n=%d)", out_117$n))]
))
scr_lab <- sprintf("CCM-screened, 117 stations (n=%d)", out_117$n)
uns_lab <- sprintf("Unscreened (n=%d)", out_all$n)
ord <- out_117$res[order(-abs(ridge_coef))]$var
plot_dt[, var := factor(var, levels = rev(ord))]
p <- ggplot(plot_dt, aes(var, ridge_coef, color = set, shape = ci_excludes_0)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_pointrange(aes(ymin = boot_ci_lo, ymax = boot_ci_hi), position = position_dodge(width = 0.5)) +
  coord_flip() +
  scale_color_manual(values = setNames(c("#C1121F", "#2166AC"), c(scr_lab, uns_lab))) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(title = "logOmega (non-dirfix) driver model: CCM-screened vs unscreened stations",
       subtitle = sprintf("screened ridge R2=%.3f | unscreened ridge R2=%.3f | filled point = bootstrap CI excludes 0\nnote: direction-fix (is_promote) not applicable to unscreened stations, excluded from both sides for fair comparison",
                           out_117$ridge_r2, out_all$ridge_r2),
       x = NULL, y = "Standardized ridge coefficient", color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", plot.subtitle = element_text(size=8, color="grey35"))
ggsave(file.path(OUT, "recovery_omega_screened_vs_unscreened_compare.png"), p, width = 9, height = 6.5, dpi = 300)
cat("\n-> recovery_omega_screened_vs_unscreened_compare.png / .csv\n")
