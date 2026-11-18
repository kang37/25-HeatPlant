#!/usr/bin/env Rscript
# =============================================================================
# recovery_phi_drivers_unscreened_allstations.R
#   目前phi(恢复力)最佳模型是178站(CCM筛选: tp 0-8内至少一个显著正向耦合)x15变量。
#   这里做同一套构造+同一套15变量+同一套三方法(OLS+VIF/RF/岭回归+bootstrap),
#   但改在"不经CCM筛选"的全量协变量完整站点(dt18, 904站)上跑一遍, 对照两者
#   差异是不是筛选本身带来的, 而不是方法差异。
#   phi(纯AR1系数)本身只用SIF自相关, 构造上不依赖VPD/CCM, 所以可以直接对全量
#   站点算, 是隔离"筛选效应"最干净的一个DV。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(car); library(glmnet); library(ranger)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 18个CORE协变量 + 15变量子集(与178站最佳模型完全一致的定义) -------------
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18  <- env$dt
VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))
cat("协变量完整站数(dt18):", nrow(dt18), "\n")

# ---- 全量站点(不做CCM筛选): SIF x VPD 全部站点, 不读stations_confirmed / negtp 文件 ---
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)[, .(meteo_stat, year, doy, date, SIF_buf1000)]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))
d   <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
setorder(d, meteo_stat, year, doy)
cat("SIF∩VPD 全量站数(未筛选):", uniqueN(d$meteo_stat), "\n")

deseason_z <- function(x, doy) {
  clim <- tapply(x, doy, mean, na.rm = TRUE); a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (a - mean(a, na.rm = TRUE)) / s
}
d[, sif_anom := deseason_z(SIF_buf1000, doy), by = meteo_stat]
d[, lag1     := shift(sif_anom, 1), by = .(meteo_stat, year)]
fit_one <- function(sub) {
  sub <- sub[!is.na(lag1) & !is.na(sif_anom)]
  if (nrow(sub) < 30) return(NULL)
  m <- tryCatch(lm(sif_anom ~ lag1, data = sub), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  co <- summary(m)$coefficients
  if (!"lag1" %in% rownames(co)) return(NULL)
  list(phi = co["lag1", "Estimate"])
}
phi_all <- d[, fit_one(.SD), by = meteo_stat]
cat("可拟合phi的全量站数(>=30个有效点):", nrow(phi_all), "\n")

dd <- merge(phi_all, dt18, by.x = "meteo_stat", by.y = "stat_id")
dd <- dd[complete.cases(dd[, VAR15, with = FALSE])]
cat("最终建模站数(全量, 未筛选CCM):", nrow(dd), "\n")

# ---- 178站CCM筛选清单(与178最佳模型完全一致的判据), 同一份phi_all里取子集 --
neg <- fread("data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv")
neg[, drho := rho - rho_min]
LC_BREAK_F <- "data_proc/output_hcsif_buf1000/landcover_break_stations.csv"
if (file.exists(LC_BREAK_F)) { lcb <- fread(LC_BREAK_F); excl <- lcb[has_break==TRUE]$stat_id; neg <- neg[!meteo_stat %in% excl] }
sig_range <- neg[tp>=0 & tp<=8]; sig_range[, pass1 := (drho>0)&(p_surr<0.05)]
sel178 <- unique(sig_range[pass1==TRUE]$meteo_stat)
dd178 <- dd[meteo_stat %in% sel178]
cat("178站筛选子集, 建模站数:", nrow(dd178), "\n")

run_battery <- function(dd, VAR15, tag) {
  f <- as.formula(paste("phi ~", paste(VAR15, collapse = "+")))
  ols <- lm(f, data = dd)
  vif_v <- sort(vif(ols), decreasing = TRUE)
  co <- as.data.table(summary(ols)$coefficients, keep.rownames = "var")
  setnames(co, c("var", "ols_coef", "se", "t", "p")); co <- co[var != "(Intercept)"]
  co[, sig := ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,".","")))]
  ols_r2 <- summary(ols)$r.squared

  set.seed(42)
  rf <- ranger(f, data = dd, num.trees = 2000, importance = "permutation", seed = 42)
  imp <- data.table(var = names(rf$variable.importance), importance = rf$variable.importance)

  X <- scale(as.matrix(dd[, VAR15, with = FALSE])); y <- dd$phi
  set.seed(42); cv <- cv.glmnet(X, y, alpha = 0, nfolds = 10, standardize = FALSE)
  ridge_coef <- as.numeric(coef(cv, s = "lambda.1se"))[-1]; names(ridge_coef) <- VAR15
  pred <- predict(cv, X, s = "lambda.1se")
  ridge_r2 <- 1 - sum((y-pred)^2)/sum((y-mean(y))^2)

  set.seed(1)
  boot_mat <- replicate(500, {
    idx <- sample(nrow(dd), replace = TRUE)
    Xb <- X[idx,,drop=FALSE]; yb <- y[idx]
    tryCatch({ cvb<-cv.glmnet(Xb,yb,alpha=0,nfolds=5,standardize=FALSE); as.numeric(coef(cvb,s="lambda.1se"))[-1] },
             error=function(e) rep(NA_real_, length(VAR15)))
  })
  rownames(boot_mat) <- VAR15
  ci_lo <- apply(boot_mat, 1, quantile, .025, na.rm = TRUE)
  ci_hi <- apply(boot_mat, 1, quantile, .975, na.rm = TRUE)

  res <- data.table(var = VAR15, ridge_coef = round(ridge_coef, 4),
                     boot_ci_lo = round(ci_lo, 4), boot_ci_hi = round(ci_hi, 4),
                     ci_excludes_0 = (ci_lo>0)|(ci_hi<0))
  res <- merge(res, co[, .(var, ols_coef, p, sig)], by = "var")
  res <- res[order(-abs(ridge_coef))]
  list(tag = tag, n = nrow(dd), ols_r2 = ols_r2, rf_r2 = rf$r.squared, ridge_r2 = ridge_r2,
       vif = vif_v, res = res, imp = imp[order(-importance)])
}

out_all <- run_battery(dd,    VAR15, "unscreened_allstations")
out_scr <- run_battery(dd178, VAR15, "ccm_screened_178")

cat("\n=== 未筛选(全量", nrow(dd), "站) ===\n")
cat("OLS R2 =", round(out_all$ols_r2,3), " | RF OOB R2 =", round(out_all$rf_r2,3),
    " | 岭回归(1se) R2 =", round(out_all$ridge_r2,3), "\n")
cat("\nVIF(降序):\n"); print(round(out_all$vif, 1))
cat("\n岭回归+bootstrap+OLS汇总:\n"); print(out_all$res)
cat("\nRF重要性:\n"); print(out_all$imp)

cat("\n=== CCM筛选(178站)重跑(seed一致, 与已存档的178最佳模型互相印证) ===\n")
cat("OLS R2 =", round(out_scr$ols_r2,3), " | RF OOB R2 =", round(out_scr$rf_r2,3),
    " | 岭回归(1se) R2 =", round(out_scr$ridge_r2,3), "\n")

fwrite(out_all$res, file.path(OUT, "recovery_phi_unscreened_allstations_inference.csv"))
fwrite(out_all$imp, file.path(OUT, "recovery_phi_unscreened_allstations_rf_importance.csv"))
fwrite(data.table(var = names(out_all$vif), VIF = out_all$vif),
       file.path(OUT, "recovery_phi_unscreened_allstations_VIF.csv"))

# ---- 筛选 vs 未筛选 对照 -----------------------------------------------------
cmp <- merge(
  out_scr$res[, .(var, ridge_coef_screened = ridge_coef, ci_excl0_screened = ci_excludes_0,
                   ols_p_screened = p, ols_sig_screened = sig)],
  out_all$res[, .(var, ridge_coef_unscreened = ridge_coef, ci_excl0_unscreened = ci_excludes_0,
                   ols_p_unscreened = p, ols_sig_unscreened = sig)],
  by = "var", all = TRUE
)
cmp[, sign_agree := sign(ridge_coef_screened) == sign(ridge_coef_unscreened)]
cmp <- cmp[order(-abs(ridge_coef_screened))]
cat("\n=== 筛选(178站) vs 未筛选(", nrow(dd), "站) 岭回归系数对照 ===\n")
print(cmp)
fwrite(cmp, file.path(OUT, "recovery_phi_screened_vs_unscreened_compare.csv"))

r2_tab <- data.table(set = c(sprintf("CCM-screened (n=%d)", nrow(dd178)), sprintf("Unscreened (n=%d)", nrow(dd))),
                      ols_r2 = c(out_scr$ols_r2, out_all$ols_r2),
                      rf_r2  = c(out_scr$rf_r2,  out_all$rf_r2),
                      ridge_r2 = c(out_scr$ridge_r2, out_all$ridge_r2))
cat("\n=== R2 对照 ===\n"); print(r2_tab)
fwrite(r2_tab, file.path(OUT, "recovery_phi_screened_vs_unscreened_R2.csv"))

# ---- 图: 系数对照 ------------------------------------------------------------
plot_dt <- rbindlist(list(
  out_all$res[, .(var, ridge_coef, boot_ci_lo, boot_ci_hi, ci_excludes_0, set = sprintf("Unscreened (n=%d)", nrow(dd)))],
  out_scr$res[, .(var, ridge_coef, boot_ci_lo, boot_ci_hi, ci_excludes_0, set = sprintf("CCM-screened, 178 stations (n=%d)", nrow(dd178)))]
))
ord <- out_scr$res[order(-abs(ridge_coef))]$var
plot_dt[, var := factor(var, levels = rev(ord))]
scr_lab <- sprintf("CCM-screened, 178 stations (n=%d)", nrow(dd178))
uns_lab <- sprintf("Unscreened (n=%d)", nrow(dd))
p <- ggplot(plot_dt, aes(var, ridge_coef, color = set, shape = ci_excludes_0)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_pointrange(aes(ymin = boot_ci_lo, ymax = boot_ci_hi), position = position_dodge(width = 0.5), fatten = 3) +
  coord_flip() +
  scale_color_manual(values = setNames(c("#C1121F", "#2166AC"), c(scr_lab, uns_lab))) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(title = "phi driver model: CCM-screened vs unscreened stations",
       subtitle = sprintf("screened ridge R2=%.3f | unscreened ridge R2=%.3f | filled point = bootstrap CI excludes 0",
                           out_scr$ridge_r2, out_all$ridge_r2),
       x = NULL, y = "Standardized ridge coefficient", color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "recovery_phi_screened_vs_unscreened_compare.png"), p, width = 9, height = 7, dpi = 300)
cat("\n-> recovery_phi_screened_vs_unscreened_compare.png / .csv\n")
