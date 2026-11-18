#!/usr/bin/env Rscript
# =============================================================================
# recovery_phi_climate_adjusted_ED.R
#   借鉴 Guo et al. 2026 (PNAS 123:e2606997123) 的两点, 重算178站的φ并检验"驱动因素结论
#   是否依赖φ的算法":
#   (1) 扣除气候强迫: Yt = a0*Y(t-1) + a1*T + a2*P + a3*S, 气候分量=a1T+a2P+a3S,
#       Y_resilience = Y - 气候分量, 再对Y_resilience算lag-1自相关。
#       版本A = 论文原版(气温/降水/辐射); 版本B = A + VPD (VPD本身是研究对象, B仅作对照)。
#   (2) TAC_ED: 遮掉扰动期再算自相关。扰动期 = 含|z|>=2的窗口、向两侧延伸到|z|<0.5
#       (同符号连续段), 按站-年在年内处理。ED_sif=用SIF距平识别(方向无关);
#       ED_vpd=用VPD去季节距平识别(只取正向, 高VPD)。
#   所有版本共用同一批"气候变量齐全"的窗口, lag只在年内相邻窗口取; 被遮掉/缺失窗口所在
#   的配对整体丢弃(不会把不相邻的窗口错配)。基线phi_raw_cc = 同样窗口、不做任何调整。
#   额外对照: drought_phys_stage1_arx_allstations.R 里联合估计的ARX φ(phi_arxA)。
#   比较: 各版本φ两两相关/噪声占比/与气候均值的相关; 同一套15变量driver模型
#   (OLS+VIF, RF, 岭回归+bootstrap, y逐版本标准化使系数可比)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(car); library(glmnet); library(ranger)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- Sys.getenv("PHI_OUT", file.path(PROJ, "data_proc/output_loose_classify")); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
NBOOT <- 300L

# ---- 178站 + 协变量 -----------------------------------------------------------
neg <- fread("data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv"); neg[, drho := rho - rho_min]
LCB <- "data_proc/output_hcsif_buf1000/landcover_break_stations.csv"
if (file.exists(LCB)) { lcb <- fread(LCB); neg <- neg[!meteo_stat %in% lcb[has_break == TRUE]$stat_id] }
sr <- neg[tp >= 0 & tp <= 8]; sr[, pass1 := (drho > 0) & (p_surr < 0.05)]
sel178 <- unique(sr[pass1 == TRUE]$meteo_stat)
sel117 <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")$meteo_stat
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]; cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18  <- env$dt
VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))
CLIM  <- c("tavg", "rh", "precip", "rsds_mean", "rsds_sd", "elev")

# ---- 8天窗口: SIF / VPD / rsds / tavg / precip (对齐方式同drought_phys_stage1_arx_allstations.R) ----
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread, select = c("meteo_stat", "year", "doy", "SIF_buf1000")), fill = TRUE)[meteo_stat %in% sel178]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[meteo_stat %in% sel178, .(meteo_stat, year, doy, vpd = vpd_mean)]
rs  <- rbindlist(lapply(2000:2022, function(y) fread(sprintf("data_raw/covariates_1km/era5_rsds_sif8d_%d.csv", y))))
rs  <- rs[stat_id %in% sel178, .(meteo_stat = stat_id, year, doy, rsds = rsds_mean)]
read_meteo_8d <- function(sid) {
  fp <- file.path("data_raw/meteo_data_1961-2023", paste0(sid, ".txt")); if (!file.exists(fp)) return(NULL)
  d <- fread(fp, skip = 1, showProgress = FALSE)
  if (!all(c("date", "tavg", "precip") %in% names(d))) return(NULL)
  d[, date := as.IDate(date)]; d <- d[year(date) %in% 2000:2022]
  d[!is.finite(tavg) | tavg < -60 | tavg > 60, tavg := NA_real_]
  d[!is.finite(precip) | precip < 0 | precip > 2000, precip := NA_real_]
  d[, yd := yday(date)]; d <- d[yd >= 124 & yd <= 275]
  d[, doy := 124L + 8L * ((yd - 124L) %/% 8L)]; d[, year := year(date)]
  d[, .(tavg = if (sum(!is.na(tavg)) >= 6) mean(tavg, na.rm = TRUE) else NA_real_,
        prcp = if (sum(!is.na(precip)) >= 6) sum(precip, na.rm = TRUE) * 8 / sum(!is.na(precip)) else NA_real_),
    by = .(year, doy)][, meteo_stat := sid][]
}
met <- rbindlist(lapply(sel178, read_meteo_8d))
d <- Reduce(function(a, b) merge(a, b, by = c("meteo_stat", "year", "doy")), list(sif, vpd, rs, met))
setorder(d, meteo_stat, year, doy)
EXY <- Sys.getenv("EXCLUDE_YEARS", "")
if (nzchar(EXY)) { d <- d[!year %in% as.integer(strsplit(EXY, ",")[[1]])]; cat("排除年份:", EXY, "\n") }
cat("178站窗口表:", nrow(d), "行, 站数", uniqueN(d$meteo_stat), "\n")

deseason_z <- function(x, doy) {
  clim <- tapply(x, doy, mean, na.rm = TRUE); a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (a - mean(a, na.rm = TRUE)) / s }
d[, `:=`(sif_z = deseason_z(SIF_buf1000, doy), vpd_z = deseason_z(vpd, doy), tavg_z = deseason_z(tavg, doy),
         rsds_z = deseason_z(rsds, doy), prcp_z = deseason_z(log1p(prcp), doy)), by = meteo_stat]
d[, cc := complete.cases(sif_z, vpd_z, tavg_z, rsds_z, prcp_z)]
cat("气候变量齐全窗口占比:", round(100 * mean(d$cc), 1), "%\n")

# ---- 扰动期遮罩(按站-年, 年内) -------------------------------------------------
mask_ed <- function(z, thr_ext = 2, thr_base = 0.5, pos_only = FALSE) {
  n <- length(z); m <- rep(FALSE, n); s <- sign(z)
  ok <- !is.na(z) & abs(z) >= thr_base & (if (pos_only) z > 0 else TRUE)
  i <- 1L
  while (i <= n) {
    if (isTRUE(ok[i])) { j <- i
      while (j < n && isTRUE(ok[j + 1L]) && s[j + 1L] == s[i]) j <- j + 1L
      if (max(abs(z[i:j])) >= thr_ext) m[i:j] <- TRUE
      i <- j + 1L } else i <- i + 1L }
  m }
d[, m_sif := mask_ed(sif_z), by = .(meteo_stat, year)]
d[, m_vpd := mask_ed(vpd_z, pos_only = TRUE), by = .(meteo_stat, year)]
cat(sprintf("遮罩占比(气候齐全窗口): ED_sif=%.1f%% | ED_vpd=%.1f%% | 并集=%.1f%%\n",
            100 * mean(d[cc == TRUE]$m_sif), 100 * mean(d[cc == TRUE]$m_vpd), 100 * mean(d[cc == TRUE]$m_sif | d[cc == TRUE]$m_vpd)))

# ---- 逐站: 气候调整(论文式) + AR1 ---------------------------------------------
ar1 <- function(y, yr) {
  lg <- data.table::shift(y, 1); same <- c(FALSE, yr[-1] == yr[-length(yr)]); lg[!same] <- NA_real_
  ok <- !is.na(y) & !is.na(lg); n <- sum(ok); if (n < 30) return(list(phi = NA_real_, se = NA_real_, n = n))
  co <- summary(lm(y[ok] ~ lg[ok]))$coefficients; list(phi = co[2, 1], se = co[2, 2], n = n) }
one_station <- function(sub) {
  sub <- sub[order(year, doy)]
  sub[, lag1 := { lg <- data.table::shift(sif_z, 1); lg[c(TRUE, year[-1] != year[-.N])] <- NA_real_; lg }]
  fit <- sub[cc == TRUE & !is.na(lag1)]
  if (nrow(fit) < 100) return(NULL)
  mA <- lm(sif_z ~ lag1 + tavg_z + prcp_z + rsds_z, data = fit)
  mB <- lm(sif_z ~ lag1 + tavg_z + prcp_z + rsds_z + vpd_z, data = fit)
  cA <- coef(mA); cB <- coef(mB)
  yA <- sub$sif_z - (cA["tavg_z"] * sub$tavg_z + cA["prcp_z"] * sub$prcp_z + cA["rsds_z"] * sub$rsds_z)
  yB <- sub$sif_z - (cB["tavg_z"] * sub$tavg_z + cB["prcp_z"] * sub$prcp_z + cB["rsds_z"] * sub$rsds_z + cB["vpd_z"] * sub$vpd_z)
  cc <- sub$cc; yr <- sub$year
  y0 <- fifelse(cc, sub$sif_z, NA_real_); yA <- fifelse(cc, yA, NA_real_); yB <- fifelse(cc, yB, NA_real_)
  msif <- sub$m_sif; mvpd <- sub$m_vpd
  r <- list(raw_cc = ar1(y0, yr), adjA = ar1(yA, yr), adjB = ar1(yB, yr),
            ED_sif = ar1(fifelse(msif, NA_real_, y0), yr), ED_vpd = ar1(fifelse(mvpd, NA_real_, y0), yr),
            adjA_ED_sif = ar1(fifelse(msif, NA_real_, yA), yr))
  out <- list(n_cc = sum(cc), g_tavg = cA["tavg_z"], g_prcp = cA["prcp_z"], g_rsds = cA["rsds_z"], r2_climate_added = summary(mA)$r.squared)
  for (nm in names(r)) { out[[paste0("phi_", nm)]] <- r[[nm]]$phi; out[[paste0("se_", nm)]] <- r[[nm]]$se; out[[paste0("n_", nm)]] <- r[[nm]]$n }
  out }
res <- d[, one_station(.SD), by = meteo_stat]

# 原始phi(与既往178最佳模型同口径, 未要求气候齐全)与联合ARX的phi
raw_all <- d[, {
  lg <- data.table::shift(sif_z, 1); ok <- !is.na(lg) & !is.na(sif_z) & c(FALSE, year[-1] == year[-.N])
  if (sum(ok) < 30) NULL else list(phi_raw_all = unname(coef(lm(sif_z[ok] ~ lg[ok]))[2])) }, by = meteo_stat]
arx <- fread(file.path(PROJ, "data_proc/output_loose_classify/drought_phys_arx_station_beta.csv"))[, .(meteo_stat, phi_arxA = phi)]
res <- Reduce(function(a, b) merge(a, b, by = "meteo_stat", all.x = TRUE), list(res, raw_all, arx))
fwrite(res, file.path(OUT, "recovery_phi_variants_178.csv"))
V <- c("raw_all", "raw_cc", "adjA", "adjB", "ED_sif", "ED_vpd", "adjA_ED_sif", "arxA")
PH <- paste0("phi_", V)
cat("\n站数:", nrow(res), "| 各版本有效站数:", paste(sapply(PH, function(p) sum(is.finite(res[[p]]))), collapse = "/"), "\n")

# ---- 版本间比较 ----------------------------------------------------------------
cat("\n=== 描述: 均值/SD/平均SE/噪声占比(mean se^2/var) ===\n")
desc <- rbindlist(lapply(V, function(v) { p <- res[[paste0("phi_", v)]]; s <- res[[paste0("se_", v)]]
  data.table(variant = v, n_stn = sum(is.finite(p)), mean = mean(p, na.rm = TRUE), sd = sd(p, na.rm = TRUE),
             mean_se = if (is.null(s)) NA_real_ else mean(s, na.rm = TRUE),
             noise_share = if (is.null(s)) NA_real_ else mean(s^2, na.rm = TRUE) / var(p, na.rm = TRUE),
             mean_pairs = if (is.null(res[[paste0("n_", v)]])) NA_real_ else mean(res[[paste0("n_", v)]], na.rm = TRUE)) }))
print(desc[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat("\n=== 相关(Spearman, 站间) ===\n")
cm <- cor(res[, ..PH], use = "pairwise.complete.obs", method = "spearman"); print(round(cm, 2))
cat("\n=== 各版本φ与站点气候均值的Spearman相关(气候依赖是否下降) ===\n")
cl <- merge(res[, c("meteo_stat", PH), with = FALSE], dt18[, c("stat_id", "tavg", "rh", "precip", "rsds_mean"), with = FALSE], by.x = "meteo_stat", by.y = "stat_id")
cor_clim <- rbindlist(lapply(V, function(v) { p <- cl[[paste0("phi_", v)]]
  data.table(variant = v, tavg = cor(p, cl$tavg, method = "spearman", use = "c"), rh = cor(p, cl$rh, method = "spearman", use = "c"),
             precip = cor(p, cl$precip, method = "spearman", use = "c"), rsds_mean = cor(p, cl$rsds_mean, method = "spearman", use = "c")) }))
print(cor_clim[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
cat(sprintf("\n气候强迫回归的R2(加入T/P/S后SIF距平解释度)中位数: %.3f | g_tavg/g_prcp/g_rsds中位数: %.3f/%.3f/%.3f\n",
            median(res$r2_climate_added), median(res$g_tavg), median(res$g_prcp), median(res$g_rsds)))
fwrite(desc, file.path(OUT, "recovery_phi_variants_desc.csv")); fwrite(as.data.table(cm, keep.rownames = "v"), file.path(OUT, "recovery_phi_variants_cor.csv"))
fwrite(cor_clim, file.path(OUT, "recovery_phi_variants_clim_cor.csv"))

# ---- 同一套driver模型(15变量; y逐版本标准化使系数可比) -------------------------
battery <- function(dd, y, VARS) {
  dd <- copy(dd); dd[, yy := as.numeric(scale(y))]
  f <- as.formula(paste("yy ~", paste(VARS, collapse = "+")))
  ols <- lm(f, data = dd); vv <- sort(vif(ols), decreasing = TRUE); ols_r2 <- summary(ols)$r.squared
  set.seed(42); rf <- ranger(f, data = dd, num.trees = 2000, importance = "permutation", seed = 42)
  X <- scale(as.matrix(dd[, VARS, with = FALSE])); yv <- dd$yy
  set.seed(42); cv <- cv.glmnet(X, yv, alpha = 0, nfolds = 10, standardize = FALSE)
  rc <- as.numeric(coef(cv, s = "lambda.1se"))[-1]; names(rc) <- VARS
  pr <- predict(cv, X, s = "lambda.1se"); rr2 <- 1 - sum((yv - pr)^2) / sum((yv - mean(yv))^2)
  set.seed(1)
  bm <- replicate(NBOOT, { ix <- sample(nrow(dd), replace = TRUE)
    tryCatch(as.numeric(coef(cv.glmnet(X[ix, , drop = FALSE], yv[ix], alpha = 0, nfolds = 5, standardize = FALSE), s = "lambda.1se"))[-1],
             error = function(e) rep(NA_real_, length(VARS))) })
  lo <- apply(bm, 1, quantile, .025, na.rm = TRUE); hi <- apply(bm, 1, quantile, .975, na.rm = TRUE)
  list(n = nrow(dd), ols_r2 = ols_r2, rf_r2 = rf$r.squared, ridge_r2 = rr2, max_vif = max(vv),
       coef = data.table(var = VARS, ridge = rc, lo = lo, hi = hi, robust = (lo > 0) | (hi < 0),
                         rf_imp = rf$variable.importance[VARS])) }
ALLV <- c("raw_all", "raw_cc", "adjA", "adjB", "ED_sif", "ED_vpd", "adjA_ED_sif", "arxA")
dd0 <- merge(res, dt18, by.x = "meteo_stat", by.y = "stat_id")
runs <- list(); sumtab <- list()
for (set_nm in c("178", "117")) {
  for (v in ALLV) {
    dd <- dd0[if (set_nm == "117") meteo_stat %in% sel117 else TRUE]
    dd <- dd[is.finite(get(paste0("phi_", v)))]; dd <- dd[complete.cases(dd[, VAR15, with = FALSE])]
    b <- battery(dd, dd[[paste0("phi_", v)]], VAR15)
    b$coef[, `:=`(variant = v, set = set_nm)]; runs[[paste(set_nm, v)]] <- b$coef
    sumtab[[paste(set_nm, v)]] <- data.table(set = set_nm, variant = v, n = b$n, ols_r2 = b$ols_r2, rf_r2 = b$rf_r2, ridge_r2 = b$ridge_r2,
      n_robust = sum(b$coef$robust), robust_vars = paste(b$coef[robust == TRUE][order(-abs(ridge))]$var, collapse = ","),
      clim_share = sum(abs(b$coef[var %in% CLIM]$ridge)) / sum(abs(b$coef$ridge)))
  }
}
st <- rbindlist(sumtab); cf <- rbindlist(runs)
cat("\n=== driver模型对比(15变量, y标准化) ===\n"); print(st[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
fwrite(st, file.path(OUT, "recovery_phi_variants_driver_summary.csv")); fwrite(cf, file.path(OUT, "recovery_phi_variants_driver_coefs.csv"))
w <- dcast(cf[set == "178"], var ~ variant, value.var = "ridge")[, c("var", ALLV), with = FALSE]
cat("\n=== 178站: 各版本岭回归系数(标准化) ===\n"); print(w[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

# ---- 图 ----------------------------------------------------------------------
lab <- c(raw_all = "raw (original rows)", raw_cc = "raw (same rows)", adjA = "climate-adj A (T,P,S)", adjB = "climate-adj B (+VPD)",
         ED_sif = "ED (mask SIF events)", ED_vpd = "ED (mask VPD events)", adjA_ED_sif = "adj A + ED (SIF)", arxA = "joint ARX phi")
pl <- rbindlist(lapply(setdiff(V, "raw_cc"), function(v) data.table(variant = lab[[v]], x = res$phi_raw_cc, y = res[[paste0("phi_", v)]])))
pl[, variant := factor(variant, levels = lab[setdiff(V, "raw_cc")])]
p1 <- ggplot(pl, aes(x, y)) + geom_abline(color = "grey60") + geom_point(alpha = 0.6, size = 1, color = "#2166AC") +
  geom_smooth(method = "lm", se = FALSE, color = "#C1121F", linewidth = 0.7, formula = y ~ x) + facet_wrap(~variant, ncol = 4) +
  labs(title = "phi variants vs. baseline phi (same rows, no adjustment), 178 stations", x = "baseline phi (raw, same rows)", y = "variant phi") + theme_bw(base_size = 10)
ggsave(file.path(OUT, "recovery_phi_variants_scatter.png"), p1, width = 11, height = 6, dpi = 300)
c178 <- cf[set == "178"]; c178[, variant := factor(lab[variant], levels = lab[ALLV])]
c178[, var := factor(var, levels = rev(w[order(-abs(raw_cc))]$var))]
p2 <- ggplot(c178, aes(variant, var, fill = ridge)) + geom_tile(color = "white") +
  geom_text(aes(label = ifelse(robust, sprintf("%.2f*", ridge), sprintf("%.2f", ridge))), size = 2.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#C1121F", midpoint = 0, name = "ridge coef\n(std. y)") +
  labs(title = "Driver model (178 stations, 15 variables): ridge coefficients by phi variant", subtitle = "* = bootstrap 95% CI excludes 0", x = NULL, y = NULL) +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT, "recovery_phi_variants_coef_heatmap.png"), p2, width = 10, height = 6.5, dpi = 300)
s178 <- st[set == "178"]; s178[, variant := factor(lab[variant], levels = lab[ALLV])]
s2 <- melt(s178[, .(variant, `RF OOB R2` = rf_r2, `Ridge R2` = ridge_r2, `Climate share of |coef|` = clim_share)], id.vars = "variant")
p3 <- ggplot(s2, aes(variant, value, fill = variant)) + geom_col(show.legend = FALSE) + facet_wrap(~variable, scales = "free_y", ncol = 3) +
  labs(title = "Driver-model fit and climate share by phi variant (178 stations)", x = NULL, y = NULL) + theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(OUT, "recovery_phi_variants_fit.png"), p3, width = 11, height = 4.8, dpi = 300)
cat("\n完成。\n")
