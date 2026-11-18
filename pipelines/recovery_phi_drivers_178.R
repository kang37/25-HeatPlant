#!/usr/bin/env Rscript
# =============================================================================
# recovery_phi_drivers_178.R
#   对照 recovery_phi_drivers_117.R: 用178站(只满足12号脚本标准(1)——tp∈[0,8]
#   内至少1个tp同时Δρ>0且p_surr<0.05——不再额外要求标准(2)optimal_tp>=0)重跑
#   phi~18变量的driver分析, 看放宽"时间顺序自洽"这条标准, 站数变多之后功效
#   是否真的提升。178站本身不是新的正式判据, 只是这次做样本量敏感性对照用,
#   不写回12号脚本、不替换117站主线。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ranger)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. 178站清单: 只用12号脚本的标准(1), 不加标准(2) ----------------------
neg <- fread("data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv")
neg[, drho := rho - rho_min]
LC_BREAK_F <- "data_proc/output_hcsif_buf1000/landcover_break_stations.csv"
if (file.exists(LC_BREAK_F)) {
  lcb <- fread(LC_BREAK_F); excl <- lcb[has_break == TRUE]$stat_id
  neg <- neg[!meteo_stat %in% excl]
}
sig_range <- neg[tp >= 0 & tp <= 8]
sig_range[, pass1 := (drho > 0) & (p_surr < 0.05)]
sel <- unique(sig_range[pass1 == TRUE]$meteo_stat)
cat("178站候选清单 n =", length(sel), "(117站是这178站里再要求optimal_tp>=0的子集)\n")

# ---- 2. HCSIF+VPD距平构造, 逐站AR1/ARX拟合(与recovery_ar1_pure.R同一套方法) --
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)[meteo_stat %in% sel, .(meteo_stat, year, doy, date, SIF_buf1000)]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[meteo_stat %in% sel]
d   <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
setorder(d, meteo_stat, year, doy)

deseason_z <- function(x, doy) {
  clim <- tapply(x, doy, mean, na.rm = TRUE)
  a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (a - mean(a, na.rm = TRUE)) / s
}
d[, sif_anom := deseason_z(SIF_buf1000, doy), by = meteo_stat]
d[, vpd_z    := as.numeric(scale(vpd_mean)), by = meteo_stat]
d[, lag1     := shift(sif_anom, 1), by = .(meteo_stat, year)]

fit_one <- function(sub) {
  sub <- sub[!is.na(lag1) & !is.na(sif_anom) & !is.na(vpd_z)]
  if (nrow(sub) < 30) return(NULL)
  m1 <- tryCatch(lm(sif_anom ~ lag1, data = sub), error = function(e) NULL)
  m2 <- tryCatch(lm(sif_anom ~ lag1 + vpd_z, data = sub), error = function(e) NULL)
  if (is.null(m1) || is.null(m2)) return(NULL)
  co1 <- summary(m1)$coefficients; co2 <- summary(m2)$coefficients
  if (!("lag1" %in% rownames(co1)) || !("lag1" %in% rownames(co2))) return(NULL)
  list(phi1_pure = co1["lag1","Estimate"], phi2_arx = co2["lag1","Estimate"],
       beta_arx = co2["vpd_z","Estimate"], n = nrow(sub))
}
phi <- d[, fit_one(.SD), by = meteo_stat]
cat("成功拟合phi的站数:", nrow(phi), "/", length(sel), "\n")

# ---- 3. 18个CORE协变量, 合并 -------------------------------------------------
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt <- env$dt; PRED <- unlist(env$CORE, use.names = FALSE)

dd <- merge(phi, dt, by.x = "meteo_stat", by.y = "stat_id")
dd <- dd[complete.cases(dd[, PRED, with = FALSE])]
cat("协变量完整、可建模的站数(178站口径):", nrow(dd), "\n")

# ---- 4. RF + OLS, 和117站结果并排对比 ---------------------------------------
set.seed(42)
f1 <- as.formula(paste("phi1_pure ~", paste(PRED, collapse = "+")))
rf1 <- ranger(f1, data = dd, num.trees = 2000, importance = "permutation", seed = 42)
imp178 <- data.table(var = names(rf1$variable.importance), importance_178 = rf1$variable.importance)
setorder(imp178, -importance_178)
cat("\n=== RF(178站口径) OOB R2 =", round(rf1$r.squared,3), "===\n")
print(imp178)

ols1 <- lm(f1, data = dd)
cat("\n=== OLS(178站口径, R2=", round(summary(ols1)$r.squared,3), ") ===\n")
co1 <- summary(ols1)$coefficients
print(round(co1[order(co1[,4]), c(1,4)], 4))

# ---- 5. 并排对比: 178 vs 117(读之前已存的重要性表) ---------------------------
imp117 <- fread(file.path(OUT, "recovery_phi_drivers_117_importance.csv"))[, .(var, importance_117 = importance)]
cmp <- merge(imp178, imp117, by = "var")
setorder(cmp, -importance_178)
cat("\n=== 178 vs 117 重要性排序对比, Spearman rho =",
    round(cor(cmp$importance_178, cmp$importance_117, method="spearman"), 3), "===\n")
print(cmp)

fwrite(dd, file.path(OUT, "recovery_phi_drivers_178_data.csv"))
fwrite(cmp, file.path(OUT, "recovery_phi_drivers_178_vs_117_importance.csv"))

cat("\n=== 汇总: 117站 vs 178站 ===\n")
cat(sprintf("117站口径: n(建模)=86, RF OOB R2=0.326\n"))
cat(sprintf("178站口径: n(建模)=%d, RF OOB R2=%.3f\n", nrow(dd), rf1$r.squared))

# ---- 6. 图 -------------------------------------------------------------------
imp178[, var := factor(var, levels = rev(var))]
p <- ggplot(imp178, aes(var, importance_178)) +
  geom_col(fill = "#C1121F") + coord_flip() +
  labs(title = "RF variable importance for phi (178-station looser criterion, no optimal_tp>=0 requirement)",
       subtitle = sprintf("n(modeled)=%d | RF OOB R2=%.2f  (vs 117-station: n=86, R2=0.33)", nrow(dd), rf1$r.squared),
       x = NULL, y = "Permutation importance") +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size = 8.5, color = "grey35"))
ggsave(file.path(OUT, "recovery_phi_drivers_178_importance.png"), p, width = 8, height = 6, dpi = 300)
cat("\n-> recovery_phi_drivers_178_importance.png\n")
