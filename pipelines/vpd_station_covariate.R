#!/usr/bin/env Rscript
# =============================================================================
# vpd_station_covariate.R
#   反复被标记的缺口: 18个CORE协变量里没有站点尺度的VPD均值/变异。补算:
#     vpd_mean_stn = 该站全部有效8天窗口(2000-2022, 5-9月)VPD均值
#     vpd_cv_stn   = 该站VPD的变异系数(sd/mean)
#   覆盖117站和178站两个候选清单(取并集算, 下游按需子集)。
# =============================================================================
suppressPackageStartupMessages(library(data.table))
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

sel117 <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")$meteo_stat
neg <- fread("data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv")
neg[, drho := rho - rho_min]
LC_BREAK_F <- "data_proc/output_hcsif_buf1000/landcover_break_stations.csv"
if (file.exists(LC_BREAK_F)) {
  lcb <- fread(LC_BREAK_F); excl <- lcb[has_break == TRUE]$stat_id
  neg <- neg[!meteo_stat %in% excl]
}
sig_range <- neg[tp >= 0 & tp <= 8]; sig_range[, pass1 := (drho > 0) & (p_surr < 0.05)]
sel178 <- unique(sig_range[pass1 == TRUE]$meteo_stat)
sel <- union(sel117, sel178)
cat("117站:", length(sel117), " | 178站:", length(sel178), " | 并集:", length(sel), "\n")

vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[meteo_stat %in% sel]
vc <- vpd[, .(vpd_mean_stn = mean(vpd_mean, na.rm = TRUE),
              vpd_cv_stn   = sd(vpd_mean, na.rm = TRUE) / mean(vpd_mean, na.rm = TRUE),
              n_windows    = .N), by = meteo_stat]
cat("\nvpd_mean_stn 分布:\n"); print(summary(vc$vpd_mean_stn))
cat("vpd_cv_stn 分布:\n"); print(summary(vc$vpd_cv_stn))

fwrite(vc, file.path(OUT, "vpd_station_covariate.csv"))
cat("\n-> vpd_station_covariate.csv (", nrow(vc), "站 )\n")
