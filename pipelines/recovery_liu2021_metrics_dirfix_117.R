#!/usr/bin/env Rscript
# =============================================================================
# recovery_liu2021_metrics_dirfix_117.R
#   修正 recovery_liu2021_metrics_117.R 的一个方法漏洞: 原版无论方向一律假设
#   "VPD高=胁迫=SIF应下探", 用 which.min(sif_anom)找"最严重日"、which.max找
#   "最佳恢复日"——这套逻辑只对 Inhibit 站成立。117站里 Promote 41 / Ambiguous
#   60 / Inhibit 16(见 classify_134.csv), 超过1/3是"VPD升高→SIF反而升高"的站,
#   对这些站硬套"找最低点"逻辑, 找到的大概率是噪声, 不是VPD事件的真实影响。
#
#   本版按 direction_2v 分别处理:
#     Inhibit:   不变 —— VPD高时段找 sif_anom 最低点=最严重, 之后找最高点=恢复
#     Promote:   翻转 —— VPD高时段找 sif_anom 最高点=峰值(助推), 之后找最低点
#                =消退(回落到常态), "严重度"用同一套|Ye-Yn|公式, 只是Ye换成峰值
#     Ambiguous: 剔除 —— 没有一致方向, "最严重/最佳恢复"这套框架概念上不适用
#   Omega/Delta 公式本身不用改(已经取绝对值), 只改 severe_i/recov_i 用
#   which.min 还是 which.max 去找。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

dirn <- fread("data_proc/smap_bivar_134/classify_134.csv")[, .(meteo_stat, direction_2v)]
cat("=== 117站方向分布 ===\n"); print(dirn[, .N, by = direction_2v])
sel <- dirn[direction_2v != "Ambiguous"]$meteo_stat
cat("剔除Ambiguous后用于本次分析的站数:", length(sel), "\n")

fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)[meteo_stat %in% sel, .(meteo_stat, year, doy, date, SIF_buf1000)]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[meteo_stat %in% sel]
d   <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
d   <- merge(d, dirn, by = "meteo_stat")
setorder(d, meteo_stat, year, doy)

clim_by_doy <- function(x, doy) tapply(x, doy, mean, na.rm = TRUE)
deseason_z <- function(x, doy) {
  clim <- clim_by_doy(x, doy)
  a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (a - mean(a, na.rm = TRUE)) / s
}
d[, sif_anom := deseason_z(SIF_buf1000, doy), by = meteo_stat]
d[, vpd_z    := as.numeric(scale(vpd_mean)), by = meteo_stat]
d[, sif_clim_doy := clim_by_doy(SIF_buf1000, doy)[as.character(doy)], by = meteo_stat]

# ---- 逐站逐年: 按方向翻转 severe_i/recov_i 的取法 ---------------------------
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
    severe_i <- seg1[which.min(yy$sif_anom[seg1])]   # VPD高->SIF压到最低=最严重
    recov_i  <- seg2[which.max(yy$sif_anom[seg2])]   # 之后回升到最高=恢复
  } else {  # Promote: 翻转
    severe_i <- seg1[which.max(yy$sif_anom[seg1])]   # VPD高->SIF被推到最高=峰值
    recov_i  <- seg2[which.min(yy$sif_anom[seg2])]   # 之后回落到最低=消退
  }

  constraint_time <- as.numeric(severe_i - onset)
  recovery_time   <- as.numeric(recov_i - end)
  recov_censored  <- as.logical(end >= n)
  if (recov_censored) recovery_time <- NA_real_

  Yn_severe <- yy$sif_clim_doy[severe_i]; Ye <- yy$SIF_buf1000[severe_i]
  Yn_recov  <- yy$sif_clim_doy[recov_i];  Ye1 <- yy$SIF_buf1000[recov_i]
  dev_e  <- Ye  - Yn_severe
  dev_e1 <- Ye1 - Yn_recov
  Omega <- as.numeric(Yn_severe / abs(dev_e))
  Delta <- as.numeric(abs(dev_e) / abs(dev_e1))

  list(onset = as.integer(onset), end = as.integer(end),
       severe_i = as.integer(severe_i), recov_i = as.integer(recov_i),
       constraint_time = constraint_time, recovery_time = recovery_time,
       recov_censored = recov_censored,
       Omega = Omega, Delta = Delta, n_windows_season = as.integer(n))
}

ev4 <- d[, calc_event(.SD, direction_2v[1]), by = .(meteo_stat, year, direction_2v)]
ev4 <- ev4[!is.null(onset)]
cat("\n站-年事件数:", nrow(ev4), " | 站数:", uniqueN(ev4$meteo_stat), "\n")
print(ev4[, .N, by = direction_2v])
cat("recovery_time 删失比例:", round(100*mean(ev4$recov_censored), 1), "%\n")

winsor <- function(x, p = 0.01) {
  q <- quantile(x, c(p, 1 - p), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}
ev4[, `:=`(Omega_w = winsor(Omega), Delta_w = winsor(Delta),
           logOmega = log(Omega), logDelta = log(Delta))]

cat("\n=== 分布, 按方向拆开看(1%缩尾) ===\n")
print(ev4[, .(constraint_time = round(mean(constraint_time, na.rm=TRUE),2),
              recovery_time   = round(mean(recovery_time, na.rm=TRUE),2),
              Omega_w = round(median(Omega_w, na.rm=TRUE),2),
              Delta_w = round(median(Delta_w, na.rm=TRUE),2)), by = direction_2v])

fwrite(ev4, file.path(OUT, "recovery_liu2021_metrics_dirfix_117_events.csv"))
cat("\n-> recovery_liu2021_metrics_dirfix_117_events.csv (", nrow(ev4), "行 )\n")

# ---- 和原版(未修正,方向硬套Inhibit逻辑)对比一下差多少 ----------------------
old <- fread(file.path(OUT, "recovery_liu2021_metrics_117_events.csv"))
cmp <- merge(ev4[, .(meteo_stat, year, direction_2v, constraint_time, recovery_time, Omega, Delta)],
             old[, .(meteo_stat, year, constraint_time_old = constraint_time,
                     recovery_time_old = recovery_time, Omega_old = Omega, Delta_old = Delta)],
             by = c("meteo_stat", "year"))
cat("\n=== 修正前后差异, 只看Promote站(受影响的那一半) ===\n")
print(cmp[direction_2v == "Promote", .(
  r_constraint = round(cor(constraint_time, constraint_time_old, use="complete.obs"), 2),
  r_recovery   = round(cor(recovery_time, recovery_time_old, use="complete.obs"), 2),
  r_Omega      = round(cor(Omega, Omega_old, use="complete.obs"), 2),
  r_Delta      = round(cor(Delta, Delta_old, use="complete.obs"), 2))])

p <- ggplot(melt(ev4[, .(meteo_stat, direction_2v, constraint_time, recovery_time, Omega_w, Delta_w)],
                  id.vars = c("meteo_stat","direction_2v")), aes(value, fill = direction_2v)) +
  geom_histogram(bins = 30, alpha = .7, position = "identity") +
  facet_wrap(~variable, scales = "free", ncol = 2) +
  scale_fill_manual(values = c(Promote="#2166AC", Inhibit="#C1121F")) +
  labs(title = "Direction-corrected recovery metrics (117 stations minus 60 Ambiguous, station-years)",
       x = NULL, y = "count", fill = NULL) + theme_bw(base_size = 11) + theme(legend.position="top")
ggsave(file.path(OUT, "recovery_liu2021_metrics_dirfix_117_hist.png"), p, width = 9, height = 7, dpi = 300)
cat("-> recovery_liu2021_metrics_dirfix_117_hist.png\n")
