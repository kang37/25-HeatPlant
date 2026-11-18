#!/usr/bin/env Rscript
# =============================================================================
# recovery_liu2021_metrics_117.R
#   在117站(CCM确认因果耦合)的事件表上补算 Liu et al. 2021 (Agric Forest
#   Meteorol) 的四个量: constraint_time, recovery_time(定义比现有Trec_windows
#   更精细), Omega(抵抗力,Isbell 2015式), Delta(恢复力,Isbell 2015式)。
#
#   和原文的适配差异(必须交代):
#   - 原文用土壤湿度(SM)阈值(WHC的50%/75%)定义干旱起止; 本项目没有实测SM,
#     改用VPD距平(vpd_z)是否>0来定义"干旱期"(含峰值窗口的连续正距平区间)。
#   - 原文"最严重干旱日"/"最佳恢复日"用的是植被指数本身的极值(不是气象变量的
#     极值)——本脚本严格照此: 在干旱期起点到季末之间找sif_anom最小值点当
#     "最严重日", 在干旱期终点到季末之间找sif_anom最大值点当"最佳恢复日"。
#   - 原文Yn/Ye/Ye+1是年尺度前后对比(2009年一次全球大干旱); 本项目是逐站逐年
#     都有一次事件, 改用该窗口所在DOY的多年气候态均值当Yn(而不是"下一年"的
#     笼统年际对比), 更贴合8天合成数据的季节结构。
#   - 原文排除温度距平>4K的事件; 本脚本暂不做这层排除(需要独立的气温8天序列,
#     不在本次范围内), 留作后续补充。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. 117站清单 + HCSIF/VPD 距平构造(同前几个脚本) -----------------------
sel <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")$meteo_stat
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)[meteo_stat %in% sel, .(meteo_stat, year, doy, date, SIF_buf1000)]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[meteo_stat %in% sel]
d   <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
setorder(d, meteo_stat, year, doy)

clim_by_doy <- function(x, doy) tapply(x, doy, mean, na.rm = TRUE)  # Yn 的来源: 逐DOY多年气候态
deseason_z <- function(x, doy) {
  clim <- clim_by_doy(x, doy)
  a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (a - mean(a, na.rm = TRUE)) / s
}
d[, sif_anom := deseason_z(SIF_buf1000, doy), by = meteo_stat]
d[, vpd_z    := as.numeric(scale(vpd_mean)), by = meteo_stat]
d[, sif_clim_doy := clim_by_doy(SIF_buf1000, doy)[as.character(doy)], by = meteo_stat]  # Yn(该DOY气候态)

# ---- 2. 逐站逐年: 干旱期(VPD正距平连续区间) -> 最严重日/最佳恢复日 ---------
calc_event <- function(sub) {
  yy <- sub[order(doy)]; n <- nrow(yy)
  i0 <- which.max(yy$vpd_z)                      # 峰值窗口(和现有Trec构造一致的锚点)
  # 干旱期: 含i0的、vpd_z>0的连续区间
  above <- yy$vpd_z > 0
  if (is.na(above[i0]) || !above[i0]) return(NULL)
  lo <- i0; while (lo > 1 && isTRUE(above[lo - 1])) lo <- lo - 1
  hi <- i0; while (hi < n && isTRUE(above[hi + 1])) hi <- hi + 1
  onset <- lo; end <- hi

  # 最严重日: onset..n 里 sif_anom 最小值点; 最佳恢复日: end..n 里最大值点
  seg1 <- onset:n
  if (all(is.na(yy$sif_anom[seg1]))) return(NULL)
  severe_i <- seg1[which.min(yy$sif_anom[seg1])]
  seg2 <- end:n
  if (all(is.na(yy$sif_anom[seg2]))) return(NULL)
  recov_i <- seg2[which.max(yy$sif_anom[seg2])]

  constraint_time <- as.numeric(severe_i - onset)
  recovery_time   <- as.numeric(recov_i - end)
  recov_censored  <- as.logical(end >= n)        # 干旱到季末都没"结束", 恢复时长不可判定
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

ev4 <- d[, calc_event(.SD), by = .(meteo_stat, year)]
ev4 <- ev4[!is.null(onset)]
cat("站-年事件数:", nrow(ev4), " | 站数:", uniqueN(ev4$meteo_stat), "\n")
cat("recovery_time 删失(干旱到季末未结束)比例:", round(100*mean(ev4$recov_censored), 1), "%\n")

# ---- 3. 极端值处理: Omega/Delta 在分母接近0时会炸到很大, 1%缩尾 ------------
winsor <- function(x, p = 0.01) {
  q <- quantile(x, c(p, 1 - p), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}
ev4[, `:=`(Omega_w = winsor(Omega), Delta_w = winsor(Delta))]

cat("\n=== 四个量的分布(1%缩尾后) ===\n")
print(ev4[, .(constraint_time = round(mean(constraint_time, na.rm=TRUE),2),
              recovery_time   = round(mean(recovery_time, na.rm=TRUE),2),
              Omega_w = round(median(Omega_w, na.rm=TRUE),2),
              Delta_w = round(median(Delta_w, na.rm=TRUE),2))])

# ---- 4. 合并方向标签(参考,不是主分析), 存表 ---------------------------------
dirn <- fread("data_proc/smap_bivar_134/classify_134.csv")[, .(meteo_stat, direction_2v)]
ev4  <- merge(ev4, dirn, by = "meteo_stat", all.x = TRUE)
fwrite(ev4, file.path(OUT, "recovery_liu2021_metrics_117_events.csv"))
cat("\n-> recovery_liu2021_metrics_117_events.csv (", nrow(ev4), "行 )\n")

# ---- 5. 四个量两两相关(含现有Trec_windows), 快速看看是不是同一个信息 -------
old <- fread("data_proc/output_loose_classify/recovery_vs_134_direction_events.csv")[
  , .(meteo_stat, year, Trec_windows)]
cmp <- merge(ev4, old, by = c("meteo_stat", "year"))
cm <- cmp[, .(constraint_time, recovery_time, Omega_w, Delta_w, Trec_windows)]
cat("\n=== 四个新量 x 旧Trec_windows 相关矩阵(pairwise complete) ===\n")
print(round(cor(cm, use = "pairwise.complete.obs"), 2))

p <- ggplot(melt(cmp[, .(meteo_stat, constraint_time, recovery_time, Omega_w, Delta_w)],
                  id.vars = "meteo_stat"), aes(value)) +
  geom_histogram(bins = 30, fill = "#2166AC") +
  facet_wrap(~variable, scales = "free", ncol = 2) +
  labs(title = "Distributions of the four Liu et al. 2021-style recovery metrics (117 stations, station-years)",
       x = NULL, y = "count") + theme_bw(base_size = 11)
ggsave(file.path(OUT, "recovery_liu2021_metrics_117_hist.png"), p, width = 9, height = 7, dpi = 300)
cat("-> recovery_liu2021_metrics_117_hist.png\n")
