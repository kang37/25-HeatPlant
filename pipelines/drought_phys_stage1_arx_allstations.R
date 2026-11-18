#!/usr/bin/env Rscript
# =============================================================================
# drought_phys_stage1_arx_allstations.R
#   两阶段法·第一阶段: 抛开CCM筛选, 对全部协变量齐全的站点(904)逐站拟合
#   De Keersmaecker et al. 2015 式 ARX:
#     模型A(主): SIF_z(t) = phi*SIF_z(t-1) + beta*VPD_z(t) + g1*Tavg_z(t) + g2*Rsds_z(t) + g3*Prcp_z(t) + e
#     模型B(对照,与此前一致): SIF_z(t) = phi*SIF_z(t-1) + beta*VPD_z(t) + e
#   全部变量=逐站去季节(按DOY气候态)后再标准化的距平, 年内滞后(不跨冬歇期)。
#   beta(模型A) = 控制气温/辐射/降水后, 1个SD的VPD距平对SIF距平(SD单位)的即期效应,
#   是"干旱-植物生理关系"的逐站度量; phi=异常持续性。窗口=HCSIF 8天合成(DOY 124..268)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]; cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
stn <- env$dt$stat_id; cat("协变量齐全站点:", length(stn), "\n")

fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread, select = c("meteo_stat", "year", "doy", "SIF_buf1000")), fill = TRUE)[meteo_stat %in% stn]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[meteo_stat %in% stn, .(meteo_stat, year, doy, vpd = vpd_mean)]
rs  <- rbindlist(lapply(2000:2022, function(y) fread(sprintf("data_raw/covariates_1km/era5_rsds_sif8d_%d.csv", y))))
rs  <- rs[stat_id %in% stn, .(meteo_stat = stat_id, year, doy, rsds = rsds_mean)]

read_meteo_8d <- function(sid) {
  fp <- file.path("data_raw/meteo_data_1961-2023", paste0(sid, ".txt"))
  if (!file.exists(fp)) return(NULL)
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
met <- rbindlist(lapply(stn, read_meteo_8d))
cat("气象8天窗口: ", nrow(met), "行, 站数", uniqueN(met$meteo_stat), "\n")

d <- Reduce(function(a, b) merge(a, b, by = c("meteo_stat", "year", "doy")), list(sif, vpd, rs, met))
setorder(d, meteo_stat, year, doy)
deseason_z <- function(x, doy) {
  clim <- tapply(x, doy, mean, na.rm = TRUE); a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (a - mean(a, na.rm = TRUE)) / s
}
d[, `:=`(sif_z = deseason_z(SIF_buf1000, doy), vpd_z = deseason_z(vpd, doy), tavg_z = deseason_z(tavg, doy),
         rsds_z = deseason_z(rsds, doy), prcp_z = deseason_z(log1p(prcp), doy)), by = meteo_stat]
d[, lag1 := shift(sif_z, 1), by = .(meteo_stat, year)]

fit_station <- function(sub) {
  sub <- sub[complete.cases(sif_z, lag1, vpd_z, tavg_z, rsds_z, prcp_z)]
  n <- nrow(sub); if (n < 100) return(NULL)
  mA <- lm(sif_z ~ lag1 + vpd_z + tavg_z + rsds_z + prcp_z, data = sub)
  mB <- lm(sif_z ~ lag1 + vpd_z, data = sub)
  cA <- summary(mA)$coefficients; cB <- summary(mB)$coefficients
  vif_vpd <- 1 / (1 - summary(lm(vpd_z ~ lag1 + tavg_z + rsds_z + prcp_z, data = sub))$r.squared)
  list(n = n, phi = cA["lag1", 1], phi_se = cA["lag1", 2],
       beta = cA["vpd_z", 1], beta_se = cA["vpd_z", 2], beta_p = cA["vpd_z", 4],
       g_tavg = cA["tavg_z", 1], g_rsds = cA["rsds_z", 1], g_prcp = cA["prcp_z", 1], r2 = summary(mA)$r.squared,
       beta_vpdonly = cB["vpd_z", 1], beta_vpdonly_se = cB["vpd_z", 2], phi_vpdonly = cB["lag1", 1], vif_vpd = vif_vpd)
}
res <- d[, fit_station(.SD), by = meteo_stat]
co  <- fread("data_raw/WS土地利用比例_1000m.csv", select = 1:3); setnames(co, c("meteo_stat", "lon", "lat"))
res <- merge(res, co, by = "meteo_stat", all.x = TRUE)
fwrite(res, file.path(OUT, "drought_phys_arx_station_beta.csv"))

cat("\n成功拟合站数:", nrow(res), " | 每站平均有效点:", round(mean(res$n)), "\n")
cat("\n=== beta(控制气温/辐射/降水后的VPD->SIF即期效应, SD单位) ===\n"); print(summary(res$beta))
cat(sprintf("beta<0: %.1f%% | beta>0: %.1f%% | p<0.05显著: %.1f%% (其中负 %.1f%% / 正 %.1f%%)\n",
            100 * mean(res$beta < 0), 100 * mean(res$beta > 0), 100 * mean(res$beta_p < 0.05),
            100 * mean(res$beta_p < 0.05 & res$beta < 0), 100 * mean(res$beta_p < 0.05 & res$beta > 0)))
cat("beta标准误中位数:", round(median(res$beta_se), 3), " | beta的SD:", round(sd(res$beta), 3),
    " | 噪声占比(mean se^2/var):", round(mean(res$beta_se^2) / var(res$beta), 2), "\n")
cat("VIF(vpd_z)中位数/最大:", round(median(res$vif_vpd), 2), "/", round(max(res$vif_vpd), 2), "\n")
cat("beta(模型A) vs beta(仅VPD) Pearson r:", round(cor(res$beta, res$beta_vpdonly), 3),
    " | beta_vpdonly<0占比:", round(100 * mean(res$beta_vpdonly < 0), 1), "%\n")
cat("phi中位数:", round(median(res$phi), 3), " | phi vs phi_vpdonly r:", round(cor(res$phi, res$phi_vpdonly), 3), "\n")

p1 <- ggplot(res, aes(beta)) + geom_histogram(bins = 50, fill = "#2166AC") + geom_vline(xintercept = 0, color = "#C1121F") +
  labs(title = sprintf("Station-level VPD -> SIF sensitivity (beta), %d stations", nrow(res)),
       subtitle = "ARX with Tavg/Rsds/Prcp anomalies as controls; SD units; red line = 0", x = "beta", y = "stations") + theme_bw(base_size = 11)
ggsave(file.path(OUT, "drought_phys_beta_hist.png"), p1, width = 7, height = 4.5, dpi = 300)
p2 <- ggplot(res[!is.na(lon)], aes(lon, lat, color = pmax(pmin(beta, 0.4), -0.4))) + geom_point(size = 1.6) +
  scale_color_gradient2(low = "#2166AC", mid = "grey90", high = "#C1121F", midpoint = 0, name = "beta\n(clipped)") + coord_fixed(1.15) +
  labs(title = "Spatial pattern of station-level beta (VPD -> SIF)", x = "Longitude", y = "Latitude") + theme_bw(base_size = 11)
ggsave(file.path(OUT, "drought_phys_beta_map.png"), p2, width = 8, height = 6, dpi = 300)
