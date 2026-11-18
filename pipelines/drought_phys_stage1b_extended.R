#!/usr/bin/env Rscript
# =============================================================================
# drought_phys_stage1b_extended.R
#   在 drought_phys_stage1_arx_allstations.R 基础上扩展, 为"投资是否缓解干旱影响"的方案网格提供
#   多种"干旱影响"度量(统一为"越大=干旱影响越严重"):
#     S_partial : -beta(VPD距平, 控制气温/辐射/降水距平)            [ARX]
#     S_vpdonly : -beta(只含VPD)                                     [ARX, 与De Keersmaecker原式一致]
#     S_spei    :  beta(SPEI-1: 取 -SPEI 作干旱方向, 控制辐射距平)   [ARX; ERA5 SPEI-1, 0.25度月尺度]
#     phi       :  AR1持续性(越大恢复越慢)
#     R_loss    :  各年VPD距平最高窗口里SIF相对气候态的损失比例 (clim-SIF)/clim 的跨年中位数
#     Z_worse   :  -(同一窗口SIF距平的跨年中位数)
#   同时保存逐窗口距平面板(含SPEI), 用于一阶段(交互项)模型。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(terra)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]; cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
stn <- env$dt$stat_id

fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread, select = c("meteo_stat", "year", "doy", "SIF_buf1000")), fill = TRUE)[meteo_stat %in% stn]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[meteo_stat %in% stn, .(meteo_stat, year, doy, vpd = vpd_mean)]
rs  <- rbindlist(lapply(2000:2022, function(y) fread(sprintf("data_raw/covariates_1km/era5_rsds_sif8d_%d.csv", y))))
rs  <- rs[stat_id %in% stn, .(meteo_stat = stat_id, year, doy, rsds = rsds_mean)]
read_meteo_8d <- function(sid) {
  fp <- file.path("data_raw/meteo_data_1961-2023", paste0(sid, ".txt")); if (!file.exists(fp)) return(NULL)
  d <- fread(fp, skip = 1, showProgress = FALSE); if (!all(c("date", "tavg", "precip") %in% names(d))) return(NULL)
  d[, date := as.IDate(date)]; d <- d[year(date) %in% 2000:2022]
  d[!is.finite(tavg) | tavg < -60 | tavg > 60, tavg := NA_real_]; d[!is.finite(precip) | precip < 0 | precip > 2000, precip := NA_real_]
  d[, yd := yday(date)]; d <- d[yd >= 124 & yd <= 275]; d[, doy := 124L + 8L * ((yd - 124L) %/% 8L)]; d[, year := year(date)]
  d[, .(tavg = if (sum(!is.na(tavg)) >= 6) mean(tavg, na.rm = TRUE) else NA_real_,
        prcp = if (sum(!is.na(precip)) >= 6) sum(precip, na.rm = TRUE) * 8 / sum(!is.na(precip)) else NA_real_), by = .(year, doy)][, meteo_stat := sid][]
}
met <- rbindlist(lapply(stn, read_meteo_8d))

# ---- SPEI-1 (月) 按站点提取, 映射到8天窗口(以窗口中点所在月) -----------------
co <- fread("data_raw/WS土地利用比例_1000m.csv", select = 1:3); setnames(co, c("meteo_stat", "lon", "lat"))
pts <- vect(co[meteo_stat %in% stn], geom = c("lon", "lat"), crs = "EPSG:4326")
spei_list <- list()
for (y in 2000:2022) for (mo in 5:9) {
  f <- list.files("data_raw/era5_drought", pattern = sprintf("^SPEI1_.*_%d%02d\\.", y, mo), full.names = TRUE)
  if (!length(f)) next
  v <- terra::extract(rast(f[1]), pts)[, 2]
  spei_list[[length(spei_list) + 1]] <- data.table(meteo_stat = co[meteo_stat %in% stn]$meteo_stat, year = y, month = mo, spei = v)
}
spei <- rbindlist(spei_list); cat("SPEI提取:", nrow(spei), "站-月, 缺失", sum(is.na(spei$spei)), "\n")

d <- Reduce(function(a, b) merge(a, b, by = c("meteo_stat", "year", "doy")), list(sif, vpd, rs, met))
d[, month := month(as.Date(paste0(year, "-01-01")) + (doy + 3L) - 1L)]
d <- merge(d, spei, by = c("meteo_stat", "year", "month"), all.x = TRUE)
setorder(d, meteo_stat, year, doy)
deseason_z <- function(x, doy) { clim <- tapply(x, doy, mean, na.rm = TRUE); a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x))); (a - mean(a, na.rm = TRUE)) / s }
d[, `:=`(sif_z = deseason_z(SIF_buf1000, doy), vpd_z = deseason_z(vpd, doy), tavg_z = deseason_z(tavg, doy),
         rsds_z = deseason_z(rsds, doy), prcp_z = deseason_z(log1p(prcp), doy)), by = meteo_stat]
d[, dry_spei := -spei]
d[, sif_clim := mean(SIF_buf1000, na.rm = TRUE), by = .(meteo_stat, doy)]
d[, lag1 := shift(sif_z, 1), by = .(meteo_stat, year)]
saveRDS(d[, .(meteo_stat, year, doy, sif_z, lag1, vpd_z, tavg_z, rsds_z, prcp_z, dry_spei)], file.path(OUT, "drought_phys_panel_anoms.rds"))

fit_station <- function(sub) {
  s1 <- sub[complete.cases(sif_z, lag1, vpd_z, tavg_z, rsds_z, prcp_z)]; n <- nrow(s1); if (n < 100) return(NULL)
  cA <- summary(lm(sif_z ~ lag1 + vpd_z + tavg_z + rsds_z + prcp_z, data = s1))$coefficients
  cB <- summary(lm(sif_z ~ lag1 + vpd_z, data = s1))$coefficients
  s2 <- sub[complete.cases(sif_z, lag1, dry_spei, rsds_z)]; n2 <- nrow(s2)
  cS <- if (n2 >= 100) summary(lm(sif_z ~ lag1 + dry_spei + rsds_z, data = s2))$coefficients else NULL
  ev <- sub[!is.na(vpd_z) & !is.na(sif_z) & sif_clim > 0][, .SD[which.max(vpd_z)], by = year]
  list(n = n, phi = cA["lag1", 1], phi_se = cA["lag1", 2], beta = cA["vpd_z", 1], beta_se = cA["vpd_z", 2],
       beta_vpdonly = cB["vpd_z", 1], beta_vpdonly_se = cB["vpd_z", 2],
       beta_spei = if (!is.null(cS)) cS["dry_spei", 1] else NA_real_, beta_spei_se = if (!is.null(cS)) cS["dry_spei", 2] else NA_real_,
       R_loss = median((ev$sif_clim - ev$SIF_buf1000) / ev$sif_clim, na.rm = TRUE), Z_worse = -median(ev$sif_z, na.rm = TRUE), n_ev = nrow(ev))
}
res <- d[, fit_station(.SD), by = meteo_stat]
res <- merge(res, co, by = "meteo_stat", all.x = TRUE)
fwrite(res, file.path(OUT, "drought_phys_station_metrics_ext.csv"))
cat("成功拟合站数:", nrow(res), "\n"); print(res[, lapply(.SD, function(x) round(median(x, na.rm = TRUE), 3)),
  .SDcols = c("phi", "beta", "beta_vpdonly", "beta_spei", "R_loss", "Z_worse")])
cat("beta_spei>0(干旱降低SIF)占比:", round(100 * mean(res$beta_spei > 0, na.rm = TRUE), 1), "%\n")
cat("各度量两两Spearman:\n"); m <- res[, .(S_partial = -beta, S_vpdonly = -beta_vpdonly, S_spei = beta_spei, phi, R_loss, Z_worse)]
print(round(cor(m, method = "spearman", use = "pairwise.complete.obs"), 2))
