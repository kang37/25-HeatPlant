#!/usr/bin/env Rscript
# =============================================================================
# panel_covariates_178.R
#   把 recovery_liu2021_panel_mixed_117.R 里的"站-年动态协变量"构造原样扩到178站
#   (CCM判据(1): tp 0-8内有 Δρ>0 且 p_surr<0.05, 排除地表覆盖断点站), 输出
#   recovery_panel_covariates_178.csv。不含urban_rate(源数据只到2008年)。
#   动态: tavg/rh/cloud/precip(5-9月均值), rsds(ERA5), ntl(log1p), invest(log1p, 仅2003-2022)
#   静态: elev + imperv/grass/crop/water/各叶型森林占比(2020快照)
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(terra)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")
YRS  <- 2000:2022
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
strip_city <- function(x) sub("(市辖区|地区|自治州|自治县|特别行政区|城区|盟|市|县|区)$", "", x)

neg <- fread("data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv")
neg[, drho := rho - rho_min]
LCB <- "data_proc/output_hcsif_buf1000/landcover_break_stations.csv"
if (file.exists(LCB)) { lcb <- fread(LCB); neg <- neg[!meteo_stat %in% lcb[has_break == TRUE]$stat_id] }
sr <- neg[tp >= 0 & tp <= 8]; sr[, pass1 := (drho > 0) & (p_surr < 0.05)]
sel <- unique(sr[pass1 == TRUE]$meteo_stat)
sel117 <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")$meteo_stat
st  <- fread("data_raw/hcsif/stations_924.csv"); setnames(st, "meteo_stat", "stat_id")
st  <- st[stat_id %in% sel]
log_msg("178站清单: ", length(sel), " | 在stations_924里: ", nrow(st), " | 其中属于117站: ", sum(st$stat_id %in% sel117))

MISS <- 999999.0
RANGE <- list(tavg = c(-60, 60), RH = c(0, 100), CF = c(0, 100), precip = c(0, 2000))
read_one_year <- function(sid) {
  fp <- file.path("data_raw/meteo_data_1961-2023", paste0(sid, ".txt"))
  if (!file.exists(fp)) return(NULL)
  d <- fread(fp, skip = 1, showProgress = FALSE)
  if (!"date" %in% names(d)) return(NULL)
  d[, date := as.IDate(date)]
  d <- d[year(date) %in% YRS & month(date) %in% 5:9]
  if (!nrow(d)) return(NULL)
  for (v in names(RANGE)) if (v %in% names(d))
    d[!is.finite(get(v)) | get(v) < RANGE[[v]][1] | get(v) > RANGE[[v]][2], (v) := NA_real_]
  d[, yr := year(date)]
  out <- d[, .(tavg = mean(tavg, na.rm = TRUE), tavg_n = mean(!is.na(tavg)),
               rh = mean(RH, na.rm = TRUE), cloud = mean(CF, na.rm = TRUE), precip = mean(precip, na.rm = TRUE)), by = yr]
  out[, stat_id := sid]; setnames(out, "yr", "year"); out
}
meteo_yr <- rbindlist(lapply(st$stat_id, read_one_year), fill = TRUE)
meteo_yr <- meteo_yr[tavg_n >= 0.5][, tavg_n := NULL]
log_msg("meteo逐年: ", nrow(meteo_yr), " 站-年, 站数 ", uniqueN(meteo_yr$stat_id))

rs <- rbindlist(lapply(YRS, function(y) fread(sprintf("data_raw/covariates_1km/era5_rsds_sif8d_%d.csv", y))))
rsds_yr <- rs[stat_id %in% sel, .(rsds = mean(rsds_mean)), by = .(stat_id, year)]
log_msg("rsds逐年: ", nrow(rsds_yr), " 站-年, 站数 ", uniqueN(rsds_yr$stat_id))

nt <- fread("data_raw/mvnl_monthly_buffer1000_2000_2024.csv", encoding = "UTF-8")
setnames(nt, names(nt), sub("^﻿", "", names(nt)))
ntl_yr <- nt[year %in% YRS & quality_flag == "OK" & stat_id %in% sel, .(ntl = log1p(mean(area_weighted_mean))), by = .(stat_id, year)]
log_msg("ntl逐年: ", nrow(ntl_yr), " 站-年, 站数 ", uniqueN(ntl_yr$stat_id), ", 年份 ", paste(range(ntl_yr$year), collapse = "-"))

cv <- vect("data_raw/china_cities/city.shp")
pt <- vect(st, geom = c("longitude", "latitude"), crs = crs(cv))
hit <- terra::relate(pt, cv, "intersects")
idx <- apply(hit, 1, function(r) { w <- which(r); if (length(w)) w[1] else NA_integer_ })
if (anyNA(idx)) idx[is.na(idx)] <- terra::nearest(pt[is.na(idx)], cv)$to_id
smap <- unique(data.table(stat_id = st$stat_id, city = strip_city(cv$ct_name[idx])), by = "stat_id")

iv <- fread("data_raw/green_invest/city_invest_data.csv", encoding = "UTF-8")
setnames(iv, 1, "city_raw"); iv[, city := strip_city(city_raw)]
ycol <- grep("^园林绿化_20(0[3-9]|1[0-9]|2[0-2])$", names(iv), value = TRUE)
iv_long <- melt(iv, id.vars = "city", measure.vars = ycol, variable.name = "yr_col", value.name = "invest_raw")
iv_long[, year := as.integer(sub("^园林绿化_", "", yr_col))]
iv_long <- iv_long[is.finite(invest_raw), .(invest = log1p(mean(invest_raw))), by = .(city, year)]
invest_yr <- merge(smap, iv_long, by = "city", allow.cartesian = TRUE)[, .(stat_id, year, invest)]
log_msg("invest逐年: ", nrow(invest_yr), " 站-年, 站数 ", uniqueN(invest_yr$stat_id), ", 年份 ", paste(range(invest_yr$year), collapse = "-"))

dem <- fread("data_raw/covariates_1km/dem_cell_1km.csv")[, .(elev = mean(elev, na.rm = TRUE)), by = stat_id]
g <- fread("data_raw/covariates_1km/glc_station_2020.csv")
L <- function(...) rowSums(as.matrix(g[, sprintf("LC%02d", c(...)), with = FALSE]))
lc <- data.table(stat_id = g$stat_id, imperv = g$LC23, grass = L(17), crop = L(0:3), water = L(27),
                 ever_needle = L(8, 9), deci_needle = L(10, 11), ever_broad = L(4, 5), deci_broad = L(6, 7), mixedleaf = L(12:13))
static <- merge(dem, lc, by = "stat_id")

grid <- CJ(stat_id = st$stat_id, year = YRS)
panel <- Reduce(function(a, b) merge(a, b, by = c("stat_id", "year"), all.x = TRUE), list(grid, meteo_yr, rsds_yr, ntl_yr, invest_yr))
panel <- merge(panel, static, by = "stat_id", all.x = TRUE)
panel[, in117 := stat_id %in% sel117]
fwrite(panel, file.path(OUT, "recovery_panel_covariates_178.csv"))
DYNVARS <- c("tavg", "rh", "cloud", "precip", "rsds", "ntl", "invest")
log_msg("面板协变量表: ", nrow(panel), " 行(站x年), 非缺失覆盖率(%):")
print(panel[, lapply(.SD, function(x) round(100 * mean(!is.na(x)), 1)), .SDcols = c(DYNVARS, "elev", "imperv", "grass")])
cat("\ninvest按年非缺失占比:\n"); print(panel[, .(share = round(mean(!is.na(invest)), 2)), by = year][order(year)][, paste0(year, ":", share, collapse = " ")])
