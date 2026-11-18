#!/usr/bin/env Rscript
# =============================================================================
# panel_covariates_all.R
#   不做CCM筛选: 对dt18里全部协变量齐全的站(904)构造逐年动态协变量
#   (构造与 panel_covariates_178.R 完全相同, 只是站点清单不同)。
#   输出 recovery_panel_covariates_all.csv: stat_id, year, tavg, rh, precip, rsds, ntl, invest
#   (静态协变量在建模脚本里直接取dt18)
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(terra)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")
YRS  <- 2000:2022
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
strip_city <- function(x) sub("(市辖区|地区|自治州|自治县|特别行政区|城区|盟|市|县|区)$", "", x)

env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]; cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
sel <- env$dt$stat_id
st  <- fread("data_raw/hcsif/stations_924.csv"); setnames(st, "meteo_stat", "stat_id")
st  <- st[stat_id %in% sel]
log_msg("dt18站数: ", length(sel), " | 在stations_924里: ", nrow(st))

RANGE <- list(tavg = c(-60, 60), RH = c(0, 100), precip = c(0, 2000))
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
  out <- d[, .(tavg = mean(tavg, na.rm = TRUE), tavg_n = mean(!is.na(tavg)), rh = mean(RH, na.rm = TRUE), precip = mean(precip, na.rm = TRUE)), by = yr]
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
log_msg("ntl逐年: ", nrow(ntl_yr), " 站-年, 站数 ", uniqueN(ntl_yr$stat_id))

cv <- vect("data_raw/china_cities/city.shp")
pt <- vect(st, geom = c("longitude", "latitude"), crs = crs(cv))
hit <- terra::relate(pt, cv, "intersects")
idx <- apply(hit, 1, function(r) { w <- which(r); if (length(w)) w[1] else NA_integer_ })
if (anyNA(idx)) idx[is.na(idx)] <- terra::nearest(pt[is.na(idx)], cv)$to_id
smap <- unique(data.table(stat_id = st$stat_id, city = strip_city(cv$ct_name[idx])), by = "stat_id")
iv <- fread("data_raw/green_invest/city_invest_data.csv", encoding = "UTF-8")
setnames(iv, 1, "city_raw"); iv[, city := strip_city(city_raw)]
ycol <- grep("^园林绿化_20(0[3-9]|1[0-9]|2[0-2])$", names(iv), value = TRUE)
iv_long <- suppressWarnings(melt(iv, id.vars = "city", measure.vars = ycol, variable.name = "yr_col", value.name = "invest_raw"))
iv_long[, year := as.integer(sub("^园林绿化_", "", yr_col))]
iv_long <- iv_long[is.finite(invest_raw), .(invest = log1p(mean(invest_raw))), by = .(city, year)]
invest_yr <- merge(smap, iv_long, by = "city", allow.cartesian = TRUE)[, .(stat_id, year, invest)]
log_msg("invest逐年: ", nrow(invest_yr), " 站-年, 站数 ", uniqueN(invest_yr$stat_id), ", 城市数 ", uniqueN(smap$city))

grid <- CJ(stat_id = st$stat_id, year = YRS)
panel <- Reduce(function(a, b) merge(a, b, by = c("stat_id", "year"), all.x = TRUE), list(grid, meteo_yr, rsds_yr, ntl_yr, invest_yr))
fwrite(panel, file.path(OUT, "recovery_panel_covariates_all.csv"))
log_msg("面板协变量表: ", nrow(panel), " 行; 非缺失覆盖率(%):")
print(panel[, lapply(.SD, function(x) round(100 * mean(!is.na(x)), 1)), .SDcols = c("tavg", "rh", "precip", "rsds", "ntl", "invest")])
