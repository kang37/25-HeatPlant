#!/usr/bin/env Rscript
# =============================================================================
# recovery_liu2021_panel_mixed_117.R
#   组装"站-年"面板协变量表(8个动态变量按年匹配 + elev + 9个仍静态的地表覆盖),
#   和 recovery_liu2021_metrics_117 的四个恢复力指标合并, 用混合效应模型
#   (随机截距=站点, nlme::lme)重新估计协变量效应——不再把117站坍缩成一个数,
#   直接在站-年层面用全部观测, 站点聚类通过随机截距处理。
#   动态变量来源(2026-09-07确认, 见对话记录):
#     tavg/rh/cloud/precip: 复用21_meteo_agg.R的原始逐日气象文件, 改成按年聚合
#     rsds:                 复用22号脚本rs中间表(已经是逐年), 不再折叠成mean/sd
#     ntl:                  data_raw/mvnl_monthly_buffer1000_2000_2024.csv 逐年
#     invest:               city_invest_data.csv 逐年列("园林绿化_2003..2022")
#     urban_rate:           中国城市数据库1990-2023.xlsx 逐年
#   静态: elev(dem) + imperv/grass/crop/water/五个叶型森林占比(glc_station_2020,
#         仍是单一年份快照, 未重做逐年提取)
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(terra); library(readxl); library(nlme)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")
YRS  <- 2000:2022
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
strip_city <- function(x) sub("(市辖区|地区|自治州|自治县|特别行政区|城区|盟|市|县|区)$", "", x)

sel <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")$meteo_stat
st  <- fread("data_raw/hcsif/stations_924.csv"); setnames(st, "meteo_stat", "stat_id")
st  <- st[stat_id %in% sel]
log_msg("目标站数: ", nrow(st))

# ---- 动态1: tavg/rh/cloud/precip 逐年(改自21_meteo_agg.R, 只改分组维度) -----
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
               rh = mean(RH, na.rm = TRUE), cloud = mean(CF, na.rm = TRUE),
               precip = mean(precip, na.rm = TRUE)), by = yr]
  out[, stat_id := sid]; setnames(out, "yr", "year")
  out
}
meteo_yr <- rbindlist(lapply(st$stat_id, read_one_year), fill = TRUE)
meteo_yr <- meteo_yr[tavg_n >= 0.5][, tavg_n := NULL]
log_msg("meteo逐年: ", nrow(meteo_yr), " 站-年")

# ---- 动态2: rsds 逐年(复用22号脚本的年度中间表, 不折叠) --------------------
rs <- rbindlist(lapply(YRS, function(y)
  fread(sprintf("data_raw/covariates_1km/era5_rsds_sif8d_%d.csv", y))))
rsds_yr <- rs[, .(rsds = mean(rsds_mean)), by = .(stat_id, year)]
log_msg("rsds逐年: ", nrow(rsds_yr), " 站-年")

# ---- 动态3: ntl 逐年 --------------------------------------------------------
nt <- fread("data_raw/mvnl_monthly_buffer1000_2000_2024.csv", encoding = "UTF-8")
setnames(nt, names(nt), sub("^﻿", "", names(nt)))
ntl_yr <- nt[year %in% YRS & quality_flag == "OK",
             .(ntl = log1p(mean(area_weighted_mean))), by = .(stat_id, year)]
log_msg("ntl逐年: ", nrow(ntl_yr), " 站-年, 覆盖年份: ", paste(range(ntl_yr$year), collapse="-"))

# ---- 站点 -> 城市映射(复用22号脚本) ----------------------------------------
cv <- vect("data_raw/china_cities/city.shp")
pt <- vect(st, geom = c("longitude", "latitude"), crs = crs(cv))
hit <- terra::relate(pt, cv, "intersects")
idx <- apply(hit, 1, function(r) { w <- which(r); if (length(w)) w[1] else NA_integer_ })
if (anyNA(idx)) idx[is.na(idx)] <- terra::nearest(pt[is.na(idx)], cv)$to_id
smap <- unique(data.table(stat_id = st$stat_id, city = strip_city(cv$ct_name[idx])), by = "stat_id")

# ---- 动态4: invest 逐年 -----------------------------------------------------
iv <- fread("data_raw/green_invest/city_invest_data.csv", encoding = "UTF-8")
setnames(iv, 1, "city_raw"); iv[, city := strip_city(city_raw)]
ycol <- grep("^园林绿化_20(0[3-9]|1[0-9]|2[0-2])$", names(iv), value = TRUE)
iv_long <- melt(iv, id.vars = "city", measure.vars = ycol, variable.name = "yr_col", value.name = "invest_raw")
iv_long[, year := as.integer(sub("^园林绿化_", "", yr_col))]
iv_long <- iv_long[is.finite(invest_raw), .(invest = log1p(mean(invest_raw))), by = .(city, year)]
invest_yr <- merge(smap, iv_long, by = "city", allow.cartesian = TRUE)[, .(stat_id, year, invest)]
log_msg("invest逐年: ", nrow(invest_yr), " 站-年, 覆盖年份: ", paste(range(invest_yr$year), collapse="-"))

# ---- 动态5: urban_rate 逐年 -------------------------------------------------
cd <- as.data.table(read_excel("data_raw/中国城市数据库1990-2023.xlsx", sheet = "原始数据", guess_max = 1e6))
num <- function(x) suppressWarnings(as.numeric(x))
cd <- cd[as.integer(num(年份)) %in% YRS,
         .(city = strip_city(城市), year = as.integer(num(年份)),
           nong = num(`非农业人口数(万人)`), hu = num(`户籍人口(万人)`))]
cd <- cd[is.finite(nong) & is.finite(hu) & hu > 0]
cd[, urban_rate := 100 * nong / hu]
urban_yr <- merge(smap, cd[, .(city, year, urban_rate)], by = "city", allow.cartesian = TRUE)[, .(stat_id, year, urban_rate)]
log_msg("urban_rate逐年: ", nrow(urban_yr), " 站-年, 覆盖年份: ", paste(range(urban_yr$year), collapse="-"))

# ---- 静态: elev + 9个地表覆盖(仍是2020单年快照) ----------------------------
dem <- fread("data_raw/covariates_1km/dem_cell_1km.csv")[, .(elev = mean(elev, na.rm = TRUE)), by = stat_id]
g <- fread("data_raw/covariates_1km/glc_station_2020.csv")
L <- function(...) rowSums(as.matrix(g[, sprintf("LC%02d", c(...)), with = FALSE]))
lc <- data.table(stat_id = g$stat_id, imperv = g$LC23, grass = L(17), crop = L(0:3), water = L(27),
                  ever_needle = L(8, 9), deci_needle = L(10, 11),
                  ever_broad = L(4, 5), deci_broad = L(6, 7), mixedleaf = L(12:13))
static <- merge(dem, lc, by = "stat_id")

# ---- 拼成"站-年"面板 --------------------------------------------------------
grid <- CJ(stat_id = st$stat_id, year = YRS)
panel <- Reduce(function(a, b) merge(a, b, by = c("stat_id", "year"), all.x = TRUE),
                 list(grid, meteo_yr, rsds_yr, ntl_yr, invest_yr, urban_yr))
panel <- merge(panel, static, by = "stat_id", all.x = TRUE)
DROP_URBAN_RATE <- TRUE  # 2026-09-07: urban_rate 数据源只到2008年, 严重压缩样本窗口, 先去掉换2003-2022覆盖
DYNVARS <- c("tavg", "rh", "cloud", "precip", "rsds", "ntl", "invest",
             if (!DROP_URBAN_RATE) "urban_rate")
STATVARS <- c("elev", "imperv", "grass", "crop", "water",
              "ever_needle", "deci_needle", "ever_broad", "deci_broad", "mixedleaf")
PRED <- c(DYNVARS, STATVARS)
fwrite(panel, file.path(OUT, "recovery_liu2021_panel_covariates.csv"))
log_msg("面板协变量表: ", nrow(panel), " 行(站x年), 覆盖率:")
print(panel[, lapply(.SD, function(x) round(100*mean(!is.na(x)),1)), .SDcols = PRED])

# ---- 合并恢复力指标(站-年), 做混合效应模型 ---------------------------------
ev <- fread(file.path(OUT, "recovery_liu2021_metrics_117_events.csv"))
ev[, `:=`(logOmega = log(Omega), logDelta = log(Delta))]
d <- merge(ev, panel, by.x = c("meteo_stat", "year"), by.y = c("stat_id", "year"))
d <- d[complete.cases(d[, PRED, with = FALSE])]
log_msg("\n可建模的站-年样本: ", nrow(d), " | 站数: ", uniqueN(d$meteo_stat))

dz <- copy(d)
for (v in PRED) dz[[v]] <- as.numeric(scale(dz[[v]]))

DVS <- c("constraint_time", "recovery_time", "logOmega", "logDelta")
run_lme <- function(dv) {
  dd <- dz[!is.na(get(dv))]
  f <- as.formula(paste(dv, "~", paste(PRED, collapse = "+")))
  m <- tryCatch(lme(f, random = ~1 | meteo_stat, data = dd, na.action = na.omit,
                     control = lmeControl(opt = "optim", maxIter = 200)),
                error = function(e) { message(dv, ": ", e$message); NULL })
  if (is.null(m)) return(NULL)
  co <- as.data.table(summary(m)$tTable, keep.rownames = "var")
  setnames(co, c("var", "est", "se", "df", "t", "p"))
  co <- co[var != "(Intercept)"][order(-abs(est))]
  co[, sig := ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(p < 0.1, ".", "")))]
  list(model = m, co = co, n = nrow(dd), n_stat = uniqueN(dd$meteo_stat))
}
res <- lapply(DVS, run_lme); names(res) <- DVS

for (dv in DVS) {
  if (is.null(res[[dv]])) { cat("\n=====", dv, ": 拟合失败 =====\n"); next }
  cat("\n=====", dv, " (混合效应, n=", res[[dv]]$n, "站-年,", res[[dv]]$n_stat, "站) =====\n")
  print(res[[dv]]$co[, .(var, std_coef = round(est, 3), p = round(p, 3), sig)])
}

coef_all <- rbindlist(lapply(DVS, function(dv)
  if (!is.null(res[[dv]])) res[[dv]]$co[, dv := dv]), fill = TRUE)
fwrite(coef_all, file.path(OUT, "recovery_liu2021_panel_mixed_coefs.csv"))
cat("\n-> recovery_liu2021_panel_mixed_coefs.csv\n")
