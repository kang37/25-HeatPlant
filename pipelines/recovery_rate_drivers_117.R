#!/usr/bin/env Rscript
# =============================================================================
# recovery_rate_drivers_117.R
#   "速率型"恢复力指标(Yi & Jackson 2021, ERL综述 Table1的Recovery rate=
#   (V_post-V_dist)/dt, 及Resilience(时间指标)=影响幅度/恢复时间同一类):
#     fall_rate  = |跌幅(dev_e)|  / constraint_time   (下跌速率)
#     climb_rate = |爬升幅度(dev_e1)| / recovery_time (回升速率, 仅未删失事件)
#   分母加1(以"窗口"为单位)避免除0——constraint_time/recovery_time=0(最严重点/
#   最佳恢复点恰好落在窗口边界)时不应该是无穷大, 而是"很快"的一个有限大数。
#   117站(CCM确认因果) x 18个CORE协变量, RF(permutation importance)+岭回归
#   两条线(OLS在这批协变量上已证实因共线性不可信, 不再单独跑)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ranger); library(glmnet)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. 117站 + HCSIF/VPD 距平构造(同前) -----------------------------------
sel <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")$meteo_stat
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)[meteo_stat %in% sel, .(meteo_stat, year, doy, date, SIF_buf1000)]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[meteo_stat %in% sel]
d   <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
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

# ---- 2. 逐站逐年事件: 同 recovery_liu2021_metrics_117.R, 这次保留 dev_e/dev_e1 --
calc_event <- function(sub) {
  yy <- sub[order(doy)]; n <- nrow(yy)
  i0 <- which.max(yy$vpd_z)
  above <- yy$vpd_z > 0
  if (is.na(above[i0]) || !above[i0]) return(NULL)
  lo <- i0; while (lo > 1 && isTRUE(above[lo - 1])) lo <- lo - 1
  hi <- i0; while (hi < n && isTRUE(above[hi + 1])) hi <- hi + 1
  onset <- lo; end <- hi
  seg1 <- onset:n
  if (all(is.na(yy$sif_anom[seg1]))) return(NULL)
  severe_i <- seg1[which.min(yy$sif_anom[seg1])]
  seg2 <- end:n
  if (all(is.na(yy$sif_anom[seg2]))) return(NULL)
  recov_i <- seg2[which.max(yy$sif_anom[seg2])]
  constraint_time <- as.numeric(severe_i - onset)
  recovery_time   <- as.numeric(recov_i - end)
  recov_censored  <- as.logical(end >= n)
  if (recov_censored) recovery_time <- NA_real_
  Yn_severe <- yy$sif_clim_doy[severe_i]; Ye <- yy$SIF_buf1000[severe_i]
  Yn_recov  <- yy$sif_clim_doy[recov_i];  Ye1 <- yy$SIF_buf1000[recov_i]
  dev_e  <- as.numeric(Ye  - Yn_severe)
  dev_e1 <- as.numeric(Ye1 - Yn_recov)
  list(constraint_time = constraint_time, recovery_time = recovery_time,
       recov_censored = recov_censored, dev_e = dev_e, dev_e1 = dev_e1)
}
ev <- d[, calc_event(.SD), by = .(meteo_stat, year)]
ev <- ev[!is.na(dev_e)]

# ---- 3. 速率: 分母 +1(窗口)避免除0 -----------------------------------------
ev[, fall_rate  := abs(dev_e)  / (constraint_time + 1)]
ev[, climb_rate := abs(dev_e1) / (recovery_time + 1)]
cat("站-年事件数:", nrow(ev), "| climb_rate 可算(未删失):", ev[!is.na(climb_rate), .N], "\n")
cat("\nfall_rate 分布: "); print(summary(ev$fall_rate))
cat("climb_rate 分布: "); print(summary(ev$climb_rate))

# 严重右偏(比值型), 沿用log处理(全部严格正)
ev[, logFallRate  := log(fall_rate)]
ev[, logClimbRate := log(climb_rate)]

fwrite(ev, file.path(OUT, "recovery_rate_117_events.csv"))

# ---- 4. 聚合到站点层面(log量取均值) + 合并18个CORE协变量 --------------------
agg <- ev[, .(logFallRate = mean(logFallRate, na.rm = TRUE),
              logClimbRate = mean(logClimbRate, na.rm = TRUE)), by = meteo_stat]

env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt <- env$dt; PRED <- unlist(env$CORE, use.names = FALSE)

dm <- merge(agg, dt, by.x = "meteo_stat", by.y = "stat_id")
dm <- dm[complete.cases(dm[, PRED, with = FALSE])]
cat("\n站点数(聚合后):", nrow(agg), " | 协变量完整可建模:", nrow(dm), "\n")
fwrite(dm, file.path(OUT, "recovery_rate_drivers_117_data.csv"))

# ---- 5. RF + 岭回归, 两个DV各一遍 -------------------------------------------
run_rf <- function(dv) {
  f <- as.formula(paste(dv, "~", paste(PRED, collapse = "+")))
  set.seed(42)
  rf <- ranger(f, data = dm, num.trees = 2000, importance = "permutation", seed = 42)
  data.table(dv = dv, var = names(rf$variable.importance), importance = rf$variable.importance, rf_r2 = rf$r.squared)
}
run_ridge <- function(dv) {
  X <- scale(as.matrix(dm[, PRED, with = FALSE])); y <- dm[[dv]]
  set.seed(42)
  cv <- cv.glmnet(X, y, alpha = 0, nfolds = 10, standardize = FALSE)
  co <- as.numeric(coef(cv, s = "lambda.1se"))[-1]
  pred <- predict(cv, X, s = "lambda.1se")
  r2 <- 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
  data.table(dv = dv, var = PRED, ridge_coef = co, ridge_r2 = r2)
}

rf_res    <- rbindlist(lapply(c("logFallRate","logClimbRate"), run_rf))
ridge_res <- rbindlist(lapply(c("logFallRate","logClimbRate"), run_ridge))
fwrite(rf_res, file.path(OUT, "recovery_rate_drivers_117_rf.csv"))
fwrite(ridge_res, file.path(OUT, "recovery_rate_drivers_117_ridge.csv"))

cat("\n=== RF OOB R2 ===\n"); print(unique(rf_res[,.(dv,rf_r2=round(rf_r2,3))]))
cat("\n=== Ridge(1se) R2 ===\n"); print(unique(ridge_res[,.(dv,ridge_r2=round(ridge_r2,3))]))

for (dv_i in c("logFallRate","logClimbRate")) {
  cat("\n--", dv_i, "RF top6 --\n")
  print(head(rf_res[dv==dv_i][order(-importance)],6)[,.(var,importance=round(importance,4))])
  cat("--", dv_i, "岭回归系数(按|系数|排序) --\n")
  print(ridge_res[dv==dv_i][order(-abs(ridge_coef))][,.(var,ridge_coef=round(ridge_coef,3))])
}

# ---- 6. 图 ------------------------------------------------------------------
rf_res[, rank_within := rank(-importance), by = dv]
p <- ggplot(rf_res[rank_within<=18], aes(reorder(var, importance), importance)) +
  geom_col(fill = "#2166AC") + coord_flip() + facet_wrap(~dv, scales = "free") +
  labs(title = "RF importance: log(fall_rate) / log(climb_rate) ~ 18 covariates (117 stations)",
       subtitle = "fall_rate=|drop magnitude|/(constraint_time+1); climb_rate=|rebound magnitude|/(recovery_time+1), Yi & Jackson 2021 (ERL) style",
       x = NULL, y = "Permutation importance") +
  theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size=8, color="grey35"))
ggsave(file.path(OUT, "recovery_rate_drivers_117_importance.png"), p, width = 10, height = 6, dpi = 300)
cat("\n-> recovery_rate_drivers_117_importance.png\n")
