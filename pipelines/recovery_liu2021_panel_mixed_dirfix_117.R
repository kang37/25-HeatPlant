#!/usr/bin/env Rscript
# =============================================================================
# recovery_liu2021_panel_mixed_dirfix_117.R
#   方向修正版事件表(57站, Promote 41/Inhibit 16, 已剔除Ambiguous) + 站-年
#   面板协变量(复用 recovery_liu2021_panel_mixed_117.R 已经组装好的
#   recovery_liu2021_panel_covariates.csv, 协变量构造和方向无关不用重算),
#   去掉 urban_rate(同上次决定, 换2003-2022更长覆盖), 混合效应模型(随机截距=站点)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(nlme)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

panel <- fread(file.path(OUT, "recovery_liu2021_panel_covariates.csv"))
ev    <- fread(file.path(OUT, "recovery_liu2021_metrics_dirfix_117_events.csv"))
ev[, `:=`(logOmega = log(Omega), logDelta = log(Delta))]

PRED <- c("tavg", "rh", "cloud", "precip", "rsds", "ntl", "invest", "elev",
          "imperv", "grass", "crop", "water",
          "ever_needle", "deci_needle", "ever_broad", "deci_broad", "mixedleaf")
d <- merge(ev, panel, by.x = c("meteo_stat", "year"), by.y = c("stat_id", "year"))
d <- d[complete.cases(d[, PRED, with = FALSE])]
cat("可建模的站-年样本:", nrow(d), " | 站数:", uniqueN(d$meteo_stat), "\n")
print(d[, .N, by = direction_2v])
cat("年份覆盖:", paste(range(d$year), collapse = "-"), "\n")

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
  list(co = co, n = nrow(dd), n_stat = uniqueN(dd$meteo_stat))
}
res <- lapply(DVS, run_lme); names(res) <- DVS

for (dv in DVS) {
  if (is.null(res[[dv]])) { cat("\n=====", dv, ": 拟合失败 =====\n"); next }
  cat("\n=====", dv, " (n=", res[[dv]]$n, "站-年,", res[[dv]]$n_stat, "站) =====\n")
  print(res[[dv]]$co[, .(var, std_coef = round(est, 3), p = round(p, 3), sig)])
}

coef_all <- rbindlist(lapply(DVS, function(dv)
  if (!is.null(res[[dv]])) res[[dv]]$co[, dv := dv]), fill = TRUE)
fwrite(coef_all, file.path(OUT, "recovery_liu2021_panel_mixed_dirfix_coefs.csv"))
cat("\n-> recovery_liu2021_panel_mixed_dirfix_coefs.csv\n")
