#!/usr/bin/env Rscript
# =============================================================================
# recovery_phi_ED_placebo.R
#   recovery_phi_climate_adjusted_ED.R 里"遮掉SIF扰动期(ED_sif)"让15变量driver模型的R²明显升高。
#   这里检验: 升高是因为"遮掉的是扰动事件", 还是仅仅"遮掉了一部分窗口"?
#   安慰剂: 对每个站-年, 把真实遮罩在年内做随机循环平移(保留遮罩窗口数和连续段长度,
#   但位置随机), 重算φ与driver模型, 重复NPL次, 与真实ED_sif和不遮罩基线比较。
#   另检查: 遮罩占比本身是否与气候/φ_raw相关(选择偏差)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(glmnet); library(ranger)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify"); NPL <- 30L
src <- readLines("pipelines/recovery_phi_climate_adjusted_ED.R", encoding = "UTF-8")
i_one <- grep("^# ---- 逐站: 气候调整", src)[1]; i_ar1 <- grep("^ar1 <- function", src)[1]; i_os <- grep("^one_station <- function", src)[1]
prefix <- src[1:(i_one - 1)]; prefix <- prefix[!grepl("^NBOOT", prefix)]
invisible(capture.output(eval(parse(text = paste(prefix, collapse = "\n")))))
eval(parse(text = paste(src[i_ar1:(i_os - 1)], collapse = "\n")))

# 真实遮罩占比 vs 气候/phi_raw
frac <- d[cc == TRUE, .(frac_sif = mean(m_sif), frac_vpd = mean(m_vpd)), by = meteo_stat]
raw  <- d[cc == TRUE, ar1(sif_z, year), by = meteo_stat][, .(meteo_stat, phi_raw = phi)]
fx <- Reduce(function(a, b) merge(a, b, by = "meteo_stat"), list(frac, raw, dt18[, .(meteo_stat = stat_id, tavg, rh, precip, rsds_mean)]))
cat("\n=== 遮罩占比(station-level)与气候均值/phi_raw的Spearman ===\n")
print(round(sapply(c("tavg", "rh", "precip", "rsds_mean", "phi_raw"), function(v) c(frac_sif = cor(fx$frac_sif, fx[[v]], method = "spearman"),
                                                                                    frac_vpd = cor(fx$frac_vpd, fx[[v]], method = "spearman"))), 2))
cat(sprintf("遮罩占比: ED_sif 站间 平均%.3f SD%.3f\n", mean(fx$frac_sif), sd(fx$frac_sif)))

fit_r2 <- function(phi_dt) {
  dd <- merge(phi_dt, dt18, by.x = "meteo_stat", by.y = "stat_id"); dd <- dd[is.finite(phi) & complete.cases(dd[, VAR15, with = FALSE])]
  f <- as.formula(paste("yy ~", paste(VAR15, collapse = "+"))); dd[, yy := as.numeric(scale(phi))]
  set.seed(42); rf <- ranger(f, data = dd, num.trees = 500, seed = 42)
  X <- scale(as.matrix(dd[, VAR15, with = FALSE])); set.seed(42); cv <- cv.glmnet(X, dd$yy, alpha = 0, nfolds = 10, standardize = FALSE)
  pr <- predict(cv, X, s = "lambda.1se"); c(n = nrow(dd), rf_r2 = rf$r.squared, ridge_r2 = 1 - sum((dd$yy - pr)^2) / sum((dd$yy - mean(dd$yy))^2)) }
phi_masked <- function(maskcol) d[, ar1(fifelse(get(maskcol) | !cc, NA_real_, sif_z), year), by = meteo_stat][, .(meteo_stat, phi)]
base <- fit_r2(d[, ar1(fifelse(cc, sif_z, NA_real_), year), by = meteo_stat][, .(meteo_stat, phi)])
real <- fit_r2(phi_masked("m_sif"))
cat(sprintf("\n基线(不遮罩): RF R2=%.3f ridge R2=%.3f | 真实ED_sif: RF=%.3f ridge=%.3f\n", base["rf_r2"], base["ridge_r2"], real["rf_r2"], real["ridge_r2"]))

set.seed(7)
pl <- rbindlist(lapply(seq_len(NPL), function(b) {
  d[, m_pl := { n <- .N; if (n > 1L && any(m_sif)) m_sif[((seq_len(n) - 1L + sample.int(n, 1L)) %% n) + 1L] else m_sif }, by = .(meteo_stat, year)]
  r <- fit_r2(phi_masked("m_pl")); data.table(b = b, rf_r2 = r["rf_r2"], ridge_r2 = r["ridge_r2"], frac = mean(d[cc == TRUE]$m_pl)) }))
cat(sprintf("\n安慰剂(随机位置、同窗口数, %d次): RF R2 均值%.3f [%.3f, %.3f] | ridge R2 均值%.3f [%.3f, %.3f] | 遮罩占比%.3f\n", NPL,
            mean(pl$rf_r2), min(pl$rf_r2), max(pl$rf_r2), mean(pl$ridge_r2), min(pl$ridge_r2), max(pl$ridge_r2), mean(pl$frac)))
cat(sprintf("真实ED_sif是否高于全部安慰剂: RF %s | ridge %s (经验p: RF %.3f, ridge %.3f)\n",
            real["rf_r2"] > max(pl$rf_r2), real["ridge_r2"] > max(pl$ridge_r2),
            (1 + sum(pl$rf_r2 >= real["rf_r2"])) / (NPL + 1), (1 + sum(pl$ridge_r2 >= real["ridge_r2"])) / (NPL + 1)))
fwrite(pl, file.path(OUT, "recovery_phi_ED_placebo.csv"))
cat("\n完成。\n")
