#!/usr/bin/env Rscript
# 用已保存的SHAP值重画依赖图(大字号), 逻辑同recovery_xgb_shap_invest.R的p3, 不重跑建模
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT <- file.path(PROJ, "data_proc/output_loose_classify")
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18 <- copy(env$dt)
cl <- fread("data_proc/smap_bivar_134/classify_134.csv")[, .(meteo_stat, direction_2v)]
dt18 <- merge(dt18, cl, by.x = "stat_id", by.y = "meteo_stat", all.x = TRUE)
dt18[, is_promote := as.integer(direction_2v == "Promote")]
impall <- fread(file.path(OUT, "recovery_xgb_shap_invest_importance.csv"))
for (dvname in c("phi", "logOmega")) {
  sh  <- fread(file.path(OUT, sprintf("recovery_xgb_shap_invest_%s_shap_values.csv", dvname)))
  imp <- impall[dv == dvname][order(-mean_abs_shap)]
  vars <- imp$var
  val <- merge(sh[, .(meteo_stat)], dt18[, c("stat_id", vars), with = FALSE], by.x = "meteo_stat", by.y = "stat_id", sort = FALSE)
  stopifnot(identical(val$meteo_stat, sh$meteo_stat))
  top6 <- unique(c("invest", vars[imp$rank <= 6]))[1:6]
  dep <- rbindlist(lapply(top6, function(v) data.table(var = v, shap = sh[[v]], value = val[[v]])))[!is.na(value)]
  dep[, var := factor(var, levels = top6)]
  p <- ggplot(dep, aes(value, shap)) + geom_hline(yintercept = 0, color = "grey70") +
    geom_point(alpha = 0.6, size = 2.4, color = "#2166AC") +
    geom_smooth(method = "loess", se = FALSE, color = "#C1121F", linewidth = 1.2, formula = y ~ x) +
    facet_wrap(~var, scales = "free", ncol = 3) +
    labs(title = sprintf("SHAP dependence: %s (invest + top features)", dvname), x = "feature value", y = "SHAP value") +
    theme_bw(base_size = 18) + theme(strip.text = element_text(size = 19), plot.title = element_text(size = 20))
  ggsave(file.path(OUT, sprintf("recovery_xgb_shap_invest_%s_dependence.png", dvname)), p, width = 13, height = 8, dpi = 300)
  cat("->", sprintf("recovery_xgb_shap_invest_%s_dependence.png", dvname), "\n")
}
