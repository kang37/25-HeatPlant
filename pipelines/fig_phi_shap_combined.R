#!/usr/bin/env Rscript
# phi / logOmega 的XGBoost+SHAP合成图: 左=重要性条形图, 右=蜂群图, 纵轴同序(按mean|SHAP|), 无标题
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
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

FS <- 20
th <- theme_bw(base_size = FS) + theme(axis.text.y = element_text(size = FS + 1, color = "black"),
                                       axis.text.x = element_text(size = FS - 2), panel.grid.minor = element_blank())

make_fig <- function(dvname) {
  sh  <- fread(file.path(OUT, sprintf("recovery_xgb_shap_invest_%s_shap_values.csv", dvname)))
  imp <- fread(file.path(OUT, "recovery_xgb_shap_invest_importance.csv"))[dv == dvname][order(-mean_abs_shap)]
  vars <- imp$var
  val <- merge(sh[, .(meteo_stat)], dt18[, c("stat_id", vars), with = FALSE], by.x = "meteo_stat", by.y = "stat_id", sort = FALSE)
  stopifnot(identical(val$meteo_stat, sh$meteo_stat))
  long <- rbindlist(lapply(vars, function(v) data.table(var = v, shap = sh[[v]], value = val[[v]])))
  long[, vs := { rg <- range(value, na.rm = TRUE); if (diff(rg) == 0) rep(0.5, .N) else (value - rg[1]) / diff(rg) }, by = var]
  long[, var := factor(var, levels = rev(vars))]
  imp[, var := factor(var, levels = rev(vars))]

  p_bar <- ggplot(imp, aes(var, mean_abs_shap, fill = var == "invest")) +
    geom_col(width = 0.75) + coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#C1121F", `FALSE` = "#2166AC"), guide = "none") +
    labs(x = NULL, y = "mean |SHAP|") + th
  p_bee <- ggplot(long, aes(shap, var, color = vs)) +
    geom_vline(xintercept = 0, color = "grey60") +
    geom_jitter(height = 0.22, size = 2.2, alpha = 0.85) +
    scale_color_gradient(low = "#2166AC", high = "#C1121F", na.value = "grey60",
                         breaks = c(0, 1), labels = c("Low", "High"), name = "Feature value") +
    labs(x = "SHAP value", y = NULL) + th +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.title = element_text(size = FS - 2), legend.text = element_text(size = FS - 2),
          legend.key.height = unit(1.4, "cm"))
  p <- p_bar + p_bee + plot_layout(widths = c(1, 1.45))
  f <- file.path(OUT, sprintf("recovery_xgb_shap_invest_%s_combined.png", dvname))
  ggsave(f, p, width = 15, height = if (length(vars) > 13) 8 else 6.8, dpi = 300)
  cat("->", basename(f), "\n")
}
for (d_ in c("phi", "logOmega")) make_fig(d_)
