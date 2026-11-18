#!/usr/bin/env Rscript
# =============================================================================
# recovery_liu2021_drivers_dirfix_117.R
#   在方向修正后的事件表(57站: Promote 41 / Inhibit 16, 已剔除Ambiguous)上
#   重跑站点层面 driver 回归(RF + OLS + 岭回归), 对照未修正版本
#   (recovery_liu2021_drivers_117.R / recovery_liu2021_ridge_117.R, 86站含全部
#   方向)。样本量比之前小(57 vs 86/117), 统计功效会更低, 结果解读要留有余地。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ranger); library(glmnet)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

ev <- fread(file.path(OUT, "recovery_liu2021_metrics_dirfix_117_events.csv"))
agg <- ev[, .(constraint_time = as.numeric(median(constraint_time, na.rm=TRUE)),
              recovery_time   = as.numeric(median(recovery_time,   na.rm=TRUE)),
              logOmega = as.numeric(mean(log(Omega), na.rm=TRUE)),
              logDelta = as.numeric(mean(log(Delta), na.rm=TRUE)),
              direction_2v = direction_2v[1], n_events = .N), by = meteo_stat]
cat("站点数(聚合后):", nrow(agg), " | 方向分布:\n"); print(agg[, .N, by = direction_2v])

env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt <- env$dt; PRED <- unlist(env$CORE, use.names = FALSE)

d <- merge(agg, dt, by.x = "meteo_stat", by.y = "stat_id")
d <- d[complete.cases(d[, PRED, with = FALSE])]
cat("协变量完整、可建模的站数:", nrow(d), "\n")
fwrite(d, file.path(OUT, "recovery_liu2021_drivers_dirfix_117_data.csv"))

DVS <- c("constraint_time", "recovery_time", "logOmega", "logDelta")

# ---- RF ----------------------------------------------------------------------
run_rf <- function(dv) {
  f <- as.formula(paste(dv, "~", paste(PRED, collapse = "+")))
  set.seed(42)
  rf <- ranger(f, data = d, num.trees = 2000, importance = "permutation", seed = 42)
  list(r2 = rf$r.squared,
       imp = data.table(dv = dv, var = names(rf$variable.importance), importance = rf$variable.importance))
}
rf_res <- lapply(DVS, run_rf); names(rf_res) <- DVS

# ---- 岭回归(处理共线性, 同之前对constraint_time/logOmega做过的) -------------
X <- scale(as.matrix(d[, PRED, with = FALSE]))
run_ridge <- function(dv) {
  y <- d[[dv]]
  set.seed(42)
  cv <- cv.glmnet(X, y, alpha = 0, nfolds = 10, standardize = FALSE)
  co_1se <- as.numeric(coef(cv, s = "lambda.1se"))[-1]
  pred <- predict(cv, X, s = "lambda.1se")
  r2 <- 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
  data.table(dv = dv, var = PRED, ridge_1se = co_1se, ridge_r2 = r2)
}
ridge_res <- rbindlist(lapply(DVS, run_ridge))

cat("\n=== 汇总: RF OOB R2 vs 岭回归(1se) R2 ===\n")
summ <- data.table(dv = DVS,
                    rf_r2 = round(sapply(rf_res, `[[`, "r2"), 3),
                    ridge_r2 = round(unique(ridge_res[, .(dv, ridge_r2)])$ridge_r2, 3))
print(summ)
fwrite(summ, file.path(OUT, "recovery_liu2021_drivers_dirfix_117_summary.csv"))

cat("\n=== 岭回归系数(|1se|排序前5), 各DV ===\n")
for (dv_i in DVS) {
  cat("\n--", dv_i, "--\n")
  print(ridge_res[dv == dv_i][order(-abs(ridge_1se))][1:5, .(var, ridge_1se = round(ridge_1se,3))])
}
fwrite(ridge_res, file.path(OUT, "recovery_liu2021_drivers_dirfix_117_ridge.csv"))

imp_all <- rbindlist(lapply(rf_res, `[[`, "imp"))
imp_all[, rank_within_dv := rank(-importance), by = dv]
p <- ggplot(imp_all[rank_within_dv <= 18], aes(reorder(var, importance), importance)) +
  geom_col(fill = "#2166AC") + coord_flip() + facet_wrap(~dv, scales = "free", ncol = 2) +
  labs(title = "RF importance, direction-corrected recovery metrics (57 stations: Promote+Inhibit only)",
       subtitle = sprintf("n=%d stations | vs uncorrected 86-station version, see docs/03_classification.md", nrow(d)),
       x = NULL, y = "Permutation importance") +
  theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size=8, color="grey35"))
ggsave(file.path(OUT, "recovery_liu2021_drivers_dirfix_117_importance.png"), p, width = 11, height = 8, dpi = 300)
cat("\n-> recovery_liu2021_drivers_dirfix_117_importance.png\n")
