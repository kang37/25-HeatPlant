#!/usr/bin/env Rscript
# =============================================================================
# recovery_liu2021_drivers_117.R
#   ① log变换 Omega/Delta(严格正值,无需offset); constraint_time/recovery_time
#      是有界计数变量,不处理。
#   ② 四个量分别聚合到站点层面(log量取均值,计数量取中位数), 和18个CORE协变量
#      做driver回归(RF permutation importance + OLS), 同一套流程跑四遍。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ranger)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. log变换 --------------------------------------------------------------
ev <- fread(file.path(OUT, "recovery_liu2021_metrics_117_events.csv"))
ev[, `:=`(logOmega = log(Omega), logDelta = log(Delta))]
cat("=== log变换后分布 ===\n")
print(ev[, .(logOmega_mean=round(mean(logOmega),2), logOmega_sd=round(sd(logOmega),2),
             logDelta_mean=round(mean(logDelta),2), logDelta_sd=round(sd(logDelta),2))])

# ---- 2. 站点层面聚合: 计数量取中位数, log量取均值(对数正态假设下更自然) -----
agg <- ev[, .(constraint_time = median(constraint_time, na.rm=TRUE),
              recovery_time   = median(recovery_time,   na.rm=TRUE),
              logOmega = mean(logOmega, na.rm=TRUE),
              logDelta = mean(logDelta, na.rm=TRUE),
              n_events = .N), by = meteo_stat]
cat("\n站点数(聚合后):", nrow(agg), "\n")

# ---- 3. 18个CORE协变量 -------------------------------------------------------
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt <- env$dt; PRED <- unlist(env$CORE, use.names = FALSE)

d <- merge(agg, dt, by.x = "meteo_stat", by.y = "stat_id")
d <- d[complete.cases(d[, PRED, with = FALSE])]
cat("协变量完整、可建模的站数:", nrow(d), "\n")
fwrite(d, file.path(OUT, "recovery_liu2021_drivers_117_data.csv"))

# ---- 4. 四个DV各跑一遍: RF(permutation importance) + OLS --------------------
DVS <- c("constraint_time", "recovery_time", "logOmega", "logDelta")
run_one <- function(dv) {
  f <- as.formula(paste(dv, "~", paste(PRED, collapse = "+")))
  set.seed(42)
  rf <- ranger(f, data = d, num.trees = 2000, importance = "permutation", seed = 42)
  imp <- data.table(dv = dv, var = names(rf$variable.importance),
                     importance = rf$variable.importance)
  ols <- lm(f, data = d)
  co  <- summary(ols)$coefficients
  sig <- rownames(co)[-1][co[-1, "Pr(>|t|)"] < 0.05]
  list(imp = imp, rf_r2 = rf$r.squared, ols_r2 = summary(ols)$r.squared,
       sig_vars = paste(sig, collapse = "; "))
}
res <- lapply(DVS, run_one)
names(res) <- DVS

cat("\n=== 四个量的driver回归汇总 ===\n")
summary_tab <- rbindlist(lapply(DVS, function(dv) data.table(
  DV = dv, RF_OOB_R2 = round(res[[dv]]$rf_r2, 3), OLS_R2 = round(res[[dv]]$ols_r2, 3),
  OLS_sig_p05 = ifelse(res[[dv]]$sig_vars == "", "(none)", res[[dv]]$sig_vars))))
print(summary_tab)
fwrite(summary_tab, file.path(OUT, "recovery_liu2021_drivers_117_summary.csv"))

imp_all <- rbindlist(lapply(res, `[[`, "imp"))
fwrite(imp_all, file.path(OUT, "recovery_liu2021_drivers_117_importance.csv"))

cat("\n=== 每个DV的RF重要性前5 ===\n")
for (dv_i in DVS) {
  cat("\n--", dv_i, "--\n")
  print(head(imp_all[dv == dv_i][order(-importance)], 5))
}

# ---- 5. 图: 四个DV的变量重要性并排 -------------------------------------------
imp_all[, rank_within_dv := rank(-importance), by = dv]
top_vars <- imp_all[rank_within_dv <= 18]
p <- ggplot(top_vars, aes(reorder(var, importance), importance)) +
  geom_col(fill = "#2166AC") + coord_flip() +
  facet_wrap(~dv, scales = "free", ncol = 2) +
  labs(title = "RF variable importance across four Liu et al. 2021-style recovery metrics (117 stations)",
       subtitle = "logOmega/logDelta = log-transformed resistance/resilience (Isbell 2015 formulas); constraint/recovery_time in 8-day windows",
       x = NULL, y = "Permutation importance") +
  theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size=8, color="grey35"))
ggsave(file.path(OUT, "recovery_liu2021_drivers_117_importance.png"), p, width = 11, height = 8, dpi = 300)
cat("\n-> recovery_liu2021_drivers_117_importance.png\n")
