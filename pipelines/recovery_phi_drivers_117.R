#!/usr/bin/env Rscript
# =============================================================================
# recovery_phi_drivers_117.R
#   因变量换成 phi(连续型AR1恢复力指标, recovery_ar1_pure.R 已算好), 而不是
#   噪声更大的 Trec(事件恢复窗口数)。自变量沿用用户指定的18个CORE协变量。
#   phi 全站样本量都是"该站全部23年数据一次性估计"(n=200-400+/站), 不像
#   Trec 中位数那样要先把每站压缩成一个高噪声的非参数汇总, 预期信噪比更好。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ranger)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. 18个CORE协变量(复用22号脚本) ----------------------------------------
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt   <- env$dt
PRED <- unlist(env$CORE, use.names = FALSE)

# ---- 2. phi(117站, 已由recovery_ar1_pure.R算好) ----------------------------
phi <- fread("data_proc/output_loose_classify/recovery_ar1_pure_vs_arx.csv")
cat("117站phi清单 n =", nrow(phi), "\n")

d <- merge(phi, dt, by.x = "meteo_stat", by.y = "stat_id")
d <- d[complete.cases(d[, PRED, with = FALSE])]
cat("协变量完整、可建模的站数:", nrow(d), "\n")

# ---- 3. RF(phi1_pure为主) + OLS(方向/显著性参考) + phi2_arx做稳健性对照 -----
set.seed(42)
f1 <- as.formula(paste("phi1_pure ~", paste(PRED, collapse = "+")))
rf1 <- ranger(f1, data = d, num.trees = 2000, importance = "permutation", seed = 42)
imp1 <- data.table(var = names(rf1$variable.importance), importance = rf1$variable.importance)
setorder(imp1, -importance)
cat("\n=== RF(phi1_pure) variable importance, OOB R2 =", round(rf1$r.squared,3), "===\n")
print(imp1)

ols1 <- lm(f1, data = d)
cat("\n=== OLS(phi1_pure, R2=", round(summary(ols1)$r.squared,3), ") ===\n")
co1 <- summary(ols1)$coefficients
print(round(co1[order(co1[,4]), c(1,4)], 4))

# 稳健性: phi2_arx(控制了VPD的版本)重要性排序是否一致
f2  <- as.formula(paste("phi2_arx ~", paste(PRED, collapse = "+")))
rf2 <- ranger(f2, data = d, num.trees = 2000, importance = "permutation", seed = 42)
imp2 <- data.table(var = names(rf2$variable.importance), importance_phi2 = rf2$variable.importance)
imp_cmp <- merge(imp1, imp2, by = "var")
setorder(imp_cmp, -importance)
cat("\n=== phi1 vs phi2 重要性排序稳健性对照, OOB R2(phi2) =", round(rf2$r.squared,3), "===\n")
print(imp_cmp)
cat("\n重要性排序 Spearman rho(phi1 vs phi2) =",
    round(cor(imp_cmp$importance, imp_cmp$importance_phi2, method="spearman"), 3), "\n")

fwrite(d, file.path(OUT, "recovery_phi_drivers_117_data.csv"))
fwrite(imp_cmp, file.path(OUT, "recovery_phi_drivers_117_importance.csv"))

# ---- 4. 图 -------------------------------------------------------------------
imp1[, var := factor(var, levels = rev(var))]
p <- ggplot(imp1, aes(var, importance)) +
  geom_col(fill = "#2166AC") + coord_flip() +
  labs(title = "RF variable importance for station-level resilience (phi, 117 CCM-confirmed stations)",
       subtitle = sprintf("DV = phi1_pure (whole-record AR1 coefficient, Forzieri-style) | RF OOB R2=%.2f",
                           rf1$r.squared),
       x = NULL, y = "Permutation importance") +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size = 8.5, color = "grey35"))
ggsave(file.path(OUT, "recovery_phi_drivers_117_importance.png"), p, width = 8, height = 6, dpi = 300)
cat("\n-> recovery_phi_drivers_117_importance.png\n")
