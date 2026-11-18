#!/usr/bin/env Rscript
# =============================================================================
# recovery_speed_drivers_117.R
#   参考 Lu et al. 2025 (Int J Digital Earth, 10.1080/17538947.2025.2507195)
#   "Quantifying global vegetation recovery speed in response to extreme
#   drought" 的设计——把逐站恢复速度当因变量，用环境协变量(随机森林)找主导
#   驱动因素——在本项目 117 个 CCM 确认因果耦合的站上做同类分析。
#
#   和原文的差异(据摘要+检索到的方法描述，全文因付费墙未能完整获取，方法细节
#   按摘要能确认的部分对齐，其余按本项目已有做法填空，需要在此说明):
#   - 原文: 全球像元尺度，2018-2021一次极端干旱事件，自变量是SM/TMP/PRE/VPD
#     这类"事件期间的动态环境条件"，驱动因素占比按"在多少%的像元上是首位重要
#     变量"统计。
#   - 本脚本: 中国117站尺度，恢复速度=该站多年(2000-2022)所有未删失站-年事件
#     的中位恢复窗口数(8天为单位)，自变量换成用户指定的18个CORE协变量——这些
#     是站点静态属性(地表覆盖/气候常年值/城市化/地形)，不是事件期动态环境条件，
#     回答的问题也相应变成"哪种类型的站点、平均而言恢复更快/更慢"，不是原文
#     "同一站点不同强度环境条件下恢复速度如何变化"，这个概念差异必须在解读时
#     交代清楚。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ranger)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. 18个CORE协变量(与已删除的44-60主线用的是同一套定义, 复用22号脚本) ----
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt   <- env$dt
PRED <- unlist(env$CORE, use.names = FALSE)
cat("18个CORE协变量:", paste(PRED, collapse=", "), "\n")

# ---- 2. 恢复速度: 117站的中位恢复窗口数(仅用未删失事件, 同Lu et al.常见简化) --
ev <- fread("data_proc/output_loose_classify/recovery_vs_134_direction_events.csv")
spd <- ev[event == 1, .(recovery_windows = median(Trec_windows), n_recovered = .N),
          by = meteo_stat]
cens <- ev[, .(n_events = .N, pct_censored = round(100*mean(censored),1)), by = meteo_stat]
spd <- merge(spd, cens, by = "meteo_stat")
cat("\n117站里可算出恢复速度的站数:", nrow(spd),
    "(", uniqueN(ev$meteo_stat) - nrow(spd), "站全部事件删失，无法计算)\n")

# ---- 3. 合并协变量, 只保留18变量完整的站 -------------------------------------
d <- merge(spd, dt, by.x = "meteo_stat", by.y = "stat_id")
d <- d[complete.cases(d[, PRED, with = FALSE])]
cat("协变量完整、可建模的站数:", nrow(d), "\n")

# ---- 4. 随机森林(仿Lu et al.用RF找主导驱动因素) + 简单OLS(方向/显著性参考) ---
set.seed(42)
f <- as.formula(paste("recovery_windows ~", paste(PRED, collapse = "+")))
rf <- ranger(f, data = d, num.trees = 2000, importance = "permutation", seed = 42)
imp <- data.table(var = names(rf$variable.importance), importance = rf$variable.importance)
setorder(imp, -importance)
imp[, share := round(100*importance/sum(pmax(importance,0)), 1)]
cat("\n=== RF variable importance(permutation), R2(OOB) =", round(rf$r.squared,3), "===\n")
print(imp)

ols <- lm(f, data = d)
cat("\n=== OLS(方向/显著性参考, R2=", round(summary(ols)$r.squared,3), ") ===\n")
print(round(summary(ols)$coefficients[,c(1,4)], 4))

fwrite(d, file.path(OUT, "recovery_speed_drivers_117_data.csv"))
fwrite(imp, file.path(OUT, "recovery_speed_drivers_117_importance.csv"))

# ---- 5. 图: 变量重要性条形图 ------------------------------------------------
imp[, var := factor(var, levels = rev(var))]
p <- ggplot(imp, aes(var, importance)) +
  geom_col(fill = "#2166AC") + coord_flip() +
  labs(title = "RF variable importance for station-level recovery speed (117 CCM-confirmed stations)",
       subtitle = sprintf("DV = median recovery windows (8-day units) across non-censored station-years | RF OOB R2=%.2f | ref: Lu et al. 2025 design, static covariates substituted for event-time SM/TMP/PRE/VPD",
                           rf$r.squared),
       x = NULL, y = "Permutation importance") +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size = 8, color = "grey35"))
ggsave(file.path(OUT, "recovery_speed_drivers_117_importance.png"), p, width = 8, height = 6, dpi = 300)
cat("\n-> recovery_speed_drivers_117_importance.png\n")
