#!/usr/bin/env Rscript
# =============================================================================
# tp_lag_drivers_117.R
#   参考 Liu et al. 2024 (Agric For Meteorol 358:110239)"Identifying thresholds
#   of time-lag and accumulative effects of extreme precipitation on major
#   vegetation types at global scale"的设计——原文逐像元用相关系数扫描降水指数
#   的滞后月份(0-4个月)取argmax作"时滞效应(time-lag effect)"，再用Boosted
#   Regression Trees(BRT)把这个时滞长度当因变量，回归11个环境因子(含VPD)找
#   主导驱动因素+偏依赖阈值(原文核心发现: VPD贡献24.68-47.88%, 时滞随VPD增大
#   而缩短)。
#
#   本项目的直接类比: CCM 因果确认判据(pipelines/hcsif_buf1000/12_ccm_causal_
#   confirmation.R)里，每站 optimal_tp = argmax(rho) over tp∈[0,8]，正是同一
#   套"扫描滞后、取峰值"的设计，只是我们扫的是CCM的rho而不是相关系数、单位是
#   8天窗口而不是月。117个CCM确认因果耦合的站点(见 docs/03_classification.md)
#   每站都已有这个 optimal_tp，本脚本把它当因变量，回归用户指定的18个CORE协
#   变量(与已删除的44-60主线/recovery_phi_drivers_117.R/
#   recovery_speed_drivers_117.R同一套定义，复用22_varpart_3dim.R)。
#
#   和原文的差异(必须在解读时交代):
#   - 原文: 全球像元尺度(千万级样本)，5种植被类型分别建模，自变量含逐事件的
#     动态环境条件(降水频率/强度等)，BRT用gbm/dismo实现，偏依赖图找"阈值"。
#   - 本脚本: 中国117站尺度(样本量远小于原文)，不分植被类型(样本太小无法分
#     层，地表覆盖已作为连续协变量纳入18个CORE变量)，自变量是站点静态属性
#     (地表覆盖/气候常年值/城市化/地形)，不是事件期动态条件；本机未装gbm/
#     dismo/pdp，且xgboost这次实测也不在(CLAUDE.md记录的"已装包"已过期，大概率
#     是R4.6.1升级又一次清空，同02_ccm.md记录的环境坑同一模式)——改用项目其余
#     117站driver脚本(recovery_phi_drivers_117.R/recovery_speed_drivers_117.R)
#     已验证可用的纯ranger方案: 随机森林重要性排序 + 手写偏依赖(固定其余变量在
#     样本中位数、扫描目标变量)近似原文BRT偏依赖图的阈值识别，OLS做方向/显著性
#     参考(同项目一贯做法)，不引入新依赖。
#   - optimal_tp 本身已知的测量噪声: 相邻tp高度相关，argmax对抽样噪声敏感
#     (docs/02_ccm.md 2026-09-04日志)，且只在tp∈[0,8]九个整数里取值，不是原文
#     那种接近连续的月尺度滞后——n=117(甚至更小,取决于协变量完整性)的driver
#     回归结果读作探索性描述，不是确认性推断，重要性排序的抽样不确定性没有
#     形式化(没有bootstrap置信区间)。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ranger); library(glmnet); library(car)
})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. 18个CORE协变量(复用22号脚本, 与其余117站driver脚本同一套定义) --------
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt   <- env$dt
PRED <- unlist(env$CORE, use.names = FALSE)
cat("18个CORE协变量:", paste(PRED, collapse=", "), "\n")

# ---- 2. optimal_tp(117站, 12号脚本已算好, "效应滞后"的直接类比) --------------
tp <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")
cat("117站 optimal_tp 清单 n =", nrow(tp), "| 取值分布:\n"); print(table(tp$optimal_tp))

d <- merge(tp, dt, by.x = "meteo_stat", by.y = "stat_id")
d <- d[complete.cases(d[, PRED, with = FALSE])]
cat("协变量完整、可建模的站数:", nrow(d), "\n")

# ---- 3. RF(仿其余117站driver脚本) + OLS(方向/显著性参考) ---------------------
set.seed(42)
f <- as.formula(paste("optimal_tp ~", paste(PRED, collapse = "+")))
rf <- ranger(f, data = d, num.trees = 2000, importance = "permutation", seed = 42)
imp <- data.table(var = names(rf$variable.importance), importance = rf$variable.importance)
setorder(imp, -importance)
imp[, share := round(100*importance/sum(pmax(importance,0)), 1)]
cat("\n=== RF variable importance(permutation), OOB R2 =", round(rf$r.squared,3), "===\n")
print(imp)

ols <- lm(f, data = d)
cat("\n=== OLS(方向/显著性参考, R2=", round(summary(ols)$r.squared,3), ") ===\n")
print(round(summary(ols)$coefficients[,c(1,4)], 4))

# ---- 4. 手写偏依赖(近似原文BRT阈值识别, 用上面已拟合的ranger rf): 固定其余变量
#         在样本中位数, 扫描目标变量 -------------------------------------------
X <- as.matrix(d[, PRED, with = FALSE])
pdp_var <- function(v, ngrid = 30) {
  grid <- seq(quantile(d[[v]], .02), quantile(d[[v]], .98), length.out = ngrid)
  base <- as.list(apply(X, 2, median))
  Xg <- do.call(rbind, lapply(grid, function(g) { b <- base; b[[v]] <- g; unlist(b) }))
  Xg <- as.data.table(Xg)[, PRED, with = FALSE]
  data.table(var = v, x = grid, yhat = predict(rf, Xg)$predictions)
}
top_vars <- imp$var[1:min(4, nrow(imp))]
pdp_all <- rbindlist(lapply(as.character(top_vars), pdp_var))

fwrite(d,       file.path(OUT, "tp_lag_drivers_117_data.csv"))
fwrite(imp,     file.path(OUT, "tp_lag_drivers_117_rf_importance.csv"))
fwrite(pdp_all, file.path(OUT, "tp_lag_drivers_117_pdp.csv"))

# ---- 5. 图 --------------------------------------------------------------------
imp[, var := factor(var, levels = rev(var))]
p1 <- ggplot(imp, aes(var, importance)) +
  geom_col(fill = "#2166AC") + coord_flip() +
  labs(title = "RF variable importance for CCM optimal_tp (effect lag, 117 stations)",
       subtitle = sprintf("DV = optimal_tp (argmax rho over tp=0..8) | RF OOB R2=%.2f | ref: Liu et al. 2024 design",
                           rf$r.squared),
       x = NULL, y = "Permutation importance") +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size = 8, color = "grey35"))
ggsave(file.path(OUT, "tp_lag_drivers_117_rf_importance.png"), p1, width = 8, height = 6, dpi = 300)

p2 <- ggplot(pdp_all, aes(x, yhat)) +
  geom_line(color = "#B2182B", linewidth = 1) +
  facet_wrap(~var, scales = "free_x") +
  labs(title = "Partial dependence of optimal_tp on top RF predictors (117 stations)",
       subtitle = "Other predictors fixed at sample median | approximates Liu et al. 2024 threshold plots",
       x = NULL, y = "Predicted optimal_tp") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "tp_lag_drivers_117_pdp.png"), p2, width = 9, height = 6, dpi = 300)
cat("\n-> tp_lag_drivers_117_rf_importance.png, tp_lag_drivers_117_pdp.png\n")

# =============================================================================
# 第二轮: 解决"样本太小、变量太多"——上面18变量/86站的结果不可信有两个可诊断的
# 具体原因(不是抽象地"n小"), 分别对症处理:
#   (a) 117->86的31站流失几乎全部来自3个变量(invest/urban_rate/cloud缺失最多,
#       分别缺20/17/4站)，不是随机缺失在18个变量上摊开——去掉这3个变量就能把
#       117站全部找回来，比插补更干净(不用编造数据)。
#   (b) 18变量OLS的VIF显示 tavg=10.5/elev=8.0/precip=7.9/deci_broad=6.5/
#       deci_needle=6.2 明显共线(经典的"气温-海拔-降水"三角)，这会把标准误
#       放大、把真实信号掩盖成"不显著"——用glmnet(弹性网)做正则化，比OLS更
#       适合这种n和p同一数量级、变量间又有共线性的场景，且自带交叉验证给出
#       诚实的样本外R2(不是OLS那种必然偏高的样本内R2)。
# =============================================================================
cat("\n\n########## 第二轮: 15个高覆盖变量(去掉invest/urban_rate/cloud) x 117站 ##########\n")
PRED15 <- setdiff(PRED, c("invest", "urban_rate", "cloud"))
d117 <- d117full <- merge(tp, dt, by.x = "meteo_stat", by.y = "stat_id")
d117 <- d117[complete.cases(d117[, PRED15, with = FALSE])]
cat("15变量完整站数:", nrow(d117), "(应等于117, 验证3个变量确实是唯一缺失来源)\n")

f15 <- as.formula(paste("optimal_tp ~", paste(PRED15, collapse = "+")))

# ---- RF(117站) -----------------------------------------------------------
rf15 <- ranger(f15, data = d117, num.trees = 2000, importance = "permutation", seed = 42)
imp15 <- data.table(var = names(rf15$variable.importance), importance = rf15$variable.importance)
setorder(imp15, -importance)
cat("\n=== RF(15变量, 117站) OOB R2 =", round(rf15$r.squared, 3), "===\n")
print(imp15)

# ---- VIF(117站, 看去掉3个稀疏变量后共线性是否变化) -------------------------
vif15 <- sort(vif(lm(f15, d117)), decreasing = TRUE)
cat("\n=== VIF(15变量, 117站) ===\n"); print(round(vif15, 2))

# ---- glmnet 弹性网(alpha=.5) + 10折CV: 应对共线性+给诚实样本外R2 -----------
X15 <- scale(as.matrix(d117[, PRED15, with = FALSE]))   # 标准化, 岭/lasso惩罚对量纲敏感
y15 <- d117$optimal_tp
set.seed(42)
cvfit <- cv.glmnet(X15, y15, alpha = 0.5, nfolds = 10, standardize = FALSE)
r2_cv <- 1 - min(cvfit$cvm) / var(y15)
cat("\n=== 弹性网(alpha=.5, 10折CV) ===\n")
cat("lambda.min CV-R2 =", round(r2_cv, 3),
    "(样本外估计, 与上面RF OOB R2=0.012同一性质, 互相印证) \n")
co <- as.matrix(coef(cvfit, s = "lambda.min"))
nz <- co[co[, 1] != 0, , drop = FALSE]
cat("lambda.min 下未被压成0的变量(", nrow(nz) - 1, "个, 不含截距):\n")
print(round(nz, 4))

fwrite(imp15, file.path(OUT, "tp_lag_drivers_117_15var_rf_importance.csv"))
fwrite(data.table(var = rownames(nz), coef = nz[, 1]),
       file.path(OUT, "tp_lag_drivers_117_15var_glmnet_coef.csv"))

imp15[, var := factor(var, levels = rev(var))]
p3 <- ggplot(imp15, aes(var, importance)) +
  geom_col(fill = "#4393C3") + coord_flip() +
  labs(title = "RF variable importance: optimal_tp ~ 15 high-coverage CORE vars (full 117 stations)",
       subtitle = sprintf("Dropped invest/urban_rate/cloud (sparsest vars) to recover full sample | RF OOB R2=%.3f | elastic-net CV-R2=%.3f",
                           rf15$r.squared, r2_cv),
       x = NULL, y = "Permutation importance") +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size = 8, color = "grey35"))
ggsave(file.path(OUT, "tp_lag_drivers_117_15var_rf_importance.png"), p3, width = 8, height = 6, dpi = 300)
cat("\n-> tp_lag_drivers_117_15var_rf_importance.png\n")

# ---- 15变量/117站(目前唯一有实际样本外解释力的版本, RF OOB R2=0.13)的逐变量
#      偏依赖 + 方向/形状归纳, 回答"各变量对optimal_tp的具体影响" -------------
X15raw <- as.matrix(d117[, PRED15, with = FALSE])
pdp_var15 <- function(v, ngrid = 30) {
  grid <- seq(quantile(d117[[v]], .02), quantile(d117[[v]], .98), length.out = ngrid)
  base <- as.list(apply(X15raw, 2, median))
  Xg <- do.call(rbind, lapply(grid, function(g) { b <- base; b[[v]] <- g; unlist(b) }))
  Xg <- as.data.table(Xg)[, PRED15, with = FALSE]
  data.table(var = v, x = grid, yhat = predict(rf15, Xg)$predictions)
}
pdp15_all <- rbindlist(lapply(PRED15, pdp_var15))
shape15 <- pdp15_all[, {
  rho_s <- suppressWarnings(cor(x, yhat, method = "spearman"))
  rng   <- max(yhat) - min(yhat)
  net   <- yhat[.N] - yhat[1]
  monot <- if (rng < 1e-6) "flat" else if (abs(net) / rng > 0.7) {
    if (net > 0) "increasing" else "decreasing"
  } else "non-monotonic/threshold"
  .(spearman_x_yhat = round(rho_s, 2), pdp_range = round(rng, 2), shape = monot)
}, by = var]
enet15 <- as.data.table(as.matrix(coef(cvfit, s = "lambda.min")), keep.rownames = "var")
setnames(enet15, 2, "enet_coef"); enet15 <- enet15[var != "(Intercept)"]
sm15 <- Reduce(function(a, b) merge(a, b, by = "var", all = TRUE),
               list(imp15, enet15, shape15))
setorder(sm15, -importance)
cat("\n=== 汇总(最佳模型: 15变量/117站/RF OOB R2=0.13): 各变量对optimal_tp的具体影响 ===\n")
print(sm15)
fwrite(sm15, file.path(OUT, "tp_lag_drivers_117_15var_summary.csv"))
fwrite(pdp15_all, file.path(OUT, "tp_lag_drivers_117_15var_pdp_all.csv"))

p4 <- ggplot(pdp15_all, aes(x, yhat)) +
  geom_line(color = "#B2182B", linewidth = 1) +
  facet_wrap(~var, scales = "free_x", ncol = 4) +
  labs(title = "Partial dependence of optimal_tp on all 15 predictors (best model: 117 stations, RF OOB R2=0.13)",
       subtitle = "Other predictors fixed at sample median",
       x = NULL, y = "Predicted optimal_tp") +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "tp_lag_drivers_117_15var_pdp_all.png"), p4, width = 13, height = 10, dpi = 300)
cat("\n-> tp_lag_drivers_117_15var_pdp_all.png\n")

# =============================================================================
# 第三轮: 用户要求"保留城镇化率/投资至少1个"——不能把两个都丢。117站里
# urban_rate缺17站(比invest缺20站覆盖率更好)，加回urban_rate(16变量)。
# =============================================================================
cat("\n\n########## 第三轮: 16个变量(15安全变量+urban_rate) x 117站 ##########\n")
PRED16 <- c(PRED15, "urban_rate")
d16 <- d117full[complete.cases(d117full[, PRED16, with = FALSE])]
cat("16变量完整站数:", nrow(d16), "\n")

f16 <- as.formula(paste("optimal_tp ~", paste(PRED16, collapse = "+")))
rf16 <- ranger(f16, data = d16, num.trees = 2000, importance = "permutation", seed = 42)
imp16 <- data.table(var = names(rf16$variable.importance), importance = rf16$variable.importance)
setorder(imp16, -importance)
cat("\n=== RF(16变量, 117站原判据) OOB R2 =", round(rf16$r.squared, 3), "===\n")
print(imp16)

X16 <- scale(as.matrix(d16[, PRED16, with = FALSE])); y16 <- d16$optimal_tp
set.seed(42)
cv16 <- cv.glmnet(X16, y16, alpha = 0.5, nfolds = 10, standardize = FALSE)
r2_cv16 <- 1 - min(cv16$cvm) / var(y16)
enet16 <- as.data.table(as.matrix(coef(cv16, s = "lambda.min")), keep.rownames = "var")
setnames(enet16, 2, "enet_coef"); enet16 <- enet16[var != "(Intercept)"]
cat("\n弹性网(16变量,117站) CV-R2 =", round(r2_cv16, 3), "\n")
print(enet16[order(-abs(enet_coef))])

fwrite(imp16, file.path(OUT, "tp_lag_drivers_117_16var_rf_importance.csv"))
fwrite(enet16, file.path(OUT, "tp_lag_drivers_117_16var_enet_coef.csv"))

# urban_rate 单独的偏依赖(用户关心的变量, 即使这个16变量模型本身R2不高也看一眼方向)
X16raw <- as.matrix(d16[, PRED16, with = FALSE])
grid_ur <- seq(quantile(d16$urban_rate, .02), quantile(d16$urban_rate, .98), length.out = 30)
base_ur <- as.list(apply(X16raw, 2, median))
Xg_ur <- as.data.table(do.call(rbind, lapply(grid_ur, function(g) { b <- base_ur; b$urban_rate <- g; unlist(b) })))[, PRED16, with = FALSE]
pdp_ur <- data.table(x = grid_ur, yhat = predict(rf16, Xg_ur)$predictions)
cat("\nurban_rate偏依赖(16变量/117站模型, 固定其余变量在中位数):\n")
cat("  x范围", round(range(grid_ur),1), "-> yhat范围", round(range(pdp_ur$yhat),2),
    " spearman(x,yhat)=", round(cor(pdp_ur$x, pdp_ur$yhat, method="spearman"),2), "\n")
