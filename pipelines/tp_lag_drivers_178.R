#!/usr/bin/env Rscript
# =============================================================================
# tp_lag_drivers_178.R
#   对照 tp_lag_drivers_117.R: 用178站(只满足12号脚本标准(1)——tp∈[0,8]内至少1
#   个tp同时Δρ>0且p_surr<0.05——不再额外要求标准(2)optimal_tp>=0, 与
#   recovery_phi_drivers_178.R的178站定义完全一致)重做 optimal_tp ~ 协变量。
#   178站本身不是新的正式判据, 只是样本量敏感性对照, 不写回12号脚本、不替换
#   117站主线(见 pipelines/hcsif_buf1000/12_ccm_causal_confirmation.R 注释)。
#
#   optimal_tp 的定义在178站口径下做了一处必要调整: 12号脚本的 optimal_tp 是
#   在 tp∈[-8,8] 全范围里 argmax(rho), 117站集合本来就是"这个全局argmax落在
#   tp>=0"的子集——178站里多出的61站, 有相当一部分全局argmax其实落在负tp(时间
#   顺序不自洽, "未来VPD预测现在SIF", 见02_ccm.md), 直接沿用负值会让因变量变成
#   不可解释的"负滞后"。本脚本改为只在 tp∈[0,8] 范围内(标准(1)本来就检验这个
#   范围)对通过显著性的tp取argmax(rho), 保证178站的optimal_tp全部非负、可解释
#   为"效应滞后"，代价是失去117站那套"用负tp扫描交叉验证时间顺序自洽"的保证——
#   这是用178站换样本量必须接受的取舍, 不是免费的。
#
#   变量集: 用户要求"保留城市化率和投资至少1个"。178站里 invest 缺29站、
#   urban_rate 缺19站、cloud 缺5站(与117站口径类似, invest/urban_rate仍是缺失
#   大户)。complete-case n对比(15个高覆盖变量为基底):
#     +urban_rate(16变量)      n=159
#     +urban_rate+cloud(17变量) n=154
#     +invest(16变量)          n=149
#     +invest+urban_rate(17变量) n=142
#   选 15安全变量+urban_rate(16变量, n=159): urban_rate覆盖率明显好于invest,
#   且上一轮117站分析里invest的RF重要性本来就接近0/负值, 丢它信息损失小;
#   不额外搭cloud, 优先保样本量。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ranger); library(glmnet); library(car)
})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. 178站清单 + optimal_tp(限定在tp∈[0,8], 理由见上) --------------------
neg <- fread("data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv")
neg[, drho := rho - rho_min]
LC_BREAK_F <- "data_proc/output_hcsif_buf1000/landcover_break_stations.csv"
if (file.exists(LC_BREAK_F)) {
  lcb <- fread(LC_BREAK_F); excl <- lcb[has_break == TRUE]$stat_id
  neg <- neg[!meteo_stat %in% excl]
}
sig_range <- neg[tp >= 0 & tp <= 8]
sig_range[, pass1 := (drho > 0) & (p_surr < 0.05)]
stations178 <- unique(sig_range[pass1 == TRUE]$meteo_stat)
cat("178站候选清单 n =", length(stations178),
    "(117站是这178站里再要求全局optimal_tp(-8..8)>=0的子集)\n")

opt178 <- sig_range[meteo_stat %in% stations178 & pass1 == TRUE][
  order(meteo_stat, -rho)][, .SD[1], by = meteo_stat]
opt178 <- opt178[, .(meteo_stat, optimal_tp = tp, optimal_rho = rho, optimal_drho = drho,
                      optimal_p_surr = p_surr)]
cat("optimal_tp(0-8限定)分布:\n"); print(table(opt178$optimal_tp))

# ---- 2. 18个CORE协变量(复用22号脚本) -----------------------------------------
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt <- env$dt; PRED_FULL <- unlist(env$CORE, use.names = FALSE)

# ---- 3. 变量集: 15安全变量 + urban_rate(保留城镇化率/投资至少1个, 见文件头说明) --
PRED <- c(setdiff(PRED_FULL, c("invest", "urban_rate", "cloud")), "urban_rate")
cat("\n本次用的16个变量:", paste(PRED, collapse=", "), "\n")

d <- merge(opt178, dt, by.x = "meteo_stat", by.y = "stat_id")
d <- d[complete.cases(d[, PRED, with = FALSE])]
cat("协变量完整、可建模的站数(178站口径):", nrow(d), "\n")

f <- as.formula(paste("optimal_tp ~", paste(PRED, collapse = "+")))

# ---- 4. RF + OOB R2 + 置换重要性 ---------------------------------------------
set.seed(42)
rf <- ranger(f, data = d, num.trees = 2000, importance = "permutation", seed = 42)
imp <- data.table(var = names(rf$variable.importance), importance = rf$variable.importance)
setorder(imp, -importance)
cat("\n=== RF(16变量, 178站口径) OOB R2 =", round(rf$r.squared, 3), "===\n")
print(imp)

# ---- 5. OLS(方向/显著性参考) + VIF -------------------------------------------
ols <- lm(f, data = d)
cat("\n=== OLS(R2=", round(summary(ols)$r.squared, 3), ") ===\n")
co <- summary(ols)$coefficients
print(round(co[order(co[, 4]), c(1, 4)], 4))
vif_v <- sort(vif(ols), decreasing = TRUE)
cat("\n=== VIF ===\n"); print(round(vif_v, 2))

# ---- 6. 弹性网(alpha=.5, 10折CV): 标准化系数, 处理共线性+给诚实样本外R2 -------
X <- scale(as.matrix(d[, PRED, with = FALSE])); y <- d$optimal_tp
set.seed(42)
cvfit <- cv.glmnet(X, y, alpha = 0.5, nfolds = 10, standardize = FALSE)
r2_cv <- 1 - min(cvfit$cvm) / var(y)
enet <- as.data.table(as.matrix(coef(cvfit, s = "lambda.min")), keep.rownames = "var")
setnames(enet, 2, "enet_coef")
enet <- enet[var != "(Intercept)"]
cat("\n=== 弹性网(alpha=.5, 10折CV), lambda.min CV-R2 =", round(r2_cv, 3), "===\n")
print(enet[order(-abs(enet_coef))])

# ---- 7. 每个变量的偏依赖(固定其余变量在中位数, 扫描目标变量): 判断方向/形状 ---
Xraw <- as.matrix(d[, PRED, with = FALSE])
pdp_var <- function(v, ngrid = 30) {
  grid <- seq(quantile(d[[v]], .02), quantile(d[[v]], .98), length.out = ngrid)
  base <- as.list(apply(Xraw, 2, median))
  Xg <- do.call(rbind, lapply(grid, function(g) { b <- base; b[[v]] <- g; unlist(b) }))
  Xg <- as.data.table(Xg)[, PRED, with = FALSE]
  data.table(var = v, x = grid, yhat = predict(rf, Xg)$predictions)
}
pdp_all <- rbindlist(lapply(PRED, pdp_var))

# 形状归纳: Spearman(x,yhat)判断单调方向, 用"净变化量占总波动的比例"判断是否单调
shape_tab <- pdp_all[, {
  rho_s <- suppressWarnings(cor(x, yhat, method = "spearman"))
  rng   <- max(yhat) - min(yhat)
  net   <- yhat[.N] - yhat[1]
  monot <- if (rng < 1e-6) "flat" else if (abs(net) / rng > 0.7) {
    if (net > 0) "increasing" else "decreasing"
  } else "non-monotonic/threshold"
  .(spearman_x_yhat = round(rho_s, 2), pdp_range = round(rng, 2), shape = monot)
}, by = var]

sm <- Reduce(function(a, b) merge(a, b, by = "var", all = TRUE),
             list(imp, enet, shape_tab))
setorder(sm, -importance)
cat("\n=== 汇总: 各变量对optimal_tp的具体影响(RF重要性排序) ===\n")
print(sm)

fwrite(d,        file.path(OUT, "tp_lag_drivers_178_data.csv"))
fwrite(sm,       file.path(OUT, "tp_lag_drivers_178_summary.csv"))
fwrite(pdp_all,  file.path(OUT, "tp_lag_drivers_178_pdp.csv"))

# ---- 8. 图 --------------------------------------------------------------------
imp[, var := factor(var, levels = rev(var))]
p1 <- ggplot(imp, aes(var, importance)) +
  geom_col(fill = "#C1121F") + coord_flip() +
  labs(title = "RF variable importance for optimal_tp (178-station looser criterion, 16 vars)",
       subtitle = sprintf("n(modeled)=%d | RF OOB R2=%.3f | elastic-net CV-R2=%.3f | urban_rate kept, invest/cloud dropped for coverage",
                           nrow(d), rf$r.squared, r2_cv),
       x = NULL, y = "Permutation importance") +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size = 8, color = "grey35"))
ggsave(file.path(OUT, "tp_lag_drivers_178_rf_importance.png"), p1, width = 8, height = 6, dpi = 300)

p2 <- ggplot(pdp_all, aes(x, yhat)) +
  geom_line(color = "#B2182B", linewidth = 1) +
  facet_wrap(~var, scales = "free_x", ncol = 4) +
  labs(title = "Partial dependence of optimal_tp on all 16 predictors (178 stations)",
       subtitle = "Other predictors fixed at sample median",
       x = NULL, y = "Predicted optimal_tp") +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "tp_lag_drivers_178_pdp_all.png"), p2, width = 13, height = 10, dpi = 300)
cat("\n-> tp_lag_drivers_178_rf_importance.png, tp_lag_drivers_178_pdp_all.png\n")
