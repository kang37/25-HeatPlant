#!/usr/bin/env Rscript
# =============================================================================
# tp_lag_var_screening_117.R
#   在对15个CORE变量做效应解读之前，先做变量筛选诊断，回答"要不要去掉一些变量"
#   三个子问题：①理论上哪些变量可以去掉(域知识判断，不依赖本次拟合结果)；
#   ②数据上哪些变量该去掉(近零方差/高度共线/持续低RF重要性，依赖本次数据)；
#   ③有没有正式方法能在解读效应之前先筛出"靠得住"的变量——用 bootstrap 稳定性
#   选择(stability selection, Meinshausen & Bühlmann 2010 的思路简化版):
#   对117站重复自助重采样(500次)，每次用弹性网(10折CV选lambda)看哪些变量系数
#   非零，统计"入选频率"——频率高的变量才值得单独解读效应，频率低的说明那个
#   变量的"重要性"很大程度是抽样噪声，换一批站就选不进模型了。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ranger); library(glmnet); library(car)
})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 0. 复用15变量/117站数据(与tp_lag_drivers_117.R第二轮同一套) -------------
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt <- env$dt; PRED_FULL <- unlist(env$CORE, use.names = FALSE)
tp <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")
d117full <- merge(tp, dt, by.x = "meteo_stat", by.y = "stat_id")
PRED15 <- setdiff(PRED_FULL, c("invest", "urban_rate", "cloud"))
d <- d117full[complete.cases(d117full[, PRED15, with = FALSE])]
cat("n =", nrow(d), "\n")

# ---- 1. 数据诊断①: 近零方差/近常数(成分数据里常见的稀疏问题) -----------------
diagn <- data.table(var = PRED15,
  pct_near0 = sapply(PRED15, function(v) round(100 * mean(d[[v]] < 0.01 * max(d[[v]], na.rm = TRUE)), 1)),
  n_unique  = sapply(PRED15, function(v) uniqueN(round(d[[v]], 4))),
  cv        = sapply(PRED15, function(v) round(sd(d[[v]]) / (abs(mean(d[[v]])) + 1e-9), 2)))
setorder(diagn, -pct_near0)
cat("\n=== 近零方差诊断(按近零比例排序, >70%视为高度稀疏) ===\n"); print(diagn)

# ---- 2. 数据诊断②: 相关矩阵 + VIF(共线性) ------------------------------------
cm <- cor(d[, PRED15, with = FALSE])
hi <- which(abs(cm) > 0.6 & abs(cm) < 1, arr.ind = TRUE); hi <- hi[hi[,1] < hi[,2], , drop = FALSE]
cat("\n=== |r|>0.6 的变量对 ===\n")
for (i in seq_len(nrow(hi))) cat(sprintf("  %-12s %-12s r=%.2f\n",
    rownames(cm)[hi[i,1]], colnames(cm)[hi[i,2]], cm[hi[i,1], hi[i,2]]))

f15 <- as.formula(paste("optimal_tp ~", paste(PRED15, collapse = "+")))
vif15 <- sort(vif(lm(f15, d)), decreasing = TRUE)
cat("\n=== VIF ===\n"); print(round(vif15, 2))

# ---- 3. 正式筛选方法: bootstrap 稳定性选择(弹性网, 500次自助重采样) -----------
set.seed(42)
NBOOT <- 500
Xall <- scale(as.matrix(d[, PRED15, with = FALSE]))
yall <- d$optimal_tp
sel_mat <- matrix(0L, NBOOT, length(PRED15), dimnames = list(NULL, PRED15))
for (b in seq_len(NBOOT)) {
  idx <- sample(nrow(d), nrow(d), replace = TRUE)
  cvb <- tryCatch(cv.glmnet(Xall[idx, ], yall[idx], alpha = 0.5, nfolds = 5, standardize = FALSE),
                   error = function(e) NULL)
  if (is.null(cvb)) next
  co <- as.matrix(coef(cvb, s = "lambda.min"))[-1, 1]
  sel_mat[b, ] <- as.integer(co != 0)
}
stab <- data.table(var = PRED15, sel_freq = round(colMeans(sel_mat), 3))
setorder(stab, -sel_freq)
cat("\n=== 稳定性选择(500次bootstrap, 弹性网非零系数入选频率) ===\n"); print(stab)
fwrite(stab, file.path(OUT, "tp_lag_var_stability_selection_117.csv"))

p <- ggplot(stab, aes(reorder(var, sel_freq), sel_freq)) +
  geom_col(fill = "#2A9D8F") + geom_hline(yintercept = 0.5, linetype = 2, color = "grey40") +
  coord_flip() +
  labs(title = "Bootstrap stability selection: how often each variable survives elastic-net (500 resamples)",
       subtitle = "optimal_tp ~ 15 CORE vars, 117 CCM-confirmed stations | dashed line = 50% selection frequency",
       x = NULL, y = "Selection frequency") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "tp_lag_var_stability_selection_117.png"), p, width = 8, height = 6, dpi = 300)
cat("\n-> tp_lag_var_stability_selection_117.png\n")

# ---- 4. 用稳定性选择推荐的精简变量集, 重跑RF+弹性网, 看是否比15变量版本更好 -----
keep <- stab[sel_freq >= 0.3]$var         # 阈值放宽到30%(500次里至少150次入选), 15个变量本来就不多
cat("\n稳定性选择(>=30%入选频率)推荐保留:", paste(keep, collapse=", "), "(", length(keep), "个)\n")
if (length(keep) >= 2) {
  fk <- as.formula(paste("optimal_tp ~", paste(keep, collapse = "+")))
  rfk <- ranger(fk, data = d, num.trees = 2000, importance = "permutation", seed = 42)
  cat("精简变量集 RF OOB R2 =", round(rfk$r.squared, 3), "(对照15变量版 R2=0.13)\n")
  impk <- data.table(var = names(rfk$variable.importance), importance = rfk$variable.importance)
  setorder(impk, -importance); print(impk)
  Xk <- scale(as.matrix(d[, keep, with = FALSE]))
  cvk <- cv.glmnet(Xk, yall, alpha = 0.5, nfolds = 10, standardize = FALSE)
  r2k <- 1 - min(cvk$cvm) / var(yall)
  cat("精简变量集 弹性网 CV-R2 =", round(r2k, 3), "(对照15变量版 CV-R2=0.03)\n")
  fwrite(impk, file.path(OUT, "tp_lag_drivers_117_stability_reduced_rf_importance.csv"))
}

# =============================================================================
# 5. 结论: bootstrap稳定性选择本身在这批数据上被"近常数变量"的数值病态干扰
#    (mixedleaf/deci_needle 91-97%近零，标准化后方差被放大，入选频率虚高但RF
#    重要性接近0甚至为负，见上面结果的矛盾)——不能直接照抄频率排序当结论。
#    改为"理论+数据诊断双重把关"手动构造精简集: 先按近零方差(>90%视为近常数)
#    丢 mixedleaf/deci_needle; 再按相关矩阵/VIF, 气候簇(tavg/rh/rsds_mean/
#    precip/elev)只保留信息量最大、历次分析最稳定的 precip+elev+rsds_sd,
#    丢冗余的 tavg/rh/rsds_mean。得到10变量精简集, 对照原15变量看R2/系数是否
#    更稳定。
# =============================================================================
cat("\n\n########## 理论+数据双重把关的精简集(10变量) ##########\n")
PRED10 <- c("imperv","grass","water","ever_needle","ever_broad","deci_broad",
            "precip","rsds_sd","elev","ntl")
cat("保留:", paste(PRED10, collapse=", "), "\n")
cat("丢弃: mixedleaf/deci_needle(近常数, >90%近零) | tavg/rh/rsds_mean(与precip/elev/rsds_sd高度共线, VIF 4.7-8.2)\n")

f10 <- as.formula(paste("optimal_tp ~", paste(PRED10, collapse = "+")))
rf10 <- ranger(f10, data = d, num.trees = 2000, importance = "permutation", seed = 42)
imp10 <- data.table(var = names(rf10$variable.importance), importance = rf10$variable.importance)
setorder(imp10, -importance)
cat("RF OOB R2 =", round(rf10$r.squared, 3), "(对照15变量版=0.13)\n"); print(imp10)

X10 <- scale(as.matrix(d[, PRED10, with = FALSE]))
cv10 <- cv.glmnet(X10, yall, alpha = 0.5, nfolds = 10, standardize = FALSE)
r2_10 <- 1 - min(cv10$cvm) / var(yall)
enet10 <- as.data.table(as.matrix(coef(cv10, s = "lambda.min")), keep.rownames = "var")
setnames(enet10, 2, "enet_coef"); enet10 <- enet10[var != "(Intercept)"]
cat("\n弹性网 CV-R2 =", round(r2_10, 3), "(对照15变量版=0.03)\n")
print(enet10[order(-abs(enet_coef))])

vif10 <- sort(vif(lm(f10, d)), decreasing = TRUE)
cat("\nVIF(10变量) 最大值 =", round(max(vif10),2), "(对照15变量版最大值=8.20)\n"); print(round(vif10,2))

fwrite(imp10,   file.path(OUT, "tp_lag_drivers_117_10var_rf_importance.csv"))
fwrite(enet10,  file.path(OUT, "tp_lag_drivers_117_10var_enet_coef.csv"))
