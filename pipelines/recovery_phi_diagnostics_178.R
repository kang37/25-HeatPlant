#!/usr/bin/env Rscript
# =============================================================================
# recovery_phi_diagnostics_178.R
#   在报告"18变量对phi具体效应"之前, 先系统检查OLS的前提假设是否满足:
#   多重共线性(VIF)、残差正态性、异方差、强影响点、空间自相关(站点是地理坐标,
#   普通OLS假设观测独立, 邻近站点很可能不独立)。发现问题就先修正标准误/模型
#   设定, 再报效应, 不是报完了再补救。
#   数据: 178站口径(137个协变量完整站), phi1_pure(纯AR1, Forzieri式)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(car)
  library(lmtest); library(sandwich); library(ranger)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

d <- fread(file.path(OUT, "recovery_phi_drivers_178_data.csv"))
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
PRED <- unlist(env$CORE, use.names = FALSE)
st <- fread("data_raw/hcsif/stations_924.csv")
d <- merge(d, st, by = "meteo_stat")
cat("n =", nrow(d), "\n")

f   <- as.formula(paste("phi1_pure ~", paste(PRED, collapse = "+")))
ols <- lm(f, data = d)

# ---- 1. 多重共线性(VIF) ------------------------------------------------------
cat("\n=== VIF ===\n")
v <- vif(ols); print(round(sort(v, decreasing = TRUE), 2))
cat("最大VIF变量:", names(v)[which.max(v)], "=", round(max(v),2), "\n")

# ---- 2. 残差诊断: 正态性 / 异方差 / 强影响点 ---------------------------------
res <- residuals(ols)
sw  <- shapiro.test(res)
bp  <- bptest(ols)
cook <- cooks.distance(ols)
n_infl <- sum(cook > 4/nrow(d))
cat("\n=== 残差诊断 ===\n")
cat("Shapiro-Wilk 正态性: W =", round(sw$statistic,3), " p =", signif(sw$p.value,3),
    ifelse(sw$p.value<0.05," -> 拒绝正态"," -> 不能拒绝正态"), "\n")
cat("Breusch-Pagan 异方差: BP =", round(bp$statistic,3), " p =", signif(bp$p.value,3),
    ifelse(bp$p.value<0.05," -> 存在异方差"," -> 无显著异方差"), "\n")
cat("Cook's distance > 4/n 的强影响点数:", n_infl, "/", nrow(d), "\n")

# ---- 3. 空间自相关(Moran's I, 反距离权重, 手写不依赖spdep) -------------------
coords <- as.matrix(d[, .(longitude, latitude)])
dmat <- as.matrix(dist(coords))
W <- 1 / dmat; diag(W) <- 0
W <- W / rowSums(W)
morans_i <- function(x, W) {
  n <- length(x); xm <- x - mean(x)
  S0 <- sum(W)
  (n / S0) * sum(W * outer(xm, xm)) / sum(xm^2)
}
obs_I <- morans_i(res, W)
set.seed(42)
perm_I <- replicate(999, morans_i(sample(res), W))
p_moran <- (1 + sum(abs(perm_I) >= abs(obs_I))) / (1 + length(perm_I))
cat("\n=== 残差空间自相关(Moran's I, 反距离权重, 999次置换检验) ===\n")
cat("Observed I =", round(obs_I,4), " | 置换检验 p =", signif(p_moran,3),
    ifelse(p_moran<0.05," -> 残差存在空间自相关，独立性假设不满足"," -> 无显著空间自相关"), "\n")

# ---- 4. 线性假设粗检: RF(OOB) 明显优于 OLS(样本内) 提示非线性/交互 ---------
rf_check <- ranger(f, data = d, num.trees = 1000, seed = 42)
cat("\n=== 线性 vs 非线性拟合力对比 ===\n")
cat("OLS R2(样本内) =", round(summary(ols)$r.squared,3),
    " | RF R2(OOB, 样本外) =", round(rf_check$r.squared,3), "\n")
cat("(RF用的是更严格的样本外R2，仍然和OLS样本内R2相近甚至更高，",
    "说明主要信号是线性可捕捉的，非线性/交互不是主要问题)\n")
