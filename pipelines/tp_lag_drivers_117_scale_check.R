#!/usr/bin/env Rscript
# =============================================================================
# tp_lag_drivers_117_scale_check.R
#   回应质疑: optimal_tp 是0-8的9档离散变量, 直接用OLS当连续变量处理是否合适?
#   用三种对"因变量是什么尺度"做不同假设的模型重跑同一套11变量/97站数据,
#   比较系数方向/显著性是否一致:
#     1) OLS(上次用的, 假设tp是区间/度量尺度, 单位间隔相等且有意义)
#     2) 累积比数有序logit(MASS::polr, 只假设tp的顺序信息, 不假设间隔相等)
#     3) 泊松/负二项GLM(假设tp是非负计数, 尊重"非负整数"这个事实,
#        用log link, 系数是"倍数效应"不是"单位效应")
#   三者对同一变量的方向若一致, 说明结论不依赖"因变量到底是什么尺度"这个假设,
#   比较稳；若不一致, 说明上次OLS的结论对尺度假设敏感, 需要谨慎对待。
#   注意: 这三个模型这次都用普通(模型自带)标准误, 没有像OLS那版一样加Conley
#   空间HAC校正——有序logit/计数模型的空间稳健标准误没有现成实现，工作量大，
#   这里的目的只是检查"方向和显著性是否对尺度假设敏感"，不是重新做完整推断，
#   这个简化必须在解读时交代清楚。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(MASS); library(car)
})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. 复用tp_lag_drivers_117_final.R同一套11变量/97站数据 ------------------
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt <- env$dt
tp <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")
d117full <- merge(tp, dt, by.x = "meteo_stat", by.y = "stat_id")

PRED <- c("imperv", "grass", "water", "ever_needle", "ever_broad", "deci_broad",
          "precip", "rsds_sd", "elev", "ntl", "invest")
d <- d117full[complete.cases(d117full[, PRED, with = FALSE])]
cat("n =", nrow(d), "\n")
cat("optimal_tp 分布(9档, 有的档观测数很少, 有序logit/计数模型对此敏感):\n")
print(table(d$optimal_tp))

f <- as.formula(paste("optimal_tp ~", paste(PRED, collapse = "+")))

# ---- 2. OLS(上次用的, 对照基准) ----------------------------------------------
ols <- lm(f, data = d)
co_ols <- summary(ols)$coefficients
cat("\n=== 1) OLS(假设tp是区间尺度) R2 =", round(summary(ols)$r.squared,3), "===\n")

# ---- 3. 有序logit(MASS::polr, 只用顺序信息不假设间隔相等) ---------------------
fo <- as.formula(paste("factor(optimal_tp, ordered=TRUE) ~", paste(PRED, collapse = "+")))
polr_fit <- polr(fo, data = d, Hess = TRUE, method = "logistic")
co_polr <- coef(summary(polr_fit))
p_polr  <- 2 * pnorm(-abs(co_polr[, "t value"]))
cat("\n=== 2) 有序logit(只假设顺序, 不假设间隔相等) ===\n")
cat("(注: polr没有伪R2的标准定义, 报似然比检验对比零模型)\n")
lr_polr <- polr_fit$deviance
cat("Residual deviance =", round(lr_polr,1), "| AIC =", round(AIC(polr_fit),1), "\n")

# ---- 4. 计数模型: 先检查是否过离散, 决定用泊松还是负二项 ----------------------
pois_fit <- glm(f, data = d, family = poisson())
disp <- sum(residuals(pois_fit, type="pearson")^2) / pois_fit$df.residual
cat("\n=== 3) 计数模型 ===\n")
cat("泊松离散度(应≈1, 越大越过离散) =", round(disp,2), "\n")
if (disp > 1.5) {
  cat("过离散明显, 改用负二项(NB)回归\n")
  count_fit <- glm.nb(f, data = d)
  count_name <- "负二项(NB)"
} else {
  count_fit <- pois_fit
  count_name <- "泊松"
}
co_count <- coef(summary(count_fit))
cat(count_name, "GLM(log link, 系数是log-count的变化) AIC =", round(AIC(count_fit),1), "\n")

# ---- 5. 三模型方向/显著性汇总对比 --------------------------------------------
dir_sign <- function(x) ifelse(x > 0, "+", ifelse(x < 0, "-", "0"))
sm <- data.table(
  var = PRED,
  ols_sign = dir_sign(co_ols[PRED, 1]), ols_p = round(co_ols[PRED, 4], 3),
  polr_sign = dir_sign(co_polr[PRED, 1]), polr_p = round(p_polr[PRED], 3),
  count_sign = dir_sign(co_count[PRED, 1]), count_p = round(co_count[PRED, 4], 3)
)
sm[, agree_all3 := (ols_sign == polr_sign) & (polr_sign == count_sign)]
sm[, any_sig := (ols_p < .1) | (polr_p < .1) | (count_p < .1)]
setorder(sm, -agree_all3, ols_p)
cat("\n=== 三种尺度假设下的方向/显著性对比(", count_name, "作计数模型) ===\n")
print(sm)
fwrite(sm, file.path(OUT, "tp_lag_drivers_117_scale_check_summary.csv"))

cat("\n三模型方向完全一致的变量数:", sum(sm$agree_all3), "/", nrow(sm), "\n")
cat("方向不一致的变量:", paste(sm[agree_all3 == FALSE]$var, collapse=", "), "\n")

# ---- 6. polr结果异常"显著"(多个p=0.000)可疑,诊断是否过参数化/准分离 ---------
cat("\n=== polr诊断: 系数量级+收敛信息(排查是否过参数化/准分离导致p值虚低) ===\n")
cat("polr参数个数(11个系数+8个切点) =", length(coef(polr_fit)) + length(polr_fit$zeta),
    " vs 观测数 n=", nrow(d), " | 每类观测数最少 =", min(table(d$optimal_tp)), "\n")
cat("收敛代码 convergence =", polr_fit$convergence, "(0=正常收敛)\n")
print(round(co_polr, 3))
