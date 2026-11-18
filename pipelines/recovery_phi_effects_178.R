#!/usr/bin/env Rscript
# =============================================================================
# recovery_phi_effects_178.R
#   诊断(recovery_phi_diagnostics_178.R)发现: 残差非正态+11个强影响点+
#   温度/降水/海拔共线(VIF 6.7-7.9)。处理方式: ①稳健回归(rlm)对照OLS,
#   检验强影响点有没有带偏结论; ②效应方向/形状以RF偏依赖图(PDP)为准,
#   不单独解读共线变量各自的OLS系数。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(MASS); library(ranger)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

d <- fread(file.path(OUT, "recovery_phi_drivers_178_data.csv"))
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
PRED <- unlist(env$CORE, use.names = FALSE)
f <- as.formula(paste("phi1_pure ~", paste(PRED, collapse = "+")))

# ---- 1. OLS vs 稳健回归(rlm), 检验强影响点有没有带偏方向/量级 ----------------
ols <- lm(f, data = d)
rob <- rlm(f, data = d, psi = psi.huber, maxit = 100)
cmp <- data.table(var = names(coef(ols)),
                   ols_coef = round(coef(ols), 4),
                   rlm_coef = round(coef(rob)[names(coef(ols))], 4))
cmp[, sign_agree := sign(ols_coef) == sign(rlm_coef)]
cat("=== OLS vs 稳健回归(rlm) 系数对比 ===\n"); print(cmp)
cat("\n符号一致的变量数:", sum(cmp$sign_agree[-1]), "/", nrow(cmp)-1, "\n")

# ---- 2. RF偏依赖图(PDP): 对重要性前6的变量, 手写PDP(固定该变量取值, 其余用 ---
# ---- 观测分布, 对全部观测取预测均值, 沿该变量的分位数网格扫描) --------------
rf <- ranger(f, data = d, num.trees = 2000, seed = 42)
imp <- fread(file.path(OUT, "recovery_phi_drivers_178_vs_117_importance.csv"))
top_vars <- head(imp[order(-importance_178)]$var, 6)
cat("\n偏依赖图变量(重要性前6):", paste(top_vars, collapse=", "), "\n")

pdp_one <- function(v) {
  grid <- quantile(d[[v]], probs = seq(0.05, 0.95, length.out = 20), na.rm = TRUE)
  pred <- sapply(grid, function(val) {
    dd <- copy(d); dd[[v]] <- val
    mean(predict(rf, data = dd)$predictions)
  })
  data.table(var = v, x = as.numeric(grid), yhat = pred)
}
pdp <- rbindlist(lapply(top_vars, pdp_one))
fwrite(pdp, file.path(OUT, "recovery_phi_pdp_178.csv"))

p <- ggplot(pdp, aes(x, yhat)) +
  geom_line(color = "#2166AC", linewidth = 1) + geom_point(color="#2166AC", size=1.5) +
  facet_wrap(~var, scales = "free_x", ncol = 3) +
  labs(title = "RF partial dependence: effect of top-6 variables on phi (resilience)",
       subtitle = "178-station criterion, n=137 | y = predicted phi (higher phi = slower return to equilibrium = lower resilience)",
       x = NULL, y = "Predicted phi") +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size=8.5, color="grey35"))
ggsave(file.path(OUT, "recovery_phi_pdp_178.png"), p, width = 11, height = 7, dpi = 300)
cat("\n-> recovery_phi_pdp_178.png\n")

# ---- 3. 每个PDP变量的净方向(首尾差)、供文字描述用 ----------------------------
trend <- pdp[, .(delta = yhat[.N] - yhat[1], range_x = paste0(round(x[1],2),"~",round(x[.N],2))), by = var]
cat("\n=== 各变量PDP首尾净变化(delta phi, 正=该变量增大时phi升高即恢复力下降) ===\n")
print(trend)
