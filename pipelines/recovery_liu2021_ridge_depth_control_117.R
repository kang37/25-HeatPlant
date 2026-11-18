#!/usr/bin/env Rscript
# =============================================================================
# recovery_liu2021_ridge_depth_control_117.R
#   给 constraint_time / recovery_time 的岭回归加一个"跌多深"的控制变量
#   (logOmega, 抵抗力), 检验18个协变量对"时长"的效应有多少是通过"深度"
#   间接起作用的(机械关系, invest 已验证过一次), 多少是独立于深度的真实效应。
#   recovery_time 额外加一版 logDelta(恢复力)当控制, 因为"回升多久"更可能被
#   "回升到多高才算数"这件事本身牵连, 逻辑上和 constraint_time-logOmega 那对
#   是同一类关系。
#   做法: 每个DV各跑两次岭回归(不控制深度 / 控制深度), 比较18个协变量的系数
#   谁掉得多(=效应主要是深度中介的, 类似invest) 谁基本不变(=独立于深度的
#   真实时长驱动因素)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(glmnet); library(ggplot2)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

d <- fread(file.path(OUT, "recovery_liu2021_drivers_117_data.csv"))
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
PRED <- unlist(env$CORE, use.names = FALSE)

fit_ridge <- function(dv, extra_pred = NULL) {
  allpred <- c(PRED, extra_pred)
  X <- scale(as.matrix(d[, allpred, with = FALSE]))
  y <- d[[dv]]
  set.seed(42)
  cv <- cv.glmnet(X, y, alpha = 0, nfolds = 10, standardize = FALSE)
  co <- as.numeric(coef(cv, s = "lambda.1se"))[-1]
  pred <- predict(cv, X, s = "lambda.1se")
  r2 <- 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
  data.table(var = allpred, coef = co, r2 = r2)
}

# ---- constraint_time: 不控深度 vs 控 logOmega ------------------------------
ct_base  <- fit_ridge("constraint_time")[, model := "no_depth_control"]
ct_ctrl  <- fit_ridge("constraint_time", "logOmega")[, model := "control_logOmega"]
ct <- rbind(ct_base, ct_ctrl)
ct_wide <- dcast(ct, var ~ model, value.var = "coef")
ct_wide[, delta := control_logOmega - no_depth_control]
ct_wide[, pct_survive := round(100 * abs(control_logOmega) / pmax(abs(no_depth_control), 1e-6), 0)]
ct_wide <- ct_wide[order(-abs(no_depth_control))]
cat("=== constraint_time: 加 logOmega 控制前后系数对比 ===\n")
cat("R2 无控制=", round(ct_base$r2[1],3), " | 有控制=", round(ct_ctrl$r2[1],3), "\n")
print(ct_wide[var != "logOmega"][, .(var, no_depth_control = round(no_depth_control,3),
              control_logOmega = round(control_logOmega,3), pct_survive)])
cat("logOmega自己的系数(控制模型里):", round(ct_ctrl[var=="logOmega", coef],3), "\n")

# ---- recovery_time: 不控深度 vs 控 logOmega vs 控 logDelta -----------------
rt_base <- fit_ridge("recovery_time")[, model := "no_depth_control"]
rt_om   <- fit_ridge("recovery_time", "logOmega")[, model := "control_logOmega"]
rt_de   <- fit_ridge("recovery_time", "logDelta")[, model := "control_logDelta"]
rt <- rbind(rt_base, rt_om, rt_de)
rt_wide <- dcast(rt, var ~ model, value.var = "coef")
rt_wide <- rt_wide[order(-abs(no_depth_control))]
cat("\n=== recovery_time: 加 logOmega/logDelta 控制前后系数对比 ===\n")
cat("R2 无控制=", round(rt_base$r2[1],3), " | 控logOmega=", round(rt_om$r2[1],3),
    " | 控logDelta=", round(rt_de$r2[1],3), "\n")
print(rt_wide[!var %in% c("logOmega","logDelta")][, .(var, no_depth_control = round(no_depth_control,3),
              control_logOmega = round(control_logOmega,3), control_logDelta = round(control_logDelta,3))])

fwrite(ct_wide, file.path(OUT, "recovery_liu2021_depth_control_constraint_time.csv"))
fwrite(rt_wide, file.path(OUT, "recovery_liu2021_depth_control_recovery_time.csv"))

# ---- 图: constraint_time 系数, 控制前 vs 控制后并排 -------------------------
plot_dt <- melt(ct_wide[var != "logOmega"], id.vars = c("var"),
                 measure.vars = c("no_depth_control", "control_logOmega"),
                 variable.name = "model", value.name = "coef")
ord <- ct_wide[var != "logOmega"][order(no_depth_control)]$var
plot_dt[, var := factor(var, levels = ord)]
p <- ggplot(plot_dt, aes(var, coef, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) + coord_flip() +
  scale_fill_manual(values = c(no_depth_control = "#C1121F", control_logOmega = "#2166AC"),
                     labels = c("no_depth_control"="not controlling for depth (logOmega)",
                                "control_logOmega"="controlling for depth (logOmega)")) +
  labs(title = "constraint_time ridge coefficients: before vs after controlling for drop depth (logOmega)",
       subtitle = "Bars that shrink a lot toward zero (blue vs red) = effect mostly mediated through depth (e.g. invest); bars that barely move = a genuine independent timing driver",
       x = NULL, y = "Ridge coefficient (standardized)", fill = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "top", plot.subtitle = element_text(size=8, color="grey35"))
ggsave(file.path(OUT, "recovery_liu2021_depth_control_constraint_time.png"), p, width = 9, height = 7, dpi = 300)
cat("\n-> recovery_liu2021_depth_control_constraint_time.png\n")
