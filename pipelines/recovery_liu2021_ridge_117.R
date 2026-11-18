#!/usr/bin/env Rscript
# =============================================================================
# recovery_liu2021_ridge_117.R
#   OLS 在 constraint_time 上因 18 个协变量高度共线(tavg/elev/precip VIF>8)
#   给不出可信的系数方向, 改用岭回归(glmnet, alpha=0)——岭回归不会像OLS那样
#   把共线变量的系数任意分摊甚至符号反转, 而是把相关变量的系数往彼此靠拢、
#   整体收缩, 更适合在共线情况下读"这一簇变量整体的方向和相对大小"。
#   四个DV都跑一遍(constraint_time是重点, 其余三个作为对照)。
#   lambda 用10折CV选(lambda.min=预测最优, lambda.1se=更保守/更抗共线，
#   高共线场景下1se通常更值得信)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(glmnet); library(ggplot2)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

d <- fread(file.path(OUT, "recovery_liu2021_drivers_117_data.csv"))
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
PRED <- unlist(env$CORE, use.names = FALSE)

Xraw <- as.matrix(d[, PRED, with = FALSE])
X <- scale(Xraw)  # 标准化, 系数=1SD变化对应的效应, 和上次OLS标准化系数可比

DVS <- c("constraint_time", "recovery_time", "logOmega", "logDelta")
run_ridge <- function(dv) {
  y <- d[[dv]]
  set.seed(42)
  cv <- cv.glmnet(X, y, alpha = 0, nfolds = 10, standardize = FALSE)
  co_min <- as.numeric(coef(cv, s = "lambda.min"))[-1]
  co_1se <- as.numeric(coef(cv, s = "lambda.1se"))[-1]
  pred <- predict(cv, X, s = "lambda.1se")
  r2_1se <- 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
  data.table(dv = dv, var = PRED, ridge_min = co_min, ridge_1se = co_1se, r2_1se = r2_1se)
}
res <- rbindlist(lapply(DVS, run_ridge))
fwrite(res, file.path(OUT, "recovery_liu2021_ridge_117.csv"))

cat("=== 岭回归(lambda.1se) R2 ===\n")
print(unique(res[, .(dv, r2_1se = round(r2_1se, 3))]))

cat("\n=== constraint_time 岭回归系数(标准化, 按|1se系数|排序) ===\n")
print(res[dv == "constraint_time"][order(-abs(ridge_1se)),
          .(var, ridge_min = round(ridge_min, 3), ridge_1se = round(ridge_1se, 3))])

cat("\n=== logOmega 岭回归系数(对照OLS, 按|1se系数|排序) ===\n")
print(res[dv == "logOmega"][order(-abs(ridge_1se)),
          .(var, ridge_min = round(ridge_min, 3), ridge_1se = round(ridge_1se, 3))])

# ---- 图: constraint_time 岭回归系数条形图, 标出方向 ------------------------
pc <- res[dv == "constraint_time"][order(ridge_1se)]
pc[, var := factor(var, levels = var)]
p1 <- ggplot(pc, aes(var, ridge_1se, fill = ridge_1se > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#C1121F"), guide = "none") +
  labs(title = "Ridge regression coefficients: constraint_time ~ 18 covariates (lambda.1se)",
       subtitle = sprintf("Standardized X | R2=%.2f | resolves the VIF>8 collinearity that made OLS coefficients unreliable",
                           pc$r2_1se[1]),
       x = NULL, y = "Ridge coefficient (standardized)") +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size=8.5, color="grey35"))
ggsave(file.path(OUT, "recovery_liu2021_ridge_constraint_time.png"), p1, width = 8, height = 6.5, dpi = 300)

pd <- res[dv == "logOmega"][order(ridge_1se)]
pd[, var := factor(var, levels = var)]
p2 <- ggplot(pd, aes(var, ridge_1se, fill = ridge_1se > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#C1121F"), guide = "none") +
  labs(title = "Ridge regression coefficients: logOmega ~ 18 covariates (lambda.1se)",
       subtitle = sprintf("Standardized X | R2=%.2f | compare against the OLS-significant variable list", pd$r2_1se[1]),
       x = NULL, y = "Ridge coefficient (standardized)") +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size=8.5, color="grey35"))
ggsave(file.path(OUT, "recovery_liu2021_ridge_logOmega.png"), p2, width = 8, height = 6.5, dpi = 300)
cat("\n-> recovery_liu2021_ridge_constraint_time.png\n-> recovery_liu2021_ridge_logOmega.png\n")
