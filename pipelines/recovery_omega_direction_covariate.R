#!/usr/bin/env Rscript
# =============================================================================
# recovery_omega_direction_covariate.R
#   logOmega(dirfix版, 已按方向正确构造) 全部57站(Promote+Inhibit)一个池子,
#   把 direction_2v 当作显式协变量(而不是像之前那样默认只看18个站点属性),
#   检验: 控制了18个协变量之后, 站点是Promote还是Inhibit本身还能不能独立
#   解释logOmega的差异? 用15var+VPD组合(phi分析里表现最好的一档)。
#   RF天然能捕捉direction x 其他协变量的交互, 不需要手写交互项; 岭回归额外
#   跑一版direction主效应, 对照"加不加direction, R2/系数变化多少"。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ranger); library(glmnet)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

ev <- fread(file.path(OUT, "recovery_liu2021_metrics_dirfix_117_events.csv"))
agg <- ev[, .(logOmega = as.numeric(mean(log(Omega), na.rm=TRUE)),
              direction_2v = direction_2v[1], n_events = .N), by = meteo_stat]
cat("站点数(57, Promote+Inhibit):\n"); print(agg[, .N, by = direction_2v])

env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18 <- env$dt
vc <- fread(file.path(OUT, "vpd_station_covariate.csv"))[, .(meteo_stat, vpd_mean_stn, vpd_cv_stn)]

d <- merge(agg, dt18, by.x = "meteo_stat", by.y = "stat_id")
d <- merge(d, vc, by = "meteo_stat")
VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))
VARSET <- c(VAR15, "vpd_mean_stn", "vpd_cv_stn")
d <- d[complete.cases(d[, VARSET, with = FALSE])]
d[, is_promote := as.integer(direction_2v == "Promote")]
cat("协变量完整可建模的站数:", nrow(d), "| Promote:", sum(d$is_promote), "Inhibit:", sum(1-d$is_promote), "\n")

# ---- RF: 不含direction vs 含direction ---------------------------------------
run_rf <- function(vars) {
  f <- as.formula(paste("logOmega ~", paste(vars, collapse = "+")))
  set.seed(42)
  rf <- ranger(f, data = d, num.trees = 2000, importance = "permutation", seed = 42)
  list(r2 = rf$r.squared, imp = data.table(var = names(rf$variable.importance), importance = rf$variable.importance))
}
rf_no  <- run_rf(VARSET)
rf_yes <- run_rf(c(VARSET, "is_promote"))
cat("\n=== RF OOB R2: 不含direction =", round(rf_no$r2,3), " | 含direction =", round(rf_yes$r2,3), "===\n")
imp_yes <- rf_yes$imp[order(-importance)]
cat("含direction模型的重要性排序(direction排第几):\n")
imp_yes[, rank := .I]
print(imp_yes)

# ---- 岭回归: direction主效应的系数和显著性方向(用bootstrap估计不确定性) ------
X_no  <- scale(as.matrix(d[, VARSET, with = FALSE]))
X_yes <- cbind(X_no, is_promote = d$is_promote)
y <- d$logOmega
run_ridge <- function(X) {
  set.seed(42); cv <- cv.glmnet(X, y, alpha = 0, nfolds = 10, standardize = FALSE)
  co <- as.numeric(coef(cv, s = "lambda.1se"))[-1]; names(co) <- colnames(X)
  pred <- predict(cv, X, s = "lambda.1se")
  list(r2 = 1 - sum((y-pred)^2)/sum((y-mean(y))^2), coef = co)
}
ridge_no  <- run_ridge(X_no)
ridge_yes <- run_ridge(X_yes)
cat("\n=== 岭回归(1se) R2: 不含direction =", round(ridge_no$r2,3), " | 含direction =", round(ridge_yes$r2,3), "===\n")
cat("含direction模型里 is_promote 的标准化系数:", round(ridge_yes$coef["is_promote"], 4),
    "(正=Promote站logOmega更高/更抗打; 0=被压缩掉,说明控制协变量后direction没有独立线性贡献)\n")

# bootstrap 500次, 给 is_promote 系数一个经验置信区间(岭回归没有解析p值)
set.seed(1)
boot_coef <- replicate(500, {
  idx <- sample(nrow(d), replace = TRUE)
  Xb <- X_yes[idx, ]; yb <- y[idx]
  if (length(unique(yb)) < 5 || sd(Xb[,"is_promote"]) == 0) return(NA_real_)
  tryCatch({
    cvb <- cv.glmnet(Xb, yb, alpha = 0, nfolds = 5, standardize = FALSE)
    as.numeric(coef(cvb, s = "lambda.1se"))[which(colnames(Xb)=="is_promote")+1]
  }, error = function(e) NA_real_)
})
ci <- quantile(boot_coef, c(.025,.5,.975), na.rm = TRUE)
cat("is_promote 系数 bootstrap 95%CI:", round(ci[1],4), "~", round(ci[3],4), "(中位数", round(ci[2],4), ")\n")

fwrite(data.table(var=names(ridge_yes$coef), coef_with_direction=ridge_yes$coef), file.path(OUT, "recovery_omega_direction_covariate_ridge.csv"))
fwrite(imp_yes, file.path(OUT, "recovery_omega_direction_covariate_rf_importance.csv"))

# ---- 图 -----------------------------------------------------------------------
imp_yes[, var := factor(var, levels = rev(var[order(-importance)]))]
imp_yes[, is_dir := var == "is_promote"]
p <- ggplot(imp_yes, aes(var, importance, fill = is_dir)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c(`TRUE`="#C1121F", `FALSE`="#2166AC"), guide = "none") +
  labs(title = "RF importance for logOmega(dirfix), direction included as an explicit covariate",
       subtitle = sprintf("57 stations, n modeled=%d | RF R2 without direction=%.3f, with direction=%.3f | red bar = is_promote",
                           nrow(d), rf_no$r2, rf_yes$r2),
       x = NULL, y = "Permutation importance") +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size=8, color="grey35"))
ggsave(file.path(OUT, "recovery_omega_direction_covariate_importance.png"), p, width = 8, height = 6.5, dpi = 300)
cat("\n-> recovery_omega_direction_covariate_importance.png\n")
