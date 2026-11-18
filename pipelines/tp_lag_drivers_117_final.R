#!/usr/bin/env Rscript
# =============================================================================
# tp_lag_drivers_117_final.R
#   在 tp_lag_var_screening_117.R 诊断基础上定稿的变量集: 10个理论+数据双重
#   把关精简变量(丢mixedleaf/deci_needle[近常数]、tavg/rh/rsds_mean[与precip/
#   elev/rsds_sd高度共线, VIF 4.7-8.2]) + 用户指定保留的社会经济变量 invest
#   (而不是上一版用的urban_rate)，共11个变量。
#
#   核心产出: 每个自变量对 optimal_tp 的效应方向 + 显著性。方向/显著性用
#   Conley 空间 HAC 标准误的 OLS(不是普通OLS)——中国气象站在空间上高度自相关
#   (邻近站点气候/地表覆盖相似)，普通OLS标准误会把这种相关性误当独立信息，
#   显著性虚高。conley_lm() 沿用已删除的45_class_logit.R(git历史里的
#   conley_glm)同一套Bartlett核+200km截断设计，改成OLS版本: bread换成
#   (X'X)^-1(不是vcov(fit)，vcov(lm)里含sigma^2会重复计入)，其余不变。
#   RF+偏依赖作非线性交叉检验(方向是否和OLS线性系数一致)，弹性网作"是否经得起
#   正则化"的稳健性对照——三者一起看才是完整结论，任何一个单独都可能误导。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ranger); library(glmnet); library(car)
})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. 数据: 11变量(10个精简+invest) x 117站 --------------------------------
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt <- env$dt
tp <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")
st <- fread("data_raw/hcsif/stations_924.csv")

d117full <- merge(tp, dt, by.x = "meteo_stat", by.y = "stat_id")
d117full <- merge(d117full, st, by = "meteo_stat")

PRED <- c("imperv", "grass", "water", "ever_needle", "ever_broad", "deci_broad",
          "precip", "rsds_sd", "elev", "ntl", "invest")
cat("最终11个变量:", paste(PRED, collapse=", "), "\n")

d <- d117full[complete.cases(d117full[, PRED, with = FALSE])]
cat("协变量完整、可建模的站数:", nrow(d), "(原117站, invest缺失导致部分流失)\n")

f <- as.formula(paste("optimal_tp ~", paste(PRED, collapse = "+")))

# ---- 2. VIF(确认精简后共线性确实降下来了) -------------------------------------
ols_raw <- lm(f, data = d)
vif_v <- sort(vif(ols_raw), decreasing = TRUE)
cat("\n=== VIF(11变量) ===\n"); print(round(vif_v, 2))

# ---- 3. Conley 空间 HAC 标准误的 OLS: 效应方向 + 显著性(核心产出) ------------
conley_lm <- function(fit, lon, lat, cut = 200) {
  X <- model.matrix(fit); u <- residuals(fit)
  S <- X * u
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  dy <- outer(lat, lat, "-") * 111.32
  K  <- pmax(0, 1 - sqrt(dx^2 + dy^2) / cut)
  dim(K) <- c(length(lon), length(lon))            # pmax()会丢dim属性,需手动补回
  bread <- solve(crossprod(X))                      # (X'X)^-1, 不是vcov(fit)(含sigma^2会重复计入)
  V <- bread %*% (crossprod(S, K %*% S)) %*% bread
  n <- nrow(X); k <- ncol(X)
  V * (n / (n - k))
}
Vc <- conley_lm(ols_raw, d$longitude, d$latitude, cut = 200)
b  <- coef(ols_raw); se_c <- sqrt(diag(Vc))
tval <- b / se_c; pval <- 2 * pt(-abs(tval), df = nrow(d) - length(b))
res <- data.table(var = names(b), coef = round(b, 4), se_conley = round(se_c, 4),
                   t = round(tval, 2), p_conley = round(pval, 4))[var != "(Intercept)"]
res[, sig := fifelse(p_conley < .01, "***", fifelse(p_conley < .05, "**",
              fifelse(p_conley < .1, "*", "")))]
res[, direction := fifelse(coef > 0, "positive(+)", "negative(-)")]
setorder(res, p_conley)
cat("\n=== OLS + Conley空间HAC(200km, Bartlett核) 标准误: 效应方向+显著性 ===\n")
cat("普通OLS R2 =", round(summary(ols_raw)$r.squared, 3), "\n")
print(res[, .(var, coef, direction, se_conley, t, p_conley, sig)])
fwrite(res, file.path(OUT, "tp_lag_drivers_117_final_conley_ols.csv"))

# ---- 4. RF + 偏依赖: 非线性交叉检验(方向是否和OLS线性系数一致) ---------------
set.seed(42)
rf <- ranger(f, data = d, num.trees = 2000, importance = "permutation", seed = 42)
imp <- data.table(var = names(rf$variable.importance), importance = rf$variable.importance)
setorder(imp, -importance)
cat("\n=== RF OOB R2 =", round(rf$r.squared, 3), "| 重要性排序 ===\n"); print(imp)

Xraw <- as.matrix(d[, PRED, with = FALSE])
pdp_var <- function(v, ngrid = 30) {
  grid <- seq(quantile(d[[v]], .02), quantile(d[[v]], .98), length.out = ngrid)
  base <- as.list(apply(Xraw, 2, median))
  Xg <- do.call(rbind, lapply(grid, function(g) { b2 <- base; b2[[v]] <- g; unlist(b2) }))
  Xg <- as.data.table(Xg)[, PRED, with = FALSE]
  data.table(var = v, x = grid, yhat = predict(rf, Xg)$predictions)
}
pdp_all <- rbindlist(lapply(PRED, pdp_var))
pdp_dir <- pdp_all[, .(pdp_spearman = round(suppressWarnings(cor(x, yhat, method = "spearman")), 2),
                        pdp_range = round(max(yhat) - min(yhat), 2)), by = var]

# ---- 5. 弹性网(10折CV): 是否经得起正则化 -------------------------------------
X <- scale(Xraw); y <- d$optimal_tp
set.seed(42)
cvfit <- cv.glmnet(X, y, alpha = 0.5, nfolds = 10, standardize = FALSE)
r2_cv <- 1 - min(cvfit$cvm) / var(y)
enet <- as.data.table(as.matrix(coef(cvfit, s = "lambda.min")), keep.rownames = "var")
setnames(enet, 2, "enet_coef"); enet <- enet[var != "(Intercept)"]
cat("\n弹性网 CV-R2 =", round(r2_cv, 3), "\n")

# ---- 6. 汇总表: 方向+显著性+非线性交叉检验+正则化稳健性, 一张表看完 ------------
sm <- Reduce(function(a, b2) merge(a, b2, by = "var", all = TRUE),
             list(res[, .(var, coef, direction, p_conley, sig)], imp, pdp_dir, enet))
setorder(sm, p_conley)
sm[, enet_survives := fifelse(enet_coef != 0, "yes", "no")]
cat("\n=== 最终汇总: 11个变量对optimal_tp的方向/显著性/稳健性 ===\n")
print(sm[, .(var, coef, direction, p_conley, sig, rf_importance = round(importance,3),
             pdp_spearman, enet_survives)])
fwrite(sm, file.path(OUT, "tp_lag_drivers_117_final_summary.csv"))
fwrite(pdp_all, file.path(OUT, "tp_lag_drivers_117_final_pdp.csv"))

# ---- 7. 图 --------------------------------------------------------------------
res[, var := factor(var, levels = rev(res$var))]
p1 <- ggplot(res, aes(var, coef, fill = sig != "")) +
  geom_col() + geom_errorbar(aes(ymin = coef - 1.96*se_conley, ymax = coef + 1.96*se_conley), width = .3) +
  coord_flip() + scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey70"), guide = "none") +
  labs(title = "OLS coefficients on optimal_tp with Conley spatial-HAC 95% CI (11 vars, 117 stations)",
       subtitle = sprintf("Red = significant at p<.1 (Conley SE, 200km Bartlett kernel) | OLS R2=%.2f",
                           summary(ols_raw)$r.squared),
       x = NULL, y = "OLS coefficient") +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size = 8.5, color = "grey35"))
ggsave(file.path(OUT, "tp_lag_drivers_117_final_conley_coef.png"), p1, width = 8, height = 6, dpi = 300)

p2 <- ggplot(pdp_all, aes(x, yhat)) +
  geom_line(color = "#B2182B", linewidth = 1) +
  facet_wrap(~var, scales = "free_x", ncol = 4) +
  labs(title = "Partial dependence of optimal_tp on all 11 final predictors (117 stations)",
       subtitle = "Other predictors fixed at sample median",
       x = NULL, y = "Predicted optimal_tp") +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "tp_lag_drivers_117_final_pdp.png"), p2, width = 13, height = 8, dpi = 300)
cat("\n-> tp_lag_drivers_117_final_conley_coef.png, tp_lag_drivers_117_final_pdp.png\n")
