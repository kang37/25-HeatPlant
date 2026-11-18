#!/usr/bin/env Rscript
# =============================================================================
# recovery_omega_10var_vpd_clean.R
#   更干净的 Omega+VPD: 用10变量精简集(已去掉tavg/rh/rsds_mean这三个和VPD
#   强共线的气候变量, 也去掉mixedleaf/deci_needle两个近零方差变量) + VPD均值/CV
#   + direction, 而不是15变量全量+VPD(上次VIF高达72, vpd_mean_stn自己VIF=49,
#   因为VPD就是用tavg+rh算出来的, 硬塞在一起是设计错误)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(car); library(glmnet); library(ranger)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18 <- env$dt
VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))
VAR10 <- setdiff(VAR15, c("mixedleaf", "deci_needle", "tavg", "rh", "rsds_mean"))
cat("10变量精简集:", paste(VAR10, collapse=", "), "\n")

ev <- fread(file.path(OUT, "recovery_liu2021_metrics_dirfix_117_events.csv"))
agg <- ev[, .(logOmega = as.numeric(mean(log(Omega), na.rm=TRUE)),
              direction_2v = direction_2v[1]), by = meteo_stat]
vc <- fread(file.path(OUT, "vpd_station_covariate.csv"))[, .(meteo_stat, vpd_mean_stn, vpd_cv_stn)]
do <- merge(agg, dt18, by.x = "meteo_stat", by.y = "stat_id")
do <- merge(do, vc, by = "meteo_stat")
ALLPRED <- c(VAR10, "vpd_mean_stn", "vpd_cv_stn", "is_promote")
do[, is_promote := as.integer(direction_2v == "Promote")]
do <- do[complete.cases(do[, c(VAR10, "vpd_mean_stn", "vpd_cv_stn"), with = FALSE])]
cat("\n可建模站数:", nrow(do), "(Promote", sum(do$is_promote), "/ Inhibit", sum(1-do$is_promote), ")\n")

# ---- OLS + VIF ---------------------------------------------------------------
ols_o <- lm(as.formula(paste("logOmega ~", paste(ALLPRED, collapse = "+"))), data = do)
vif_o <- sort(vif(ols_o), decreasing = TRUE)
cat("\n=== VIF(降序) ===\n"); print(round(vif_o, 1))
co_o <- as.data.table(summary(ols_o)$coefficients, keep.rownames = "var")
setnames(co_o, c("var","ols_coef","se","t","p")); co_o <- co_o[var != "(Intercept)"]
co_o[, sig := ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,".","")))]
cat("\n=== OLS(R2=", round(summary(ols_o)$r.squared,3), ") ===\n")
print(co_o[order(p)][, .(var, ols_coef=round(ols_coef,4), p=round(p,3), sig)])

# ---- RF ------------------------------------------------------------------------
set.seed(42)
rf_o <- ranger(as.formula(paste("logOmega ~", paste(ALLPRED, collapse = "+"))),
                data = do, num.trees = 2000, importance = "permutation", seed = 42)
imp_o <- data.table(var = names(rf_o$variable.importance), importance = rf_o$variable.importance)
cat("\n=== RF OOB R2 =", round(rf_o$r.squared,3), "===\n"); print(imp_o[order(-importance)])

# ---- 岭回归 + bootstrap 500次 -------------------------------------------------
Xo <- scale(as.matrix(do[, ALLPRED, with = FALSE])); yo <- do$logOmega
set.seed(42); cvo <- cv.glmnet(Xo, yo, alpha = 0, nfolds = 10, standardize = FALSE)
ridge_o <- as.numeric(coef(cvo, s = "lambda.1se"))[-1]; names(ridge_o) <- ALLPRED
predo <- predict(cvo, Xo, s = "lambda.1se")
ridge_r2_o <- 1 - sum((yo-predo)^2)/sum((yo-mean(yo))^2)
cat("\n岭回归(1se) R2 =", round(ridge_r2_o,3), "\n")

set.seed(1)
boot_o <- replicate(500, {
  idx <- sample(nrow(do), replace = TRUE)
  Xb <- Xo[idx,,drop=FALSE]; yb <- yo[idx]
  tryCatch({ cvb<-cv.glmnet(Xb,yb,alpha=0,nfolds=5,standardize=FALSE); as.numeric(coef(cvb,s="lambda.1se"))[-1] },
           error=function(e) rep(NA_real_, length(ALLPRED)))
})
rownames(boot_o) <- ALLPRED
ci_lo_o <- apply(boot_o,1,quantile,.025,na.rm=TRUE); ci_hi_o <- apply(boot_o,1,quantile,.975,na.rm=TRUE)
res_o <- data.table(var=ALLPRED, ridge_coef=round(ridge_o,4), boot_ci_lo=round(ci_lo_o,4),
                     boot_ci_hi=round(ci_hi_o,4), ci_excludes_0=(ci_lo_o>0)|(ci_hi_o<0))
res_o <- merge(res_o, co_o[,.(var,ols_coef,p,sig)], by="var")
res_o <- res_o[order(-abs(ridge_coef))]
cat("\n=== 汇总: 岭回归+bootstrap+OLS ===\n"); print(res_o)
fwrite(res_o, file.path(OUT, "recovery_omega_10var_vpd_clean_inference.csv"))
fwrite(imp_o, file.path(OUT, "recovery_omega_10var_vpd_clean_rf_importance.csv"))

# ---- 图 -----------------------------------------------------------------------
imp_o[, var := factor(var, levels = rev(var[order(importance)]))]
p1 <- ggplot(imp_o, aes(var, importance, fill = var=="is_promote")) +
  geom_col() + coord_flip() +
  scale_fill_manual(values=c(`TRUE`="#C1121F",`FALSE`="#2166AC"), guide="none") +
  labs(title = "logOmega RF importance: 10var(去共线) + VPD + direction",
       subtitle = sprintf("n=%d | RF OOB R2=%.3f | red = is_promote", nrow(do), rf_o$r.squared),
       x=NULL, y="Permutation importance") + theme_bw(base_size=12)
ggsave(file.path(OUT, "recovery_omega_10var_vpd_clean_importance.png"), p1, width=7.5, height=6, dpi=300)

res_o[, var := factor(var, levels = rev(var))]
p2 <- ggplot(res_o, aes(var, ridge_coef, color = ci_excludes_0)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_pointrange(aes(ymin = boot_ci_lo, ymax = boot_ci_hi)) + coord_flip() +
  scale_color_manual(values = c(`TRUE`="#C1121F", `FALSE`="grey50"), guide = "none") +
  labs(title = "logOmega (10var+VPD+direction): ridge coefficients + bootstrap 95% CI",
       subtitle = sprintf("n=%d | ridge R2=%.3f | red = 95%% CI excludes 0", nrow(do), ridge_r2_o),
       x = NULL, y = "Standardized ridge coefficient") + theme_bw(base_size = 12)
ggsave(file.path(OUT, "recovery_omega_10var_vpd_clean_coef_ci.png"), p2, width = 8, height = 6, dpi = 300)

vif_dt <- data.table(var = names(vif_o), VIF = vif_o)
vif_dt[, var := factor(var, levels = rev(var))]
p3 <- ggplot(vif_dt, aes(var, VIF, fill = VIF > 5)) +
  geom_col() + coord_flip() + geom_hline(yintercept = 5, linetype = "dashed") +
  scale_fill_manual(values = c(`TRUE`="#C1121F", `FALSE`="#2166AC"), guide = "none") +
  labs(title = "VIF: logOmega 10var(decorrelated)+VPD+direction model", x=NULL, y="VIF") + theme_bw(base_size=12)
ggsave(file.path(OUT, "recovery_omega_10var_vpd_clean_VIF.png"), p3, width=7, height=5.5, dpi=300)
cat("\n-> 3张图已出\n")
