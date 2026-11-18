#!/usr/bin/env Rscript
# =============================================================================
# viz_phi178_and_omega_inference.R
#   1) phi(178站x15变量)的单独重要性图 + VIF条形图(回答"共线性怎么解决")
#   2) Omega(57站,15var+VPD+direction)补齐OLS+VIF+岭回归+bootstrap全变量推断,
#      和phi那边同一套方法, 图表可比
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(car); library(glmnet); library(ranger)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# =============================================================================
# PART 1: phi(178x15var) 重要性单图 + VIF图
# =============================================================================
imp_all <- fread(file.path(OUT, "recovery_phi_drivers_varsets_importance.csv"))
imp178 <- imp_all[station_set == "178" & var_set == "15var"][order(-importance)]
imp178[, var := factor(var, levels = rev(var))]
p_imp <- ggplot(imp178, aes(var, importance)) +
  geom_col(fill = "#2166AC") + coord_flip() +
  labs(title = "phi RF importance: 178 stations x 15 variables",
       subtitle = "RF OOB R2 = 0.342 (n=178)",
       x = NULL, y = "Permutation importance") +
  theme_bw(base_size = 12)
ggsave(file.path(OUT, "recovery_phi_178_15var_rf_importance_solo.png"), p_imp, width = 7.5, height = 6, dpi = 300)

vif_tab <- fread(file.path(OUT, "recovery_phi_178_15var_inference.csv"))
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))

# 重新拟合一次OLS只为拿VIF画图(和上次推断脚本同一份数据/公式)
neg <- fread("data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv")
neg[, drho := rho - rho_min]
LC_BREAK_F <- "data_proc/output_hcsif_buf1000/landcover_break_stations.csv"
if (file.exists(LC_BREAK_F)) { lcb <- fread(LC_BREAK_F); excl <- lcb[has_break==TRUE]$stat_id; neg <- neg[!meteo_stat %in% excl] }
sig_range <- neg[tp>=0 & tp<=8]; sig_range[, pass1 := (drho>0)&(p_surr<0.05)]
sel178 <- unique(sig_range[pass1==TRUE]$meteo_stat)
fs <- list.files("data_raw/hcsif/station_v3", pattern="^hcsif_station_[0-9]{4}\\.csv$", full.names=TRUE)
sif <- rbindlist(lapply(fs, fread), fill=TRUE)[meteo_stat %in% sel178, .(meteo_stat, year, doy, date, SIF_buf1000)]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[meteo_stat %in% sel178]
d <- merge(sif, vpd, by=c("meteo_stat","year","doy")); setorder(d, meteo_stat, year, doy)
deseason_z <- function(x, doy) { clim<-tapply(x,doy,mean,na.rm=TRUE); a<-x-clim[as.character(doy)]
  s<-sd(a,na.rm=TRUE); if(!is.finite(s)||s==0) return(rep(NA_real_,length(x))); (a-mean(a,na.rm=TRUE))/s }
d[, sif_anom := deseason_z(SIF_buf1000, doy), by=meteo_stat]
d[, lag1 := shift(sif_anom,1), by=.(meteo_stat,year)]
fit_one <- function(sub){ sub<-sub[!is.na(lag1)&!is.na(sif_anom)]; if(nrow(sub)<30) return(NULL)
  m<-tryCatch(lm(sif_anom~lag1,data=sub),error=function(e)NULL); if(is.null(m)) return(NULL)
  co<-summary(m)$coefficients; if(!"lag1"%in%rownames(co)) return(NULL); list(phi=co["lag1","Estimate"]) }
phi <- d[, fit_one(.SD), by=meteo_stat]
dt18 <- env$dt
dd <- merge(phi, dt18, by.x="meteo_stat", by.y="stat_id")
dd <- dd[complete.cases(dd[, VAR15, with=FALSE])]
ols <- lm(as.formula(paste("phi ~", paste(VAR15, collapse="+"))), data=dd)
vif_v <- sort(vif(ols), decreasing = TRUE)
vif_dt <- data.table(var = names(vif_v), VIF = vif_v)
vif_dt[, var := factor(var, levels = rev(var))]
p_vif <- ggplot(vif_dt, aes(var, VIF, fill = VIF > 5)) +
  geom_col() + coord_flip() +
  geom_hline(yintercept = 5, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c(`TRUE`="#C1121F", `FALSE`="#2166AC"), guide = "none") +
  labs(title = "VIF, 178-station x 15-variable phi model",
       subtitle = "dashed line = VIF=5 conventional concern threshold; red = above threshold",
       x = NULL, y = "VIF") + theme_bw(base_size = 12)
ggsave(file.path(OUT, "recovery_phi_178_15var_VIF.png"), p_vif, width = 7, height = 6, dpi = 300)
cat("-> phi 178x15var: 重要性图 + VIF图 已出\n")

# =============================================================================
# PART 2: Omega(57站, 15var+VPD+direction) 补全OLS+VIF+岭回归+bootstrap
# =============================================================================
ev <- fread(file.path(OUT, "recovery_liu2021_metrics_dirfix_117_events.csv"))
agg <- ev[, .(logOmega = as.numeric(mean(log(Omega), na.rm=TRUE)),
              direction_2v = direction_2v[1]), by = meteo_stat]
vc <- fread(file.path(OUT, "vpd_station_covariate.csv"))[, .(meteo_stat, vpd_mean_stn, vpd_cv_stn)]
do <- merge(agg, dt18, by.x = "meteo_stat", by.y = "stat_id")
do <- merge(do, vc, by = "meteo_stat")
VARSET_O <- c(VAR15, "vpd_mean_stn", "vpd_cv_stn")
do <- do[complete.cases(do[, VARSET_O, with = FALSE])]
do[, is_promote := as.integer(direction_2v == "Promote")]
ALLPRED <- c(VARSET_O, "is_promote")
cat("\nOmega模型可建模站数:", nrow(do), "\n")

ols_o <- lm(as.formula(paste("logOmega ~", paste(ALLPRED, collapse = "+"))), data = do)
vif_o <- tryCatch(sort(vif(ols_o), decreasing = TRUE), error = function(e) NULL)
cat("\n=== Omega模型 VIF ===\n"); print(round(vif_o, 1))
co_o <- as.data.table(summary(ols_o)$coefficients, keep.rownames = "var")
setnames(co_o, c("var","ols_coef","se","t","p")); co_o <- co_o[var != "(Intercept)"]
co_o[, sig := ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,".","")))]
cat("\n=== Omega OLS(R2=", round(summary(ols_o)$r.squared,3), ") ===\n")
print(co_o[order(p)][, .(var, ols_coef=round(ols_coef,4), p=round(p,3), sig)])

Xo <- scale(as.matrix(do[, ALLPRED, with = FALSE])); yo <- do$logOmega
set.seed(42); cvo <- cv.glmnet(Xo, yo, alpha = 0, nfolds = 10, standardize = FALSE)
ridge_o <- as.numeric(coef(cvo, s = "lambda.1se"))[-1]; names(ridge_o) <- ALLPRED
predo <- predict(cvo, Xo, s = "lambda.1se")
ridge_r2_o <- 1 - sum((yo-predo)^2)/sum((yo-mean(yo))^2)
cat("\nOmega岭回归(1se) R2 =", round(ridge_r2_o,3), "\n")

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
cat("\n=== Omega 汇总: 岭回归+bootstrap+OLS ===\n"); print(res_o)
fwrite(res_o, file.path(OUT, "recovery_omega_full_inference.csv"))

# RF重要性(已有, 复用)
imp_omega <- fread(file.path(OUT, "recovery_omega_direction_covariate_rf_importance.csv"))
imp_omega[, var := factor(var, levels = rev(var[order(importance)]))]
p_imp_o <- ggplot(imp_omega, aes(var, importance, fill = var=="is_promote")) +
  geom_col() + coord_flip() +
  scale_fill_manual(values=c(`TRUE`="#C1121F",`FALSE`="#2166AC"), guide="none") +
  labs(title = "logOmega RF importance: 57 stations, 15var+VPD+direction",
       subtitle = "RF OOB R2 = 0.512 | red = is_promote", x=NULL, y="Permutation importance") +
  theme_bw(base_size=12)
ggsave(file.path(OUT, "recovery_omega_rf_importance_solo.png"), p_imp_o, width=7.5, height=6.5, dpi=300)

res_o[, var := factor(var, levels = rev(var))]
p_coef_o <- ggplot(res_o, aes(var, ridge_coef, color = ci_excludes_0)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_pointrange(aes(ymin = boot_ci_lo, ymax = boot_ci_hi)) +
  coord_flip() +
  scale_color_manual(values = c(`TRUE`="#C1121F", `FALSE`="grey50"), guide = "none") +
  labs(title = "logOmega: ridge coefficients + bootstrap 95% CI",
       subtitle = sprintf("n=%d | ridge R2=%.3f | red = 95%% CI excludes 0", nrow(do), ridge_r2_o),
       x = NULL, y = "Standardized ridge coefficient") +
  theme_bw(base_size = 12)
ggsave(file.path(OUT, "recovery_omega_coef_ci.png"), p_coef_o, width = 8, height = 6.5, dpi = 300)
cat("\n-> Omega: 重要性图 + 系数CI图 已出\n")
