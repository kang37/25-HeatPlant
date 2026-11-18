#!/usr/bin/env Rscript
# =============================================================================
# tp_lag_drivers_summary_all.R
#   汇总本轮"optimal_tp(效应滞后) ~ 协变量"驱动因素分析里试过的全部变量集/
#   站点集组合，统一用同一套方法(OLS + Conley空间HAC标准误, 200km Bartlett核)
#   重新计算显著变量数，做成一张可比的总结表——之前每一轮用的显著性检验方法
#   不完全一致(有的只做了RF+弹性网、没做正式OLS显著性检验；有的用了普通OLS
#   标准误没做空间校正；只有最后"定稿版"和"尺度稳健性检查"用了Conley HAC)，
#   直接把历史数字拼到一张表里不可比，所以这里全部按统一标准补算一遍。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(car)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/2025-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

conley_lm <- function(fit, lon, lat, cut = 200) {
  X <- model.matrix(fit); u <- residuals(fit)
  S <- X * u
  dx <- outer(lon, lon, "-") * cos(mean(lat) * pi / 180) * 111.32
  dy <- outer(lat, lat, "-") * 111.32
  K  <- pmax(0, 1 - sqrt(dx^2 + dy^2) / cut); dim(K) <- c(length(lon), length(lon))
  bread <- solve(crossprod(X))
  V <- bread %*% (crossprod(S, K %*% S)) %*% bread
  n <- nrow(X); k <- ncol(X); V * (n / (n - k))
}
fit_conley <- function(f, d) {
  m <- lm(f, data = d)
  Vc <- conley_lm(m, d$longitude, d$latitude, cut = 200)
  b <- coef(m); se <- sqrt(diag(Vc)); tt <- b/se; pv <- 2*pt(-abs(tt), df=nrow(d)-length(b))
  pv <- pv[names(pv) != "(Intercept)"]
  list(r2 = summary(m)$r.squared, n = nrow(d), p = pv,
       n_sig05 = sum(pv < .05, na.rm=TRUE), n_sig10 = sum(pv < .1, na.rm=TRUE),
       sig_vars05 = names(pv)[pv < .05], sig_vars10 = names(pv)[pv < .1 & pv >= .05])
}

# ---- 数据基座: 18个CORE变量 + 117站/178站两套optimal_tp --------------------
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt <- env$dt
st <- fread("data_raw/hcsif/stations_924.csv")

tp117 <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")
d117 <- merge(merge(tp117, dt, by.x="meteo_stat", by.y="stat_id"), st, by="meteo_stat")

neg <- fread("data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv")
neg[, drho := rho - rho_min]
lcb <- fread("data_proc/output_hcsif_buf1000/landcover_break_stations.csv")
neg <- neg[!meteo_stat %in% lcb[has_break==TRUE]$stat_id]
sig_range <- neg[tp>=0 & tp<=8]; sig_range[, pass1 := (drho>0) & (p_surr<0.05)]
opt178 <- sig_range[meteo_stat %in% unique(sig_range[pass1==TRUE]$meteo_stat) & pass1==TRUE][
  order(meteo_stat,-rho)][, .SD[1], by=meteo_stat][, .(meteo_stat, optimal_tp=tp)]
d178 <- merge(merge(opt178, dt, by.x="meteo_stat", by.y="stat_id"), st, by="meteo_stat")

CORE18 <- unlist(env$CORE, use.names=FALSE)
S15 <- setdiff(CORE18, c("invest","urban_rate","cloud"))
S16u <- c(S15, "urban_rate")
S10 <- c("imperv","grass","water","ever_needle","ever_broad","deci_broad","precip","rsds_sd","elev","ntl")
S11i <- c(S10, "invest")

configs <- list(
  list(name="A. 18变量(含invest+urban_rate+cloud)", data=d117, pred=CORE18),
  list(name="B. 15变量(去掉3个稀疏变量)",           data=d117, pred=S15),
  list(name="C. 16变量(15+urban_rate)",             data=d117, pred=S16u),
  list(name="D. 16变量(15+urban_rate), 178站松判据", data=d178, pred=S16u),
  list(name="E. 10变量(理论+数据精简)",              data=d117, pred=S10),
  list(name="F. 11变量(10精简+invest, 定稿版)",       data=d117, pred=S11i)
)

res <- rbindlist(lapply(configs, function(cfg) {
  d <- cfg$data[complete.cases(cfg$data[, cfg$pred, with=FALSE])]
  f <- as.formula(paste("optimal_tp ~", paste(cfg$pred, collapse="+")))
  fit <- fit_conley(f, d)
  data.table(配置 = cfg$name, 因变量 = "optimal_tp",
             自变量数 = length(cfg$pred), 建模站数 = fit$n,
             OLS_R2 = round(fit$r2, 3),
             显著p05 = fit$n_sig05, 显著p10含p05 = fit$n_sig10 + fit$n_sig05,
             显著变量p05 = paste(fit$sig_vars05, collapse=";"),
             边缘显著p10 = paste(fit$sig_vars10, collapse=";"))
}))
cat("=== 全部配置统一用OLS+Conley空间HAC(200km)重算, 可比总结表 ===\n")
print(res)
fwrite(res, file.path(OUT, "tp_lag_drivers_summary_all_configs.csv"))
