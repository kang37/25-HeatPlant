#!/usr/bin/env Rscript
# =============================================================================
# recovery_phi_178_15var_inference.R
#   178站x15变量是178站口径里表现最好的phi(恢复力)模型(RF R2=0.342)。RF只给
#   重要性排序,不给方向/显著性,这里补OLS(+VIF诊断)和岭回归(+bootstrap CI)
#   两条线,回答"每个变量往哪个方向、有多大把握"。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(car); library(glmnet)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 178站清单(与recovery_phi_drivers_varsets.R一致的判据) ------------------
neg <- fread("data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv")
neg[, drho := rho - rho_min]
LC_BREAK_F <- "data_proc/output_hcsif_buf1000/landcover_break_stations.csv"
if (file.exists(LC_BREAK_F)) {
  lcb <- fread(LC_BREAK_F); excl <- lcb[has_break == TRUE]$stat_id
  neg <- neg[!meteo_stat %in% excl]
}
sig_range <- neg[tp >= 0 & tp <= 8]; sig_range[, pass1 := (drho > 0) & (p_surr < 0.05)]
sel178 <- unique(sig_range[pass1 == TRUE]$meteo_stat)

# ---- phi(纯AR1)逐站重拟合, 178站全量 ----------------------------------------
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)[meteo_stat %in% sel178, .(meteo_stat, year, doy, date, SIF_buf1000)]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[meteo_stat %in% sel178]
d   <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
setorder(d, meteo_stat, year, doy)
deseason_z <- function(x, doy) {
  clim <- tapply(x, doy, mean, na.rm = TRUE); a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (a - mean(a, na.rm = TRUE)) / s
}
d[, sif_anom := deseason_z(SIF_buf1000, doy), by = meteo_stat]
d[, lag1     := shift(sif_anom, 1), by = .(meteo_stat, year)]
fit_one <- function(sub) {
  sub <- sub[!is.na(lag1) & !is.na(sif_anom)]
  if (nrow(sub) < 30) return(NULL)
  m <- tryCatch(lm(sif_anom ~ lag1, data = sub), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  co <- summary(m)$coefficients
  if (!"lag1" %in% rownames(co)) return(NULL)
  list(phi = co["lag1", "Estimate"])
}
phi <- d[, fit_one(.SD), by = meteo_stat]

# ---- 18个CORE协变量, 定义15变量子集 -----------------------------------------
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18 <- env$dt
VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))

dd <- merge(phi, dt18, by.x = "meteo_stat", by.y = "stat_id")
dd <- dd[complete.cases(dd[, VAR15, with = FALSE])]
cat("178站x15变量, 可建模站数:", nrow(dd), "\n")

# ---- OLS + VIF ---------------------------------------------------------------
f <- as.formula(paste("phi ~", paste(VAR15, collapse = "+")))
ols <- lm(f, data = dd)
vif_v <- sort(vif(ols), decreasing = TRUE)
cat("\n=== VIF(降序) ===\n"); print(round(vif_v, 1))

co <- as.data.table(summary(ols)$coefficients, keep.rownames = "var")
setnames(co, c("var", "ols_coef", "se", "t", "p"))
co <- co[var != "(Intercept)"]
co[, sig := ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(p < 0.1, ".", "")))]
cat("\n=== OLS(R2=", round(summary(ols)$r.squared,3), ") 系数(未标准化) ===\n")
print(co[order(p)][, .(var, ols_coef = round(ols_coef,4), p = round(p,3), sig)])

# ---- 岭回归(标准化, 处理VIF>5的共线簇) + bootstrap 500次 CI -----------------
X <- scale(as.matrix(dd[, VAR15, with = FALSE])); y <- dd$phi
set.seed(42); cv <- cv.glmnet(X, y, alpha = 0, nfolds = 10, standardize = FALSE)
ridge_coef <- as.numeric(coef(cv, s = "lambda.1se"))[-1]; names(ridge_coef) <- VAR15
pred <- predict(cv, X, s = "lambda.1se")
ridge_r2 <- 1 - sum((y-pred)^2)/sum((y-mean(y))^2)
cat("\n=== 岭回归(1se) R2 =", round(ridge_r2,3), "===\n")

set.seed(1)
boot_mat <- replicate(500, {
  idx <- sample(nrow(dd), replace = TRUE)
  Xb <- X[idx, , drop=FALSE]; yb <- y[idx]
  tryCatch({
    cvb <- cv.glmnet(Xb, yb, alpha = 0, nfolds = 5, standardize = FALSE)
    as.numeric(coef(cvb, s = "lambda.1se"))[-1]
  }, error = function(e) rep(NA_real_, length(VAR15)))
})
rownames(boot_mat) <- VAR15
ci_lo <- apply(boot_mat, 1, quantile, .025, na.rm = TRUE)
ci_hi <- apply(boot_mat, 1, quantile, .975, na.rm = TRUE)
excl0 <- (ci_lo > 0) | (ci_hi < 0)

res <- data.table(var = VAR15, ridge_coef = round(ridge_coef, 4),
                   boot_ci_lo = round(ci_lo, 4), boot_ci_hi = round(ci_hi, 4),
                   ci_excludes_0 = excl0)
res <- merge(res, co[, .(var, ols_coef, p, sig)], by = "var")
res[, direction_agree := sign(ridge_coef) == sign(ols_coef)]
res <- res[order(-abs(ridge_coef))]
cat("\n=== 汇总: 岭回归系数 + bootstrap 95%CI + OLS方向/p值对照 ===\n")
print(res)

fwrite(res, file.path(OUT, "recovery_phi_178_15var_inference.csv"))

# ---- 图 -----------------------------------------------------------------------
res[, var := factor(var, levels = rev(var))]
p <- ggplot(res, aes(var, ridge_coef, color = ci_excludes_0)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_pointrange(aes(ymin = boot_ci_lo, ymax = boot_ci_hi), fatten = 3) +
  coord_flip() +
  scale_color_manual(values = c(`TRUE`="#C1121F", `FALSE`="grey50"), guide = "none") +
  labs(title = "178-station x 15-variable phi model: ridge coefficients + bootstrap 95% CI",
       subtitle = sprintf("n=%d | ridge R2=%.3f | red = 95%% CI excludes 0", nrow(dd), ridge_r2),
       x = NULL, y = "Standardized ridge coefficient") +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size=8, color="grey35"))
ggsave(file.path(OUT, "recovery_phi_178_15var_inference.png"), p, width = 8, height = 6.5, dpi = 300)
cat("\n-> recovery_phi_178_15var_inference.png\n")
