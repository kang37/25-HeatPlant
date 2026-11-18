#!/usr/bin/env Rscript
# =============================================================================
# recovery_phi_variants_invest_check.R
#   recovery_phi_climate_adjusted_ED.R 的15变量里没有invest(覆盖不全, 2026-09起被剔除)。
#   这里在同一批φ版本(raw_cc / adjA / ED_vpd / ED_sif)上把invest加回去, 只用invest非缺失
#   的站, 比较"15变量"与"15变量+invest": CV R²/RF OOB R²的变化, invest的岭回归系数+bootstrap CI,
#   OLS系数, 边际与偏相关。另报告invest与ntl/imperv/气候的相关(它是不是"管理"的独立信息)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(glmnet); library(ranger)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify"); NBOOT <- 500L
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]; cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18 <- env$dt; VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))
ph <- fread(file.path(OUT, "recovery_phi_variants_178.csv"))
cat("178站φ版本站数:", nrow(ph), "\n")
dd0 <- merge(ph, dt18[, c("stat_id", VAR15, "invest", "urban_rate"), with = FALSE], by.x = "meteo_stat", by.y = "stat_id")
cat("invest非缺失:", sum(!is.na(dd0$invest)), "/", nrow(dd0), " | urban_rate非缺失:", sum(!is.na(dd0$urban_rate)), "/", nrow(dd0), "\n")
cat("\ninvest与其它变量的Spearman(站点层, 非缺失站):\n")
for (ex in c("invest", "urban_rate")) { xx <- dd0[!is.na(get(ex)) & complete.cases(dd0[, VAR15, with = FALSE])]
  cat("\n", ex, "(n=", nrow(xx), ") 与其它变量的Spearman:\n"); print(round(sapply(c("ntl", "imperv", "tavg", "rh", "precip", "elev", "grass", "ever_broad", if (ex == "urban_rate") "invest"), function(v) cor(xx[[ex]], xx[[v]], method = "spearman", use = "pairwise")), 2)) }

cvr2 <- function(X, y, reps = 20) mean(sapply(seq_len(reps), function(r) { set.seed(r); cv <- cv.glmnet(X, y, alpha = 0, nfolds = 10, standardize = FALSE)
  1 - cv$cvm[cv$lambda == cv$lambda.1se] / var(y) * (length(y) / (length(y) - 1)) }))
run_one <- function(v, ex) {
  x <- dd0[!is.na(get(ex)) & complete.cases(dd0[, VAR15, with = FALSE])]; setnames(x, ex, "EXV")
  dd <- copy(x); dd[, yy := as.numeric(scale(get(paste0("phi_", v))))]; dd <- dd[is.finite(yy)]
  n <- nrow(dd); res <- list()
  for (nm in c("without", "with")) {
    VARS <- if (nm == "with") c(VAR15, "EXV") else VAR15
    f <- as.formula(paste("yy ~", paste(VARS, collapse = "+")))
    set.seed(42); rf <- ranger(f, data = dd, num.trees = 2000, seed = 42)
    X <- scale(as.matrix(dd[, VARS, with = FALSE]))
    res[[nm]] <- c(rf_oob = rf$r.squared, ridge_cv = cvr2(X, dd$yy))
  }
  VARS <- c(VAR15, "EXV"); X <- scale(as.matrix(dd[, VARS, with = FALSE])); yv <- dd$yy
  set.seed(42); cv <- cv.glmnet(X, yv, alpha = 0, nfolds = 10, standardize = FALSE)
  rc <- as.numeric(coef(cv, s = "lambda.1se"))[-1]; names(rc) <- VARS
  set.seed(1); bm <- replicate(NBOOT, { ix <- sample(n, replace = TRUE)
    tryCatch(as.numeric(coef(cv.glmnet(X[ix, , drop = FALSE], yv[ix], alpha = 0, nfolds = 5, standardize = FALSE), s = "lambda.1se"))[-1], error = function(e) rep(NA_real_, length(VARS))) })
  lo <- quantile(bm[which(VARS == "EXV"), ], .025, na.rm = TRUE); hi <- quantile(bm[which(VARS == "EXV"), ], .975, na.rm = TRUE)
  ols <- summary(lm(as.formula(paste("yy ~", paste(VARS, collapse = "+"))), data = dd))$coefficients["EXV", ]
  ctrl <- c("tavg", "rh", "precip", "ntl", "imperv")
  r1 <- resid(lm(as.formula(paste("yy ~", paste(ctrl, collapse = "+"))), data = dd)); r2 <- resid(lm(as.formula(paste("EXV ~", paste(ctrl, collapse = "+"))), data = dd))
  data.table(extra = ex, variant = v, n = n, marginal_rho = cor(dd$yy, dd$EXV, method = "spearman"), partial_r = cor(r1, r2),
             rf_without = res$without["rf_oob"], rf_with = res$with["rf_oob"], ridgecv_without = res$without["ridge_cv"], ridgecv_with = res$with["ridge_cv"],
             extra_ridge = rc["EXV"], ci_lo = lo, ci_hi = hi, ci_excl0 = (lo > 0) | (hi < 0), ols_coef = ols[1], ols_p = ols[4])
}
out <- rbindlist(lapply(c("invest", "urban_rate"), function(ex) rbindlist(lapply(c("raw_cc", "adjA", "ED_vpd", "ED_sif"), run_one, ex = ex))))
out[, `:=`(d_rf = rf_with - rf_without, d_ridgecv = ridgecv_with - ridgecv_without)]
cat("\n=== 加入额外变量(invest / urban_rate)前后(同一批站, y标准化) ===\n")
print(out[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])
fwrite(out, file.path(OUT, "recovery_phi_variants_invest_check.csv"))
cat("\n完成。\n")
