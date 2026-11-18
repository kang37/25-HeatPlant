#!/usr/bin/env Rscript
# =============================================================================
# recovery_phi_drivers_varsets.R
#   phi(AR1恢复力) driver 分析的系统对比: 2个站集(117严格/178宽松) x 3个自变量
#   集(15var去稀疏 / 10var去共线精简 / 10var+VPD站级均值CV) = 6个RF模型。
#   phi 现场重新拟合(不复用之前已按18变量完整案例过滤过的data.csv, 那样会在
#   测小变量集时白白损失本可以救回来的站)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ranger)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

# ---- 1. 两个站集清单 --------------------------------------------------------
sel117 <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")$meteo_stat
neg <- fread("data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv")
neg[, drho := rho - rho_min]
LC_BREAK_F <- "data_proc/output_hcsif_buf1000/landcover_break_stations.csv"
if (file.exists(LC_BREAK_F)) {
  lcb <- fread(LC_BREAK_F); excl <- lcb[has_break == TRUE]$stat_id
  neg <- neg[!meteo_stat %in% excl]
}
sig_range <- neg[tp >= 0 & tp <= 8]; sig_range[, pass1 := (drho > 0) & (p_surr < 0.05)]
sel178 <- unique(sig_range[pass1 == TRUE]$meteo_stat)
cat("117站:", length(sel117), " | 178站:", length(sel178), "\n")

# ---- 2. 逐站现场重拟合phi(纯AR1, Forzieri式), 覆盖178站全量 -----------------
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
  list(phi = co["lag1", "Estimate"], n = nrow(sub))
}
phi <- d[, fit_one(.SD), by = meteo_stat]
cat("成功拟合phi的站数(178口径全量):", nrow(phi), "\n")

# ---- 3. 协变量: 18个CORE + VPD站级均值/CV -----------------------------------
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]
cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18 <- env$dt
vc <- fread(file.path(OUT, "vpd_station_covariate.csv"))[, .(meteo_stat, vpd_mean_stn, vpd_cv_stn)]

d_all <- merge(phi, dt18, by.x = "meteo_stat", by.y = "stat_id")
d_all <- merge(d_all, vc, by = "meteo_stat")

VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))
VAR10 <- setdiff(VAR15, c("mixedleaf", "deci_needle", "tavg", "rh", "rsds_mean"))
VAR10VPD <- c(VAR10, "vpd_mean_stn", "vpd_cv_stn")
cat("15var:", length(VAR15), "| 10var:", length(VAR10), "| 10var+VPD:", length(VAR10VPD), "\n")

# ---- 4. 6个组合: 2站集 x 3变量集 --------------------------------------------
run_one <- function(stset_name, stset_ids, varset_name, varset) {
  dd <- d_all[meteo_stat %in% stset_ids]
  dd <- dd[complete.cases(dd[, varset, with = FALSE])]
  if (nrow(dd) < 20) return(list(summary = data.table(station_set=stset_name, var_set=varset_name, n=nrow(dd), rf_r2=NA_real_), imp=NULL))
  f <- as.formula(paste("phi ~", paste(varset, collapse = "+")))
  set.seed(42)
  rf <- ranger(f, data = dd, num.trees = 2000, importance = "permutation", seed = 42)
  imp <- data.table(station_set=stset_name, var_set=varset_name,
                     var = names(rf$variable.importance), importance = rf$variable.importance)
  list(summary = data.table(station_set=stset_name, var_set=varset_name, n=nrow(dd), rf_r2=rf$r.squared), imp = imp)
}

combos <- list(
  list("117","15var"), list("117","10var"), list("117","10var+VPD"),
  list("178","15var"), list("178","10var"), list("178","10var+VPD")
)
varsets <- list(`15var`=VAR15, `10var`=VAR10, `10var+VPD`=VAR10VPD)
stsets  <- list(`117`=sel117, `178`=sel178)

res <- lapply(combos, function(cb) run_one(cb[[1]], stsets[[cb[[1]]]], cb[[2]], varsets[[cb[[2]]]]))
summ <- rbindlist(lapply(res, `[[`, "summary"))
cat("\n=== 6个组合汇总 ===\n"); print(summ[order(-rf_r2)])
fwrite(summ, file.path(OUT, "recovery_phi_drivers_varsets_summary.csv"))

imp_all <- rbindlist(lapply(res, `[[`, "imp"), fill = TRUE)
fwrite(imp_all, file.path(OUT, "recovery_phi_drivers_varsets_importance.csv"))

best <- summ[order(-rf_r2)][1]
cat("\n=== 最优组合:", best$station_set, "x", best$var_set, "(R2=", round(best$rf_r2,3), ") 重要性前8 ===\n")
print(imp_all[station_set==best$station_set & var_set==best$var_set][order(-importance)][1:8])

# ---- 5. 图 -------------------------------------------------------------------
imp_all[, rank_within := rank(-importance), by = .(station_set, var_set)]
p <- ggplot(imp_all[rank_within<=12], aes(reorder(var, importance), importance)) +
  geom_col(fill="#2166AC") + coord_flip() +
  facet_grid(station_set ~ var_set, scales = "free_y") +
  labs(title = "phi RF importance across 2 station sets x 3 variable sets",
       subtitle = paste(sprintf("%s/%s: n=%d R2=%.3f", summ$station_set, summ$var_set, summ$n, summ$rf_r2), collapse=" | "),
       x=NULL, y="Permutation importance") +
  theme_bw(base_size=9) + theme(plot.subtitle=element_text(size=6.5,color="grey35"))
ggsave(file.path(OUT, "recovery_phi_drivers_varsets_importance.png"), p, width=11, height=7, dpi=300)
cat("\n-> recovery_phi_drivers_varsets_importance.png\n")
