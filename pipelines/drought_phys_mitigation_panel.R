#!/usr/bin/env Rscript
# =============================================================================
# drought_phys_mitigation_panel.R
#   一阶段(交叉层/站内交互)检验"投资是否缓解干旱对SIF的影响":
#     SIF_z(t) = a*SIF_z(t-1) + b*D(t) + 控制距平 + 年份固定效应
#                + c_between * [D(t) x 投资(站间: 站点均值)] + c_within * [D(t) x 投资(站内: 当年-站点均值)]
#                + D(t) x (precip, vpd均值/CV, ntl, imperv 标准化调节变量) + e
#   D = VPD距平 或 -SPEI(越大越干); 城市聚类稳健SE。
#   c>0: 投资越高, 干旱斜率越不负(干旱对SIF的抑制越弱) = 缓解。
#   c_within利用同一站点不同年份的投资变化(年份FE控制共同冲击); 另做"上一年投资"版本。
#   数据: drought_phys_stage1b_extended.R 保存的逐窗口距平面板。
# =============================================================================
suppressPackageStartupMessages({library(data.table)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")
strip_city <- function(x) sub("(市辖区|地区|自治州|自治县|特别行政区|城区|盟|市|县|区)$", "", x)
env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]; cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18 <- env$dt; smap <- env$smap

P <- readRDS(file.path(OUT, "drought_phys_panel_anoms.rds"))
vp <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))
vc <- vp[, .(vpd_mean_stn = mean(vpd_mean, na.rm = TRUE), vpd_cv_stn = sd(vpd_mean, na.rm = TRUE) / mean(vpd_mean, na.rm = TRUE)), by = meteo_stat]
iv <- fread("data_raw/green_invest/city_invest_data.csv", encoding = "UTF-8"); setnames(iv, 1, "city_raw"); iv[, city := strip_city(city_raw)]
ycol <- grep("^园林绿化_20(0[3-9]|1[0-9]|2[0-2])$", names(iv), value = TRUE)
ivl <- melt(iv, id.vars = "city", measure.vars = ycol, variable.name = "yr_col", value.name = "raw")
ivl[, year := as.integer(sub("^园林绿化_", "", yr_col))]; ivl <- ivl[is.finite(raw), .(inv_y = log1p(mean(raw))), by = .(city, year)]
sm <- smap[, .(stat_id, city)]; sm[, city := as.character(city)]
inv_yr <- merge(sm, ivl, by = "city", allow.cartesian = TRUE)[, .(meteo_stat = stat_id, year, inv_y, city)]
inv_yr[, inv_between_raw := mean(inv_y), by = meteo_stat]; inv_yr[, inv_within_raw := inv_y - inv_between_raw]
setorder(inv_yr, meteo_stat, year); inv_yr[, inv_within_lag_raw := shift(inv_within_raw, 1), by = meteo_stat]
inv_yr[year == 2018, inv_within_lag_raw := NA_real_]   # 2017年无投资数据, 2018年的"上一年"不可用

mods <- dt18[, .(meteo_stat = stat_id, precip, ntl, imperv)]; mods <- merge(mods, vc, by = "meteo_stat")
st <- unique(inv_yr[, .(meteo_stat, city, inv_between_raw)])
d <- merge(P, mods, by = "meteo_stat"); d <- merge(d, inv_yr[, .(meteo_stat, year, inv_within_raw, inv_within_lag_raw)], by = c("meteo_stat", "year"), all.x = TRUE)
d <- merge(d, st, by = "meteo_stat", all.x = TRUE)
d <- d[!is.na(inv_between_raw) & !is.na(inv_within_raw)]
zs <- function(x) as.numeric(scale(x))
stn_tab <- unique(d[, .(meteo_stat, precip, vpd_mean_stn, vpd_cv_stn, ntl, imperv, inv_between_raw)])
for (v in c("precip", "vpd_mean_stn", "vpd_cv_stn", "ntl", "imperv", "inv_between_raw")) stn_tab[[paste0(v, "_z")]] <- zs(stn_tab[[v]])
d <- merge(d, stn_tab[, .SD, .SDcols = c("meteo_stat", grep("_z$", names(stn_tab), value = TRUE))], by = "meteo_stat")
sdw <- sd(d$inv_within_raw, na.rm = TRUE); d[, inv_within_z := inv_within_raw / sdw]; d[, inv_within_lag_z := inv_within_lag_raw / sdw]
MED_P <- median(stn_tab$precip)
cat(sprintf("面板: %d窗口 | %d站 | %d城市 | 年份%d-%d | invest站内SD(原尺度)=%.3f, 站间SD=%.3f\n", nrow(d), uniqueN(d$meteo_stat), uniqueN(d$city),
            min(d$year), max(d$year), sdw, sd(stn_tab$inv_between_raw)))

cl_rob <- function(X1, y, g) {
  fit <- lm.fit(X1, y); e <- fit$residuals; XtXi <- solve(crossprod(X1)); S_g <- rowsum(X1 * e, g); meat <- crossprod(S_g)
  G <- nrow(S_g); n <- nrow(X1); kk <- ncol(X1); V <- XtXi %*% meat %*% XtXi * (G / (G - 1)) * ((n - 1) / (n - kk))
  list(coef = fit$coefficients, se = sqrt(diag(V)), G = G, n = n) }
fit_spec <- function(dd, drought, ctrl, within_var, label, sub) {
  dd <- copy(dd); dd[, D := get(drought)]
  need <- c("sif_z", "lag1", "D", ctrl, within_var); dd <- dd[complete.cases(dd[, ..need])]
  if (sub == "dry") dd <- dd[precip <= MED_P]
  yrs <- model.matrix(~ factor(year), dd)[, -1, drop = FALSE]
  inter <- cbind(D_x_between = dd$D * dd$inv_between_raw_z, D_x_within = dd$D * dd[[within_var]],
                 D_x_precip = dd$D * dd$precip_z, D_x_vpdmean = dd$D * dd$vpd_mean_stn_z, D_x_vpdcv = dd$D * dd$vpd_cv_stn_z,
                 D_x_ntl = dd$D * dd$ntl_z, D_x_imperv = dd$D * dd$imperv_z)
  cm <- as.matrix(dd[, ..ctrl]); colnames(cm) <- ctrl
  X1 <- cbind(Intercept = 1, lag1 = dd$lag1, D = dd$D, cm, within_main = dd[[within_var]], inter, yrs)
  rb <- cl_rob(X1, dd$sif_z, dd$city)
  cf <- function(nm) { i <- match(nm, colnames(X1)); c(est = unname(rb$coef[i]), se = unname(rb$se[i]), p = unname(2 * pnorm(-abs(rb$coef[i] / rb$se[i])))) }
  b <- cf("D"); cb <- cf("D_x_between"); cw <- cf("D_x_within")
  data.table(spec = label, drought = drought, subset = sub, n = rb$n, G = rb$G, base_slope = b[["est"]],
             between_est = cb[["est"]], between_se = cb[["se"]], between_p = cb[["p"]],
             within_est = cw[["est"]], within_se = cw[["se"]], within_p = cw[["p"]]) }
res <- list()
for (drought in c("vpd_z", "dry_spei")) {
  ctrl <- if (drought == "vpd_z") c("tavg_z", "rsds_z", "prcp_z") else c("rsds_z")
  for (sub in c("all", "dry")) for (wv in c("inv_within_z", "inv_within_lag_z")) {
    lab <- if (wv == "inv_within_z") "within(same-year invest)" else "within(previous-year invest)"
    res[[length(res) + 1]] <- tryCatch(fit_spec(d, drought, ctrl, wv, lab, sub), error = function(e) { message(e$message); NULL }) } }
res <- rbindlist(res)
res[, drought := ifelse(drought == "vpd_z", "VPD", "-SPEI")]
cat("\n===== 一阶段交互模型(系数>0=投资越高, 干旱对SIF的抑制越弱=缓解) =====\n基线斜率b: SIF_z对1个单位干旱距平的响应(同一模型内)\n")
print(res[, .(drought, subset, spec, n, G, base_slope = round(base_slope, 3),
              between = sprintf("%.4f (se %.4f, p=%.3f)", between_est, between_se, between_p),
              within = sprintf("%.4f (se %.4f, p=%.3f)", within_est, within_se, within_p))])
allp <- rbind(res[, .(drought, subset, spec, term = "between", est = between_est, p = between_p)], res[, .(drought, subset, spec, term = "within", est = within_est, p = within_p)])
allp[, p_holm := p.adjust(p, "holm")]; allp[, direction := ifelse(est > 0, "mitigating", "opposite")]
cat("\n全部", nrow(allp), "个交互系数: 缓解方向", sum(allp$est > 0), "个; p<0.05:", sum(allp$p < 0.05), "个 (其中缓解方向", sum(allp$p < 0.05 & allp$est > 0), "个); Holm校正后p<0.05:", sum(allp$p_holm < 0.05), "个\n")
print(allp[order(p)][1:6, .(drought, subset, spec, term, est = round(est, 4), p = round(p, 4), p_holm = round(p_holm, 3), direction)])
fwrite(res, file.path(OUT, "drought_phys_mitigation_panel.csv")); fwrite(allp, file.path(OUT, "drought_phys_mitigation_panel_terms.csv"))
cat("\n完成。\n")
