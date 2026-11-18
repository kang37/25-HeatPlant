#!/usr/bin/env Rscript
# =============================================================================
# drought_phys_mitigation_spec_curve.R
#   问题: 投资是否缓解干旱对植物(SIF)的影响? 不挑方案——预设一个网格, 全部跑完并公开:
#     6种"干旱影响"度量(统一为越大=影响越严重) x 3种投资度量 x 4套控制变量 x 3个站点子集 = 216方案
#     每个方案: 站点水平WLS(权重=第一阶段误差权重), 城市聚类稳健SE, DV与invest均标准化。
#     "缓解效应" = -coef(invest)  (>0 表示投资越高, 干旱影响越轻)
#   校正"从很多方案里挑最好的"带来的偏差:
#     (a) 置换检验: 在降水三分位内把整座城市的投资值打乱(保留与气候的关联), 重跑全网格, 得到
#         "最大缓解t""显著缓解方案个数"的零分布;
#     (b) 留出验证: 随机把城市分成两半, 在一半里选最优方案, 在另一半里检验它。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(readxl); library(parallel)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")
NCORE <- max(1L, min(6L, detectCores() - 1L)); NPERM <- 500L; NSPLIT <- 300L
strip_city <- function(x) sub("(市辖区|地区|自治州|自治县|特别行政区|城区|盟|市|县|区)$", "", x)

env <- new.env(); lines <- readLines("22_varpart_3dim.R", encoding = "UTF-8")
lines <- lines[!grepl("^setwd", lines)]; cut <- grep("^CORE <- lapply", lines)[1]
invisible(capture.output(eval(parse(text = paste(lines[1:cut], collapse = "\n")), envir = env)))
dt18 <- env$dt; smap <- env$smap
VAR15 <- setdiff(unlist(env$CORE, use.names = FALSE), c("invest", "urban_rate", "cloud"))
VAR10 <- setdiff(VAR15, c("mixedleaf", "deci_needle", "tavg", "rh", "rsds_mean"))
m <- fread(file.path(OUT, "drought_phys_station_metrics_ext.csv"))
vp <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))
vc <- vp[, .(vpd_mean_stn = mean(vpd_mean, na.rm = TRUE), vpd_cv_stn = sd(vpd_mean, na.rm = TRUE) / mean(vpd_mean, na.rm = TRUE)), by = meteo_stat]
g <- as.data.table(read_excel("data_raw/城市绿地面积数据_2003-2023.xlsx", sheet = 1)); setnames(g, 1, "city_raw"); g[, city := strip_city(city_raw)]
gc <- grep("_绿地面积$", names(g), value = TRUE); g[, green := rowMeans(sapply(gc, function(cn) suppressWarnings(as.numeric(g[[cn]]))), na.rm = TRUE)]
g <- g[is.finite(green) & green > 0, .(green = mean(green)), by = city]

d <- merge(m[, .(meteo_stat, phi, phi_se, beta, beta_se, beta_vpdonly, beta_vpdonly_se, beta_spei, beta_spei_se, R_loss, Z_worse)],
           dt18[, c("stat_id", unique(c(VAR10, "tavg", "rh", "invest"))), with = FALSE], by.x = "meteo_stat", by.y = "stat_id")
d <- merge(d, vc, by = "meteo_stat"); d <- merge(d, smap[, .(stat_id, city)], by.x = "meteo_stat", by.y = "stat_id", all.x = TRUE)
d[is.na(city), city := paste0("stn", meteo_stat)]; d <- merge(d, g, by = "city", all.x = TRUE)
d <- d[complete.cases(d[, c(VAR10, "vpd_mean_stn", "vpd_cv_stn", "tavg", "rh"), with = FALSE])]
tau_w <- function(est, se) { ok <- is.finite(est) & is.finite(se); v <- var(est[ok]); t2 <- max(v - mean(se[ok]^2), 0.25 * v); w <- 1 / (se^2 + t2); w / mean(w[ok]) }
d[, `:=`(S_partial = -beta, S_vpdonly = -beta_vpdonly, S_spei = -beta_spei)]
d[, `:=`(w_S_partial = tau_w(S_partial, beta_se), w_S_vpdonly = tau_w(S_vpdonly, beta_vpdonly_se), w_S_spei = tau_w(S_spei, beta_spei_se),
         w_phi = tau_w(phi, phi_se), w_R_loss = 1, w_Z_worse = 1)]
d[, inv_mean := invest]; d[, inv_green := ifelse(!is.na(invest) & !is.na(green), log(expm1(invest) / green), NA_real_)]
make_resid <- function(dd) { ok <- !is.na(dd$inv_mean); f <- lm(inv_mean ~ precip + tavg + rh + elev + vpd_mean_stn + vpd_cv_stn, data = dd[ok])
  r <- rep(NA_real_, nrow(dd)); r[ok] <- resid(f); r }
d[, inv_resid := make_resid(d)]
MED_P <- median(d$precip); MED_I <- median(d$imperv)
cat(sprintf("站点=%d | 有invest=%d | 城市=%d | inv_green可用=%d\n", nrow(d), sum(!is.na(d$inv_mean)), uniqueN(d[!is.na(inv_mean)]$city), sum(!is.na(d$inv_green))))

DVS <- c("S_partial", "S_vpdonly", "S_spei", "phi", "R_loss", "Z_worse"); INV <- c("inv_mean", "inv_resid", "inv_green"); SS <- c("all", "dry", "urban")
C1 <- c("precip", "vpd_mean_stn", "vpd_cv_stn")
CS <- list(C0 = character(0), C1 = C1, C2 = c(C1, "ntl", "imperv"), C3 = c(VAR10, "vpd_mean_stn", "vpd_cv_stn"))
GRID <- as.data.table(expand.grid(dv = DVS, inv = INV, cs = names(CS), ss = SS, stringsAsFactors = FALSE))

cl_rob <- function(X1, y, w, g) {
  fit <- lm.wfit(X1, y, w); e <- fit$residuals; bread <- solve(crossprod(X1 * sqrt(w)))
  S_g <- rowsum(X1 * (w * e), g); meat <- crossprod(S_g); G <- nrow(S_g); n <- nrow(X1); kk <- ncol(X1)
  V <- bread %*% meat %*% bread * (G / (G - 1)) * ((n - 1) / (n - kk)); list(coef = fit$coefficients, se = sqrt(diag(V)), G = G) }
run_spec <- function(dd, dv, inv, cs, ss) {
  x <- dd
  if (ss == "dry") x <- x[precip <= MED_P] else if (ss == "urban") x <- x[imperv > MED_I]
  x <- x[!is.na(get(inv)) & is.finite(get(dv))]
  if (nrow(x) < 60 || uniqueN(x$city) < 15) return(NULL)
  ctrl <- CS[[cs]]
  Xm <- cbind(1, as.numeric(scale(x[[inv]])), if (length(ctrl)) scale(as.matrix(x[, ctrl, with = FALSE])))
  if (anyNA(Xm)) return(NULL)
  y <- as.numeric(scale(x[[dv]])); w <- x[[paste0("w_", dv)]]; w[!is.finite(w)] <- 1
  rb <- tryCatch(cl_rob(Xm, y, w, x$city), error = function(e) NULL); if (is.null(rb)) return(NULL)
  list(n = nrow(x), G = rb$G, mit = unname(-rb$coef[2]), se = unname(rb$se[2]), t = unname(-rb$coef[2] / rb$se[2]))
}
run_grid <- function(dd) {
  r <- lapply(seq_len(nrow(GRID)), function(i) run_spec(dd, GRID$dv[i], GRID$inv[i], GRID$cs[i], GRID$ss[i]))
  out <- copy(GRID); out[, `:=`(n = NA_integer_, G = NA_integer_, mit = NA_real_, se = NA_real_, t = NA_real_)]
  for (i in seq_along(r)) if (!is.null(r[[i]])) set(out, i, c("n", "G", "mit", "se", "t"), r[[i]])
  out }

obs <- run_grid(d); obs[, `:=`(p2 = 2 * pnorm(-abs(t)), lo = mit - 1.96 * se, hi = mit + 1.96 * se)]
fwrite(obs, file.path(OUT, "drought_phys_mitigation_spec_curve.csv"))
cat(sprintf("\n== 观测网格: %d个方案 ==\n", sum(!is.na(obs$t))))
cat(sprintf("缓解方向(mit>0): %d (%.0f%%) | 显著缓解(t>1.96): %d | 反方向显著(t<-1.96): %d | 效应中位数 %.3f SD (范围 %.3f ~ %.3f)\n",
            sum(obs$mit > 0, na.rm = TRUE), 100 * mean(obs$mit > 0, na.rm = TRUE), sum(obs$t > 1.96, na.rm = TRUE), sum(obs$t < -1.96, na.rm = TRUE),
            median(obs$mit, na.rm = TRUE), min(obs$mit, na.rm = TRUE), max(obs$mit, na.rm = TRUE)))
cat("\n按干旱影响度量汇总:\n"); print(obs[, .(n_spec = sum(!is.na(t)), share_mitigating = round(mean(mit > 0, na.rm = TRUE), 2), n_sig_mit = sum(t > 1.96, na.rm = TRUE),
  n_sig_opp = sum(t < -1.96, na.rm = TRUE), median_mit = round(median(mit, na.rm = TRUE), 3), median_t = round(median(t, na.rm = TRUE), 2)), by = dv])
cat("\n按站点子集汇总:\n"); print(obs[, .(n_spec = sum(!is.na(t)), share_mitigating = round(mean(mit > 0, na.rm = TRUE), 2), n_sig_mit = sum(t > 1.96, na.rm = TRUE),
  n_sig_opp = sum(t < -1.96, na.rm = TRUE), median_mit = round(median(mit, na.rm = TRUE), 3)), by = ss])
cat("\n按投资度量汇总:\n"); print(obs[, .(n_spec = sum(!is.na(t)), share_mitigating = round(mean(mit > 0, na.rm = TRUE), 2), n_sig_mit = sum(t > 1.96, na.rm = TRUE),
  n_sig_opp = sum(t < -1.96, na.rm = TRUE), median_mit = round(median(mit, na.rm = TRUE), 3)), by = inv])
cat("\n按控制变量汇总:\n"); print(obs[, .(n_spec = sum(!is.na(t)), share_mitigating = round(mean(mit > 0, na.rm = TRUE), 2), n_sig_mit = sum(t > 1.96, na.rm = TRUE),
  n_sig_opp = sum(t < -1.96, na.rm = TRUE), median_mit = round(median(mit, na.rm = TRUE), 3)), by = cs])
cat("\n最支持'缓解'的前8个方案:\n"); print(obs[order(-t)][1:8, .(dv, inv, cs, ss, n, G, mit = round(mit, 3), t = round(t, 2), p2 = round(p2, 4))])
cat("\n最反对'缓解'的前5个方案:\n"); print(obs[order(t)][1:5, .(dv, inv, cs, ss, n, G, mit = round(mit, 3), t = round(t, 2), p2 = round(p2, 4))])

# ---- (a) 置换检验: 降水三分位内整城打乱 ----------------------------------------
city_tab <- d[!is.na(inv_mean), .(precip_med = median(precip), inv_mean = inv_mean[1], inv_green = inv_green[1]), by = city]
city_tab[, tert := cut(precip_med, quantile(precip_med, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = FALSE)]
perm_once <- function(seed) {
  set.seed(seed); ct <- copy(city_tab); ct[, src := city[sample(.N)], by = tert]
  mp <- ct[, .(city, inv_mean_p = inv_mean[match(src, city)], inv_green_p = inv_green[match(src, city)])]
  dd <- copy(d); dd <- merge(dd, mp, by = "city", all.x = TRUE); dd[, `:=`(inv_mean = inv_mean_p, inv_green = inv_green_p)]
  dd[, inv_resid := make_resid(dd)]; g2 <- run_grid(dd)
  c(max_t = max(g2$t, na.rm = TRUE), n_sig = sum(g2$t > 1.96, na.rm = TRUE), share_pos = mean(g2$mit > 0, na.rm = TRUE), median_t = median(g2$t, na.rm = TRUE)) }
cat("\n运行置换检验(", NPERM, "次, 每次重跑全网格)...\n")
PERM_F <- file.path(OUT, "drought_phys_mitigation_perm_null.csv")
pm <- if (file.exists(PERM_F)) fread(PERM_F) else as.data.table(do.call(rbind, mclapply(seq_len(NPERM), function(i) perm_once(7000 + i), mc.cores = NCORE)))
o <- c(max_t = max(obs$t, na.rm = TRUE), n_sig = sum(obs$t > 1.96, na.rm = TRUE), share_pos = mean(obs$mit > 0, na.rm = TRUE), median_t = median(obs$t, na.rm = TRUE))
cat(sprintf("观测: 最大t=%.2f | 显著缓解方案数=%d | 缓解方向占比=%.2f | t中位数=%.2f\n", o["max_t"], o["n_sig"], o["share_pos"], o["median_t"]))
cat(sprintf("置换零分布: 最大t 均值=%.2f (95%%分位=%.2f) | 显著缓解方案数 均值=%.1f (95%%分位=%.0f) | 缓解方向占比 均值=%.2f | t中位数 均值=%.2f\n",
            mean(pm$max_t), quantile(pm$max_t, .95), mean(pm$n_sig), quantile(pm$n_sig, .95), mean(pm$share_pos), mean(pm$median_t)))
cat(sprintf("置换p值: 最大t p=%.3f | 显著缓解方案数 p=%.3f | t中位数(整体方向) p=%.3f\n",
            (1 + sum(pm$max_t >= o["max_t"])) / (NPERM + 1), (1 + sum(pm$n_sig >= o["n_sig"])) / (NPERM + 1), (1 + sum(pm$median_t >= o["median_t"])) / (NPERM + 1)))
fwrite(pm, file.path(OUT, "drought_phys_mitigation_perm_null.csv"))

# ---- (b) 留出验证: 一半城市选最优, 另一半检验 -----------------------------------
inv_cities <- unique(d[!is.na(inv_mean)]$city)
split_once <- function(seed) {
  set.seed(seed); a <- sample(inv_cities, floor(length(inv_cities) / 2)); A <- d[city %in% a | is.na(inv_mean)]; B <- d[!(city %in% a) | is.na(inv_mean)]
  gA <- run_grid(A); i <- which.max(gA$t); if (!length(i) || is.na(gA$t[i])) return(c(tA = NA, tB = NA, mitB = NA))
  rB <- run_spec(B, GRID$dv[i], GRID$inv[i], GRID$cs[i], GRID$ss[i]); c(tA = gA$t[i], tB = if (is.null(rB)) NA else rB$t, mitB = if (is.null(rB)) NA else rB$mit) }
cat("\n运行留出验证(", NSPLIT, "次随机分半)...\n")
sp <- as.data.table(do.call(rbind, mclapply(seq_len(NSPLIT), function(i) split_once(9000 + i), mc.cores = NCORE))); sp <- sp[!is.na(tB)]
cat(sprintf("发现集最优方案的t均值=%.2f; 同一方案在留出集的t均值=%.2f (若真有效应应接近发现集)\n留出集: t>1.96占%.1f%% | t<-1.96占%.1f%% | 缓解方向占%.1f%% | 平均缓解效应=%.3f SD\n",
            mean(sp$tA), mean(sp$tB), 100 * mean(sp$tB > 1.96), 100 * mean(sp$tB < -1.96), 100 * mean(sp$mitB > 0), mean(sp$mitB)))
fwrite(sp, file.path(OUT, "drought_phys_mitigation_holdout.csv"))

# ---- 图 --------------------------------------------------------------------------
ob <- obs[!is.na(t)][order(mit)]; ob[, rk := .I]
ggsave(file.path(OUT, "drought_phys_mitigation_spec_curve.png"),
  ggplot(ob, aes(rk, mit, color = dv)) + geom_hline(yintercept = 0, color = "grey60") +
    geom_linerange(aes(ymin = lo, ymax = hi), alpha = 0.35) + geom_point(size = 1.2) +
    labs(title = sprintf("Specification curve: does invest mitigate drought impact? (%d specs)", nrow(ob)),
         subtitle = "mitigation effect = -coef(invest) in SD of 'drought impact' per SD invest; >0 = invest associated with milder impact; bars = 95% city-cluster CI",
         x = "Specifications ranked by effect", y = "Mitigation effect (SD)", color = "Impact metric") + theme_bw(base_size = 10),
  width = 10, height = 5.5, dpi = 300)
ggsave(file.path(OUT, "drought_phys_mitigation_heatmap.png"),
  ggplot(obs[!is.na(t)], aes(paste(inv, cs, sep = " | "), dv, fill = pmax(pmin(t, 4), -4))) + geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f", t)), size = 2.3) + facet_wrap(~ss, ncol = 1, labeller = as_labeller(c(all = "all stations", dry = "dry half", urban = "high-imperv half"))) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#C1121F", midpoint = 0, name = "t (mitigation)") +
    labs(title = "Mitigation t-statistics across the specification grid (red = invest associated with milder drought impact)", x = "invest measure | control set", y = NULL) +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1)), width = 11, height = 9, dpi = 300)
cat("\n完成。\n")
