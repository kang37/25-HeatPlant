#!/usr/bin/env Rscript
# =============================================================================
# csd_rolling_ar1.R  —— 临界慢化(CSD)式滑动窗口分析
#   Verbesselt 2016 / Forzieri 2022 的做法: 对去季节(去趋势)的植被距平序列,
#   用滑动窗口逐窗算一阶自相关AR(1)(及方差), 看AR(1)是否随时间上升(恢复力下降)。
#   和项目里静态的phi(全年合并一个值)不同, 这里关心的是"同一个站的AR(1)有没有
#   随时间变化"。
#
#   指标: AR(1)=窗口内生长季内相邻8天窗口对的OLS斜率(跨年的断点不配对);
#         rel_sd=窗口内距平标准差/窗口内SIF均值(相对波动, 避免"SIF整体变大
#         导致方差机械变大")。
#   趋势: Kendall tau + Sen斜率(每十年)。窗口重叠, 标准p值虚高, 故用
#         "整年打乱顺序"的置换检验(所有站共用同一置换, 保留空间相关), 给出
#         单站p值与"站集合整体"的场显著性。
#   主设定: 8年窗口、步长1年、去线性趋势; 敏感性: 6年/10年窗口、不去趋势。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(parallel)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")
NCORE <- max(1L, min(6L, detectCores() - 1L)); B <- 1000L; MIN_PAIRS <- 60L

# ---- 站集 ---------------------------------------------------------------------
sel117 <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")$meteo_stat
neg <- fread("data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv")
neg[, drho := rho - rho_min]
LCB <- "data_proc/output_hcsif_buf1000/landcover_break_stations.csv"
if (file.exists(LCB)) { lcb <- fread(LCB); neg <- neg[!meteo_stat %in% lcb[has_break == TRUE]$stat_id] }
sr <- neg[tp >= 0 & tp <= 8]; sr[, pass1 := (drho > 0) & (p_surr < 0.05)]
sel178 <- unique(sr[pass1 == TRUE]$meteo_stat)
dirn <- fread("data_proc/smap_bivar_134/classify_134.csv")[, .(meteo_stat, direction_2v)]

# ---- 数据: 与phi分析同一份SIF∩VPD(生长季8天窗口) ------------------------------
fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)[, .(meteo_stat, year, doy, SIF_buf1000)]
vpd <- as.data.table(readRDS("data_raw/hcsif/vpd_8day_cache.rds"))[, .(meteo_stat, year, doy)]
d   <- merge(sif, vpd, by = c("meteo_stat", "year", "doy")); setorder(d, meteo_stat, year, doy)
YRS <- sort(unique(d$year)); NY <- length(YRS)
stns <- sort(unique(d$meteo_stat)); NS <- length(stns)
cat(sprintf("站数 %d | 年份 %d-%d (%d年) | 117站 %d | 178站 %d\n", NS, min(YRS), max(YRS), NY, length(sel117), length(sel178)))

detr <- function(a, tt) { ok <- !is.na(a); if (sum(ok) < 30) return(rep(NA_real_, length(a)))
  b <- coef(lm(a[ok] ~ tt[ok])); a - (b[1] + b[2] * tt) }
prep <- function(drop) {
  x <- if (length(drop)) d[!year %in% drop] else copy(d)
  x[, clim := mean(SIF_buf1000, na.rm = TRUE), by = .(meteo_stat, doy)]
  x[, anom := SIF_buf1000 - clim]; x[, tt := year + doy / 366]
  x[, anom_dt := detr(anom, tt), by = meteo_stat]; x
}

# ---- 站x年 充分统计量 --------------------------------------------------------------
build_stats <- function(avar, x0) {
  x <- copy(x0); x[, a := get(avar)]
  x[, `:=`(a_lag = shift(a, 1), doy_lag = shift(doy, 1)), by = .(meteo_stat, year)]
  x[, pair_ok := !is.na(a) & !is.na(a_lag) & (doy - doy_lag == 8)]
  ps <- x[pair_ok == TRUE, .(np = .N, sx = sum(a_lag), sy = sum(a), sxx = sum(a_lag^2), sxy = sum(a_lag * a)),
          by = .(meteo_stat, year)]
  ls <- x[!is.na(a) & !is.na(SIF_buf1000), .(na = .N, sa = sum(a), saa = sum(a^2), sS = sum(SIF_buf1000)), by = .(meteo_stat, year)]
  st <- merge(ps, ls, by = c("meteo_stat", "year"), all = TRUE)
  arr <- array(0, dim = c(NS, NY, 9)); K <- c("np", "sx", "sy", "sxx", "sxy", "na", "sa", "saa", "sS")
  si <- match(st$meteo_stat, stns); yi <- match(st$year, YRS)
  for (k in seq_along(K)) { v <- st[[K[k]]]; v[is.na(v)] <- 0; arr[cbind(si, yi, k)] <- v }
  arr
}

# ---- 窗口指标 + 置换 ---------------------------------------------------------------
win_matrix <- function(W) { nw <- NY - W + 1; M <- matrix(0, NY, nw); for (s in seq_len(nw)) M[s:(s + W - 1), s] <- 1; M }
indicators <- function(arr, ord, Wm) {
  g <- function(k) arr[, ord, k] %*% Wm
  np <- g(1); sx <- g(2); sy <- g(3); sxx <- g(4); sxy <- g(5); na <- g(6); sa <- g(7); saa <- g(8); sS <- g(9)
  den <- np * sxx - sx^2
  ar1 <- (np * sxy - sx * sy) / den; ar1[np < MIN_PAIRS | den <= 0] <- NA
  v <- saa / na - (sa / na)^2; rsd <- sqrt(pmax(v, 0)) / (sS / na); rsd[na < MIN_PAIRS] <- NA
  list(ar1 = ar1, rsd = rsd)
}
tau_rows <- function(V) {
  nw <- ncol(V); ij <- which(upper.tri(matrix(0, nw, nw)), arr.ind = TRUE)
  S <- sign(V[, ij[, 2], drop = FALSE] - V[, ij[, 1], drop = FALSE])
  tau <- rowMeans(S, na.rm = TRUE); nvalid <- rowSums(!is.na(V)); tau[nvalid < max(6, ceiling(nw * 0.75))] <- NA
  tau
}
sen_rows <- function(V) {
  nw <- ncol(V); ij <- which(upper.tri(matrix(0, nw, nw)), arr.ind = TRUE)
  Sl <- (V[, ij[, 2], drop = FALSE] - V[, ij[, 1], drop = FALSE]) / matrix(rep(ij[, 2] - ij[, 1], each = nrow(V)), nrow(V))
  apply(Sl, 1, median, na.rm = TRUE) * 10
}

run_config <- function(label, avar, W, seed, drop = integer(0)) {
  cat(sprintf("\n[%s] 窗口=%d年, 变量=%s, 排除年份=%s ... ", label, W, avar, paste(drop, collapse = "&"))); t0 <- Sys.time()
  arr <- build_stats(avar, prep(drop)); Wm <- win_matrix(W); nw <- ncol(Wm)
  obs <- indicators(arr, seq_len(NY), Wm)
  tau_obs <- list(ar1 = tau_rows(obs$ar1), rsd = tau_rows(obs$rsd))
  sen_obs <- list(ar1 = sen_rows(obs$ar1), rsd = sen_rows(obs$rsd))
  set.seed(seed)
  perm <- mclapply(seq_len(B), function(b) { o <- sample(NY); ind <- indicators(arr, o, Wm)
    list(ar1 = tau_rows(ind$ar1), rsd = tau_rows(ind$rsd)) }, mc.cores = NCORE)
  tp <- list(ar1 = do.call(rbind, lapply(perm, `[[`, "ar1")), rsd = do.call(rbind, lapply(perm, `[[`, "rsd")))
  cat(sprintf("%.0fs\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  list(label = label, W = W, avar = avar, nw = nw, obs = obs, tau_obs = tau_obs, sen_obs = sen_obs, tau_perm = tp)
}

CFGS <- list(list("主设定(8年,去趋势)", "anom_dt", 8L, 1L), list("6年窗口", "anom_dt", 6L, 2L),
             list("10年窗口", "anom_dt", 10L, 3L), list("8年窗口,不去趋势", "anom", 8L, 4L),
             list("8年,去2018", "anom_dt", 8L, 5L, 2018L),
             list("8年,去2016&2018", "anom_dt", 8L, 6L, c(2016L, 2018L)),
             list("6年,去2016&2018", "anom_dt", 6L, 7L, c(2016L, 2018L)),
             list("10年,去2016&2018", "anom_dt", 10L, 8L, c(2016L, 2018L)))
res <- lapply(CFGS, function(cf) do.call(run_config, c(cf[1:4], if (length(cf) > 4) list(drop = cf[[5]]))))

# ---- 汇总 ----------------------------------------------------------------------------
dd <- merge(data.table(meteo_stat = stns), dirn, by = "meteo_stat", all.x = TRUE)
SETS <- list(`全部站` = rep(TRUE, NS), `178站(CCM筛选)` = stns %in% sel178, `117站(CCM严格)` = stns %in% sel117,
             `117站-Promote` = stns %in% sel117 & dd$direction_2v %in% "Promote",
             `117站-Inhibit` = stns %in% sel117 & dd$direction_2v %in% "Inhibit",
             `117站-Ambiguous` = stns %in% sel117 & dd$direction_2v %in% "Ambiguous")
p2 <- function(obs, perm) (1 + colSums(abs(perm) >= matrix(abs(obs), nrow(perm), length(obs), byrow = TRUE), na.rm = TRUE)) / (nrow(perm) + 1)

stn_out <- list(); set_out <- list()
for (r in res) for (ind in c("ar1", "rsd")) {
  tau <- r$tau_obs[[ind]]; tp <- r$tau_perm[[ind]]; sen <- r$sen_obs[[ind]]
  pv <- p2(tau, tp); hi <- apply(tp, 2, quantile, .975, na.rm = TRUE); lo <- apply(tp, 2, quantile, .025, na.rm = TRUE)
  stn_out[[length(stn_out) + 1]] <- data.table(config = r$label, indicator = ind, meteo_stat = stns, tau = tau, p_perm = pv,
                                               sen_per_decade = sen, first_win = r$obs[[ind]][, 1], last_win = r$obs[[ind]][, ncol(r$obs[[ind]])])
  for (sn in names(SETS)) { m <- SETS[[sn]] & !is.na(tau); if (sum(m) < 5) next
    med_obs <- median(tau[m]); med_null <- apply(tp[, m, drop = FALSE], 1, median, na.rm = TRUE)
    fh_obs <- mean(tau[m] > hi[m]); fl_obs <- mean(tau[m] < lo[m])
    fh_null <- rowMeans(tp[, m, drop = FALSE] > matrix(hi[m], nrow(tp), sum(m), byrow = TRUE), na.rm = TRUE)
    fl_null <- rowMeans(tp[, m, drop = FALSE] < matrix(lo[m], nrow(tp), sum(m), byrow = TRUE), na.rm = TRUE)
    set_out[[length(set_out) + 1]] <- data.table(config = r$label, indicator = ind, set = sn, n = sum(m),
      median_tau = med_obs, pct_tau_pos = 100 * mean(tau[m] > 0), median_sen_per_decade = median(sen[m]),
      pct_sig_up = 100 * fh_obs, pct_sig_down = 100 * fl_obs,
      p_field_median = (1 + sum(abs(med_null) >= abs(med_obs))) / (B + 1),
      p_field_up = (1 + sum(fh_null >= fh_obs)) / (B + 1), p_field_down = (1 + sum(fl_null >= fl_obs)) / (B + 1)) }
}
stn_tab <- rbindlist(stn_out); set_tab <- rbindlist(set_out)
stn_tab[, `:=`(in117 = meteo_stat %in% sel117, in178 = meteo_stat %in% sel178)]
fwrite(stn_tab, file.path(OUT, "csd_rolling_ar1_station_trends.csv"))
fwrite(set_tab, file.path(OUT, "csd_rolling_ar1_set_summary.csv"))
options(width = 200)
RAW <- "主设定(8年,去趋势)"; CLN <- "8年,去2016&2018"
show <- function(cfg) print(set_tab[config == cfg, .(indicator, set, n, median_tau = round(median_tau, 3), pct_tau_pos = round(pct_tau_pos, 1),
  sen = round(median_sen_per_decade, 4), pct_sig_up = round(pct_sig_up, 1), pct_sig_down = round(pct_sig_down, 1),
  p_med = round(p_field_median, 3), p_up = round(p_field_up, 3), p_down = round(p_field_down, 3))])
cat("\n=============== 原始数据(含2016/2018) ===============\n"); show(RAW)
cat("\n=============== 排除2018 ===============\n"); show("8年,去2018")
cat("\n=============== 清洗版: 排除2016&2018 ===============\n"); show(CLN)
cat("\n=============== 敏感性(全部站 / 117站) ===============\n")
print(set_tab[set %in% c("全部站", "117站(CCM严格)"), .(config, indicator, set, n, median_tau = round(median_tau, 3), pct_tau_pos = round(pct_tau_pos, 1),
  sen = round(median_sen_per_decade, 4), pct_sig_up = round(pct_sig_up, 1), p_up = round(p_field_up, 3), p_med = round(p_field_median, 3))][order(indicator, set, config)])

# ---- 轨迹: 原始 vs 清洗 ----------------------------------------------------------------
mk_traj <- function(r, ver) { ctr <- YRS[1] + (seq_len(r$nw) - 1) + (r$W - 1) / 2
  rbindlist(lapply(c("ar1", "rsd"), function(ind) rbindlist(lapply(names(SETS)[1:3], function(sn) {
    V <- r$obs[[ind]][SETS[[sn]], , drop = FALSE]
    data.table(version = ver, indicator = ind, set = sn, center_year = ctr, med = apply(V, 2, median, na.rm = TRUE),
               q25 = apply(V, 2, quantile, .25, na.rm = TRUE), q75 = apply(V, 2, quantile, .75, na.rm = TRUE)) })))) }
getr <- function(lbl) res[[which(sapply(res, `[[`, "label") == lbl)]]
traj <- rbind(mk_traj(getr(RAW), "Raw data"), mk_traj(getr(CLN), "Excluding 2016 & 2018"))
fwrite(traj, file.path(OUT, "csd_rolling_ar1_trajectory.csv"))
cat("\n轨迹(全部站, rsd 中位数):\n"); print(dcast(traj[indicator == "rsd" & set == "全部站"], center_year ~ version, value.var = "med")[, lapply(.SD, round, 4)])

# ---- 图 ------------------------------------------------------------------------------
Sys.setlocale("LC_ALL", "en_US.UTF-8")
en <- c(`全部站` = "All stations", `178站(CCM筛选)` = "CCM-screened (178)", `117站(CCM严格)` = "CCM-strict (117)")
traj[, set_en := factor(en[set], levels = en)]
traj[, ind_en := factor(ifelse(indicator == "ar1", "Lag-1 autocorrelation (AR1)", "Relative SD (SD / mean SIF)"))]
traj[, version := factor(version, levels = c("Raw data", "Excluding 2016 & 2018"))]
FS <- 16; COL <- c("#555555", "#2166AC", "#C1121F")
p1 <- ggplot(traj, aes(center_year, med, color = set_en, fill = set_en)) +
  geom_ribbon(data = traj[set == "全部站"], aes(ymin = q25, ymax = q75), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.1) + geom_point(size = 1.8) +
  facet_grid(ind_en ~ version, scales = "free_y") +
  scale_color_manual(values = COL) + scale_fill_manual(values = COL) +
  labs(x = "Window centre year (8-year window)", y = "Median across stations", color = NULL, fill = NULL) +
  theme_bw(base_size = FS) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "csd_rolling_ar1_trajectory.png"), p1, width = 12, height = 8, dpi = 300)

tt <- stn_tab[config %in% c(RAW, CLN) & !is.na(tau)]
tt[, `:=`(ind_en = factor(ifelse(indicator == "ar1", "AR1", "Relative SD"), levels = c("AR1", "Relative SD")),
          version = factor(ifelse(config == RAW, "Raw data", "Excluding 2016 & 2018"), levels = c("Raw data", "Excluding 2016 & 2018")))]
tt_l <- rbindlist(list(tt[, .(ind_en, version, tau, set_en = "All stations")], tt[in178 == TRUE, .(ind_en, version, tau, set_en = "CCM-screened (178)")],
                       tt[in117 == TRUE, .(ind_en, version, tau, set_en = "CCM-strict (117)")]))
tt_l[, set_en := factor(set_en, levels = c("All stations", "CCM-screened (178)", "CCM-strict (117)"))]
p2g <- ggplot(tt_l, aes(tau, fill = set_en)) + geom_vline(xintercept = 0, color = "grey40") +
  geom_density(alpha = 0.35, color = NA) + facet_grid(version ~ ind_en) +
  scale_fill_manual(values = COL) +
  labs(x = "Kendall tau of window indicator vs. time (per station)", y = "Density", fill = NULL) +
  theme_bw(base_size = FS) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "csd_rolling_ar1_tau_distribution.png"), p2g, width = 12, height = 8, dpi = 300)

ss <- set_tab[set %in% c("全部站", "117站(CCM严格)")]
lv <- c("主设定(8年,去趋势)", "8年,去2018", "8年,去2016&2018", "6年窗口", "6年,去2016&2018", "10年窗口", "10年,去2016&2018", "8年窗口,不去趋势")
lb <- c("8y raw", "8y excl 2018", "8y excl 2016&18", "6y raw", "6y excl 2016&18", "10y raw", "10y excl 2016&18", "8y raw, no detrend")
ss[, `:=`(set_en = ifelse(set == "全部站", "All stations", "CCM-strict (117)"), ind_en = ifelse(indicator == "ar1", "AR1", "Relative SD"),
          cfg = factor(config, levels = lv, labels = lb))]
p3 <- ggplot(ss, aes(cfg, median_tau, fill = set_en)) + geom_hline(yintercept = 0, color = "grey40") +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = ifelse(p_field_median < 0.05, "*", "")), position = position_dodge(0.8), vjust = -0.3, size = 8) +
  facet_wrap(~ind_en, ncol = 2) + scale_fill_manual(values = c("#555555", "#C1121F")) +
  labs(x = NULL, y = "Median Kendall tau across stations", fill = NULL, caption = "* field-significance p < 0.05 (shared year-permutation)") +
  theme_bw(base_size = FS) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(OUT, "csd_rolling_ar1_sensitivity.png"), p3, width = 13, height = 6.5, dpi = 300)
cat("\n完成。\n")
