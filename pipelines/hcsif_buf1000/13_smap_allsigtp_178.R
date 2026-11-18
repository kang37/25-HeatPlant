#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 13_smap_allsigtp_178.R
#
# 12_ccm_causal_confirmation.R 的标准(1)("tp∈[0,8]内至少1个tp同时满足
# drho>0且p_surr<0.05")本身就是一个独立可用的站点集合(178站, 见funnel.csv
# 第4行"1_drho_and_psurr_same_tp_in_0to8")，比标准(1)+(2)合起来的117站更宽
# ——(2)是"全局optimal_tp>=0"的时间顺序检验，和"p值/收敛性检验"是两件事。
#
# 用户要画一张热力图: 横轴tp(0..8), 纵轴站点, 格子在该(站,tp)本身显著时上色
# (不是只在06/08脚本用的单个optimal_tp上算)，颜色=S-map方向。这需要给178站
# 的全部324个显著(站,tp)格子逐个算二变量S-map方向，06脚本只在optimal_tp
# 一个点上算，这里的范围更大，所以另开一个脚本，不改06。
#
# 方法与06/10完全一致: 状态空间[SIF(t),VPD(t-tp)], E=2, theta=2, embedded=TRUE,
# 中位数系数+主导符号占比>=75%阈值判Promote/Inhibit/Ambiguous。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(rEDM)
})

PROJ     <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"
setwd(PROJ)
BUF_R    <- 1000L
STAT_DIR <- file.path(PROJ, "data_raw/hcsif/station_v3")
CACHE    <- file.path(PROJ, "data_raw/hcsif/vpd_8day_cache.rds")
SIF_COL  <- sprintf("SIF_buf%d", BUF_R)
OUT_DIR  <- file.path(PROJ, "data_proc/smap_bivar_134")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_PTS      <- 40L
MIN_DAYS_WIN <- 5L
SEED_BASE    <- 20260905L
THETA        <- 2
DRHO_MIN     <- 0
PSURR_MAX    <- 0.05
TP_SIG_LO    <- 0; TP_SIG_HI <- 8

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# ===========================================================================
# 1. 重建178站x显著tp格子清单(照抄12号脚本标准(1)的逻辑, 不改标准, 只是不在
#    这里额外加标准(2)的 optimal_tp>=0 过滤)
# ===========================================================================
SRC <- "data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv"
d0  <- fread(SRC)
d0[, drho := rho - rho_min]

LC_BREAK_F <- file.path(PROJ, "data_proc/output_hcsif_buf1000/landcover_break_stations.csv")
lcb  <- fread(LC_BREAK_F)
excl <- lcb[has_break == TRUE]$stat_id
d0   <- d0[!meteo_stat %in% excl]

sig_range <- d0[tp >= TP_SIG_LO & tp <= TP_SIG_HI]
sig_range[, pass1 := (drho > DRHO_MIN) & (p_surr < PSURR_MAX)]
cells <- sig_range[pass1 == TRUE, .(meteo_stat, tp, rho, drho, p_surr, E)]
st178 <- unique(cells$meteo_stat)
log_msg("178站标准(仅p值+收敛性, 不含optimal_tp方向检验): ", length(st178), " 站, ",
        nrow(cells), " 个显著(站,tp)格子")
stopifnot(length(st178) == 178L)

# ===========================================================================
# 2. SIF + VPD 合并, 归一化(与06/10一致)
# ===========================================================================
fs  <- list.files(STAT_DIR, pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
sif <- rbindlist(lapply(fs, fread), fill = TRUE)
sif <- sif[meteo_stat %in% st178, c("meteo_stat", "year", "doy", "date", SIF_COL), with = FALSE]

vpd <- readRDS(CACHE)[n_days >= MIN_DAYS_WIN & meteo_stat %in% st178]

d <- merge(sif, vpd, by = c("meteo_stat", "year", "doy"))
setorder(d, meteo_stat, year, doy)
d[, idx8 := as.integer(round(as.numeric(
  as.Date(as.character(date), format = "%Y%m%d") - as.Date("2000-01-01")) / 8))]

znorm <- function(x) {
  if (sum(!is.na(x)) < 3) return(rep(NA_real_, length(x)))
  s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
for (v in c(SIF_COL, "vpd_mean")) d[, paste0(v, "_dt") := znorm(get(v)), by = meteo_stat]
SIF_V <- paste0(SIF_COL, "_dt"); VPD_V <- "vpd_mean_dt"

log_msg("SIF-VPD 匹配后: ", nrow(d), " 条, ", uniqueN(d$meteo_stat), " 站")

# ===========================================================================
# 3. 逐(站,tp)二变量 S-map (与06一致的函数, 只是这里对每个显著格子都算一次)
# ===========================================================================
bivar_smap <- function(sid, tp_x, theta = THETA) {
  dd <- d[meteo_stat == sid]
  dd[, x_lag := shift(get(VPD_V), n = tp_x, type = "lag"), by = year]
  dd <- dd[!is.na(idx8)]
  tmp  <- dd[, .(idx8, y = get(SIF_V), x = x_lag)]
  full <- tmp[data.table(idx8 = seq.int(min(dd$idx8), max(dd$idx8))), on = "idx8"]
  setorder(full, idx8)
  n_valid <- full[!is.na(y) & !is.na(x), .N]
  if (n_valid < MIN_PTS) return(NULL)

  df <- data.frame(time = seq_len(nrow(full)), sif = full$y, vpd = full$x)
  n <- nrow(df)
  set.seed(SEED_BASE + as.integer(sid) + tp_x * 1000L)

  sm <- tryCatch(rEDM::SMap(dataFrame = df, E = 2, theta = theta,
                             lib = paste("1", n), pred = paste("1", n),
                             columns = c("sif", "vpd"), target = "sif", embedded = TRUE),
                 error = function(e) NULL)
  if (is.null(sm)) return(NULL)
  co <- as.data.table(sm$coefficients)
  cc <- which(grepl("vpd", colnames(co), ignore.case = TRUE))[1]
  if (is.na(cc)) return(NULL)
  cv <- co[[cc]]
  ok <- !is.na(cv) & !is.nan(cv)
  cv <- cv[ok]
  if (!length(cv)) return(NULL)
  list(coef = cv, n_valid = n_valid)
}

direction_call <- function(cv, thresh = 0.75) {
  fp <- mean(cv > 0)
  if (fp >= thresh) "Promote" else if (fp <= (1 - thresh)) "Inhibit" else "Ambiguous"
}

log_msg("开始逐(站,tp)二变量 S-map, 共 ", nrow(cells), " 个格子 ...")
res <- rbindlist(lapply(seq_len(nrow(cells)), function(i) {
  r <- cells[i]
  out <- tryCatch(bivar_smap(r$meteo_stat, r$tp), error = function(e) NULL)
  if (i %% 50 == 0) log_msg("  ", i, "/", nrow(cells))
  if (is.null(out)) {
    return(data.table(meteo_stat = r$meteo_stat, tp = r$tp, n_coef_pts = 0L,
                       median_coef_2v = NA_real_, frac_positive_2v = NA_real_,
                       direction_2v = NA_character_))
  }
  data.table(meteo_stat = r$meteo_stat, tp = r$tp, n_coef_pts = length(out$coef),
             median_coef_2v = median(out$coef), frac_positive_2v = mean(out$coef > 0),
             direction_2v = direction_call(out$coef))
}), fill = TRUE)

cells <- merge(cells, res, by = c("meteo_stat", "tp"))
log_msg("完成: ", cells[!is.na(direction_2v), .N], "/", nrow(cells), " 格子算出方向")
print(cells[, .N, by = direction_2v])

saveRDS(cells, file.path(OUT_DIR, "smap_allsigtp_178.rds"))
fwrite(cells, file.path(OUT_DIR, "smap_allsigtp_178.csv"))
log_msg("已写出 ", file.path(OUT_DIR, "smap_allsigtp_178.{rds,csv}"))
