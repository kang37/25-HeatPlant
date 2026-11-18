#!/usr/bin/env Rscript
# =============================================================================
# recovery_phi_block_pooling.R
#   phi(AR1持续性)按单年估计噪声太大。这里把k个相邻年份的窗口合并成一个"站-块"
#   估一个phi, 比较不同块长k下: 估计噪声(SE)、站间信号占比(ICC)、相邻块一致性、
#   可用行数, 用于决定面板版phi的块长。
#   块长k: 互不重叠的连续k年, 不足k年的尾块丢弃。lag1只在同一年内取(与原phi构造一致)。
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(nlme)})
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"; setwd(PROJ)
OUT  <- file.path(PROJ, "data_proc/output_loose_classify")

sel117 <- fread("data_proc/ccm_hcsif_buf1000_causal_confirmed/stations_confirmed.csv")$meteo_stat
neg <- fread("data_proc/ccm_hcsif_buf1000_fwd_negtp_full/fwd_negtp_full_20260904_1126.csv")
neg[, drho := rho - rho_min]
LCB <- "data_proc/output_hcsif_buf1000/landcover_break_stations.csv"
if (file.exists(LCB)) { lcb <- fread(LCB); neg <- neg[!meteo_stat %in% lcb[has_break == TRUE]$stat_id] }
sr <- neg[tp >= 0 & tp <= 8]; sr[, pass1 := (drho > 0) & (p_surr < 0.05)]
sel178 <- unique(sr[pass1 == TRUE]$meteo_stat)

fs  <- list.files("data_raw/hcsif/station_v3", pattern = "^hcsif_station_[0-9]{4}\\.csv$", full.names = TRUE)
d   <- rbindlist(lapply(fs, fread), fill = TRUE)[meteo_stat %in% sel178, .(meteo_stat, year, doy, SIF_buf1000)]
setorder(d, meteo_stat, year, doy)
deseason_z <- function(x, doy) { clim <- tapply(x, doy, mean, na.rm = TRUE); a <- x - clim[as.character(doy)]
  s <- sd(a, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x))); (a - mean(a, na.rm = TRUE)) / s }
d[, sif_anom := deseason_z(SIF_buf1000, doy), by = meteo_stat]
d[, lag1 := shift(sif_anom, 1), by = .(meteo_stat, year)]
d <- d[!is.na(lag1) & !is.na(sif_anom)]
YEARS <- sort(unique(d$year)); cat("年份:", min(YEARS), "-", max(YEARS), "| 每站每年平均有效lag对:",
  round(d[, .N, by = .(meteo_stat, year)][, mean(N)], 1), "\n")

fit_blk <- function(sub) {
  if (nrow(sub) < 15) return(NULL)
  m <- lm(sif_anom ~ lag1, data = sub); co <- summary(m)$coefficients
  list(phi = co["lag1", "Estimate"], se = co["lag1", "Std. Error"], n = nrow(sub))
}
KS <- c(1, 2, 3, 4, 5, 6, 8, 11)
res <- list(); blocks_all <- list()
for (k in KS) {
  nb <- floor(length(YEARS) / k)
  map <- data.table(year = YEARS[seq_len(nb * k)], blk = rep(seq_len(nb), each = k))
  dk <- merge(d, map, by = "year")
  pb <- dk[, fit_blk(.SD), by = .(meteo_stat, blk)]
  pb[, k := k]; blocks_all[[as.character(k)]] <- pb
  for (set in c("178", "117")) {
    q <- if (set == "117") pb[meteo_stat %in% sel117] else pb
    q <- q[is.finite(phi)]
    vc <- tryCatch(VarCorr(lme(phi ~ 1, random = ~ 1 | meteo_stat, data = as.data.frame(q))), error = function(e) NULL)
    icc <- if (is.null(vc)) NA_real_ else { v <- as.numeric(vc[, "Variance"]); v[1] / sum(v) }
    adj <- merge(q[, .(meteo_stat, blk, phi)], q[, .(meteo_stat, blk = blk - 1L, phi_next = phi)], by = c("meteo_stat", "blk"))
    res[[length(res) + 1]] <- data.table(set = set, k = k, n_blocks_per_stn = nb, rows = nrow(q),
      pairs_per_block = mean(q$n), mean_se = mean(q$se), sd_phi = sd(q$phi),
      noise_share = mean(q$se^2) / var(q$phi),
      icc = icc, adj_block_cor = if (nrow(adj) > 10) cor(adj$phi, adj$phi_next) else NA_real_)
  }
}
res <- rbindlist(res)
res[, c("mean_se", "sd_phi", "noise_share", "icc", "adj_block_cor", "pairs_per_block") :=
      lapply(.SD, round, 3), .SDcols = c("mean_se", "sd_phi", "noise_share", "icc", "adj_block_cor", "pairs_per_block")]
print(res[set == "178"]); cat("\n--- 117站 ---\n"); print(res[set == "117"])
fwrite(res, file.path(OUT, "recovery_phi_block_pooling_summary.csv"))
fwrite(rbindlist(blocks_all), file.path(OUT, "recovery_phi_block_pooling_blocks.csv"))

pl <- melt(res[set == "178", .(k, `Mean SE of phi` = mean_se, `Noise share of variance` = noise_share,
                               `ICC (between-station share)` = icc, `Adjacent-block correlation` = adj_block_cor,
                               `Rows (station-blocks)` = rows)], id.vars = "k")
p <- ggplot(pl, aes(k, value)) + geom_line(color = "#2166AC") + geom_point(color = "#2166AC", size = 2) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) + scale_x_continuous(breaks = KS) +
  labs(title = "phi block pooling: choosing the block length k (years)", subtitle = "178 stations, non-overlapping blocks",
       x = "Block length k (years)", y = NULL) + theme_bw(base_size = 11)
ggsave(file.path(OUT, "recovery_phi_block_pooling.png"), p, width = 10, height = 5.5, dpi = 300)
