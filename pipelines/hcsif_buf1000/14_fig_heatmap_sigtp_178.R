#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 14_fig_heatmap_sigtp_178.R
#
# 178站(12号脚本标准(1): p值+收敛性检验, 不含optimal_tp方向检验)的
# 站x tp 显著性+方向热力图, 两个版本(用户2026-09-15要求):
#   fig5: 主导符号占比法, 阈值从0.75改成0.60(frac_positive>=0.6记Promote,
#         <=0.4记Inhibit, 否则Ambiguous)
#   fig6: 中位数符号法(median_coef_2v>0记Promote, 否则Inhibit; 没有Ambiguous
#         类别——中位数总能给出一个确定符号, 这是它和占比法的本质区别)
#
# 排序(用户2026-09-15改进要求): 纵轴先按显著tp个数分组(1个显著tp的站排最前),
# 组内按显著tp列表数值升序排(沿用上一版), 新增第三级: 同一"tp列表"组内先画
# 促进、再画抑制(fig5还有Ambiguous, 排最后)。
#
# 纵轴标签左边新增各nsig组的促进/抑制站点数标注(同一ggplot对象里加一个"info"
# 伪tp列, 只在每组的中间那一行写文字, 保证和热力图行完全对齐); 站点数太少的
# nsig组自动合并标注(尾部若干组累计凑到>=MIN_LABEL_N才单独成一条标注)。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(showtext); library(sysfonts)
})
Sys.setlocale("LC_ALL", "en_US.UTF-8")
PROJ <- "/Users/Kang/Library/CloudStorage/Dropbox/RCloud/25-HeatPlant"
setwd(PROJ)
font_add("cjk", "/System/Library/Fonts/STHeiti Medium.ttc")
showtext_auto(); showtext_opts(dpi = 300)

OUT_DIR <- file.path(PROJ, "data_proc/smap_bivar_134")
cells0  <- readRDS(file.path(OUT_DIR, "smap_allsigtp_178.rds"))
dir_colors3 <- c(Promote = "#C2703A", Inhibit = "#2E7D74", Ambiguous = "grey65")
dir_colors2 <- c(Promote = "#C2703A", Inhibit = "#2E7D74")
TP_ALL <- 0:8
MIN_LABEL_N <- 10L   # nsig组标注: 尾部小组累计站点数达到这个数才单独成一条

direction_call_thresh <- function(fp, thresh = 0.6) {
  fifelse(fp >= thresh, "Promote", fifelse(fp <= (1 - thresh), "Inhibit", "Ambiguous"))
}

# ===========================================================================
# 通用函数: 给定一份带 direction 列的 cells 表, 算排序 + 画图 + 存标注表
# ===========================================================================
build_and_plot <- function(cells, dir_levels, dir_colors, title, subtitle, out_file, legend_name) {
  cells <- copy(cells)
  cells[, direction := factor(direction, levels = dir_levels)]

  # ---- 每站显著tp列表(升序) + 显著个数 + 第一显著tp处的方向(排序/计数用) ----
  by_st <- cells[, .(tp_list = list(sort(tp))), by = meteo_stat]
  by_st[, nsig := lengths(tp_list)]
  max_nsig <- max(by_st$nsig)
  for (k in seq_len(max_nsig)) {
    by_st[, paste0("tpk_", k) := sapply(tp_list, function(v) if (length(v) >= k) v[k] else Inf)]
  }
  first_tp <- sapply(by_st$tp_list, `[`, 1)
  dir_first <- cells[match(paste(by_st$meteo_stat, first_tp), paste(meteo_stat, tp)), direction]
  by_st[, dir_rank := as.integer(factor(dir_first, levels = dir_levels))]

  sort_cols <- c("nsig", paste0("tpk_", seq_len(max_nsig)), "dir_rank")
  setorderv(by_st, sort_cols)
  st_order <- by_st$meteo_stat
  n_st <- length(st_order)

  # ---- nsig组标注: 促进/抑制(/不定)站点数, 小组自动合并 --------------------
  grp_n <- by_st[, .N, by = nsig][order(nsig)]
  grp_n[, merge_id := nsig]
  buf <- 0L; buf_ids <- integer(0)
  for (i in seq_len(nrow(grp_n))) {
    if (grp_n$N[i] >= MIN_LABEL_N && buf == 0L) next   # 自己够大, 单独一组
    buf <- buf + grp_n$N[i]; buf_ids <- c(buf_ids, grp_n$nsig[i])
    if (buf >= MIN_LABEL_N || i == nrow(grp_n)) {
      grp_n[nsig %in% buf_ids, merge_id := min(buf_ids)]
      buf <- 0L; buf_ids <- integer(0)
    }
  }
  by_st <- merge(by_st, grp_n[, .(nsig, merge_id)], by = "nsig")
  setorderv(by_st, sort_cols)   # merge 后重新确保顺序不变

  by_st[, dir_count := factor(dir_first, levels = dir_levels)]
  dir_abbr <- c(Promote = "促", Inhibit = "抑", Ambiguous = "不定")
  lab_tab <- by_st[, {
    tb <- table(dir_count)
    nsig_rng <- range(nsig)
    lab_head <- if (nsig_rng[1] == nsig_rng[2]) sprintf("nsig=%d", nsig_rng[1])
                else sprintf("nsig=%d-%d", nsig_rng[1], nsig_rng[2])
    lab_body <- paste(sprintf("%s%d", dir_abbr[names(tb)], as.integer(tb)), collapse = " ")
    .(n_total = .N, mid_row = meteo_stat[ceiling(.N / 2)],
      label = paste0(lab_head, "(", .N, ")  ", lab_body))
  }, by = merge_id]

  # ---- 拼完整 站 x tp(0..8) 网格(数值x轴), 左侧留一块空白区放标注文字 --------
  grid <- CJ(meteo_stat = st_order, tp = TP_ALL)
  plot_df <- merge(grid, cells[, .(meteo_stat, tp, direction)], by = c("meteo_stat", "tp"), all.x = TRUE)
  plot_df[, meteo_stat_f := factor(meteo_stat, levels = rev(st_order))]

  lab_df <- merge(lab_tab, data.table(meteo_stat = st_order), by.x = "mid_row", by.y = "meteo_stat")
  lab_df[, meteo_stat_f := factor(mid_row, levels = rev(st_order))]

  fig_height <- max(6, 1.4 + n_st * 0.075)
  LABEL_SIZE <- 2.15 * 3   # 用户要求: 分组数量文字放大到原来3倍
  gap_x    <- -0.65   # 标注区和 tp=0 列之间的竖线
  x_anchor <- gap_x - 0.15   # 文字右对齐锚点, 紧贴竖线(用户要求: 拉到图旁边, 不要留大片空白)
  x_lab    <- gap_x - 7.2   # 面板左边界, 刚好够放下最长的标注文字(字体放大3倍后需要的宽度)

  p <- ggplot(plot_df, aes(x = tp, y = meteo_stat_f)) +
    geom_tile(aes(fill = direction), color = "grey92", linewidth = 0.15, width = 0.92, height = 0.92) +
    geom_text(data = lab_df, aes(x = x_anchor, label = label), size = LABEL_SIZE, family = "cjk",
              hjust = 1, lineheight = 0.85) +
    geom_vline(xintercept = gap_x, linewidth = 0.4, color = "grey50") +
    scale_fill_manual(values = dir_colors, na.value = "white", name = legend_name,
                       na.translate = FALSE, breaks = dir_levels) +
    scale_x_continuous(breaks = TP_ALL, labels = as.character(TP_ALL),
                        limits = c(x_lab, max(TP_ALL) + 0.6), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(title = title, subtitle = subtitle, x = "tp (滞后, 每步8天)", y = NULL) +
    theme_minimal(base_size = 11, base_family = "cjk") +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),     # 用户要求: 去掉纵轴站点ID
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8.5),
          plot.margin = margin(5.5, 5.5, 5.5, 5.5),
          legend.position = "right")

  ggsave(file.path(OUT_DIR, out_file), p, width = 13, height = fig_height, dpi = 300, limitsize = FALSE)
  cat("已保存 ", out_file, " (height=", round(fig_height, 1), "in, 站数=", n_st, ")\n", sep = "")
  invisible(list(by_st = by_st, lab_tab = lab_tab))
}

# ===========================================================================
# fig5: 主导符号占比法, 阈值0.60
# ===========================================================================
cells_thresh <- copy(cells0)
cells_thresh[, direction := direction_call_thresh(frac_positive_2v, thresh = 0.6)]
cat("=== fig5(占比法, 阈值0.60) 方向分布 ===\n")
print(cells_thresh[, .N, by = direction])
r5 <- build_and_plot(
  cells_thresh, dir_levels = c("Promote", "Inhibit", "Ambiguous"), dir_colors = dir_colors3,
  title = "178站 tp x 站点 因果显著性+方向热力图(主导符号占比法, 阈值60%)",
  subtitle = paste0("格子上色条件: 该(站,tp)同时满足 drho>0 且 p_surr<0.05(12号脚本标准(1));\n",
                     "方向=系数为正的时间点占比>=60%记Promote/<=40%记Inhibit/否则Ambiguous;\n",
                     "纵轴: 先按显著tp个数分组, 组内按显著tp列表升序, 同列表内先促进后抑制(Ambiguous最后)"),
  out_file = "fig5_heatmap_sigtp_178.png", legend_name = "S-map方向\n(占比法,60%)")

# ===========================================================================
# fig6: 中位数符号法(二分类, 无Ambiguous)
# ===========================================================================
cells_med <- copy(cells0)
cells_med[, direction := fifelse(median_coef_2v > 0, "Promote", "Inhibit")]
cat("\n=== fig6(中位数法) 方向分布 ===\n")
print(cells_med[, .N, by = direction])
r6 <- build_and_plot(
  cells_med, dir_levels = c("Promote", "Inhibit"), dir_colors = dir_colors2,
  title = "178站 tp x 站点 因果显著性+方向热力图(中位数符号法)",
  subtitle = paste0("格子上色条件: 该(站,tp)同时满足 drho>0 且 p_surr<0.05(12号脚本标准(1));\n",
                     "方向=该(站,tp)的S-map系数中位数(median_coef_2v)符号, 正记Promote/负记Inhibit(无Ambiguous类别);\n",
                     "纵轴: 先按显著tp个数分组, 组内按显著tp列表升序, 同列表内先促进后抑制"),
  out_file = "fig6_heatmap_sigtp_178_median.png", legend_name = "S-map方向\n(中位数法)")

# ---- 两版方向标签的一致率(写文档用) -----------------------------------------
cmp <- merge(cells_thresh[, .(meteo_stat, tp, dir_thresh = direction)],
             cells_med[, .(meteo_stat, tp, dir_med = direction)], by = c("meteo_stat", "tp"))
cmp[, agree := dir_thresh == dir_med]
cat("\n=== 两版方向一致率 ===\n")
cat("全部324格子: 一致", cmp[, sum(agree)], "/", nrow(cmp),
    "(", round(mean(cmp$agree) * 100, 1), "%)\n")
cat("仅看占比法非Ambiguous的格子(即占比法本身给出确定方向的格子): ")
sub <- cmp[dir_thresh != "Ambiguous"]
cat(sub[, sum(agree)], "/", nrow(sub), "(", round(mean(sub$agree) * 100, 1), "%)\n")

fwrite(cmp, file.path(OUT_DIR, "sigtp_178_thresh_vs_median.csv"))
