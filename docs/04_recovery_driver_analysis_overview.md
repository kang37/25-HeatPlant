# 恢复力/抵抗力驱动因素分析——项目总览（截至 2026-09-16）

本文档是一份**横向总览快照**，汇总"CCM确认因果耦合站点 → 各类恢复力/因果属性指标 →
driver回归"这条探索线上目前尝试过的全部因变量、自变量组合与分析方法。不按单个
模块的"状态速览+更新日志"格式写（内容横跨模块二/三，且是快照式总结，不是持续
追加的日志），过时后需要人工核对更新，不会随分析进展自动同步。

具体每个模块的持续性状态记录仍以 [02_ccm.md](02_ccm.md)、[03_classification.md](03_classification.md)
为准；本文档只做"全局导航图"用途。

## 一、恢复力/抵抗力相关因变量（7个）

同一批117/178站的SIF/VPD原始序列，用不同文献的构造方式衍生出7个因变量：

| 因变量 | 构造方式来源 | 主要脚本 |
|---|---|---|
| Trec_windows | Schwalm et al. 2017式：跌破基线到回正所需窗口数 | `recovery_schwalm_style_demo.R`、`recovery_vs_134_direction.R` |
| φ(phi) | Forzieri et al. 2022式：纯AR1持续性系数 | `recovery_ar1_pure.R`、`recovery_phi_drivers_117/178.R`、`recovery_phi_drivers_varsets.R`、`recovery_phi_178_15var_inference.R` |
| constraint_time | Liu et al. 2021式：受抑制时长 | `recovery_liu2021_metrics_117.R` 及 `_dirfix` 版 |
| recovery_time | Liu et al. 2021式：恢复时长 | 同上 |
| logOmega | Isbell et al. 2015式：抵抗力Ω=Yn/\|Ye-Yn\|，取对数 | `recovery_liu2021_drivers_117.R`、`_ridge_117.R`、`_ridge_depth_control_117.R`、`_drivers_dirfix_117.R`、`_panel_mixed_117/_dirfix.R`、`recovery_omega_direction_covariate.R`、`recovery_omega_10var_vpd_clean.R` |
| logDelta | Isbell et al. 2015式：恢复力Δ=\|Ye-Yn\|/\|Ye+1-Yn\|，取对数 | 同logOmega系列脚本 |
| fall_rate / climb_rate | Yi & Jackson 2021式：跌幅/爬升幅度 除以 对应时长 | `recovery_rate_drivers_117.R` |

**自变量组合(7种)**：18变量全量CORE；15变量(去invest/urban_rate/cloud)；16变量(15+urban_rate)；
10变量精简(去共线簇+近零方差变量)；10变量+VPD站级均值/CV(`vpd_station_covariate.R`)；
10变量+VPD+方向哑变量(is_promote)；动态面板(8个逐年动态变量+9个静态地表覆盖+elev,
由 `recovery_liu2021_panel_mixed_117.R` 组装)。

**分析方法(6种)**：OLS(+VIF诊断)；岭回归(glmnet alpha=0)；岭回归+bootstrap(500次重抽样)
经验置信区间；随机森林(ranger)+置换重要性；混合效应模型(nlme::lme,站点随机截距,
用于站-年面板)；生存分析(Cox回归+Kaplan-Meier,用于处理Trec的右删失)。

**站点集(3种)**：117(严格,通过CCM统一判据的全部两条标准)；178(宽松,只通过标准(1),
不要求optimal_tp≥0)；57(117里方向明确的Promote+Inhibit站,剔除60个Ambiguous)。

**约35-40次独立模型拟合**（非完整7×7×6×3全排列，是按每轮发现的问题针对性调整）。

**结果要点**：
- fall_rate/climb_rate：**全部无信号**(RF R²为负)，最早排除的方向。
- Trec/constraint_time/recovery_time：信号弱到中等；发现过一次实质性方法漏洞——
  未按促进/抑制方向翻转"最严重/最佳恢复"取法，混用方向会把噪声算进抵抗力，
  修正(`_dirfix`版)后logOmega解释力从R²≈0.33升到0.57。
- **logOmega 是这组里最强信号**：干净版(57站,10变量精简+VPD+direction)
  RF R²≈0.49-0.51、岭回归R²≈0.49-0.54。稳定驱动因素：precip/tavg/rh/rsds_sd等
  气候变量(logOmega=log(Yn/|Ye-Yn|)越大=抵抗力越强，这些变量系数为正=越暖湿抵抗力反而越强)+ is_promote(独立负向,Promote站抵抗力更差,
  与生存分析结论互相印证)。VPD站级均值/CV在此因变量上**不提供独立信息**
  （与tavg/rh共线，去共线后仍不显著）。
- **logDelta：从头到尾测不出信号**，是本项目里唯一稳定复现"无解释力"的因变量。
- **φ：仅次于logOmega的第二强信号**，178站×15变量最优配置RF R²=0.372(全样本，
  非18变量版打折后的86站)。VPD变异系数(非均值)被证实是有效补充自变量，
  在φ模型里进入重要性前3，但在logOmega模型里不显著——两个因变量对VPD统计量
  的响应不同。

产出目录：`data_proc/output_loose_classify/recovery_*`、`vpd_station_covariate.csv`。
配套Word报告：`note/20260915_phi178_omega_report.docx`（φ178站×15变量 + Omega清洁版
两个最终模型的完整可视化+推断）。

### logOmega 促进/抑制分开建模（sign(median)重分类，2026-09-16）

**动机**：上面的57站logOmega模型把Promote(41)+Inhibit(16)两组用`is_promote`哑变量
一起回归，隐含"两组共享同一套斜率、方向只做加性偏移"的假设；同时75%占优符号
阈值把178个CCM显著站里104个（58%）划成Ambiguous弃用，样本浪费严重。

**方案**：不用75%占优符号阈值，改用每站最佳(最大rho)显著tp处`median_coef_2v`的
**符号**直接二分（不留Ambiguous）；置信度权重=`|frac_positive_2v-0.5|`（越接近0.5
=越像抛硬币=权重越低）。诊断显示：旧Ambiguous组的`frac_positive_2v`中位数=0.4888，
几乎精确卡在抛硬币中点，104站按新规则几乎对半分(50→Promote/54→Inhibit)——说明
样本增益主要来自"以前因为不够一致被弃用、现在符号一确定就能用、但可信度天然较低"
的站点，因此必须配权重使用，不是免费午餐。

**结果**（脚本`recovery_omega_178_bydirection.R`，10变量精简+VPD协变量，加权RF+
加权岭回归+加权bootstrap 500次，VIF全部<5）：

| 组 | n | RF OOB R² | 岭回归R² | 稳健正向驱动(bootstrap CI不含0) | 稳健负向驱动 |
|---|---|---|---|---|---|
| Promote | 90 | 0.359 | 0.41 | precip(+0.082), rsds_sd(+0.039), ever_broad(+0.020), water(+0.020) | 无（vpd_mean_stn负向但CI含0） |
| Inhibit | 88 | 0.289 | 0.14 | rsds_sd(+0.017), vpd_cv_stn(+0.012), precip(+0.011), ever_needle(+0.011) | imperv(-0.007) |

**解读**：
- 站点数从41/16翻倍到90/88，两组量级终于可比，且都远大于57站合并版。
- Promote组信号明显强于Inhibit组（岭回归R²0.41 vs 0.14）——分开后才看出两组
  可解释程度本身不对称，这一点在57站合并模型里被平均掉、看不出来。
- precip、rsds_sd在两组里都是稳健正向驱动，方向一致，与合并模型"暖湿气候不利于
  抵抗力"的结论**不冲突**（这里的正向是指对logOmega/抵抗力的提升，需注意与旧
  15var+VPD混杂模型的符号约定核对，不要跨脚本直接比系数值，只比方向）。
- **vpd_cv_stn只在Inhibit组里稳健**（+0.012，CI不含0），Promote组不显著——
  与φ模型"VPD变异系数是有效因子"的发现部分印证，但揭示了方向依赖性：VPD
  变异对VPD-促进型站点的抵抗力没有独立解释力，只对VPD-抑制型站点有效。
- Inhibit组多个bootstrap CI下界卡在0.0000（而非负值），提示该组加权后bootstrap
  分布明显偏态（多数重抽样样本权重集中、少数高权重站点主导），岭回归R²也低
  很多，结论要比Promote组更谨慎对待。
- 产出：`data_proc/output_loose_classify/recovery_omega_178_bydirection_{data,inference,importance}.csv`、
  `recovery_omega_178_bydirection_coef_ci.png`。

## 二、optimal_tp（CCM效应滞后）driver分析

**因变量**：optimal_tp——CCM里"效应滞后几个tp(每个8天)达到峰值"，1个。参照
Liu et al. 2024(*Agric For Meteorol* 358:110239)"扫描滞后取峰值"设计。

**自变量组合(5种)**：18变量全量→15变量(去稀疏)→16变量(15+urban_rate)→
10变量(理论+数据双重精简)→11变量(10精简+invest,用户拍板定稿版)，逐步迭代。

**分析方法(7种)**：OLS；OLS+Conley空间HAC标准误(处理站点空间自相关)；
随机森林+置换重要性；弹性网(10折CV)；bootstrap稳定性选择(500次重抽样筛变量,
Meinshausen & Bühlmann 2010思路)；有序logit(MASS::polr)；负二项GLM(处理0-8档
计数型因变量的过离散问题)。

**站点集(2种)**：117、178。

**约20+次独立拟合(6轮迭代)**。**结果要点**：**本项目所有因变量里信号最弱、最不
稳定的一个**——18变量/86站 RF R²=0.012(几乎无解释力)；精简到10变量后RF R²
仅升到0.13；178站(放宽站集)信号不升反降(0.027)；6个变量集配置里没有一个变量
能跨配置稳定复现p<0.05显著性；有序logit因9档分类准完全分离拟合失败作废。
**结论**：现有18个站点静态协变量解释不了optimal_tp的站间差异，不是方法问题。
**已知未填补的缺口**：这条线始终没有把站级VPD均值/变异加进自变量表（后来在
φ/Omega分析里补了，但optimal_tp这条线还没有补做），若信号仍弱，需另外考虑
是否本质上就缺乏站点层面可读的驱动因素。

脚本：`tp_lag_drivers_117.R`、`_117_final.R`、`_117_scale_check.R`、`_178.R`、
`_summary_all.R`、`tp_lag_var_screening_117.R`。
产出目录：`data_proc/output_loose_classify/tp_lag_drivers_*`。

## 三、历史支线（仅供参考，不同站点口径，未在当前体系重做）

`pipelines/loose_classify_drivers.R`——在切换到现在117/178统一站点清单**之前**，
基于旧的383站"宽松分类"口径，对`rho_peak`(CCM因果强度)和`peak_tp`做过一次
RF+OLS driver分析，发现干旱/高辐射是强度的头号预测因子。这是不同站点口径下的
旧结果，未在当前117/178体系里重新验证，如需要可补一版。

## 四、总体统计

| 维度 | 数量 |
|---|---|
| 已测试因变量 | 9个（7个恢复力/抵抗力类 + optimal_tp + 历史rho_peak/peak_tp） |
| 已测试自变量组合 | 约8种 |
| 已测试分析方法 | 约9种 |
| 已测试站点集口径 | 3种（117/178/57） |
| 累计独立模型拟合次数(估算) | 约60次 |

## 五、目前最值得信赖的结论（跨方法/跨口径复现过的）

1. **暖湿气候(降水/温度/湿度)对"恢复持续性"和"抵抗力"的方向相反**——φ为正(越暖湿，
   异常持续越久=恢复越慢)，logOmega也为正(越暖湿，事件期偏离越小=抵抗力越强)，
   两者各自在多套自变量集、OLS与岭回归两种方法下复现。(2026-10-06更正：此前误把
   logOmega方向当成"越大越差"，写成了"两个指标方向一致地不利"。)
2. **Promote方向站点独立地抵抗力更差**（logOmega模型的is_promote系数），与
   生存分析"Promote组恢复更慢"的结论方向一致，是跨方法复现次数最多的发现。
3. **VPD站级统计量对φ有用、对logOmega和optimal_tp暂未证实有用**——尚未在
   optimal_tp上测试VPD协变量，这是明确的下一步。
4. **离散型/计数型指标(optimal_tp)普遍比连续型指标(φ、logOmega)更难测出稳定
   信号**——不确定这是数据本身的限制还是这类指标信噪比结构性更差，值得作为
   方法论讨论点写入论文。

## 六、抛开CCM的全站点分析：干旱-植物生理关系及其调节因素（2026-10-07）

**设计**（参考De Keersmaecker 2015 ARX + Jiao 2019两阶段）：不做CCM筛选，用协变量齐全的904站(883站有≥100个有效点)。
- 第一阶段(`pipelines/drought_phys_stage1_arx_allstations.R`)：逐站ARX，SIF距平 = φ·上期距平 + β·VPD距平 + 气温/辐射/降水距平，
  全部为去季节后标准化(SD单位)，年内滞后。β=控制气温/辐射/降水后VPD对SIF的即期效应。结果：β<0占90%，p<0.05显著为负55.5%、
  显著为正0.6%；中位β=-0.17；β标准误中位0.076，噪声占比0.36；VPD的VIF中位3.2(最大13.9)。**只控制VPD(旧做法)时β<0仅48%，
  与控制气温/辐射后的β相关仅0.39**，提示此前"约35%站点VPD促进SIF"的方向分类很可能与气温/辐射混杂有关(仅为线索，未验证)。
  空间上北方/西北干旱半干旱区最敏感，南方湿润区接近0。
- 第二阶段(`pipelines/drought_phys_stage2_moderators.R`)：S=-β(越大=VPD抑制越强)对VAR10+VPD均值/CV+invest回归，权重1/(se²+τ²)，
  城市聚类(259城)。WLS+城市聚类稳健SE、岭回归+城市整簇bootstrap、XGBoost+SHAP(空间分块CV，invest缺失保留)。
- **结果**：XGBoost空间分块CV加权R²=0.266(RF 0.281，岭回归0.216，置换p=0.048最小值)，显著高于此前φ/面板模型。
  稳定调节因素：precip(负，SHAP占37.6%，bootstrap 100%第1，约在5-6处阈值式下降)、vpd_cv_stn(正)、vpd_mean_stn(正)、rsds_sd(正)；
  岭回归CI不跨0的还有ntl(负)、ever_broad(负)、ever_needle(负)。
- **invest**：WLS系数-0.0017(聚类SE 0.0064, p=0.79)，岭回归+0.0026(CI -0.006~0.011)，Spearman(S,invest)=0.066，
  按三分位S均值0.152/0.188/0.171(非单调)，XGBoost名次第11(SHAP占2.5%)，去掉invest后R²反而+0.008，**无调节作用**；
  CI上限说明每增加1个SD的invest，S的变化不超过约0.01(≈0.08个S的SD)。
产出：`drought_phys_arx_station_beta.csv`、`drought_phys_stage2_{wls_cluster,ridge_boot,xgb_cv,xgb_importance}.csv`、
`drought_phys_{beta_hist,beta_map}.png`、`drought_phys_stage2_{ridge_ci,shap_beeswarm,shap_dependence}.png`。

## 七、"投资是否缓解干旱对植物的影响"：全方案检验（2026-10-07）

**做法**：不挑方案，预设网格全部公开，并用置换检验和留出验证校正"挑最好"的偏差。
- 站点水平网格(`pipelines/drought_phys_mitigation_spec_curve.R`，216个方案)：6种"干旱影响"度量(S_partial=-β[VPD,控制气温/辐射/降水]、
  S_vpdonly、S_spei[ERA5 SPEI-1]、φ、R_loss、Z_worse；统一为越大越严重，数据来自`drought_phys_stage1b_extended.R`)×3种投资度量
  (绝对投资inv_mean、对气候残差化inv_resid、单位绿地面积投资inv_green)×4套控制变量×3个站点子集(全部/干旱半区/高不透水面半区)；
  WLS+城市聚类SE，缓解效应=-coef(invest)(SD单位，>0=缓解)。
- 一阶段面板(`drought_phys_mitigation_panel.R`)：SIF距平~干旱距平×投资(站间+站内年度变化)，年份FE，城市聚类；
  干旱=VPD或-SPEI，子集=全部/干旱半区，投资=当年/前一年。
**结果**：
- 216方案：缓解方向仅21个(10%)，显著缓解(t>1.96)**0个**，显著反方向38个，效应中位数-0.052 SD(范围-0.214~+0.103)。
  最支持缓解的方案=S_spei×inv_green×干旱半区(n=393站/132城, +0.103 SD, t=1.55, p=0.12，不显著)。
- 置换检验(500次，降水三分位内整城打乱投资)：零分布下"最大t"均值2.34(95%分位3.64)，观测最大t=1.55，p=0.87；零分布平均有8.9个方案显著缓解，
  观测0个；t中位数观测-1.04 vs 零分布0.14(p=0.99，整体偏向反方向)。
- 留出验证(300次，一半城市选最优方案，另一半检验)：发现集t均值1.95→留出集-0.33；留出集显著缓解仅1%，显著反方向5%，缓解方向39%。
- 面板16个交互系数：缓解方向9个，p<0.05共5个(缓解2、反方向3)，Holm校正后仅1个显著且为**反方向**(站间、VPD、全部站点：-0.0143，p=0.001，Holm 0.019；
  即高投资城市的VPD→SIF斜率更负)。名义上最支持缓解的是站内、VPD、干旱半区、当年投资：+0.0087(p=0.021，Holm 0.29；基线斜率-0.243，约削弱3.6%)，
  但用上一年投资(+0.0007, p=0.87)和用SPEI(≈0)均不复现。
- 注意：R_loss/Z_worse(最高VPD窗口的SIF损失)中位数≈0，VPD最高窗口往往也是气温/辐射最高、光合旺盛的时候，这两个指标不能干净度量干旱影响。
**结论**：没有任何一类方案能为"投资缓解干旱影响"提供可信支持；最强的显著结果反而是反方向(关联而非因果，可能与城市规模/气候/反向因果混杂)。
产出：`drought_phys_mitigation_{spec_curve.csv,spec_curve.png,heatmap.png,perm_null.csv,holdout.csv,panel.csv,panel_terms.csv}`。
