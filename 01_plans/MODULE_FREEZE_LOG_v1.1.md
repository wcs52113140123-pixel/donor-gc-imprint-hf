# 模块冻结记录 v1.1（2026-10-03T05:36:57Z）

**文件**：`01_config/modules/modules_v1.1.tsv`，sha256 `3a2abdad2b609002752b4d4833d8219114fd305b2339c7ddf406da86448591be`（冻结后设为只读）
**代码**：`03_code/04_models/pbtools_v11.py`（修正 B1、B6）；签名推导脚本 `03_code/04_models/sig_v11.py`
**依据**：审计报告 `00_admin/审计报告_既有分析管线_2026-10-03.md`；叙事方向 A+C（用户 2026-10-03 确认）

## 使用的数据（发现集，仅此两项）
- SCP3342：HFpEF 15 例 vs 对照 24 例，巨噬细胞伪 bulk（CellBender），基因层面效应
- GSE290367 成人左室髓系：DCM 67 例 vs ND 60 例；HCM 15 例 vs ND 16 例（HCM 只用来报告一致性）

## 冻结时尚未打开的验证集（v1.1 签名在其中从未计算过）
Hahn 2021 心内膜活检、MAGNet（GSE141910）、PXD060431、E-MTAB-12975 Visium、Kuppe 2022 Visium、PERIHEART。
**需要披露**：GTEx v11 和 Age-OME 在 v1.0 阶段已经用 `E_efferocytosis_receptors`（与 RES_EFF 有 7 个基因重叠）做过年龄分析，所以 GTEx 和 Age-OME 中的衰老检验**不是完全盲法**。

## 新模块
| 模块 | 层级 | 基因 | 规则 |
|---|---|---|---|
| RES_EFF_core | 主终点 | C1QA C1QB C1QC CD163 F13A1 LRP1 LYVE1 MARCO MERTK  | 来自文献候选池；在 SCP3342 HFpEF 和 GSE290367 DCM 中均下调，p<0.05，且 SCP3342 中 logCPM≥3（规则在运行前写在脚本头部） |
| MHCII_core | 共同主终点 | HLA-DRA HLA-DRB1 HLA-DPA1 HLA-DPB1 HLA-DQA1 HLA-DQB1 HLA-DMA HLA-DMB CD74 CIITA  | 固定的文献列表，不做筛选 |
| RES_EFF_lit | 敏感性分析 | 候选池全部 32 个基因 | 不做数据筛选 |
| ISG_IFN | 探索性 | CXCL9 CXCL10 GBP1 GBP2 GBP5 STAT1 IRF1 IDO1  | 在发现集中方向不一致（事后从原 MHCII_IFN 列表中按基因类别拆出，已披露） |
| MAC_PAN_ID | 调整用 | AIF1 CD14 CD68 CYBB FCGR1A MS4A6A SPI1 TYROBP  | 规则：在所有发现对比中 \|β\|≤1；用于 bulk 反卷积和身份校正。CSF1R、LY86、FCER1G 因随疾病变化被剔除 |

**复合指标（次要终点）**：EFF–MHCII 状态指数 = z(RES_EFF_core) − z(MHCII_core)，z 以各数据集内的参照组标准化。

## 发现集中值得注意的现象（不作为验证证据）
- 常驻身份基因 FOLR2、TIMD4、SIGLEC1 无变化；MAF、SLC40A1、SELENOP、CD163L1 上调。这提示变化是"状态"而不是常驻巨噬细胞整体丢失，需要在亚状态分析中检验。
- 本签名在 SCP3342 和 GSE290367 DCM 中的效应是**循环的**（由这两组数据选出），不能当作验证结果报告。

## 预先指定的验证假设
- H1 衰竭：在独立衰竭数据（Hahn HFpEF/HFrEF、MAGNet DCM/HCM、PXD060431 蛋白）中，RES_EFF_core↓、MHCII_core↑，并在校正巨噬细胞含量（MAC_PAN_ID 或反卷积）后仍然成立。
- H2 衰老：在非衰竭/健康供体中，RES_EFF_core 和 MHCII_core 的年龄斜率接近 0；以等效性界限 ±0.1 SD/10 年为标准，用 TOST 检验（数据集：MAGNet NF、GSE290367 ND、GTEx〔非盲〕）。
- H3 组成 vs 状态：在 GSE290367 Resident 和 MonoDerived 亚群内分别检验疾病效应。
- 每个假设一个检验族，族内 BH-FDR。

## 代码修正（相对 v1.0）
- B1：HC3 在杠杆值 h≥0.999 时改用 HC2 并打标记。
- B2：阈值按"供体×细胞类型"设定；间质为 ≥50 个核、≥10k UMI，髓系在两个数据集中统一为 ≥20 个核、≥5k UMI。
- B6：竞争性检验的随机基因无放回抽样，且排除模块基因。
- v1.0 的文件和结果保持不变。
