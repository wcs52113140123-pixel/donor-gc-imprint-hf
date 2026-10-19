# 预注册分析计划 v1.3：替代解释检验、影响审计与校正工具

**冻结时间**：见文末 UTC 时间戳；sha256 写在服务器上的只读副本旁边。
**上游文件**：v1.2 计划（sha256 26db5075…02b2）及其结果（`00_admin/v1.2_预注册检验结果_2026-10-03.md`）、`modules_v1.1.tsv`（3a2abdad…）、体外 GC 参考集 `gc_gene_classification.csv` 与 `meta_all_genes.csv.gz`（4 个通过 QC 的数据集）。
**本计划写作前已经查看的数据**：
- 心脏数据：v1.2 的全部结果；SCP3342 各区室 HFpEF 全基因差异结果（作为 GC 富集初查，已写入《选题充分性评估》报告）；文献中报告的基因。
- 体外参考集：昼夜钟基因和 cAMP 基因的体外分类结果，仅用于本计划第 1 节的基因筛选。
- **尚未查看**：任何心脏数据中的钟基因或 cAMP 基因评分；校正后的全基因差异结果；阳性对照基因在校正前后的结果。

## 1. 固定的评分定义（只用体外数据筛选，未参考心脏数据）
统一筛选规则：基因必须**未被**体外参考集判为 GC 诱导或 GC 抑制，且合并 |log2FC| < 0.5；体外参考集中未测到的基因可以纳入。
- **GC_specific3** = FKBP5、TSC22D3、ZBTB16。这 3 个是体外诱导最强、也最常被引用的 GC 靶基因；PER1、KLF9、DDIT4 因同时受钟或 cAMP 调控而排除。
- **CLOCK_PHASE**（昼夜相位指数）= mean z(晨峰臂) − mean z(夜峰臂)。
  - 晨峰臂：NR1D2、DBP、TEF、PER2、BHLHE41、CIART。
  - 夜峰臂：ARNTL、NPAS2、CLOCK、CRY1。
  - 已排除：PER1（GC 诱导）、PER3（GC 抑制）、BHLHE40（GC 抑制）、NFIL3（GC 诱导）、NR1D1（|log2FC| = 1.05）、CRY2（相位归属不明确）。
- **CAMP_IEG**（cAMP/儿茶酚胺即早基因）= NR4A1、NR4A2、CREM、FOSB、JUNB、PDE4B、PDE4D。
  - 已排除：FOS（GC 诱导）、NR4A3（|log2FC| = 2.45）、ATF3 和 EGR1（GC 抑制）、SIK1（|log2FC| = 0.67）。
- 打分方式与 v1.1/v1.2 相同：以参照组为基准的 z 分数（pbtools.ref_z / module_score），逐细胞类型的伪 bulk，阈值同 v1.2。
- 数据集与对比同 v1.2：SCP3342 HFpEF vs 对照（pool 作区组）；GSE290367 DCM vs ND（去重、只用含 DCM 的研究、研究作固定效应、校正年龄和性别）；MAGNet 剔除重叠样本后的 DCM vs NF（校正年龄和性别）。

## 2. 检验 C：昼夜节律能否解释 GC 印记（判定）
- **C1（描述）**：在每个数据集、每个细胞类型中，报告 CLOCK_PHASE 的病例−对照差异。
- **C2（判定）**：在 SCP3342 和 GSE290367 DCM 的 5 个细胞类型中，GC_specific3 的"对照−病例"差异分别用两个模型估计：(a) 原模型；(b) 原模型 + CLOCK_PHASE。
  - 计算每个细胞类型的保留比例 = β(b) / β(a)，只在 β(a) > 0 的细胞类型中计算。
  - 同时检验巨噬细胞 G4（GC 诱导子集，v1.2 定义）在 GSE290367 Resident DCM 和 MAGNet 剔除重叠 DCM 中的效应，模型为原模型 + CLOCK_PHASE。
  - **通过**（昼夜节律不足以解释）需同时满足：两个数据集各自的保留比例中位数 ≥ 0.5；G4 在两个对比中都是 β < 0 且 P < 0.05。
  - **失败的后果**：在稿件中把"取材时间/昼夜相位"列为与 GC 并列、无法区分的解释，并相应弱化所有暴露相关表述。
- **C3（描述）**：在对照组内，报告 GC_specific3 与 CLOCK_PHASE 的 Spearman ρ。

## 3. 检验 K：儿茶酚胺/cAMP 能否解释 GC 印记（判定）
- **K1（判定）**：与 C2 相同，只是协变量改为 CAMP_IEG，通过标准也与 C2 相同。
- **K2（描述）**：在每个细胞类型中，比较 GC_specific3 与 CAMP_IEG 区分对照和病例的 AUC。
- **CK（敏感性）**：同时加入 CLOCK_PHASE 和 CAMP_IEG，报告保留比例。
- **失败的后果**：同 C2，解释改为"应激激素（GC/儿茶酚胺）混合印记"。

## 4. 检验 A：对已有结论的影响审计（描述性，预先指定指标）
- **对比**：
  - SCP3342 的 Cardiomyocyte、Fibroblast、Endothelial、Pericyte、Myeloid、VascStromal；
  - GSE290367 DCM 的同 5 个 L1 细胞类型；
  - MAGNet 剔除重叠 DCM vs NF；
  - Hahn 2021 RVS HFpEF vs 对照。
- **逐基因模型**：在参照组标准化的基因 z 值上做 OLS（设计同上），BH-FDR。
- **A1**：
  - FDR < 0.05 的下调基因中 GC 诱导基因的富集（Fisher OR），上调基因中 GC 抑制基因的富集；
  - Top-100 下调和 Top-100 上调基因中 GC 基因的比例。
- **A2 校正**：在设计中加入 GC_specific3。报告校正前后 FDR < 0.05 基因数、两组基因的 Jaccard 指数，以及校正后的 A1 指标。校正分数本身的 3 个基因不计入统计。
- **A3 阳性对照（判定校正工具是否"不伤真信号"）**：在 MAGNet 剔除重叠 DCM 中检验 NPPA↑、NPPB↑、POSTN↑、ASPN↑、SFRP4↑、FRZB↑、MYH6↓。
  - **通过**：校正后 ≥ 6/7 个基因保持预期方向且 FDR < 0.05。
  - 同时报告校正前的结果，以及 GSE290367 心肌细胞（NPPA、NPPB、MYH6）和成纤维细胞（POSTN、ASPN、SFRP4、FRZB）中的结果（描述）。
- **A4（描述）**：把文献中被解读为心衰或恢复生物学的基因（FKBP5、ZBTB16、TSC22D3、DDIT4、KLF9、SGK1、PDK4、FOXO3、PTPN13 等）与体外 GC 分类逐一对照，列成表。

## 5. 校正工具 T
- **内容**：
  - (i) 体外 GC 参考集的全基因分类表（诱导、抑制、不敏感、不确定）；
  - (ii) `gc_score()`：计算 GC_specific3 和 GC_response 评分；
  - (iii) `gc_adjusted_de()`：在差异分析设计中加入 GC 协变量；
  - (iv) `gc_flag()`：标记结果中的 GC 敏感基因。
- **验证**：以 A3 的通过为准。A3 失败时，工具只作为"标记"用途发布，不推荐作协变量校正。

## 6. 多重比较与偏离规则
- C2 和 K1 是两个独立的判定检验，各自按上面的规则判定，不再做跨检验校正。
- A 部分全部是描述性的，报告 BH-FDR。
- 不允许在看到结果后修改基因列表、阈值或判定规则。
- 任何偏离都写进结果报告和稿件方法部分。
- LVAD 前后配对数据（GSE226314）**不在本计划内**，属于后续探索性分析。

---
冻结时间（UTC）：2026-10-03T15:52:44Z
