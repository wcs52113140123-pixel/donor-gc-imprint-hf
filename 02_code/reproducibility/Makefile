# =====================================================================================================================
# v1.2 reproducibility master pipeline  (03_code/07_v12/reproducibility/Makefile)
# Regenerates, in dependency order, every result table and figure used by manuscript v0.1 from the processed
# pseudobulk / raw files, by calling the ORIGINAL scripts in 03_code/{04_models,05_validation,06_synthesis} unmodified
# wherever their paths derive from PROJ_ROOT/argv/cwd, and path-only wrapped copies (wrapped/, see WRAP_DIFFS.patch)
# where a script hard-codes the project root. Outputs go to 04_results/rerun_v12/ only.
# Interpreters: v1.0/v1.1/synthesis scripts -> system python3 (as originally run); proteomics and SCP3342 substates -> venv
# (harmonypy 2.0.2, verified to reproduce the archived label transfer byte-identically); synthesis UMAP -> system python3 +
# harmonypy 0.0.10 in $(RR)/pylib (as originally run).
#   make all          v1.0 + v1.1 + synthesis + figures + v1.0-with-v1.1-fixes + comparison tables
#   make pseudobulk_check   rebuild both pseudobulks from raw h5ad and compare with 02_data/processed (slow, I/O bound)
# Threads: each recipe runs with $(NT) BLAS threads; run with  make -j$(J)  so that J*NT <= 32 (shared machine policy).
# =====================================================================================================================
ROOT ?= /home/weichangsheng/projects/cardiac_ageing_target_clearance
HERE := $(ROOT)/03_code/07_v12/reproducibility
RR   ?= $(ROOT)/04_results/rerun_v12
SH   := $(RR)/root
SHF  := $(RR)/root_v11fix
SHD  := $(RR)/root_v11fix_dedup
S    := $(RR)/stamps
L    := $(RR)/logs
W    := $(RR)/work
PY   ?= /usr/bin/python3
VPY  ?= /home/weichangsheng/envs/cardiac_ageing/bin/python
NT   ?= 8
THR  := OMP_NUM_THREADS=$(NT) OPENBLAS_NUM_THREADS=$(NT) MKL_NUM_THREADS=$(NT) NUMBA_NUM_THREADS=$(NT)
C    := $(ROOT)/03_code
R4   := $(SH)/04_results
V11  := $(R4)/validation_v1.1
OUTD ?= $(ROOT)/04_results/validation_v1.2/reproducibility
SHELL := /bin/bash
.SHELLFLAGS := -eo pipefail -c

.PHONY: all v10 v11 synthesis derived figures v11fix compare pseudobulk_check clean_stamps
all: v10 v11 synthesis derived figures v11fix compare
v10: $(S)/v10_fig2 $(S)/v10_drift $(S)/v10_fig4 $(S)/v10_hfpef $(S)/v10_audit $(S)/v11_sig
v11: $(S)/v11_gc $(S)/v11_bulk $(S)/v11_bulkgc $(S)/v11_gse_sub $(S)/v11_scp_sub $(S)/v11_overlap $(S)/v11_prot $(S)/v11_visium
synthesis: $(S)/syn_main $(S)/syn_umap
derived: $(S)/der_substates $(S)/der_spatial $(S)/der_prot $(S)/der_syn
figures: $(S)/figures
v11fix: $(S)/fix_fig2 $(S)/fix_hfpef $(S)/fix_fig4 $(S)/fixd_fig2 $(S)/fixd_hfpef $(S)/fixd_fig4

$(S)/setup:
	mkdir -p $(S) $(L) $(W); ROOT=$(ROOT) RR=$(RR) bash $(HERE)/setup_shadow.sh > $(L)/setup.log 2>&1; touch $@

# ------------------------------------------------ v1.0 layer (original code, original pbtools) ------------------------
$(S)/v10_fig2: $(S)/setup
	cd $(SH) && PROJ_ROOT=$(SH) $(THR) $(PY) $(C)/04_models/fig2_normative_coupling.py > $(L)/v10_fig2.log 2>&1; touch $@
$(S)/v10_drift: $(S)/setup
	cd $(SH) && PROJ_ROOT=$(SH) $(THR) $(PY) $(C)/04_models/drift_diagnosis.py > $(L)/v10_drift.log 2>&1; touch $@
$(S)/v10_fig4: $(S)/setup
	cd $(SH) && PROJ_ROOT=$(SH) $(THR) $(PY) $(C)/04_models/fig4_disease_specificity.py > $(L)/v10_fig4.log 2>&1; touch $@
$(S)/v10_hfpef: $(S)/v10_fig2
	cd $(SH) && PROJ_ROOT=$(SH) $(THR) $(PY) $(C)/04_models/hfpef_scp3342.py $(RR)/frozen/gene_age_effects_GSE290367.csv.gz > $(L)/v10_hfpef.log 2>&1; touch $@
# audit / signature derivation were run from a staged cwd: audit/ = project mirror, pb/G pb/S = pseudobulk dirs
$(S)/v10_audit: $(S)/setup
	rm -rf $(W)/audit && mkdir -p $(W)/audit/pb && cd $(W)/audit && ln -s $(SH) audit && ln -s $(ROOT)/02_data/processed/GSE290367_pseudobulk pb/G && ln -s $(ROOT)/02_data/processed/SCP3342_pseudobulk pb/S && \
	$(THR) $(PY) $(C)/04_models/audit_controls_2026-10-03.py > $(L)/v10_audit.log 2>&1 && mkdir -p $(R4)/audit_2026-10-03 && cp audit_controls_results.json audit_out_*.csv $(R4)/audit_2026-10-03/; touch $@
$(S)/v11_sig: $(S)/v10_hfpef
	rm -rf $(W)/sig && mkdir -p $(W)/sig/pb && cd $(W)/sig && ln -s $(SH) audit && ln -s $(ROOT)/02_data/processed/GSE290367_pseudobulk pb/G && ln -s $(ROOT)/02_data/processed/SCP3342_pseudobulk pb/S && \
	$(THR) $(PY) $(C)/04_models/sig_v11.py > $(L)/v11_sig.log 2>&1 && mkdir -p $(R4)/audit_2026-10-03 && cp sig_v11_*.csv $(R4)/audit_2026-10-03/; touch $@

# ------------------------------------------------ v1.1 validation layer ---------------------------------------------
$(S)/v11_gc: $(S)/setup
	rm -rf $(W)/gc && mkdir -p $(W)/gc/pb && cd $(W)/gc && ln -s $(SH) audit && ln -s $(ROOT)/02_data/processed/GSE290367_pseudobulk pb/G && ln -s $(ROOT)/02_data/processed/SCP3342_pseudobulk pb/S && \
	cp $(ROOT)/01_config/modules/modules_v1.1.tsv $(V11)/donor_overlap/gse290367_duplicate_individuals.csv . && \
	$(THR) $(PY) $(C)/05_validation/gc/gc_discovery.py > $(L)/v11_gc.log 2>&1 && mkdir -p $(V11)/gc_discovery && cp gc_discovery_results.json gc_*.csv $(V11)/gc_discovery/; touch $@
$(S)/v11_bulk: $(S)/setup
	cd $(SH) && PROJ_ROOT=$(SH) INP=$(C)/05_validation/bulk $(THR) $(PY) $(C)/05_validation/bulk/bulk_validation_v11.py $(V11)/bulk > $(L)/v11_bulk.log 2>&1; touch $@
$(S)/v11_bulkgc: $(S)/setup
	cd $(SH) && PROJ_ROOT=$(SH) INP=$(C)/05_validation/bulk $(THR) $(PY) $(C)/05_validation/bulk/bulk_gc_diagnostic.py $(V11)/bulk/posthoc_gc > $(L)/v11_bulkgc.log 2>&1; touch $@
$(S)/v11_gse_sub: $(S)/setup
	cd $(SH) && PROJ_ROOT=$(SH) $(THR) $(PY) $(HERE)/wrapped/gse_substates.py > $(L)/v11_gse_sub.log 2>&1; touch $@
$(S)/v11_scp_sub: $(S)/setup
	cd $(SH) && PROJ_ROOT=$(SH) $(THR) $(VPY) $(HERE)/wrapped/scp_substates.py > $(L)/v11_scp_sub.log 2>&1; touch $@
$(S)/v11_overlap: $(S)/setup
	rm -rf $(W)/ov && mkdir -p $(W)/ov && cd $(W)/ov && ln -s $(R4) 04_results && ln -s $(ROOT)/02_data 02_data && cp $(V11)/donor_overlap/MAGNet_phenoData.csv . && \
	$(PY) $(C)/05_validation/donor_overlap/overlap_matching.py > $(L)/v11_overlap.log 2>&1 && cp donor_overlap_tuple_matching*.csv donor_overlap_shared_ids.csv donor_overlap_scp3342_id_null.csv $(V11)/donor_overlap/; touch $@
$(S)/v11_prot: $(S)/setup
	mkdir -p $(V11)/protein_spatial/proteomics && cd $(W) && $(THR) $(VPY) $(C)/05_validation/protein_spatial/proteomics_v11.py $(ROOT)/02_data/raw/PXD060431/Full_dataset-3_FOR_PRIDE_1-31-2025.csv \
	$(ROOT)/02_data/raw/PXD060431/PXD060431_community_annotated.sdrf.tsv $(ROOT)/01_config/modules/modules_v1.1.tsv $(V11)/protein_spatial/proteomics > $(L)/v11_prot.log 2>&1 && \
	cp $(RR)/frozen/prot_grep.txt $(V11)/protein_spatial/proteomics/prot_grep.txt && \
	cp $(V11)/protein_spatial/proteomics/proteomics_detection.csv $(V11)/protein_spatial/proteomics_detection.csv; touch $@
$(S)/v11_visium: $(S)/setup
	mkdir -p $(V11)/protein_spatial/spatial && cd $(W) && $(THR) $(PY) $(C)/05_validation/protein_spatial/visium_v11.py $(SH) $(ROOT)/01_config/modules/modules_v1.1.tsv $(V11)/protein_spatial/spatial > $(L)/v11_visium.log 2>&1; touch $@

# ------------------------------------------------ sandbox post-processing (code recovered from artifact lineage) -------
$(S)/der_substates: $(S)/v11_gse_sub $(S)/v11_scp_sub
	mkdir -p $(V11)/substates/final && cd $(V11)/substates/final && for f in substates_GSE290367_effects GSE290367_decomposition_shares substates_SCP3342_effects; do \
	IN_DIR=$(V11)/substates OUT_DIR=$(V11)/substates/final $(PY) $(HERE)/wrapped/sandbox_derived/derive_$$f.py >> $(L)/der_substates.log 2>&1; done; touch $@
$(S)/der_spatial: $(S)/v11_visium
	for f in spatial_niche spatial_niche_by_section spatial_scores_by_sample spatial_zone_tests; do \
	IN_DIR=$(V11)/protein_spatial/spatial OUT_DIR=$(V11)/protein_spatial $(PY) $(HERE)/wrapped/sandbox_derived/derive_$$f.py >> $(L)/der_spatial.log 2>&1; done; touch $@
$(S)/der_prot: $(S)/v11_prot
	IN_DIR=$(V11)/protein_spatial/proteomics OUT_DIR=$(V11)/protein_spatial $(PY) $(HERE)/wrapped/sandbox_derived/derive_proteomics_effects_final.py > $(L)/der_prot.log 2>&1; touch $@

# ------------------------------------------------ synthesis --------------------------------------------------------
$(S)/syn_main: $(S)/v11_scp_sub
	cd $(SH) && PROJ_ROOT=$(SH) $(THR) $(PY) $(HERE)/wrapped/synthesis_reanalysis.py > $(L)/syn_main.log 2>&1; touch $@
$(S)/syn_umap: $(S)/syn_main
	cd $(SH) && PROJ_ROOT=$(SH) PYTHONPATH=$(RR)/pylib $(THR) $(PY) $(HERE)/wrapped/synthesis_scp_umap.py > $(L)/syn_umap.log 2>&1; touch $@

# sandbox post-processing of synthesis outputs (per-study GC table + combined wide results file)
$(S)/der_syn: $(S)/syn_umap
	cd $(R4)/synthesis_v1.1 && IN_DIR=$(R4)/synthesis_v1.1 OUT_DIR=$(R4)/synthesis_v1.1 $(PY) $(HERE)/wrapped/sandbox_derived/derive_A3_GSE290367_DCM_by_study_GC.py > $(L)/der_syn.log 2>&1 && \
	IN_DIR=$(R4)/synthesis_v1.1 OUT_DIR=$(R4)/synthesis_v1.1 $(PY) $(HERE)/wrapped/sandbox_derived/derive_synthesis_results_combined.py >> $(L)/der_syn.log 2>&1; touch $@

# ------------------------------------------------ figures (original fig*.py via shim) --------------------------------
V10S := $(S)/v10_fig2 $(S)/v10_drift $(S)/v10_fig4 $(S)/v10_hfpef $(S)/v10_audit $(S)/v11_sig
V11S := $(S)/v11_gc $(S)/v11_bulk $(S)/v11_bulkgc $(S)/v11_gse_sub $(S)/v11_scp_sub $(S)/v11_overlap $(S)/v11_prot $(S)/v11_visium
FIXS := $(S)/fix_fig2 $(S)/fix_hfpef $(S)/fix_fig4 $(S)/fixd_fig2 $(S)/fixd_hfpef $(S)/fixd_fig4
$(S)/figures: $(S)/der_syn $(S)/der_substates $(S)/der_spatial $(S)/der_prot $(V10S) $(V11S)
	mkdir -p $(R4)/synthesis_v1.1/figures && for n in 1 2 3 4 5 6; do \
	MPLBACKEND=Agg $(PY) $(HERE)/wrapped/figures/figrun.py $(C)/06_synthesis/figures/fig$$n.py $(R4) $(ROOT)/01_config $(W)/fig$$n >> $(L)/figures.log 2>&1 && \
	cp $(W)/fig$$n/figures/* $(R4)/synthesis_v1.1/figures/; done; touch $@

# ------------------------------------------------ v1.0 analyses with v1.1 fixes (B1 HC2, B6 no-replacement, B2 thresholds)
define FIXRULES
$(S)/$(1)_fig2: $(S)/setup
	cd $(2) && PROJ_ROOT=$(2) $(3) $(THR) $(PY) $(HERE)/wrapped/v11fix/fig2_normative_coupling_v11fix.py > $(L)/$(1)_fig2.log 2>&1; touch $$@
$(S)/$(1)_hfpef: $(S)/$(1)_fig2
	cd $(2) && PROJ_ROOT=$(2) $(3) $(THR) $(PY) $(HERE)/wrapped/v11fix/hfpef_scp3342_v11fix.py $(RR)/frozen/gene_age_effects_GSE290367.csv.gz > $(L)/$(1)_hfpef.log 2>&1; touch $$@
$(S)/$(1)_fig4: $(S)/setup
	cd $(2) && PROJ_ROOT=$(2) $(3) $(THR) $(PY) $(HERE)/wrapped/v11fix/fig4_disease_specificity_v11fix.py > $(L)/$(1)_fig4.log 2>&1; touch $$@
endef
$(eval $(call FIXRULES,fix,$(SHF),V12_DEDUP=0))
$(eval $(call FIXRULES,fixd,$(SHD),V12_DEDUP=1 V12_DUPFILE=$(V11)/donor_overlap/gse290367_duplicate_individuals.csv))

# ------------------------------------------------ comparison -------------------------------------------------------
$(S)/compare: $(S)/figures $(FIXS)
	mkdir -p $(OUTD) && $(PY) $(HERE)/compare/compare_outputs.py $(ROOT)/04_results $(R4) $(OUTD)/rerun_diff_table.csv > $(L)/compare.log 2>&1 && \
	$(PY) $(HERE)/compare/v10_vs_v11fix.py $(ROOT)/04_results $(SHF)/04_results $(SHD)/04_results $(OUTD)/v10_vs_v11fix_endpoints.csv >> $(L)/compare.log 2>&1 && \
	$(PY) $(HERE)/compare/annotate_diff.py $(OUTD)/rerun_diff_table.csv >> $(L)/compare.log 2>&1 && \
	$(PY) $(HERE)/compare/key_endpoints.py $(OUTD)/v10_vs_v11fix_endpoints.csv $(OUTD)/v10_vs_v11fix_key_endpoints.csv >> $(L)/compare.log 2>&1; touch $@
compare: $(S)/compare

# ------------------------------------------------ pseudobulk rebuild (optional, ~100 GB read) ----------------------------
pseudobulk_check: $(S)/pb_check
$(S)/pb_check: $(S)/setup
	mkdir -p $(RR)/pseudobulk && cd $(ROOT) && \
	$(PY) $(C)/03_pseudobulk/make_pseudobulk_scp3342.py 02_data/raw/SCP3342/HFpEF_snRNAseq_single_cell_portal_withRawCounts_09.24.2026.h5ad $(RR)/pseudobulk/SCP3342_pseudobulk > $(L)/pb_scp.log 2>&1 && \
	$(PY) $(C)/02_pseudobulk/make_pseudobulk.py 02_data/interim/GSE290367_integrated_RNA_adata.h5ad $(RR)/pseudobulk/GSE290367_pseudobulk > $(L)/pb_gse.log 2>&1 && \
	$(PY) $(HERE)/compare/compare_pseudobulk.py $(ROOT)/02_data/processed $(RR)/pseudobulk $(OUTD)/pseudobulk_rebuild_check.csv > $(L)/pb_check.log 2>&1; touch $@

clean_stamps:
	rm -f $(S)/*
