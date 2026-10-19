#!/usr/bin/env bash
# One-command rerun of every manuscript result table/figure into 04_results/rerun_v12/ (see Makefile header).
#   bash 03_code/07_v12/reproducibility/run_all.sh            # main pipeline, 4 parallel recipes x 8 BLAS threads = 32
#   PB_CHECK=1 bash .../run_all.sh                             # additionally rebuild pseudobulk from raw h5ad (~100 GB read)
set -euo pipefail
ROOT=${ROOT:-/home/weichangsheng/projects/cardiac_ageing_target_clearance}
HERE=$ROOT/03_code/07_v12/reproducibility
J=${J:-4}; NT=${NT:-8}
[ $((J*NT)) -le 32 ] || { echo "J*NT must be <= 32 on the shared machine"; exit 2; }
cd "$HERE"
mkdir -p "$ROOT/04_results/rerun_v12"
make -k -j"$J" NT="$NT" ROOT="$ROOT" all > "$ROOT/04_results/rerun_v12/make_all.log" 2>&1 || echo "some targets failed: see make_all.log"
if [ "${PB_CHECK:-0}" = 1 ]; then make NT="$NT" ROOT="$ROOT" pseudobulk_check; fi
ls "$ROOT/04_results/rerun_v12/stamps"
