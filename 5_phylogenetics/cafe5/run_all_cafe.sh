#!/bin/bash
# CAFE5 model suite driver. All models on filtered counts (max<100 copies/species).
set -e
CAFE=/hpcfs/users/a1864358/miniconda/envs/v2r_phylo/bin/cafe5
CNT=inputs/cafe_counts_filtered.tsv
TREE=trees/species_tree_ultrametric.nwk
LT=trees/lambda_tree_marine.txt
C=32
cd /hpcfs/users/a1864358/sanders_lab/phylogenomics/cafe
log(){ echo "[$(date +%H:%M:%S)] $*"; }

log "M_final gamma k=4 + error model"
$CAFE -i $CNT -t $TREE -p -k 4 -e$EM -c $C -o runs/gamma_k4_error > logs/gamma_k4_error.log 2>&1
    # MODEL SELECTED ^^^^

log "M2 error-model + base"
$CAFE -i $CNT -t $TREE -p -e -c $C -o runs/base_error         > logs/base_error.log 2>&1

log "M3 gamma k=2"
$CAFE -i $CNT -t $TREE -p -k 2 -c $C -o runs/gamma_k2         > logs/gamma_k2.log 2>&1

log "M4 gamma k=3"
$CAFE -i $CNT -t $TREE -p -k 3 -c $C -o runs/gamma_k3         > logs/gamma_k3.log 2>&1

log "M5 multi-lambda (marine vs background)"
$CAFE -i $CNT -t $TREE -p -y $LT -c $C -o runs/multilambda    > logs/multilambda.log 2>&1

log "M6 gamma k=3 + error model"
EM=runs/base_error/Base_error_model.txt
$CAFE -i $CNT -t $TREE -p -k 3 -e$EM -c $C -o runs/gamma_k3_error > logs/gamma_k3_error.log 2>&1

log "M7 multi-lambda + error model"
$CAFE -i $CNT -t $TREE -p -y $LT -e$EM -c $C -o runs/multilambda_error > logs/multilambda_error.log 2>&1

log "ALL DONE"
