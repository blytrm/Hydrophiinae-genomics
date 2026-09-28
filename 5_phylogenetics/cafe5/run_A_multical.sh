#!/bin/bash
cd /hpcfs/users/a1864358/sanders_lab/phylogenomics/cafe

C5=/hpcfs/users/a1864358/miniconda/envs/v2r_phylo/bin/cafe5
IQ=/hpcfs/users/a1864358/miniconda/envs/v2r_phylo/bin/iqtree2
TREE=trees/species_tree_ultrametric_multical.nwk

log(){ echo "[$(date +%H:%M:%S)] $*"; }

# 1) CAFE base + multi-lambda on multical tree
log "CAFE base_multical"
$C5 -i inputs/cafe_counts_filtered.tsv -t $TREE -p -c 16 -o runs/base_multical > logs/base_multical.log 2>&1
log "CAFE multilambda_multical"
$C5 -i inputs/cafe_counts_filtered.tsv -t $TREE -p -y trees/lambda_tree_marine.txt -c 16 -o runs/multilambda_multical > logs/multilambda_multical.log 2>&1

# 2) drop-one cross-validation dating
cd trees

declare -A CALS=( [caspera_hcy]="caspera,hcy -100:-90" [vberus_hcy]="vberus,hcy -55:-50" [nhel_hcy]="nhel,hcy -38:-30" [hcure_hmaj]="hcure,hmaj -5:-3" )
mkdir -p dating/crossval

for drop in caspera_hcy vberus_hcy nhel_hcy hcure_hmaj; do
  df=dating/crossval/drop_${drop}.txt
  { for t in caspera hcure hcurw hcy hmaj horn nhel nscut ptext tele vberus vkomodo; do echo "$t 0"; done
    for k in caspera_hcy vberus_hcy nhel_hcy hcure_hmaj; do [ "$k" != "$drop" ] && echo "${CALS[$k]}"; done
  } > $df
  $IQ -s SpeciesTreeAlignment.named.fa -te dated2.treefile -m LG+G4 --date $df --date-root -170 --dating LSD --date-options "-u 0" -pre dating/crossval/drop_${drop} -nt 8 -redo > dating/crossval/drop_${drop}.run.log 2>&1
done


cd /hpcfs/users/a1864358/sanders_lab/phylogenomics/cafe
log "A_MULTICAL_DONE"
