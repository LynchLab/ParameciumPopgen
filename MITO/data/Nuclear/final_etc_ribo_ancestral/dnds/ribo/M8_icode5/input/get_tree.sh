#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# FULL PATHS
# ============================================================

BASE="/data/LynchLabCME/Paramecium/popGenomics_Farhan_Ravina/Ravina/nuclear_etc_ribo/burs_blast/diploid_out/dnds_analysis/input"

CSV="${BASE}/oma_species.csv"
MASTER="${BASE}/master.nwk"

OUTDIR="${BASE}/pruned_trees"
OUTTSV="${BASE}/all_pruned_trees.tsv"

mkdir -p "${OUTDIR}"
: > "${OUTTSV}"
echo -e "oma_id\tnewick" >> "${OUTTSV}"

master_tree="$(tr -d '\n' < "${MASTER}")"

# master tree tips
ALL_TIPS=(PTET PBUR PMMN PCAU PSEX PBIA PPRIM)

remove_tip () {
  # remove tip "TIP:branchlen" from Newick, handling commas
  local tree="$1"
  local tip="$2"
  tree="$(echo "$tree" | sed -E "s/(,)?${tip}:[0-9.eE+-]+//g")"
  tree="$(echo "$tree" | sed -E "s/${tip}:[0-9.eE+-]+,?//g")"
  echo "$tree"
}

clean_newick () {
  sed -E \
    -e 's/,\)/)/g' \
    -e 's/\(,/\(/g' \
    -e 's/,,+/,/g' \
    -e 's/;+$//'
}

# ============================================================
# CSV format:
# oma_id,species
# OMA01173,PBIA,PBUR,PCAU,PMMN,PPRIM,PSEX,PTET
# ============================================================

while IFS=',' read -r oma rest; do
  [[ "${oma}" == "oma_id" || -z "${oma}" ]] && continue

  IFS=',' read -r -a cols <<< "$(echo "${oma},${rest}")"
  oma_id="${cols[0]}"

  keep_tips=()
  for ((i=1; i<${#cols[@]}; i++)); do
    t="$(echo "${cols[i]}" | tr -d '[:space:]')"
    [[ -n "$t" ]] && keep_tips+=("$t")
  done

  # remove tips = ALL_TIPS minus keep_tips
  remove_tips=()
  for tip in "${ALL_TIPS[@]}"; do
    found=0
    for k in "${keep_tips[@]}"; do
      [[ "$tip" == "$k" ]] && found=1 && break
    done
    [[ $found -eq 0 ]] && remove_tips+=("$tip")
  done

  tree="$master_tree"

  for r in "${remove_tips[@]}"; do
    tree="$(remove_tip "$tree" "$r")"
  done

  tree="$(echo "$tree" | clean_newick)"
  tree="$(echo "$tree" | clean_newick)"

  tree="${tree%;}"
  final="${tree};"

  echo "$final" > "${OUTDIR}/${oma_id}.tree"
  echo -e "${oma_id}\t${final}" >> "${OUTTSV}"
done < "${CSV}"

echo "✅ Done"
echo "📂 Trees: ${OUTDIR}/OMA*.tree"
echo "📄 Combined TSV: ${OUTTSV}"

