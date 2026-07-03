#!/usr/bin/env bash
set -euo pipefail

CSV="oma_species.csv"
MASTER="master.nwk"
REN="rename.sed"

mkdir -p pruned_trees
: > all_98_pruned_trees.tsv
echo -e "oma_id\tnewick" >> all_98_pruned_trees.tsv

while IFS=',' read -r oma s1 s2 s3 s4 s5 s6 s7; do
  [[ "$oma" == "oma_id" || -z "$oma" ]] && continue

  # species present for this OMA (alignment labels)
  keep="$s1 $s2 $s3 $s4 $s5 $s6 $s7"

  # start with all PB taxa, then remove those present
  remove="PTET PPRIM PBIA PMMN PCAU PBUR PSEX"

  for k in $keep; do
    case "$k" in
      tetr) remove=${remove//PTET/} ;;
      prim) remove=${remove//PPRIM/} ;;
      biau) remove=${remove//PBIA/} ;;
      mult) remove=${remove//PMMN/} ;;
      caud) remove=${remove//PCAU/} ;;
      burs) remove=${remove//PBUR/} ;;
      sexa) remove=${remove//PSEX/} ;;
      *) : ;;
    esac
  done

  tree=$(cat "$MASTER")

  # remove unwanted tips (also handles leading commas)
  for r in $remove; do
    tree=$(echo "$tree" | sed -E "s/(,)?$r:[0-9.eE+-]+//g")
  done

  # clean syntax artifacts
  tree=$(echo "$tree" | sed -E 's/,\)/)/g; s/\(,/\(/g; s/,,/,/g')

  # rename tips to alignment names
  tree=$(echo "$tree" | sed -f "$REN")

  # remove trailing ';' then wrap to ensure FULLY BRACKETED ROOT WITH 0.0
  core="${tree%;}"
  final="((${core}):0.0);"

  # write per OMA tree
  echo "$final" > "pruned_trees/${oma}.tree"

  # write combined file
  echo -e "${oma}\t${final}" >> all_98_pruned_trees.tsv

done < "$CSV"

echo "✅ Done"
echo "📂 Trees: pruned_trees/OMA*.tree"
echo "📄 Combined: all_98_pruned_trees.tsv"

