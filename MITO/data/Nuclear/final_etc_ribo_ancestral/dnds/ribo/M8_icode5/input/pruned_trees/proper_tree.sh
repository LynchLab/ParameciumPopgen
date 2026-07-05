#!/usr/bin/env bash
set -euo pipefail

for f in OMA*.tree; do
  sed -E '
    # collapse ((X):a):b  → (X):(a+b)
    :again
    s/\(\(([^()]+)\):([0-9.eE+-]+)\):([0-9.eE+-]+)/(\1):\2+\3/;
    t again

    # clean leftover syntax artifacts
    s/,\)/)/g
    s/\(,/\(/g
    s/,,+/,/g
  ' "$f" \
  | python3 - <<'PY'
import sys,re
s=sys.stdin.read().strip()
# evaluate summed branch lengths like ":0.12+0.34"
s=re.sub(r':([0-9.eE+-]+)\+([0-9.eE+-]+)',
         lambda m: f":{float(m.group(1))+float(m.group(2)):.6f}", s)
print(s if s.endswith(";") else s+";")
PY
done

