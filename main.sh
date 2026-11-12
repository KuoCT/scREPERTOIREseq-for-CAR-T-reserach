#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

Rscript scripts/GEX_analysis.R
Rscript scripts/REP_analysis.R

rm -f Rplots.pdf