#!/usr/bin/env bash
# Idempotent environment setup for a fresh cloud container (about 3-5 min the first time, ~1 s afterwards).
# Installs R from Ubuntu apt, binary R packages from apt, and robustmatrix + cellWise from the
# GitHub CRAN mirror (CRAN itself may be blocked by the network policy).
# Usage: bash tools/overnight/setup.sh   (from repo root)
set -euo pipefail
LOG=/tmp/overnight_setup.log
if Rscript -e 'quit(status = !all(sapply(c("robustmatrix","cellWise","MASS"), requireNamespace, quietly=TRUE)))' >/dev/null 2>&1; then
  echo "setup: already done ($(Rscript -e 'cat(R.version.string, "| robustmatrix", as.character(packageVersion("robustmatrix")))'))"; exit 0
fi
export DEBIAN_FRONTEND=noninteractive
echo "setup: installing R + apt packages (log: $LOG)"
apt-get update -qq >>"$LOG" 2>&1 || true
apt-get install -y -qq r-base-core r-base-dev \
  r-cran-rcpp r-cran-rcpparmadillo r-cran-rcppprogress r-cran-rdpack r-cran-mass \
  r-cran-ggplot2 r-cran-matrixstats r-cran-rrcov r-cran-robustbase r-cran-reshape2 \
  r-cran-scales r-cran-gridextra r-cran-shape >>"$LOG" 2>&1
SRC=$(mktemp -d)
for pkg in svd robustmatrix cellWise; do
  if ! Rscript -e "quit(status = !requireNamespace('$pkg', quietly=TRUE))" >/dev/null 2>&1; then
    echo "setup: building $pkg from github.com/cran/$pkg"
    git clone -q --depth 1 "https://github.com/cran/$pkg.git" "$SRC/$pkg" >>"$LOG" 2>&1
    MAKEFLAGS="-j$(nproc)" R CMD INSTALL --no-docs --no-multiarch "$SRC/$pkg" >>"$LOG" 2>&1
  fi
done
Rscript -e 'for (p in c("robustmatrix","cellWise","MASS","Rcpp","RcppArmadillo")) cat(p, as.character(packageVersion(p)), "\n")'
echo "setup: done"
