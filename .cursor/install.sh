#!/usr/bin/env bash
# Idempotent setup for the Function bioinformatics notes repository.
# Installs R and the commonly used CRAN packages that the repository's
# runnable snippets rely on (base stats, ggplot2 plotting, tidyverse/igraph
# data wrangling). Uses Ubuntu's prebuilt r-cran-* binaries so the install is
# fast and reproducible. Safe to run repeatedly.
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive

APT_PACKAGES=(
  r-base
  r-base-dev
  r-cran-ggplot2
  r-cran-igraph
  r-cran-tidyverse
  r-cran-dplyr
  r-cran-tibble
  r-cran-readr
  r-cran-purrr
  r-cran-stringr
  r-cran-data.table
  r-cran-matrix
  r-cran-reshape2
  r-cran-rcolorbrewer
  r-cran-patchwork
  r-cran-circlize
  r-cran-ggsci
  r-cran-gtools
  r-cran-hmisc
)

echo "==> Updating apt package lists"
sudo apt-get update -qq

echo "==> Installing R and CRAN binary packages"
sudo apt-get install -y --no-install-recommends "${APT_PACKAGES[@]}"

echo "==> Verifying R installation"
Rscript -e 'cat("R", as.character(getRversion()), "ready\n")'

echo "==> install.sh complete"
