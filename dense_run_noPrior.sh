#!/usr/bin/env bash
# Align the human and yeast BioGRID interactomes in dense mode with a uniform
# prior, that is, on topology alone with no sequence similarity.
set -euo pipefail
make -s
./tri-match -t smat \
    -G input/BioGRID_full_human_net.smat \
    -H input/BioGRID_full_yeast_net.smat \
    -x uniform -Y 0 -a 0.15 -b 10 --iter 3 --post_iter 3 -o output
