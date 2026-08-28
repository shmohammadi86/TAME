#!/usr/bin/env bash
# Align the human and yeast BioGRID interactomes in sparse mode (Y=2),
# seeded with BLAST sequence similarity.
set -euo pipefail
make -s
./tri-match -t smat \
    -G input/BioGRID_full_human_net.smat \
    -H input/BioGRID_full_yeast_net.smat \
    -S input/BioGRID_full_human-yeast.smat \
    -x seqsim -Y 2 -a 0.85 -b 0.1 --iter 3 --post_iter 3 -o output
