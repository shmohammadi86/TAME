#!/usr/bin/env bash
# Align the human and yeast BioGRID interactomes in dense mode (Y=0),
# seeded with BLAST sequence similarity.
set -euo pipefail
make -s
./tri-match -t smat \
    -G input/BioGRID_full_human_net.smat \
    -H input/BioGRID_full_yeast_net.smat \
    -S input/BioGRID_full_human-yeast.smat \
    -x seqsim -Y 0 -a 0.15 -b 10 --iter 3 --post_iter 3 -o output
