./tri-match -t smat -H input/BioGRID/yeast_net.smat -G input/BioGRID/human_net.smat -S input/BioGRID/human-yeast.smat --iter 0 -x seqsim -a 0.15 -b 10 -C output/X_dense.mat
./tri-match -t smat -H input/BioGRID/yeast_net.smat -G input/BioGRID/human_net.smat -S input/BioGRID/human-yeast.smat --iter 0 -x seqsim -a 0.5 -b 10 -C output/X_dense.mat
time ./tri-match -t smat -H input/BioGRID/yeast_net.smat -G input/BioGRID/human_net.smat -S input/BioGRID/human-yeast.smat --iter 0 -x seqsim -a 0.85 -b 10 -C output/X_dense.mat
