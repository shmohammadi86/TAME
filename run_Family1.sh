G_net=input/NAPA/Family_1/A.smat
H_net=input/NAPA/Family_1/B.smat
SeqSim=input/NAPA/Family_1/A-B.smat

./tri-match -t smat -G $G_net -H $H_net -S $SeqSim -Y 2 -a 0.85 -b 0.1 --iter 3 -x seqsim --post_iter 0

