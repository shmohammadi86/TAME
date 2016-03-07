alpha=('0')
beta=('1')
epsilon=1e-5

NDIM=1
MAX_IT=15
INPUT_PATH=input/NAPA

echo > tri_table_sparse.tsv
for instance_id in {1..10}; do		
	instance=Family_$instance_id
	G_net=$INPUT_PATH/$instance/A.smat
	H_net=$INPUT_PATH/$instance/B.smat
	SeqSim=$INPUT_PATH/$instance/A-B.smat

	rm output/*
	./tri-match -t smat -H $H_net -G $G_net -S $SeqSim -Y 0 -b 1 -a 0.1 --iter 3 -x seqsim
	fname=`ls output/ -1t | head -1`
	outpath=/home/shahin/Dropbox/Projects/HONEA/experiment/output/NAPA/sorted/Family_$instance_id/TAME_alpha=0.1_it=3
	mkdir -p $outpath
	mv output/$fname $outpath/X.smat
done


