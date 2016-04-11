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

	mkdir Family$i_output
	./tri-match -t smat -H $H_net -G $G_net -S $SeqSim -Y 0 -a 0.5 -b 0.1 --iter 3 -x seqsim -o Family$i_output &

done


