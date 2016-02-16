alpha=('0' '1e-3' '1e-2' '1e-1' '1' '1e1' '1e2' '1e3')
beta=('1' '0.85' '0.5' '0.15')
epsilon=1e-5

NDIM=1
MAX_IT=15
INPUT_PATH=../../input/NAPA/pairwise/CG_set

echo > tri_table_sparse.tsv
for instance_id in {1..10}; do		
	instance=Family_$instance_id
	echo $instance >>  tri_table_sparse.tsv
	for i in "${!alpha[@]}"; do
		G_net=$INPUT_PATH/$instance/A.smat
		H_net=$INPUT_PATH/$instance/B.smat
		SeqSim=$INPUT_PATH/$instance/A-B.smat
		for k in "${!beta[@]}"; do
			outPath=./output/NAPA/$instance/TAME_alpha=${alpha[$i]}\_beta=${beta[$k]}\_it=$MAX_IT\_type=sparse
			mkdir -p $outPath
			Name=$instance\_beta=${beta[$k]}_sparse
			echo logfile $Name.log > $Name\_screenrc.conf
			screen  -c $Name\_screenrc.conf -dmLS $Name ./tri-match -t smat -G $G_net -H $H_net -S $SeqSim --dim $NDIM --iter $MAX_IT -x seqsim -s combined -a ${alpha[$i]} -b ${beta[$k]} -e $epsilon --output $outPath -Y 2
		MAX=`cat $outPath/max.txt`
		printf "%s," "$MAX" >>  tri_table_sparse.tsv
		done
	printf "\n"  >>  tri_table_sparse.tsv
	done
	echo
done


