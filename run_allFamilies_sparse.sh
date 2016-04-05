INPUT_PATH=input/NAPA
alpha=1

echo > tri_table_sparse.tsv
for instance_id in {1..1}; do		
	instance=Family_$instance_id
	G_net=$INPUT_PATH/$instance/A.smat
	H_net=$INPUT_PATH/$instance/B.smat
	SeqSim=$INPUT_PATH/$instance/A-B.smat

	rm output/*
	echo ./tri-match -t smat -H $H_net -G $G_net -S $SeqSim -Y 2 -a $alpha -b 0 --iter 3 -x seqsim
	#fname=`ls output/ -1t | head -1`
	#outpath=/home/shahin/Dropbox/Projects/HONEA/experiment/output/NAPA/sorted/Family_$instance_id/cTAME_alpha=0.5
	#mkdir -p $outpath
	#mv output/$fname $outpath/X.smat
done


