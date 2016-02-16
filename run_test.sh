alpha=('0' '1e-2' '1e-1' '1')
beta=('0.85' '0.5' '0.15')

epsilon=1e-5

NDIM=1
MAX_IT=20
INPUT_PATH=../../input
OUTPUT_PATH=../../output
prefix="BioGRID_full"

for i in "${!alpha[@]}"; do
	G_net=$INPUT_PATH/$prefix\_human_net.smat
	H_net=$INPUT_PATH/$prefix\_yeast_net.smat
	SeqSim=$INPUT_PATH/$prefix\_human-yeast.smat
	for k in "${!beta[@]}"; do
		Name=TAME\_alpha=${alpha[$i]}\_beta=${beta[$k]}

		echo ./tri-match -t smat -G $G_net -H $H_net -S $SeqSim --output $outPath --dim $NDIM --iter $MAX_IT -x seqsim -s combined -a ${alpha[$i]} -b ${beta[$k]} -e $epsilon
	done
done


