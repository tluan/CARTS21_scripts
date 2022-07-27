#!/bin/bash





echo "Input reads as fastq - FAST5 not required: " $1
reads=$1
echo "Input contigs as fasta: " $2
ctg=$2
echo "Evaluation sequence as fasta: " $3
eval=$3
echo "Sample name: " $4
samplename=$4
echo "Poished reads outputs: " $5
results_dir=$5





mkdir ${samplename}_long_eval
cd ${samplename}_long_eval
cp $ctg ${samplename}_assembly.fasta
ctg=${samplename}_assembly.fasta

echo echo "Contig Name Changed to: " $ctg

outdir='racon_polished_medaka'
mkdir $outdir

tmpovl='tmp.ovl.paf'




maxiter=4

threads=8

cp $ctg tmp.0.fasta


echo '#################################################'
echo 'Starting Minimap and Racon'
echo '#################################################'



for ((i=1; i<=$maxiter; i++))
do
	echo '#################################################'
	echo 'Minimap round: ' $i
	echo '#################################################'
	h=$(expr $i - 1)
	minimap2 -x map-ont -t $threads tmp.$h.fasta  "$reads" > $tmpovl



	echo '#################################################'
	echo 'Racon round: ' $i
	echo '#################################################'

	racon -t $threads "$reads" $tmpovl tmp.$h.fasta > tmp.$i.fasta
	
done


racon_only='racon_polish_only'
mkdir $racon_only

mv tmp.$maxiter.fasta racon_out_ctg.fasta
mv tmp.1.fasta racon_out_1_iter_ctg.fasta

cp racon_out_ctg.fasta $racon_only
cp racon_out_1_iter_ctg.fasta $racon_only


echo '#################################################'
echo 'Finished Minimap and Racon'
echo '#################################################'




echo '#################################################'
echo 'Starting Medaka'
echo '#################################################'

echo "Info: Run Minimap for Medaka"
tmpsam="medakatmp.sam"
minimap2 -a -x map-ont -t $threads racon_out_ctg.fasta  "$reads" > $tmpsam

echo "Converting SAM to BAM"
x=$tmpsam
samtools view -@ $threads -bhS $x > $x.bam
samtools sort -@ $threads $x.bam -o $x.s.bam
samtools index $x.s.bam
rm $x.bam


conda activate medaka



nopolish='nopolish_medaka'
mkdir $nopolish

medaka_consensus -i "$reads" -d racon_out_ctg.fasta -o $outdir -t $threads
medaka_consensus -i "$reads" -d $ctg -o $nopolish -t $threads



conda deactivate
echo '#################################################'
echo 'Finished Medaka'
echo '#################################################'



medaka_racon='medaka_racon'
mkdir $medaka_racon

rm tmp.?.fasta
cp $nopolish/consensus.fasta tmp.0.fasta



echo '#################################################'
echo 'Starting Minimap and Racon'
echo '#################################################'



for ((k=1; k<=$maxiter; k++))
do
	echo '#################################################'
	echo 'Minimap round: ' $k
	echo '#################################################'
	h_=$(expr $k - 1)
	minimap2 -x map-ont -t $threads tmp.$h_.fasta  "$reads" > $tmpovl



	echo '#################################################'
	echo 'Racon round: ' $k
	echo '#################################################'

	racon -t $threads "$reads" $tmpovl tmp.$h_.fasta > tmp.$k.fasta
	
done



cp tmp.$maxiter.fasta $medaka_racon/racon_out_ctg.fasta




racon_only='racon_polish_only'
mkdir $racon_only

mv tmp.$maxiter.fasta racon_out_ctg.fasta
mv tmp.1.fasta racon_out_1_iter_ctg.fasta

cp racon_out_ctg.fasta $racon_only
cp racon_out_1_iter_ctg.fasta $racon_only


echo '#################################################'
echo 'Finished Minimap and Racon'
echo '#################################################'


















rm tmp.ovl.paf
rm $tmpsam
rm tmp.?.fasta

echo '#################################################'
echo 'Everything has finished'
echo '#################################################'


eval_dir='eval_dir'
eval_dir=${samplename}_$eval_dir
mkdir $eval_dir

mkdir $results_dir
mv $medaka_racon/racon_out_ctg.fasta $results_dir/${samplename}_medaka_racon_4x.fasta
mv $outdir/consensus.fasta $results_dir/${samplename}_racon_4x_medaka_consensus.fasta
mv $nopolish/consensus.fasta $results_dir/${samplename}_only_medaka_consensus.fasta
mv $racon_only/racon_out_ctg.fasta $results_dir/${samplename}_only_racon_4x_consensus.fasta 
mv $racon_only/racon_out_1_iter_ctg.fasta $results_dir/${samplename}_only_racon_1x_consensus.fasta

conda activate quast-env 
quast $results_dir/${samplename}_only_medaka_consensus.fasta -o $eval_dir/eval_medaka_only -r $eval
quast $results_dir/${samplename}_racon_4x_medaka_consensus.fasta -o $eval_dir/eval_racon_medaka_iter -r $eval
quast $results_dir/${samplename}_only_racon_4x_consensus.fasta -o $eval_dir/eval_racon_iter_only -r $eval
quast $results_dir/${samplename}_only_racon_1x_consensus.fasta -o $eval_dir/eval_racon_single_only -r $eval
quast $ctg  -o $eval_dir/eval_no_polish -r $eval
quast $results_dir/${samplename}_medaka_racon_4x.fasta  -o $eval_dir/eval_medaka_racon_iter -r $eval
conda deactivate



cd $eval_dir



multiqc -n ${samplename}_eval_longread_polish .
