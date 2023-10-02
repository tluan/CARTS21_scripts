#!/bin/sh
#SBATCH -J long_polish_pipeline # Job name
#SBATCH -o long_polish_pipeline.o%j # Name of Output File
#SBATCH -e long_polish_pipeline.e%j # Name of Error File

#SBATCH --mail-type=all
#SBATCH --time=7-2:00:00
#SBATCH --qos=large
#SBATCH --mem=30gb
#SBATCH --ntasks=16


module load racon
module load blast



echo "Input reads as fastq - FAST5 not required: " $1
reads=$1
echo "Input contigs as fasta: " $2
ctg=$2
echo "Sample name: " $3
samplename=$3
echo "Poished reads outputs: " $4
results_dir=$4
echo  "Method: " $5
method=$5




mkdir ${samplename}_long_eval
cd ${samplename}_long_eval
cp $ctg ${samplename}_assembly.fasta
ctg=${samplename}_assembly.fasta

echo echo "Contig Name Changed to: " $ctg

outdir='racon_polished_medaka'
mkdir $outdir

tmpovl='tmp.ovl.paf'




maxiter=4

threads=16

cp $ctg tmp.0.fasta


echo '#################################################'
echo 'Starting Minimap and Racon'
echo '#################################################'



for ((i=1; i<=$maxiter; i++))
do
	echo '#################################################'
	echo 'Minimap round: ' $i
	echo '#################################################'
	echo 'Minimap round: ' $i >&2
	h=$(expr $i - 1)
	time minimap2 -x map-ont -t $threads tmp.$h.fasta  "$reads" > $tmpovl



	echo '#################################################'
	echo 'Racon round: ' $i
	echo '#################################################'
	echo 'Racon round: ' $i  >&2
	time racon -t $threads "$reads" $tmpovl tmp.$h.fasta > tmp.$i.fasta

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
echo 'Starting Medaka minimap'>&2
echo "Info: Run Minimap for Medaka"
tmpsam="medakatmp.sam"
time minimap2 -a -x map-ont -t $threads racon_out_ctg.fasta  "$reads" > $tmpsam

echo "Converting SAM to BAM"
echo "Converting SAM to BAM" >&2
x=$tmpsam
time samtools view -@ $threads -bhS $x > $x.bam
time samtools sort -@ $threads $x.bam -o $x.s.bam
time samtools index $x.s.bam
rm $x.bam

module unload racon
conda activate medaka



nopolish='nopolish_medaka'
mkdir $nopolish
echo 'Starting racon_Medaka'>&2
time medaka_consensus -i "$reads" -d racon_out_ctg.fasta -o $outdir -t $threads
echo 'Starting only_Medaka'>&2
time medaka_consensus -i "$reads" -d $ctg -o $nopolish -t $threads



conda deactivate
echo '#################################################'
echo 'Finished Medaka'
echo '#################################################'


module load racon
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
	echo 'Minimap round: ' $k >&2
	echo '#################################################'
	h_=$(expr $k - 1)
	time minimap2 -x map-ont -t $threads tmp.$h_.fasta  "$reads" > $tmpovl



	echo '#################################################'
	echo 'Racon round: ' $k
	echo 'Racon round: ' $k  >&2
	echo '#################################################'

	time racon -t $threads "$reads" $tmpovl tmp.$h_.fasta > tmp.$k.fasta

done



cp tmp.$maxiter.fasta $medaka_racon/racon_out_ctg.fasta





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
mv $medaka_racon/racon_out_ctg.fasta $results_dir/${samplename}_${method}_medaka_1x_racon_4x.fasta
mv $outdir/consensus.fasta $results_dir/${samplename}_${method}_racon_4x_medaka_1x.fasta
mv $nopolish/consensus.fasta $results_dir/${samplename}_${method}_medaka_1x.fasta
mv $racon_only/racon_out_ctg.fasta $results_dir/${samplename}_${method}_racon_4x.fasta
mv $racon_only/racon_out_1_iter_ctg.fasta $results_dir/${samplename}_${method}_racon_1x.fasta

