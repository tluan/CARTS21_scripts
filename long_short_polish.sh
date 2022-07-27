#!/bin/bash




echo "Input long reads as fastq: " $1
long_reads=$1
echo "Input contigs as fasta: " $2
ctg=$2
echo "Input short reads 1 as fastq: " $3
short_reads_1=$3
echo "Input short reads 2 as fastq: " $4
short_reads_2=$4
echo "Sample name: " $5
samplename=$5
echo "Poished reads outputs directory: " $6
results_dir=$6
echo "Number of Threads: " $7
threads=$7




mkdir ${samplename}_long_eval
cd ${samplename}_long_eval
cp $ctg ${samplename}_assembly.fasta
ctg=${samplename}_assembly.fasta

echo echo "Contig Name Changed to: " $ctg

outdir='racon_polished_medaka'
mkdir $outdir

tmpovl='tmp.ovl.paf'




maxiter=4


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





mv tmp.$maxiter.fasta racon_out_ctg.fasta



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
medaka_consensus -i "$reads" -d racon_out_ctg.fasta -o $outdir -t $threads
conda deactivate

echo '#################################################'
echo 'Finished Medaka'
echo '#################################################'




mv $outdir/consensus.fasta $outdir/${samplename}_racon_4x_medaka_consensus.fasta



























