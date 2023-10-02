#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --time=18:00:00
#SBATCH --mem=36gb
#SBATCH --ntasks=8




echo "Sample name: " $1
samplename=$1

echo "Input short read 1 as fastq: " $2
short_reads_1=$2

echo "Input short read 2 as fastq: " $3
short_reads_2=$3

echo "Input long read as fastq: " $4
long_read=$4

echo "Poished reads outputs: " $5
results_dir=$5

echo "Number of Threads: " $6
threads=$6

unicycler -1 $short_reads_1 -2 $short_reads_2 -l $long_read  -o ${results_dir}/${samplename}_unicycler  -t $threads
