#!/bin/bash
#SBATCH -J flye # Job name
#SBATCH -o flye.o # Name of Output File
#SBATCH -e flye.e # Name of Error File
o
#SBATCH --mail-type=all
#SBATCH --time=18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --ntasks=8


flye --nano-raw /fs/cbcb-scratch/tluan/iso_assembler_eval/exp_samples/filtered_reads/fltd_CFSAN110838.fastq  --threads 8 --asm-coverage 200 --genome-size 50m --out-dir /fs/cbcb-scratch/tluan/iso_assembler_eval/exp_samples/CFSAN110838_fltd/flye_CFSAN110838_fltd_cov_50
