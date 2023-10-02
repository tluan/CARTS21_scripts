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

long_in=$1
out=$2
flye --nano-raw ${long_in}  --threads 8 --asm-coverage 200 --genome-size 10m --out-dir $out
