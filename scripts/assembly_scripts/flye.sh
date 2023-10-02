#!/bin/bash
#SBATCH -J flye # Job name
#SBATCH -o flye.o # Name of Output File
#SBATCH -e flye.e # Name of Error File
#SBATCH --mail-type=all
#SBATCH --time=18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --ntasks=8

longread=$1
out=$2
flye --nano-raw ${longread}  --threads 8 --out-dir ${out}
