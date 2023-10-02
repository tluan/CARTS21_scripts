#!/bin/bash
#SBATCH -J short_read_polish # Job name
#SBATCH -o short_read_polish.o%j # Name of Output File
#SBATCH -e short_read_polish.e%j # Name of Error File
#SBATCH --mail-type=all
#SBATCH --time=18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --ntasks=8





echo "Input contigs as fasta: " $1
ctg=$1

echo "Sample name: " $2
samplename=$2

echo "Input short read 1 as fastq: " $3
short_reads_1=$3

echo "Input short read 2 as fastq: " $4
short_reads_2=$4


echo "Poished reads outputs: " $5
results_dir=$5

echo  "Method: " $6
method=$6




threads=8
master=$results_dir
mkdir $master
cd $master
maxiter=4
threads=8
mkdir all_polished_${samplename}



pilon_dir='pilon'
mkdir $pilon_dir
cd $pilon_dir

consensus='illumina_mapping_consensus'
pilon_polished='pilon_polished'
mkdir $consensus
mkdir $pilon_polished
pilon_maxiter=$maxiter
input=$ctg
cp $ctg pilon_round_0.fasta
for ((j=1; j<=$pilon_maxiter; j++))
do
	echo '#################################################'
	echo 'Start - Pilon round: ' $j
	echo 'Start - Pilon round: ' $j >&2
	echo '#################################################'



    b=$(expr $j - 1)

    time {
    echo pilon_round_$b.fasta
    bwa index pilon_round_$b.fasta
    bwa mem -t $threads pilon_round_$b.fasta "$short_reads_1" "$short_reads_2" | samtools view -Sb  | samtools sort -@$threads -o mapping_pilon.$b.sorted.bam
    samtools index mapping_pilon.$b.sorted.bam
    pilon  --genome pilon_round_$b.fasta  --fix all --changes --frags mapping_pilon.$b.sorted.bam -Xmx3G --threads $threads --output pilon_round_$j
    }

	echo '#################################################'
	echo 'Done - Pilon round: ' $j
	echo 'Done - Pilon round: ' $j >&2
	echo '#################################################'



done

cd ..

mv $pilon_dir/pilon_round_1.fasta all_polished_${samplename}/${method}_pilon_1x.fasta
mv $pilon_dir/pilon_round_${pilon_maxiter}.fasta all_polished_${samplename}/${method}_pilon_4x.fasta


polypolish='polypolish'
mkdir $polypolish
cd $polypolish
polypolish_input=$ctg
poly_maxiter=$maxiter

for ((ppc=1; ppc<=$poly_maxiter; ppc++))
do
	echo '#################################################'
	echo 'Start - Polypolish round: ' $ppc

	echo 'Start - Polypolish round: ' $ppc >&2
	echo '#################################################'


    time {
    bwa index $polypolish_input
    bwa mem -t 16 -a $polypolish_input $short_reads_1 > alignments_1_$ppc.sam
    bwa mem -t 16 -a $polypolish_input $short_reads_2 > alignments_2_$ppc.sam
    polypolish_insert_filter.py --in1 alignments_1_$ppc.sam --in2 alignments_2_$ppc.sam --out1 filtered_1_$ppc.sam --out2 filtered_2_$ppc.sam
    polypolish $polypolish_input filtered_1_$ppc.sam filtered_2_$ppc.sam > polished_iter_$ppc.fasta
    }
	echo '#################################################'
	echo 'Done - Polypolish round: ' $ppc
	echo 'Done - Polypolish round: ' $ppc >&2
	echo '#################################################'
    polypolish_input=polished_iter_$ppc.fasta


done

cd ..

mv $polypolish/polished_iter_1.fasta all_polished_${samplename}/${method}_polypolish_1x.fasta
mv $polypolish/polished_iter_${poly_maxiter}.fasta all_polished_${samplename}/${method}_polypolish_${poly_maxiter}x.fasta



ntEdit='ntEdit'
ntEdit_maxiter=$maxiter
mkdir $ntEdit
cd $ntEdit
cp $ctg ntEdit_0_edited.fa
ntEditinput=ntEdit_0_edited.fa
nthits --outbloom -p solidBF -b 36 -k 25 -t $threads ${short_reads_1} ${short_reads_2}
for ((ntc=1; ntc<=${ntEdit_maxiter};ntc++));
do

    echo '#################################################'
    echo 'Start - ntEdit round: ' $ntc
    echo 'Start - ntEdit round: ' $ntc  >&2
    echo '#################################################'

    time ntedit -f $ntEditinput -r solidBF_k25.bf -b ntEdit_$ntc -t $threads

    echo '#################################################'
    echo 'Done - ntEdit round: ' $ntc
    echo 'Done - ntEdit round: ' $ntc >&2
    echo '#################################################'
    ntEditinput=ntEdit_${ntc}_edited.fa
done

cd ..
mv $ntEdit/ntEdit_1_edited.fa all_polished_${samplename}/${method}_ntEdit_1x.fa
mv $ntEdit/ntEdit_${ntEdit_maxiter}_edited.fa all_polished_${samplename}/${method}_ntEdit_${ntEdit_maxiter}x.fa




polca_maxiter=$maxiter
polca='polca'
mkdir $polca
cd $polca
cp $ctg polca_0_edited.fa
polca_input=polca_0_edited.fa

for ((pc=1; pc<=${polca_maxiter};pc++));
do



    echo '#################################################'
    echo 'Start - polca round: ' $pc
    echo 'Start - polca round: ' $pc >&2
    echo '#################################################'
    time {
    x="'$short_reads_1 $short_reads_2'"
    echo /fs/cbcb-scratch/tluan/conda_2/envs/polca_env/bin/polca.sh -a $polca_input -r $x -t $threads |bash
    }
    echo '#################################################'
    echo 'Done - polca round: ' $pc
    echo 'Done - polca round: ' $pc >&2
    echo '#################################################'
    mv "${polca_input}.PolcaCorrected.fa" "polca_${pc}_edited.PolcaCorrected.fa"
    polca_input=polca_${pc}_edited.PolcaCorrected.fa




done

cd ..
mv $polca/polca_1_edited.PolcaCorrected.fa all_polished_${samplename}/${method}_Polca_1x.fa
mv $polca/polca_${polca_maxiter}_edited.PolcaCorrected.fa all_polished_${samplename}/${method}_Polca_${polca_maxiter}x.fa


NextPolish='NextPolish'
mkdir $NextPolish
cd $NextPolish
NextPolishmaxiter=$maxiter
NextPolishinput=$ctg

for ((npc=1; npc<=${NextPolishmaxiter};npc++));
do
   echo '#################################################'
   echo 'Start - NextPolish round: ' $npc >&2
   echo 'Start - NextPolish round: ' $npc
   echo '#################################################'
   time {
   bwa index ${NextPolishinput};
   bwa mem -t ${threads} ${NextPolishinput} ${short_reads_1} ${short_reads_2}|samtools view --threads $threads -F 0x4 -b -|samtools fixmate -m --threads $threads  - -|samtools sort --threads $threads -|samtools markdup --threads $threads -r - sgs.sort.bam
   samtools index -@ ${threads} sgs.sort.bam;
   samtools faidx ${NextPolishinput};
   python /fs/cbcb-scratch/tluan/iso_assembler_eval/NextPolish/lib/nextpolish1.py -g ${NextPolishinput} -t 1 -p ${threads} -s sgs.sort.bam > genome.polishtemp.fa;
   NextPolishinput=genome.polishtemp.fa;
   bwa index ${NextPolishinput};
   bwa mem -t ${threads} ${NextPolishinput} ${short_reads_1} ${short_reads_2}|samtools view --threads 3 -F 0x4 -b -|samtools fixmate -m --threads 3  - -|samtools sort -m 2g --threads 5 -|samtools markdup --threads 5 -r - sgs.sort.bam
   samtools index -@ ${threads} sgs.sort.bam;
   samtools faidx ${NextPolishinput};
   python /fs/cbcb-scratch/tluan/iso_assembler_eval/NextPolish/lib/nextpolish1.py -g ${NextPolishinput} -t 2 -p ${threads} -s sgs.sort.bam > NextPolish_iter_$npc.fasta;
  }
   echo '#################################################'
   echo 'Done - NextPolish round: ' $npc >&2
   echo 'Done - NextPolish round: ' $npc
   echo '#################################################'
   NextPolishinput=NextPolish_iter_$npc.fasta
done;


cd ..


mv $NextPolish/NextPolish_iter_1.fasta all_polished_${samplename}/${method}_NextPolish_1x.fasta
mv $NextPolish/NextPolish_iter_${NextPolishmaxiter}.fasta all_polished_${samplename}/${method}_NextPolish_${NextPolishmaxiter}x.fasta









cd ..
