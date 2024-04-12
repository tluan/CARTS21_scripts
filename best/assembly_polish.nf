#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info "assembly_polish Pipeline starting..."

// Define parameters for easy adjustment
params.forward_reads = ""
params.reverse_reads = ""
params.long_reads = ""
params.output = ""
params.help = false
params.threads = 8
params.maxiter = 4
params.sample_id = ""
params.nextpolish_path=""
// Helper function for displaying usage information
def usage(status) {
    log.info "Usage: \nnextflow run assembly_polish.nf \n" +
            "[--forward_reads /path/to/forwardReads --reverse_reads /path/to/reverseReads --long_reads /path/to/unpairedReads] \n" +
            "--output /path/to/outputDir\n"
    log.info "Required:\n"
    log.info " --forward_reads        Path to forward paired-end read."
    log.info " --reverse_reads        Path to reverse paired-end read."
    log.info " --long_reads           Path to unpaired read fasta file(s)."
    log.info " --output               Path to the output folder."
    log.info " --nextpolish_path      Path to your nextpolish executable." 
    log.info " --help                 Print help message."

    System.exit(status)
}

// Check for help request
if (params.help) {
    usage(0)
}

// Check for required parameters
if (!params.forward_reads || !params.reverse_reads || !params.long_reads || !params.output|| !params.nextpolish_path) {
    log.error "Error: Missing required parameter(s)."
    usage(1)
}

// Process for assembling genomes using Unicycler
process UnicyclerAssembly{
    publishDir "${params.output}/unicycler", mode: 'copy'

    input:
    path forward_reads
    path reverse_reads
    path long_reads

    output:
    val true

    script:
    """
    mkdir -p ${params.output}/unicycler
    unicycler -1 ${forward_reads} -2 ${reverse_reads} -l ${long_reads} -o ${params.output}/unicycler -t ${params.threads}
    echo "test"
    """
}

// Process for polishing assembled genomes using Racon
process RaconPolishing {
    publishDir "${params.output}/racon", mode: 'copy', pattern: "${sample_id}_racon*.fasta"
    maxForks 4

    input:
    val ready
    path long_reads

    output:
    val true

    script:
    """
 
    polished=${params.output}/unicycler/assembly.fasta
    for i in {1..${params.maxiter}}
    do
     
       minimap2 -x map-ont -t ${params.threads} \${polished} ${long_reads} > mappings.paf

       racon -t ${params.threads} ${long_reads} mappings.paf \${polished} > racon_round\${i}.fasta
       polished=racon_round\${i}.fasta
    done
    mv racon_round${params.maxiter}.fasta ${params.output}/racon_round${params.maxiter}.fasta
    """
}

process MedakaPolishing {
    publishDir "${params.output}/medaka", mode: 'copy'

    input:
    val ready
    path long_reads

    output:
    val true

    script:
    """
 
    medaka_consensus -i ${long_reads} -d ${params.output}/racon_round${params.maxiter}.fasta  -o . -t ${params.threads}
    mv consensus.fasta ${params.output}/medaka.fasta
    """
}

process NextPolishPolishing {
    publishDir "${params.outdir}/nextpolish", mode: 'copy'

    input:
    val ready
    path forward_reads
    path reverse_reads

    output:
    val true

    script:
    """
   
    NextPolishinput=${params.output}/medaka.fasta
    for ((i=1; i<=${params.maxiter};i++)); do

        bwa index \${NextPolishinput};
        bwa mem -t ${params.threads} \${NextPolishinput} ${forward_reads} ${reverse_reads}|samtools view --threads 3 -F 0x4 -b -|samtools fixmate -m --threads 3  - -|samtools sort -m 2g --threads 5 -|samtools markdup --threads 5 -r - sgs.sort.bam

        samtools index -@ ${params.threads} sgs.sort.bam;
        samtools faidx \${NextPolishinput};

        python ${params.nextpolish_path}/NextPolish/lib/nextpolish1.py -g \${NextPolishinput} -t 1 -p ${params.threads} -s sgs.sort.bam > genome.polishtemp.fa;
        NextPolishinput=genome.polishtemp.fa;

        bwa index \${NextPolishinput};
        bwa mem -t ${params.threads} \${NextPolishinput} ${forward_reads} ${reverse_reads}|samtools view --threads 3 -F 0x4 -b -|samtools fixmate -m --threads 3  - -|samtools sort -m 2g --threads 5 -|samtools markdup --threads 5 -r - sgs.sort.bam
        samtools index -@ ${params.threads} sgs.sort.bam;
        samtools faidx \${NextPolishinput};
        python ${params.nextpolish_path}/NextPolish/lib/nextpolish1.py -g \${NextPolishinput} -t 2 -p ${params.threads} -s sgs.sort.bam > genome.nextpolish.fa;
        NextPolishinput=genome.nextpolish.fa;

    done;
    cp genome.nextpolish.fa ${params.output}/genome.nextpolish.fa
    """
}

// Define the workflow logic
workflow {
    forward_reads = file(params.forward_reads)
    reverse_reads = file(params.reverse_reads)
    long_reads = file(params.long_reads)

    // Define the process sequence
    UnicyclerAssembly(forward_reads, reverse_reads, long_reads)
    RaconPolishing(UnicyclerAssembly.out,long_reads)
    MedakaPolishing(RaconPolishing.out, long_reads)
    NextPolishPolishing(MedakaPolishing.out, forward_reads, reverse_reads)
}
