

---

# Assembly Polish Pipeline Achieving Near-perfect Accuracy 

## Overview

This Assembly Polish Pipeline is a robust Nextflow-based workflow designed to streamline genome assembly and polishing processes. It leverages Unicycler for assembly, along with Racon and Medaka for long-read polishing, and NextPolish for short-read polishing. This pipeline facilitates the assembly of genomes from paired-end and long-read sequencing data, followed by iterative polishing to significantly improve assembly quality.

## Requirements

To be able to reproduce our results, we recommend using the versions listed below. However, later versions are also compatible.

### Software
- **Nextflow**: Version 20.04.1 or later. [Get started with Nextflow](https://www.nextflow.io/docs/latest/getstarted.html).
- **Java**: Version 8 or later.

### Bioinformatics Tools
- **Unicycler**: Version 0.5.0. [GitHub repository](https://github.com/rrwick/Unicycler).
- **Racon**: Version 1.4.12. [GitHub repository](https://github.com/isovic/racon).
- **Minimap2**: Version 2.26. [GitHub repository](https://github.com/lh3/minimap2).
- **Medaka**: Version 1.7.2. [GitHub repository](https://github.com/nanoporetech/medaka).
- **BWA**: Version 0.7.17. [Official website](http://bio-bwa.sourceforge.net/).
- **Samtools**: Version 1.6. [Official website](http://www.htslib.org/).
- **Python**: Version 3.x for specific scripts.
- **NextPolish**: Version 1.4.0 [Official website](https://github.com/Nextomics/NextPolish). This tool must be installed in the same directory as the Nextflow script for proper functionality.

---

This version enhances readability by organizing requirements into categories and providing direct links for more information. If there are additional sections you want to revise or add, such as installation instructions or usage examples, feel free to let me know!
## Installation

Before running the pipeline, ensure all the required tools are installed and accessible in your environment. This may include setting up Conda environments or installing software locally.

## Configuration

Parameters can be adjusted in the script or passed as command-line arguments:

- `--forward_reads`: Path to forward paired-end reads FASTQ file.
- `--reverse_reads`: Path to reverse paired-end reads FASTQ file.
- `--long_reads`: Path to long-read FASTQ file(s).
- `--output`: Path to the output directory where results will be saved.
- `--threads`: Number of threads to use (default: 8).
- `--maxiter`: Maximum number of polishing iterations (default: 4).
- `--nextpolish_path`:  Specifies the path to the NextPolish executable if you wish to use a different version of NextPolish. The default path is the current directory.


## Usage

To run the pipeline, use the following command:

```bash
nextflow run assembly_polish.nf \
  --forward_reads /path/to/forwardReads \
  --reverse_reads /path/to/reverseReads \
  --long_reads /path/to/unpairedReads \
  --output /path/to/outputDir
```

### Help

To view help information, use:

```bash
nextflow run assembly_polish.nf --help
```

## Workflow Processes

1. **UnicyclerAssembly**: Assembles genomes using Unicycler combining both short-read and long-read data. There are no additional parameters set for this process beyond the long and short reads. 
2. **RaconPolishing**: Polishes the assembled genome using Racon through multiple iterations. There are no additional parameters set for this process beyond the long reads and the baseline assembly. 
3. **MedakaPolishing**: Uses Medaka for further polishing of the genome assembled and polished by previous steps. Here the default model for medaka Version 1.7.2 is used. There are no additional parameters set for this process beyond the long reads and the racon polished assembly. 
4. **NextPolishPolishing**: Final polishing using NextPolish to refine assembly based on error correction using paired reads. There are no additional parameters set for this process beyond the short reads and the medaka polished assembly. 

## Outputs

The pipeline will generate polished genome assemblies in the specified output directory, organized into subdirectories for each tool used:

- `unicycler/assembly.fasta`: The initial assembly results from Unicycler.
- `racon_round4.fasta`:  The Racon polished genome files.
- `medaka.fasta`: The medaka polished genome files.
- `genome.nextpolish.fa`: The final polished genome files from NextPolish.


## Troubleshooting

Common issues include path errors, environment misconfigurations, or missing dependencies. Ensure all paths are correct and all required software is properly installed and accessible. Please put an issue ticket if you encounter any unresolvable issues. 

## Contributing

Contributions to the pipeline are welcome. Please fork the repository, make your changes, and submit a pull request.

This README provides a comprehensive overview for running and understanding the Assembly Polish Pipeline. For detailed information on each tool's parameters and capabilities, consult the respective tool documentation.
