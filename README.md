# RiboBA

RiboBA is an R package for bias-aware ORF detection from ribosome profiling (Ribo-seq) data.

## Citation

If you use RiboBA, please cite:

Bai J, Yang R. RiboBA: a bias-aware probabilistic framework for robust ORF identification across diverse ribosome profiling protocols. bioRxiv [Preprint]. 2026.

## Requirements

- R >= 4.1.0
- Dependencies listed in `DESCRIPTION` (`Imports` / `Suggests`)
- For `R/prepare_bam.sh`, install:
  - `bowtie` / `bowtie-build` (Bowtie 1.3)
  - `samtools` (recommended >= 1.10)
  - `curl` (recommended >= 7.0)

## Installation

```r
install.packages("devtools")
devtools::install_github("Bai-JunYu/RiboBA")
```

## Example (Human GRCh38)

Use placeholder paths and replace them with your own paths.

Step 1: Build transcript annotation resources from genome FASTA and GTF/GFF.

```r
library(RiboBA)

build_annotation_index(
  output_dir = "/path/to/output_annotation_dir",
  genome_file = "/path/to/reference_and_input_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa", # reference genome FASTA
  gff_file = "/path/to/reference_and_input_data/Homo_sapiens.GRCh38.109.gtf"                      # annotation GTF/GFF
)
```

Step 2: Prepare transcriptome-aligned BAM inputs from adapter-trimmed FASTQ files.

```bash
bash R/prepare_bam.sh \
  --workdir "/path/to/workdir_for_prepare_bam" \                                      # working directory for indices/logs/BAM
  --fastq-dir "/path/to/adapter_trimmed_fastq_dir" \                                  # input FASTQ directory (adapter-trimmed)
  --annotate-dir "/path/to/output_annotation_dir" \                                   # output_dir used in build_annotation_index()
  --genome-file "/path/to/reference_and_input_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \  # same genome FASTA
  --species "Homo sapiens" \                          # species name used by the script
  --rrna-fa "/path/to/reference_and_input_data/tr_rna.fasta" \                        # rRNA sequences FASTA for depletion step
  --threads 8
```

Step 3: Run the RiboBA pipeline with generated transcript metadata and BAM directory.

```r
res <- run_riboba_pipeline(
  transcript.file = "/path/to/output_annotation_dir/transcript/tx_information.RData", # from build_annotation_index()
  bam_dir = "/path/to/transcriptome_bam_dir"                                         # BAM directory from prepare_bam.sh
)
```
