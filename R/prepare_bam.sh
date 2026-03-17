#!/usr/bin/env bash
set -euo pipefail

# No circRNA branch in this pipeline.
# Steps:
# 1) Build Bowtie indices (rRNA, CDS, lncRNA)
# 2) Remove rRNA reads from adapter-trimmed FASTQ
# 3) Map remaining reads to CDS (keep unmapped as no-CDS reads)
# 4) Map no-CDS reads to lncRNA
# 5) Convert SAM to sorted BAM and build BAM index
#
# Usage:
#   bash R/prepare_bam.sh \
#     --config /path/to/prepare_bam.yaml \
#     --workdir /path/to/work \
#     --fastq-dir /path/to/adapter_trimmed_fastq \
#     --annotate-dir /path/to/annotate_out \
#     --genome-file /path/to/genome.fa \
#     --species "Homo sapiens" \
#     --rrna-fa /path/to/rrna.fa \
#     --cds-fa /path/to/longest_txs.fa \
#     --lncrna-fa /path/to/lncRNA_txs.fa \
#     --threads 12

CONFIG_FILE=""
WORKDIR=""
FASTQ_DIR=""
ANNOTATE_DIR=""
GENOME_FILE=""
SPECIES=""
RRNA_FA=""
CDS_FA=""
LNCRNA_FA=""
THREADS="12"

read_config_value() {
  local key="$1"
  local file="$2"
  awk -v key="$key" '
    BEGIN { FS=":" }
    /^[[:space:]]*#/ { next }
    /^[[:space:]]*$/ { next }
    {
      raw=$0
      sub(/^[[:space:]]*/, "", raw)
      if (index(raw, key ":") != 1) next
      val=substr(raw, length(key)+2)
      sub(/^[[:space:]]*/, "", val)
      sub(/[[:space:]]*$/, "", val)
      gsub(/^["'"'"']|["'"'"']$/, "", val)
      print val
      exit
    }
  ' "$file"
}

load_config() {
  local cfg="$1"
  if [[ ! -f "$cfg" ]]; then
    echo "Config file not found: $cfg" >&2
    exit 1
  fi
  [[ -z "$WORKDIR" ]] && WORKDIR="$(read_config_value workdir "$cfg")"
  [[ -z "$FASTQ_DIR" ]] && FASTQ_DIR="$(read_config_value fastq_dir "$cfg")"
  [[ -z "$ANNOTATE_DIR" ]] && ANNOTATE_DIR="$(read_config_value annotate_dir "$cfg")"
  [[ -z "$GENOME_FILE" ]] && GENOME_FILE="$(read_config_value genome_file "$cfg")"
  [[ -z "$SPECIES" ]] && SPECIES="$(read_config_value species "$cfg")"
  [[ -z "$RRNA_FA" ]] && RRNA_FA="$(read_config_value rrna_fa "$cfg")"
  [[ -z "$CDS_FA" ]] && CDS_FA="$(read_config_value cds_fa "$cfg")"
  [[ -z "$LNCRNA_FA" ]] && LNCRNA_FA="$(read_config_value lncrna_fa "$cfg")"
  local cfg_threads
  cfg_threads="$(read_config_value threads "$cfg")"
  [[ -n "$cfg_threads" && "$THREADS" == "12" ]] && THREADS="$cfg_threads"
}

infer_transcript_fastas() {
  if [[ -n "$ANNOTATE_DIR" ]]; then
    [[ -z "$CDS_FA" ]] && CDS_FA="$ANNOTATE_DIR/transcript/longest_txs.fa"
    [[ -z "$LNCRNA_FA" ]] && LNCRNA_FA="$ANNOTATE_DIR/transcript/lncRNA_txs.fa"
  fi
}

infer_species_from_genome_file() {
  local genome_file="$1"
  local base no_ext match
  [[ -z "$genome_file" ]] && return 1
  base="$(basename "$genome_file")"
  no_ext="$(echo "$base" | sed -E 's/\.(fa|fasta|fa\.gz|fasta\.gz|2bit)$//I')"
  match="$(echo "$no_ext" | sed -nE 's/^([A-Z][a-z]+)[._]([a-z]+).*/\1 \2/p' | head -n 1)"
  if [[ -n "$match" ]]; then
    echo "$match"
    return 0
  fi
  return 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG_FILE="$2"; shift 2 ;;
    --workdir) WORKDIR="$2"; shift 2 ;;
    --fastq-dir) FASTQ_DIR="$2"; shift 2 ;;
    --annotate-dir) ANNOTATE_DIR="$2"; shift 2 ;;
    --genome-file) GENOME_FILE="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --rrna-fa) RRNA_FA="$2"; shift 2 ;;
    --cds-fa) CDS_FA="$2"; shift 2 ;;
    --lncrna-fa) LNCRNA_FA="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

if [[ -n "$CONFIG_FILE" ]]; then
  load_config "$CONFIG_FILE"
fi

infer_transcript_fastas
if [[ -z "$SPECIES" ]]; then
  SPECIES="$(infer_species_from_genome_file "$GENOME_FILE" || true)"
  if [[ -n "$SPECIES" ]]; then
    echo "Inferred species from genome_file: $SPECIES"
  fi
fi

if [[ -z "$WORKDIR" || -z "$FASTQ_DIR" || -z "$RRNA_FA" || -z "$CDS_FA" || -z "$LNCRNA_FA" ]]; then
  echo "Missing required arguments. Need: --workdir --fastq-dir --rrna-fa and CDS/lncRNA FASTA (or --annotate-dir)." >&2
  exit 1
fi
for f in "$CDS_FA" "$LNCRNA_FA"; do
  if [[ ! -f "$f" ]]; then
    echo "Input FASTA not found: $f" >&2
    exit 1
  fi
done

check_dependencies() {
  local missing=()
  local cmd
  for cmd in bowtie bowtie-build samtools curl; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
      missing+=("$cmd")
    fi
  done

  if [[ ${#missing[@]} -gt 0 ]]; then
    echo "[Dependency error] Missing required command(s): ${missing[*]}" >&2
    echo "Please install Bowtie(1), Samtools, and curl, then re-run." >&2
    echo "Examples:" >&2
    echo "  conda install -c bioconda bowtie samtools" >&2
    echo "  # or use your system package manager if preferred" >&2
    exit 127
  fi

  echo "[Dependency check] OK"
}

check_dependencies

mkdir -p "$WORKDIR"/{index/tx_rRNA,index/tx_cds,index/tx_lncRNA,no_rrna,map_cds,map_lncrna,bam}

download_rrna_fasta_if_missing() {
  local rrna_fa="$1"
  local species="$2"
  local out_dir="$3"
  local esearch_xml ids_csv

  if [[ -f "$rrna_fa" ]]; then
    return 0
  fi
  if [[ -z "$species" ]]; then
    echo "rrna_fa not found: $rrna_fa" >&2
    echo "Please provide --species (or --genome-file for auto-inference) to download rRNA FASTA from NCBI automatically." >&2
    exit 1
  fi

  echo "[rRNA] rrna_fa not found. Downloading from NCBI for species: $species"
  mkdir -p "$out_dir"
  mkdir -p "$(dirname "$rrna_fa")"
  esearch_xml="$out_dir/rrna_esearch.xml"

  curl -sG "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi" \
    --data-urlencode "db=nuccore" \
    --data-urlencode "retmode=xml" \
    --data-urlencode "retmax=500" \
    --data-urlencode "term=${species}[Organism] AND biomol_rrna[PROP]" \
    > "$esearch_xml"

  ids_csv="$(
    awk -F'[<>]' '/<Id>[0-9]+<\/Id>/{print $3}' "$esearch_xml" | paste -sd, -
  )"
  if [[ -z "$ids_csv" ]]; then
    echo "Failed to find NCBI rRNA records for species: $species" >&2
    exit 1
  fi

  curl -sG "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" \
    --data-urlencode "db=nuccore" \
    --data-urlencode "rettype=fasta" \
    --data-urlencode "retmode=text" \
    --data-urlencode "id=${ids_csv}" \
    > "$rrna_fa"

  if [[ ! -s "$rrna_fa" ]]; then
    echo "Downloaded rRNA FASTA is empty: $rrna_fa" >&2
    exit 1
  fi
  echo "[rRNA] Downloaded: $rrna_fa"
}

download_rrna_fasta_if_missing "$RRNA_FA" "$SPECIES" "$WORKDIR/index/tx_rRNA"

index_complete() {
  local prefix="$1"
  if [[ -f "${prefix}.1.ebwt" && -f "${prefix}.2.ebwt" && -f "${prefix}.rev.1.ebwt" && -f "${prefix}.rev.2.ebwt" ]]; then
    return 0
  fi
  if [[ -f "${prefix}.1.ebwtl" && -f "${prefix}.2.ebwtl" && -f "${prefix}.rev.1.ebwtl" && -f "${prefix}.rev.2.ebwtl" ]]; then
    return 0
  fi
  return 1
}

index_any_exists() {
  local prefix="$1"
  compgen -G "${prefix}"'*.ebwt*' >/dev/null 2>&1
}

build_index_if_missing() {
  local fa="$1"
  local prefix="$2"
  local label="$3"
  local log_file="${prefix}.build.log"
  if index_complete "$prefix"; then
    echo "  - Reusing existing ${label} index: ${prefix}"
  elif index_any_exists "$prefix"; then
    echo "  - Found existing but incomplete ${label} index files for prefix: ${prefix}" >&2
    echo "    To avoid repeated rebuilding, this script will not run bowtie-build automatically." >&2
    echo "    Please clean ${prefix}*.ebwt* manually if you want to rebuild." >&2
    exit 1
  else
    echo "  - Building ${label} index: ${prefix} (log: ${log_file})"
    if ! bowtie-build "$fa" "$prefix" >"$log_file" 2>&1; then
      echo "    Failed building ${label} index. Check log: ${log_file}" >&2
      exit 1
    fi
  fi
}

echo "[Input] FASTQ directory should contain adapter-trimmed reads only: $FASTQ_DIR"
shopt -s nullglob
FASTQ_FILES=( "$FASTQ_DIR"/*.fq "$FASTQ_DIR"/*.fastq "$FASTQ_DIR"/*.fq.gz "$FASTQ_DIR"/*.fastq.gz )
if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
  echo "No FASTQ files found in: $FASTQ_DIR" >&2
  exit 1
fi
echo "[Start] prepare_bam.sh with ${#FASTQ_FILES[@]} FASTQ file(s)"

echo "[1/5] Building indices..."
build_index_if_missing "$RRNA_FA" "$WORKDIR/index/tx_rRNA/rrna" "rRNA"
build_index_if_missing "$CDS_FA" "$WORKDIR/index/tx_cds/cds" "CDS"
build_index_if_missing "$LNCRNA_FA" "$WORKDIR/index/tx_lncRNA/lncrna" "lncRNA"

echo "[2/5] Removing rRNA reads..."
for fq in "${FASTQ_FILES[@]}"; do
  base="$(basename "$fq")"
  sample="${base%%.*}"
  echo "  [rRNA] ${sample}"
  bowtie -p "$THREADS" -v 2 \
    --un "$WORKDIR/no_rrna/${sample}_norrna.fq" \
    -x "$WORKDIR/index/tx_rRNA/rrna" \
    "$fq" \
    > /dev/null 2>> "$WORKDIR/no_rrna/rrna.log"
done

echo "[3/5] Mapping to CDS..."
for fq in "$WORKDIR/no_rrna/"*_norrna.fq; do
  [[ -e "$fq" ]] || continue
  base="$(basename "$fq")"
  sample="${base%_norrna.fq}"
  echo "  [CDS] ${sample}"
  bowtie -p "$THREADS" -v 1 -a -m 10 --best --strata -S \
    --un "$WORKDIR/map_cds/${sample}_nocds.fq" \
    --no-unal \
    -x "$WORKDIR/index/tx_cds/cds" \
    "$fq" \
    > "$WORKDIR/map_cds/${sample}.sam" 2>> "$WORKDIR/map_cds/cds.log"
done

echo "[4/5] Mapping no-CDS reads to lncRNA..."
for fq in "$WORKDIR/map_cds/"*_nocds.fq; do
  [[ -e "$fq" ]] || continue
  base="$(basename "$fq")"
  sample="${base%_nocds.fq}"
  echo "  [lncRNA] ${sample}"
  bowtie -p "$THREADS" -v 1 -a -m 10 --best --strata -S \
    --no-unal \
    -x "$WORKDIR/index/tx_lncRNA/lncrna" \
    "$fq" \
    > "$WORKDIR/map_lncrna/${sample}.sam" 2>> "$WORKDIR/map_lncrna/lncrna.log"
done

echo "[5/5] Converting SAM to sorted BAM + index..."
for sam in "$WORKDIR/map_cds/"*.sam; do
  [[ -e "$sam" ]] || continue
  base="$(basename "$sam" .sam)"
  out="$WORKDIR/bam/${base}_cds_txsorted.bam"
  echo "  [BAM:CDS] ${base}"
  samtools view -@ "$THREADS" -bS "$sam" | samtools sort -@ "$THREADS" -o "$out" -
  samtools index -@ "$THREADS" "$out"
done

for sam in "$WORKDIR/map_lncrna/"*.sam; do
  [[ -e "$sam" ]] || continue
  base="$(basename "$sam" .sam)"
  out="$WORKDIR/bam/${base}_lncrna_txsorted.bam"
  echo "  [BAM:lncRNA] ${base}"
  samtools view -@ "$THREADS" -bS "$sam" | samtools sort -@ "$THREADS" -o "$out" -
  samtools index -@ "$THREADS" "$out"
done

echo "Done. Outputs:"
echo "  CDS SAM:      $WORKDIR/map_cds/*.sam"
echo "  lncRNA SAM:   $WORKDIR/map_lncrna/*.sam"
echo "  sorted BAM:   $WORKDIR/bam/*_txsorted.bam"
