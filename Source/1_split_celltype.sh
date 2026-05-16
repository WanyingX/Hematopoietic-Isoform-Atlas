#!/bin/bash
#SBATCH --job-name=split_celltype
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=8:00:00
#SBATCH --array=0-25
#SBATCH --output=logs/split_%A_%a.out
#SBATCH --error=logs/split_%A_%a.err

set -euo pipefail

# ══════════════════════════════════════════════════════════════════
# USER CONFIG — modify these paths before running
# ══════════════════════════════════════════════════════════════════
RAWBAM_DIR=/path/to/raw_bam              # directory containing per-donor BAM folders
OUTDIR=/path/to/output/celltype_bam      # output directory
META=/path/to/metadata.txt               # cell metadata with celltype annotation
SCRIPT_DIR=/path/to/scripts              # directory containing helper Python scripts
# ══════════════════════════════════════════════════════════════════

TMPDIR=$OUTDIR/tmp
mkdir -p "$OUTDIR" "$TMPDIR" logs

# ── Sample list ───────────────────────────────────────────────────
# Edit this list to match your donor/sample names
SAMPLES=(
  sample1
  sample2
  sample3
  # add more samples here
)

f=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
inbam="$RAWBAM_DIR/$f/gex_possorted_bam.bam"

echo "=== START: $f at $(date) ==="

if [ ! -f "$inbam" ]; then
  echo "ERROR: $inbam not found"
  exit 1
fi

# ── Step 1: Remove PCR duplicates ────────────────────────────────
# Removes reads flagged as PCR duplicates by Cell Ranger (-F 1024)
if [ ! -f "$TMPDIR/$f.dedup.bam" ]; then
  echo "[$f] Removing duplicates..."
  samtools view -F 1024 -b -@ 4 "$inbam" > "$TMPDIR/$f.dedup.bam"
  samtools index "$TMPDIR/$f.dedup.bam"
else
  echo "[$f] dedup.bam exists, skipping."
fi

# ── Step 2: Add donor-specific barcode prefix ─────────────────────
# Prevents barcode collision when merging across donors
if [ ! -f "$TMPDIR/$f.CBpref.bam" ]; then
  echo "[$f] Adding barcode prefix..."
  python "$SCRIPT_DIR/add_prefix_to_CB.py" \
    --in_bam  "$TMPDIR/$f.dedup.bam" \
    --out_bam "$TMPDIR/$f.CBpref.bam" \
    --prefix  "${f}_" \
    --only_if_missing
  samtools index "$TMPDIR/$f.CBpref.bam"
else
  echo "[$f] CBpref.bam exists, skipping."
fi

# Clean up dedup BAM to save disk space
rm -f "$TMPDIR/$f.dedup.bam" "$TMPDIR/$f.dedup.bam.bai"

# ── Step 3: Split BAM by cell type ───────────────────────────────
# Assigns reads to cell types based on metadata
# Metadata must contain a cell barcode column and a celltype column
if [ ! -d "$TMPDIR/$f" ] || [ -z "$(ls -A $TMPDIR/$f 2>/dev/null)" ]; then
  echo "[$f] Splitting by cell type..."
  mkdir -p "$TMPDIR/$f"
  python "$SCRIPT_DIR/split_bam_by_celltype_from_meta_v3.py" \
    --bam          "$TMPDIR/$f.CBpref.bam" \
    --meta         "$META" \
    --outdir       "$TMPDIR/$f" \
    --celltype_col celltype \
    --sep          $'\t' \
    --drop_no_cb \
    --write_unknown
else
  echo "[$f] Split already done, skipping."
fi

echo "=== DONE: $f at $(date) ==="
