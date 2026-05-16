#!/bin/bash
#SBATCH --job-name=merge_stringtie
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=logs/merge_%j.out
#SBATCH --error=logs/merge_%j.err

set -euo pipefail

# ══════════════════════════════════════════════════════════════════
# USER CONFIG — modify these paths before running
# ══════════════════════════════════════════════════════════════════
OUTDIR=/path/to/output/celltype_bam      # same as in 1_split_celltype.sh
REF_GTF=/path/to/genome_annotation.gtf   # reference genome annotation
RAWBAM_DIR=/path/to/raw_bam              # used to detect cell type list
# ══════════════════════════════════════════════════════════════════

TMPDIR=$OUTDIR/tmp
mkdir -p "$OUTDIR" logs

# ── Step 4: Merge BAMs per cell type across donors ────────────────
echo "=== Merging by cell type at $(date) ==="

first_sample=$(ls "$TMPDIR" | grep -v "\.bam" | head -1)

for ct_bam in "$TMPDIR/$first_sample"/*.bam; do
  ct=$(basename "$ct_bam" .bam)
  
  bam_list=$(ls "$TMPDIR"/*/"${ct}.bam" 2>/dev/null || true)
  
  if [ -z "$bam_list" ]; then
    echo "No BAM found for cell type: $ct, skipping."
    continue
  fi
  
  if [ ! -f "$OUTDIR/${ct}.bam" ]; then
    echo "Merging: $ct"
    samtools merge -f -@ 4 "$OUTDIR/${ct}.bam" $bam_list
    samtools index "$OUTDIR/${ct}.bam"
  else
    echo "$ct.bam exists, skipping."
  fi
done

# ── Step 5: StringTie transcript quantification ───────────────────
echo "=== Running StringTie at $(date) ==="

for bam in "$OUTDIR"/*.bam; do
  name=$(basename "$bam" .bam)
  
  if [ ! -f "$OUTDIR/$name.stringtie.gtf" ]; then
    echo "StringTie: $name"
    stringtie "$bam" \
      -G "$REF_GTF" \
      -p 4 \
      -o "$OUTDIR/$name.stringtie.gtf" \
      -A "$OUTDIR/$name.abundance.tsv"
  else
    echo "$name.stringtie.gtf exists, skipping."
  fi
done

# ── Step 6: Clean up temporary files ─────────────────────────────
rm -rf "$TMPDIR"
echo "=== ALL DONE at $(date) ==="