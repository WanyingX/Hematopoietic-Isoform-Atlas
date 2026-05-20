# Hematopoietic Isoform Atlas

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Cell-type resolved transcript quantification across 36 hematopoietic cell types derived from single-cell multiomics (10x Genomics Multiome, paired scATAC/scRNA-seq).

---

## Overview

This pipeline processes scRNA-seq BAM files from 10x Genomics Multiome to generate cell-type-specific transcript expression profiles. Reads are assigned to individual cell types using single-cell metadata, merged across donors, and quantified using StringTie.

**Live demo:** `https://WanyingX.github.io/Hematopoietic-Isoform-Atlas/`

---


## Pipeline

```
scRNA-seq BAM (per donor)
        ↓
1. PCR deduplication (-F 1024)
        ↓
2. Donor-specific barcode prefixing
        ↓
3. Cell-type splitting (via metadata)
        ↓
4. Cross-donor merging (per cell type)
        ↓
5. StringTie transcript quantification
        ↓
Cell-type resolved isoform expression matrix
```

---

## Scripts

| Script | Description |
|--------|-------------|
| `1_split_celltype.sh` | SLURM array job: dedup, prefix, and split BAM by cell type |
| `2_merge_stringtie.sh` | Merge BAMs per cell type and run StringTie |
| `add_prefix_to_CB.py` | Append donor-specific prefix to cell barcode (CB) tag |
| `split_bam_by_celltype_from_meta_v3.py` | Assign reads to cell types using metadata table |

---

## Requirements

- `samtools` ≥ 1.10
- `StringTie` ≥ 3.0
- `Python` ≥ 3.8
- SLURM workload manager

---

## Cell Types

36 hematopoietic cell types spanning HSC, MPP, erythroid, myeloid, and lymphoid lineages.

| Lineage | Cell Types |
|---------|------------|
| Stem & Progenitor | HSC, MPP_MkEry, MPP_MyLy, LMPP, MLP, MDP, MKP |
| Erythroid | MEP, BFU_E, CFU_E, ProE, BasoE, EoBasoP, OrthoE, PolyE |
| Myeloid | GMP, Early_GMP, ProMono, CD14Mono, CD16Mono, cDC, pDC |
| Lymphoid | CLP, CyclingP, PreProB, VDJ_ProB, Cycling_ProB, ImmatureB, Large_PreB, Small_PreB, MatureB, Plasma, CD4T, CD8T, NK |
| Other | Stromal |

---

## Usage

```bash
# Step 1: Dedup and split by cell type (parallel across donors)
jid=$(sbatch --parsable 1_split_celltype.sh)

# Step 2: Merge and quantify (runs after Step 1 completes)
sbatch --dependency=afterok:$jid 2_merge_stringtie.sh
```

---

## Output

```
celltype_bams/
├── HSC.bam
├── HSC.stringtie.gtf
├── HSC.abundance.tsv
├── ProE.bam
├── ProE.stringtie.gtf
└── ...
```

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contact

Wanying Xu · Sankaran Lab · Boston Children's Hospital  
GitHub: [@WanyingX](https://github.com/WanyingX)
