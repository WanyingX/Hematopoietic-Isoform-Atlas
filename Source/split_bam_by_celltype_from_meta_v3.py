#!/usr/bin/env python3
import argparse
import os
import re
import pysam

def sanitize(name: str) -> str:
    name = name.strip()
    name = re.sub(r"[^\w\.\-]+", "_", name)
    return name if name else "NA"

def split_tsv(line: str, sep: str):
    return line.rstrip("\n").split(sep)

def load_map(meta_path: str, celltype_col: str, sep: str = "\t"):
    """
    Handles two common exports:

    Case 1 (no rownames column):
      header: gex_barcode ... source celltype
      row:    young1t1bmmc_AAAC...-1 ... redeem CD4T
      -> barcode = column 'gex_barcode'

    Case 2 (rownames column exists with NO header field):
      header: gex_barcode ... source celltype        (N columns)
      row:    young1t1bmmc_AAAC...-1  AAAC...-1 ... redeem CD4T   (N+1 columns)
            ^rownames(barcode)       ^gex_barcode col
      -> barcode = FIRST field (rownames)
      -> all named columns are shifted by +1 in the data rows
    """

    with open(meta_path) as f:
        header_line = f.readline()
        if not header_line:
            raise RuntimeError("Empty metadata file")

        header = split_tsv(header_line, sep)
        if celltype_col not in header:
            raise RuntimeError(f"celltype_col='{celltype_col}' not found in header: {header}")
        ci_header = header.index(celltype_col)

        # Peek first non-empty data line to detect extra rownames column
        first_data = None
        pos = f.tell()
        for line in f:
            if line.strip():
                first_data = line
                break
        if first_data is None:
            raise RuntimeError("Metadata has header but no data rows")
        f.seek(pos)

        parts0 = split_tsv(first_data, sep)
        extra_rownames = (len(parts0) == len(header) + 1)

        # Indices in data rows
        if extra_rownames:
            barcode_idx = 0                     # rownames barcode
            celltype_idx = ci_header + 1        # shift by +1
            mode = "EXTRA_ROWNAME_COL"
        else:
            # barcode expected in 'gex_barcode' column
            if "gex_barcode" not in header:
                raise RuntimeError("No 'gex_barcode' column and no extra rowname col detected.")
            barcode_idx = header.index("gex_barcode")
            celltype_idx = ci_header
            mode = "NORMAL"

        bc2ct = {}
        examples = []
        bad = 0

        for line in f:
            if not line.strip():
                continue
            parts = split_tsv(line, sep)
            if len(parts) <= max(barcode_idx, celltype_idx):
                bad += 1
                continue
            bc = parts[barcode_idx].strip()
            ct = parts[celltype_idx].strip()
            if bc and ct:
                sct = sanitize(ct)
                bc2ct[bc] = sct
                if len(examples) < 10:
                    examples.append(sct)

    return bc2ct, examples, bad, mode

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--celltype_col", default="celltype")
    ap.add_argument("--sep", default="\t")
    ap.add_argument("--cb_tag", default="CB")
    ap.add_argument("--drop_no_cb", action="store_true")
    ap.add_argument("--write_unknown", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    bc2ct, examples, bad, mode = load_map(args.meta, args.celltype_col, args.sep)
    print(f"[info] meta_mode={mode} loaded_barcodes={len(bc2ct)} bad_lines_skipped={bad}")
    print(f"[info] example celltypes: {examples}")

    ibam = pysam.AlignmentFile(args.bam, "rb")

    out_handles = {}
    counts = {}
    total = written = no_cb = unknown = 0

    def get_handle(ct: str):
        if ct not in out_handles:
            out_path = os.path.join(args.outdir, f"{ct}.bam")
            out_handles[ct] = pysam.AlignmentFile(out_path, "wb", template=ibam)
            counts[ct] = 0
        return out_handles[ct]

    for r in ibam.fetch(until_eof=True):
        total += 1
        if not r.has_tag(args.cb_tag):
            no_cb += 1
            if args.drop_no_cb:
                continue
            else:
                continue

        cb = r.get_tag(args.cb_tag)
        ct = bc2ct.get(cb)

        if ct is None:
            unknown += 1
            if args.write_unknown:
                h = get_handle("UNKNOWN")
                h.write(r)
                counts["UNKNOWN"] += 1
                written += 1
            continue

        h = get_handle(ct)
        h.write(r)
        counts[ct] += 1
        written += 1

    ibam.close()
    for h in out_handles.values():
        h.close()

    print(f"[done] total_reads={total} written={written} no_CB={no_cb} unknown_CB={unknown}")
    top = sorted(counts.items(), key=lambda x: x[1], reverse=True)[:15]
    print("[top celltypes by reads]")
    for ct, n in top:
        print(ct, n)

if __name__ == "__main__":
    main()
