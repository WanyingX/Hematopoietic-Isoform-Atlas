#!/usr/bin/env python3
import argparse
import pysam

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_bam", required=True)
    ap.add_argument("--out_bam", required=True)
    ap.add_argument("--prefix", required=True, help='e.g. "aged1bmmc_"')
    ap.add_argument("--only_if_missing", action="store_true",
                    help="only add prefix if CB doesn't already start with prefix")
    args = ap.parse_args()

    ibam = pysam.AlignmentFile(args.in_bam, "rb")
    obam = pysam.AlignmentFile(args.out_bam, "wb", template=ibam)

    total = kept = changed = 0
    for r in ibam.fetch(until_eof=True):
        total += 1
        if r.has_tag("CB"):
            cb = r.get_tag("CB")
            if (not args.only_if_missing) or (not cb.startswith(args.prefix)):
                r.set_tag("CB", args.prefix + cb, value_type="Z")
                changed += 1
        obam.write(r)
        kept += 1

    ibam.close()
    obam.close()
    pysam.index(args.out_bam)
    print(f"[done] total={total} written={kept} CB_changed={changed} out={args.out_bam}")

if __name__ == "__main__":
    main()
