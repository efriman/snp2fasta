#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pysam
from snp2fasta.snp2fasta_functions import *
import argparse
import logging
import warnings
logging.basicConfig(format="%(message)s", level="INFO")


def parse_args_overlap_peaks():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input_table", 
        type=str, 
        help="Tab-delimited input table containing at least columns chrom, start, ref, alt"
    )
    parser.add_argument(
        "--flank",
        "--flanks",
        type=int,
        required=True,
        help="""The amount of flanking reference genome you want to append surrounding the SNPs""",
    )
    parser.add_argument(
        "--outname",
        type=str,
        required=True,
        help="""Prefix for output files. This should NOT be the path to a directory (set with --outdir)""",
    )
    parser.add_argument(
        "--fasta",
        type=str,
        required=True,
        help="""Path to fasta file""",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=".",
        required=False,
        help="""Directory to save in. Defaults to current""",
    )
    parser.add_argument(
        "--ignore_char",
        type=str,
        default="-",
        required=False,
        help="""Character used to signify deletions""",
    )
    return parser
    
def main():
    parser = parse_args_overlap_peaks()
    args = parser.parse_args()

    logging.debug(args)

    data = pd.read_table(args.input_table)
    
    flank = 85
    
    header_cols = data.columns
    
    snp = preprocess_snp_table(data, flank=flank)
    
    fasta = pysam.FastaFile(args.fasta)
    
    snp["seq"] = snp.apply(lambda x: fetch_fa_from_bed_series(x, fasta), axis=1)
    
    snp = process_snp_fasta_table(snp)
    
    snp_mismatch = snp.loc[snp["ref_mismatch"]].reset_index(drop=True)
    
    print(f"{snp_mismatch.shape[0]} mismatched reference alleles compared to genome sequence")
    if not snp_mismatch.empty:
        snp_mismatch.to_csv(f"{args.outname}_mismatched_refs.txt", sep="\t", index=False)
        print(f"Saved as {args.outname}_mismatched.txt")
    
    snp_match = snp.loc[~snp["ref_mismatch"]].reset_index(drop=True)
    
    snp_match['seq_alt'] = snp_match.apply(lambda x: replace_str(x["seq"], 
                                                                 x['snp_pos'], 
                                                                 x["ref_length"], 
                                                                 x["alt"], 
                                                                 ignore_char=args.ignore_char), axis=1)
    
    snp_match[['seq', 'seq_alt']] = snp_match.apply(lambda x: trim_strings(x["seq"], x['seq_alt'], flank*2), axis=1, result_type="expand")
    
    snp_match["fasta"] = snp_match.apply(lambda x: series_to_fasta(x, header_cols=header_cols), axis=1)
    
    fasta_out = "".join(snp_match["fasta"])
    
    text_file = open(f"{args.outname}_matched.fa", "w")
    text_file.write(fasta_out)
    text_file.close()
    
    print(f"Saved {snp_match.shape[0]} entries as {args.outname}_matched.fa")

if __name__ == "__main__":
    main()