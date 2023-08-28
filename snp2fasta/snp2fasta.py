#!/usr/bin/env python3

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
        "--combinations",
        type=int,
        default=1,
        required=False,
        help="""Specify >1 to generate up to this number of allele combinations when SNPs are within a certain distance""",
    )
    parser.add_argument(
        "--maxdist",
        type=int,
        required=False,
        help="""When using combinations, maximum distance between SNPs to combine. By default, uses flank""",
    )
    parser.add_argument(
        "--ignore_char",
        type=str,
        default="-",
        required=False,
        help="""Character used to signify deletions""",
    )
    parser.add_argument(
        "--no_trim",
        action="store_true",
        default=False,
        required=False,
        help="""Set to not trim the ends to maintain a maximum size of 2*flank""",
    )
    return parser
    
def main():
    parser = parse_args_overlap_peaks()
    args = parser.parse_args()

    logging.debug(args)

    data = pd.read_table(args.input_table)
    
    flank = args.flank
    
    header_cols = data.columns
    
    snp = preprocess_snp_table(data, flank=flank)
    
    fasta = pysam.FastaFile(args.fasta)
    
    snp["seq"] = snp.apply(lambda x: fetch_fa_from_bed_series(x, fasta), axis=1)
    
    snp = process_snp_fasta_table(snp, ignore_char=args.ignore_char)
    
    snp_mismatch = snp.loc[snp["ref_mismatch"]].reset_index(drop=True)
    
    logging.info(f"{snp_mismatch.shape[0]} mismatched reference alleles compared to genome sequence")
    if not snp_mismatch.empty:
        snp_mismatch.to_csv(f"{args.outdir}/{args.outname}_mismatched_refs.txt", sep="\t", index=False)
        logging.info(f"Saved as {args.outdir}/{args.outname}_mismatched.txt")
    
    snp_match = snp.loc[~snp["ref_mismatch"]].reset_index(drop=True)

    snp_match['seq'] = snp_match.apply(lambda x: replace_str(x["seq"], 
                                                             x['snp_pos'], 
                                                             x["ref_length"], 
                                                             x["ref"], 
                                                             ignore_char=args.ignore_char), axis=1)
    
    snp_match['seq_alt'] = snp_match.apply(lambda x: replace_str(x["seq"], 
                                                                 x['snp_pos'], 
                                                                 x["ref_length"], 
                                                                 x["alt"], 
                                                                 ignore_char=args.ignore_char), axis=1)

    if not args.no_trim:    
        snp_match[['seq', 'seq_alt']] = snp_match.apply(lambda x: trim_strings(x["seq"], x['seq_alt'], flank*2), axis=1, result_type="expand")
    
    snp_match["fasta"] = snp_match.apply(lambda x: series_to_fasta(x, header_cols=header_cols), axis=1)
    
    fasta_out = "".join(snp_match["fasta"])
    
    text_file = open(f"{args.outdir}/{args.outname}_matched.fa", "w")
    text_file.write(fasta_out)
    text_file.close()
    
    logging.info(f"Saved {snp_match.shape[0]} entries as {args.outdir}/{args.outname}_matched.fa")

    if args.combinations > 1:
        if args.maxdist:
            if args.maxdist > 2*args.flank:
                raise ValueError(f"maxdist can't be more than 2*flank")
            else:
                maxdist = args.maxdist
        else:
            maxdist = args.flank
        closest = snp_match[header_cols]
        closest = extract_closest(closest, flank=args.flank, maxdist=maxdist, k=(args.combinations))
        closest["seq"] = closest.rename(columns=lambda x: re.sub('1','',x)).apply(lambda x: fetch_fa_from_bed_series(x, fasta), axis=1)
        closest["seq"] = closest["seq"].str.lower()
        closest["id"] = closest["chrom1"] + "_" + closest["start1"].astype(str) + "_" + closest["alt1"]
        ncomb = len(closest["id"].unique())
        logging.info(f"{ncomb} instances of up to {args.combinations} combinations within {maxdist} bp")
        fasta_comb = ""
        for id in closest["id"].unique():
            fa_id = extract_combinations_fasta(closest[closest["id"] == id], flank=args.flank, trim=not args.no_trim)
            for fa_entry in fa_id.split("\n"):
                if fa_entry not in fasta_comb:
                    fasta_comb = f"{fasta_comb}{fa_entry}\n"
    
        text_file = open(f"{args.outdir}/{args.outname}_combinations.fa", "w")
        text_file.write(fasta_comb)
        text_file.close()
    
        logging.info(f"Saved combinations as {args.outdir}/{args.outname}_combinations.fa")

if __name__ == "__main__":
    main()