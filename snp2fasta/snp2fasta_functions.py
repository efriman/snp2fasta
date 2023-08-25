import numpy as np
import pandas as pd
import bioframe
import re
import pysam
import itertools

def preprocess_snp_table(snp, flank):
    if not set(["chrom", "start", "ref", "alt"]).issubset(snp.columns):
        raise ValueError("""file needs to be tab-delimited and contain 
        at least columns chrom, start, ref, alt""")
    snp["ref"] = snp["ref"].str.upper()
    snp["alt"] = snp["alt"].str.upper()
    snp["flank_start"] = snp["start"].astype(int)-(flank)
    snp["flank_end"] = snp["start"].astype(int)+(flank)
    snp["snp_pos"] = (snp["start"] - snp["flank_start"]) - 1
    snp["ref_length"] = snp.apply(lambda x: len(x["ref"]), axis=1)
    return snp

def fetch_fa_from_bed_series(series, pysam_fa):
    fa = pysam_fa.fetch(start=series["flank_start"], 
                        end=series["flank_end"], 
                        region=series["chrom"])
    return fa

def slice_str(string, start, end):
    assert isinstance(string, str)
    assert isinstance(start, int)
    assert isinstance(end, int)
    stringout = string[start: end]
    return stringout

def process_snp_fasta_table(data, ignore_char="-"):
    data["seq"] = data["seq"].str.lower()
    data['ref_pos'] = data.apply(lambda x: slice_str(x["seq"], 
                                                     x['snp_pos'], 
                                                     x['snp_pos']+x["ref_length"]), 
                                 axis=1)
    data["alt"] = np.where(data["ref"] == ignore_char,
                             data['ref_pos'].str.upper() + data['alt'],
                             data['alt'])
    data["ref"] = np.where(data["ref"] == ignore_char,
                             data['ref_pos'].str.upper(),
                             data['ref'])
    data["ref_mismatch"] = np.where(data["ref"].str.lower() != data["ref_pos"], True, False) 
    return data

def replace_str(string, pos, length, replacewith, ignore_char="-"):
    assert isinstance(string, str)
    assert isinstance(pos, int)
    assert isinstance(length, int)
    assert isinstance(replacewith, str)
    replacewith = replacewith.replace(ignore_char, "")
    stringout = string[:pos] + replacewith + string[(pos + length):]
    return stringout

def trim_strings(string1, string2, length):
    assert isinstance(string1, str)
    assert isinstance(string2, str)
    assert isinstance(length, int)
    trim = 0
    if len(string1) > length:
        trim = len(string1) - length
    elif len(string2) > length:
        trim = len(string2) - length
    trim = int(np.ceil(trim / 2))
    string1 = string1[trim:len(string1)-trim]
    string2 = string2[trim:len(string2)-trim]
    return string1, string2

def series_to_fasta(series, header_cols):
    trim_cols = [col for col in header_cols if col not in ["ref", "alt"]]
    header = ">" + "_".join(series[trim_cols].astype(str))
    fasta = ""
    for seq in ["seq", "seq_alt"]:
        refalt = "ref" if seq == "seq" else "alt"
        fasta = fasta + f"{header}_{series[refalt]}\n{series[seq]}\n"
    return fasta

def extract_closest(df, flank, k):
    closest = bioframe.closest(df, df.copy(), k=k, suffixes=["1", "2"], ignore_upstream=True)
    closest = closest[~((closest["distance"] == 0) & 
                      (closest["start1"] == closest["start2"]) &
                      (closest["alt1"] == closest["alt2"])) &
                      (closest["distance"] < flank)].reset_index(drop=True)
    closest["maxstart"] = closest.groupby(["chrom1", "start1", "ref1"])["start2"].transform("max")
    closest["center"] = np.ceil((closest["start1"] + closest["maxstart"])/2).astype(int)
    closest["flank_start"] = closest["center"] - flank
    closest["flank_end"] = closest["center"] + flank
    return closest

def extract_combinations_fasta(df, flank, trim=True):
    if len(df["seq"].unique()) > 1:
        raise ValueError("Non-unique sequences in DataFrame")
    df = df.reset_index(drop=True)
    flank_start = df.loc[0,"flank_start"]
    sequence = df.loc[0, "seq"]
    df2 = pd.concat([df[["chrom1","start1", "end1", "ref1", "alt1", "rsID1", "seq"]].drop_duplicates().rename(columns=lambda x: re.sub('1','',x)), 
                     df[["chrom2","start2", "end2", "ref2", "alt2", "rsID2", "seq"]].drop_duplicates().rename(columns=lambda x: re.sub('2','',x))]).reset_index(drop=True)
    df2["pos"] = (df2["start"] - flank_start) - 1
    df2["ref_length"] = df2.apply(lambda x: len(x["ref"]), axis=1)
    positions = list(df2["pos"].astype(int))
    length = df2.shape[0]
    allele_combinations = itertools.product(*[df2[["ref", "alt"]].iloc[i].tolist() for i in range(length)])
    reflengths = list(df2["ref_length"])
    df2["coord"] = df2["chrom"] + "_" + df2["start"].astype(str)
    coords = list(df2["coord"])
    sequences = ""
    insertion = 0
    if trim:
        df2["alt_length"] = df2.apply(lambda x: len(x["alt"]), axis=1)
        insertions = np.array(df2["alt_length"] - df2["ref_length"])
        insertion = int(np.ceil(insertions[insertions>0].sum()/2))
    for rep in allele_combinations:
        sequence_mod = sequence
        header = ">"
        for i in range(length):
            pos = positions[i] + (len(sequence_mod)-2*flank)
            sequence_mod = replace_str(sequence_mod, pos, reflengths[i], rep[i])
            header = header + "_" + coords[i] + "_" + rep[i]
        header = header.replace(">_", ">")
        sequence_mod = sequence_mod[insertion:len(sequence_mod)-insertion]
        sequences =  sequences + f"{header}\n{sequence_mod}\n"
    return sequences
