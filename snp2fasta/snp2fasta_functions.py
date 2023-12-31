import numpy as np
import pandas as pd
import bioframe
import re
import pysam
import itertools
import logging

def preprocess_snp_table(snp, flank):
    if not set(["chrom", "start", "ref", "alt"]).issubset(snp.columns):
        raise ValueError("""file needs to be tab-delimited and contain 
        at least columns chrom, start, ref, alt""")
    snp = snp.sort_values(["chrom", "start", "ref", "alt"])
    snp["ref"] = snp["ref"].str.upper()
    snp["alt"] = snp["alt"].str.upper()
    snp["flank_start"] = snp["start"].astype(int)-(flank) - 1
    snp["snp_pos"] = (snp["start"] - snp["flank_start"]) - 1
    snp["ref_length"] = snp.apply(lambda x: len(x["ref"]), axis=1)
    snp["flank_end"] = snp["start"].astype(int)+(flank) + snp["ref_length"] - 1
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
    data['seq_pos'] = data.apply(lambda x: slice_str(x["seq"], 
                                                     x['snp_pos'], 
                                                     x['snp_pos']+x["ref_length"]), 
                                 axis=1)
    data["alt"] = np.where(data["ref"] == ignore_char,
                           data['seq_pos'].str.upper() + data['alt'],
                           data['alt'])
    data["ref"] = np.where(data["ref"] == ignore_char,
                           data['seq_pos'].str.upper(),
                           data['ref'])
    data["ref_mismatch"] = np.where(data["ref"].str.lower() != data["seq_pos"], True, False) 
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

def extract_closest(df, flank, maxdist, k, fasta):
    closest = bioframe.closest(df, df.copy(), k=k, suffixes=["1", "2"], ignore_upstream=True)
    closest = closest[~((closest["distance"] == 0) & 
                      (closest["start1"] >= closest["start2"])) &
                      (closest["distance"] < maxdist)].reset_index(drop=True)
    overlapping = closest.loc[((closest["start1"] + closest["ref_length1"]) > closest["start2"]) |
                              (closest["start1"] == closest["start2"]) , ["chrom1", "start1", "start2"]]
    closest = closest.iloc[[idx for idx in list(closest.index) if idx not in list(overlapping.index)]]
    overlapping = overlapping.drop_duplicates().reset_index(drop=True)
    closest["maxstart"] = closest.groupby(["chrom1", "start1", "ref1", "alt1"])["start2"].transform("max")
    closest["center"] = np.ceil((closest["start1"] + closest["maxstart"])/2).astype(int)
    closest["flank_start"] = closest["center"] - flank
    closest["flank_end"] = closest["center"] + flank
    closest["seq"] = closest.rename(columns=lambda x: re.sub('1','',x)).apply(lambda x: fetch_fa_from_bed_series(x, fasta), axis=1)
    closest["seq"] = closest["seq"].str.lower()
    closest["id"] = closest["chrom1"] + "_" + closest["start1"].astype(str) + "_" + closest["ref1"].astype(str) + "_" + closest["alt1"].astype(str)
    return closest, overlapping

def concat_overlaps(df):
    if len(df["seq"].unique()) > 1:
        raise ValueError("Non-unique sequences in DataFrame")
    df = df.reset_index(drop=True)
    flank_start = df.loc[0,"flank_start"]
    df2 = pd.concat([df[["chrom1", "start1", "ref1", "alt1", "ref_length1", "seq"]].rename(columns=lambda x: re.sub('1','',x)), 
                     df[["chrom2", "start2", "ref2", "alt2", "ref_length2", "seq"]].rename(columns=lambda x: re.sub('2','',x)),
                    ]
                   ).sort_values(["chrom", "start"]).drop_duplicates().reset_index(drop=True)
    df2["pos"] = (df2["start"] - flank_start) - 1
    df2["previous_end"] = df2['pos'].shift(1) + df2['ref_length'].shift(1) - 1
    df2["coord"] = df2["chrom"] + "_" + df2["start"].astype(str)
    return df2

def generate_combination_seq(df, flank, fasta):
    assert set(['chrom','start', 'ref', 'alt']).issubset(df.columns)
    assert df.shape[0] > 1
    df = df.sort_values(["chrom", "start", "ref", "alt"]).reset_index(drop=True)
    df["previous_end"] = df['pos'].shift(1) + df['ref_length'].shift(1) - 1
    assert df[(df["pos"] <= df["previous_end"].fillna(0))].shape[0] == 0
    df["flank_start"] = np.ceil((min(df["start"]) + max(df["start"]))/2).astype(int) - flank
    df["flank_end"] = np.ceil((min(df["start"]) + max(df["start"]))/2).astype(int) + flank
    df["seq"] = df.apply(lambda x: fetch_fa_from_bed_series(x, fasta), axis=1)
    df["seq"] = df["seq"].str.lower()
    return df

def check_if_entry_part_of_list(checkfor, list_of_lists):
    assert isinstance(checkfor, list)
    assert isinstance(list_of_lists, list)
    for list_entry in list_of_lists:
        for comb in itertools.combinations(list_entry, len(checkfor)):
            if tuple(checkfor) == comb:
                return True
    return False

def generate_nonoverlapping_indices(df):
    index_combs = list()
    for i in reversed(range(2, len(df.index))):
        for combs in itertools.combinations(df.index, r=i):
            df2 = df.iloc[list(combs)].copy()
            df2["previous_end"] = df2["start"].shift(1) + df2["ref_length"].shift(1) - 1
            df2["exclude"] = np.where(df2["start"] <= df2["previous_end"].fillna(0),
                                     True,
                                     False)
            if not df2.exclude.any():
                if not check_if_entry_part_of_list(list(df2.index), index_combs):
                    index_combs = index_combs + [list(tuple(df2.index))]
    return index_combs

def extract_combinations_fasta(df, flank, trim=True):
    if len(df["seq"].unique()) > 1:
        raise ValueError("Non-unique sequences in DataFrame")
    sequence = df.loc[0, "seq"]
    positions = list(df["pos"].astype(int))
    reflengths = list(df["ref_length"])
    coords = list(df["coord"])  
    insertion = 0
    if trim:
        df["alt_length"] = df.apply(lambda x: len(x["alt"]), axis=1)
        indels = np.array(df["alt_length"] - df["ref_length"])
        insertion = int(np.ceil(indels[indels>=0].sum()/2))
    entries = df.shape[0]
    allele_combinations = itertools.product(*[df[["ref", "alt"]].iloc[i].tolist() for i in range(entries)])
    sequences = ""
    for alleles in allele_combinations:
        sequence_mod = sequence
        header = ">"
        for i in range(entries):
            pos = positions[i]
            if (i > 0) and (pos < positions[(i-1)]+reflengths[(i-1)]):
                raise ValueError("Overlapping SNPs should not occur at this stage")
            pos = pos + (len(sequence_mod)-2*flank)
            sequence_mod = replace_str(sequence_mod, pos, reflengths[i], alleles[i])
            header = header + "_" + coords[i] + "_" + alleles[i]
        header = header.replace(">_", ">")
        sequence_mod = sequence_mod[insertion:len(sequence_mod)-insertion]
        sequences =  sequences + f"{header}\n{sequence_mod}\n"
    return sequences
