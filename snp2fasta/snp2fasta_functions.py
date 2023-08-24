import numpy as np

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

def process_snp_fasta_table(data):
    data["seq"] = data["seq"].str.lower()
    data['ref_pos'] = data.apply(lambda x: slice_str(x["seq"], 
                                                     x['snp_pos'], 
                                                     x['snp_pos']+x["ref_length"]), 
                                 axis=1)
    data["alt"] = np.where(data["ref"] == "-",
                             data['ref_pos'].str.upper() + data['alt'],
                             data['alt'])
    data["ref"] = np.where(data["ref"] == "-",
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
