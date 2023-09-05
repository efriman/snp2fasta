# snp2fasta
Take list of SNPs and generate fasta files with flanks

## Installation
`git clone https://github.com/efriman/snp2fasta.git`

`cd snp2fasta`

`pip install .`

Requires pysam, which in some cases may need to be installed separately using `conda install -c bioconda pysam`

## Usage
`snp2fasta input_table.txt --flank [INTEGER] --fasta reference_genome.fa --outname output [OPTIONS]`

By default, character "-" is treated as deletion. Can be changed with `--ignore_char`

By default, both ends of ref and alt are trimmed for insertions so that the length does not exceed 2*flank. Can be changed by specifying `--no_trim`

Combinations of SNPs can be generated using `--combinations k` where `k` is the maximum number of combinations allowed (has to be at least 2). This will flank around the center of the SNPs within the maximum distance and generate up to `k` combinations of alleles in a separate file. Set maximum distance between SNPs using `--maxdist`.

## Example inputs and output
`snp2fasta input_table.txt --flank 5 --fasta reference_genome.fa --outname test --combinations 2 --maxdist 10`

input_table.bed
| chrom  | start | ref | alt |
| ------------- | ------------- | ------------- |  ------------- |
| chr1  | 6  | A | G |
| chr2  | 11  | C | T |

test_matched.fa

\>chr1_6_A

atgatAtagcc

\>chr1_6_G

atgatGtagcc

\>chr2_11_C

atagcCgtacg

\>chr2_11_T

atagcTgtacg


test_combinations.fa

\>chr1_6_A_chr1_11_C

atAtagcCgt

\>chr1_6_G_chr1_11_C

atGtagcCgt

\>chr1_6_A_chr1_11_T

atAtagcTgt

\>chr1_6_G_chr1_11_T

atGtagcTgt
