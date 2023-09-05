# snp2fasta
Take list of SNPs and generate fasta files with flanks

## Installation
`git clone https://github.com/efriman/snp2fasta.git`

`cd snp2fasta`

`pip install .`

Requires pysam, which in some cases may need to be installed separately using `conda install -c bioconda pysam`

## Usage

`snp2fasta input_table.txt --flank 5 --fasta reference_genome.fa --outname test`

By default, character "-" is treated as deletion. Can be changed with `--ignore_char`

By default, both ends of ref and alt are trimmed for insertions so that the length does not exceed 2*flank. Can be unset with `--no_trim`

Experimental function (still yields some missed combinations):

Combinations of SNPs within the flanking distance can be generated using `--combinations k` where `k` is the maximum number of combinations allowed. This will flank around the center of all nearby SNPs and generate all combinations in a separate file. Set maximum distance using `--maxdist`.

## Example inputs and output
input_table.bed
| chrom  | start | ref | alt |
| ------------- | ------------- | ------------- |  ------------- |
| chr1  | 1  | A | G |
| chr2  | 10  | C | T |

test_matched.fa

\>chr1_1_A

atgatAatgat

\>chr1_1_G

atgatGatgat

\>chr2_10_C

ctagcCgtacg

\>chr2_10_T

ctagcTgtacg