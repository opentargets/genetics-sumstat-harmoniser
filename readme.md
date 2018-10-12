GWAS summary statistics harmoniser
==================================

Scripts for harmonising GWAS summary stats against a reference VCF sitelist. See [flowchart](flowchart_v3.svg) for implementation details.

#### Requirements

- python3
- [HTSlib](http://www.htslib.org/download/) (for tabix)

```
# Download
git clone https://github.com/opentargets/sumstat_harmoniser.git
cd sumstat_harmoniser

# Install dependencies into isolated environment
conda env create -n sumstat_harmoniser --file environment.yaml

# Activate environment
source activate sumstat_harmoniser

# Run tests
cd tests
./run_tests.py
```

#### Usage

```
$ bin/sumstat_harmoniser --help
usage: main.py [-h] --sumstats <file> --vcf <file> [--hm_sumstats <file>]
               [--hm_statfile <file>] [--strand_counts <file>]
               [--palin_mode [infer|forward|reverse|drop]]
               [--af_vcf_field <str>] [--infer_maf_threshold <float>]
               --chrom_col <str> --pos_col <str> --effAl_col <str>
               --otherAl_col <str> [--beta_col <str>] [--or_col <str>]
               [--or_col_lower <str>] [--or_col_upper <str>] [--eaf_col <str>]
               [--rsid_col <str>] [--only_chrom <str>] [--in_sep <str>]
               [--out_sep <str>] [--na_rep_in <str>] [--na_rep_out <str>]
               [--chrom_map <str> [<str> ...]]

Summary statistc harmoniser

optional arguments:
  -h, --help            show this help message and exit

Input files:
  --sumstats <file>     GWAS summary statistics file
  --vcf <file>          Reference VCF file. Use # as chromosome wildcard.

Output files:
  --hm_sumstats <file>  Harmonised sumstat output file (use 'gz' extension to
                        gzip)
  --hm_statfile <file>  Statistics from harmonisation process output file.
                        Should only be used in conjunction with --hm_sumstats.
  --strand_counts <file>
                        Output file showing number of variants that are
                        forward/reverse/palindromic

Harmonisation mode:
  --palin_mode [infer|forward|reverse|drop]
                        Mode to use for palindromic variants: (a) infer strand
                        from effect-allele freq, (b) assume forward strand,
                        (c) assume reverse strand, (d) drop palindromic
                        variants

Strand inference options:
  Options that are specific to strand inference (--palin_mode infer)

  --af_vcf_field <str>  VCF info field containing alt allele freq (default:
                        AF_NFE)
  --infer_maf_threshold <float>
                        Max MAF that will be used to infer palindrome strand
                        (default: 0.42)

Input column names:
  --chrom_col <str>     Chromosome column
  --pos_col <str>       Position column
  --effAl_col <str>     Effect allele column
  --otherAl_col <str>   Other allele column
  --beta_col <str>      beta column
  --or_col <str>        Odds ratio column
  --or_col_lower <str>  Odds ratio lower CI column
  --or_col_upper <str>  Odds ratio upper CI column
  --eaf_col <str>       Effect allele frequency column
  --rsid_col <str>      rsID column in the summary stat file

Other args:
  --only_chrom <str>    Only process this chromosome
  --in_sep <str>        Input file column separator (default: tab)
  --out_sep <str>       Output file column separator (default: tab)
  --na_rep_in <str>     How NA are represented in the input file (default: "")
  --na_rep_out <str>    How to represent NA values in output (default: "")
  --chrom_map <str> [<str> ...]
                        Map summary stat chromosome names, e.g. `--chrom_map
                        23=X 24=Y`
```

#### Examples

```
# Count number of variants that are forward strand/reverse strand/palindromic,
# without harmonisation

bin/sumstat_harmoniser \
  --sumstats {input sumstats} \
  --vcf {reference vcf} \
  --chrom_col {chrom name column} \
  --pos_col {bp position column} \
  --effAl_col {effect allele column} \
  --otherAl_col {non effect allele column} \
  --stand_counts {output count file}

# Harmonise variants whilst inferring strand of palindromic variants (requires
# known effect allele and VCF alt allele frequencies)

bin/sumstat_harmoniser \
  --sumstats {input sumstats} \
  --vcf {reference vcf} \
  --chrom_col {chrom name column} \
  --pos_col {bp position column} \
  --effAl_col {effect allele column} \
  --otherAl_col {non effect allele column} \
  --beta_col {beta column} \
  --palin_mode infer \
  --eaf_col {effect allele freq column} \
  --af_vcf_field {alt allele freq VCF info field} \
  --infer_maf_threshold {freq threshold for inference} \
  --hm_sumstats {output harmonised sumstats} \
  --hm_statfile {output outcome code stats file}

# Harmonise variants whilst assuming palindromic variants are on forward strand

bin/sumstat_harmoniser \
  --sumstats {input sumstats} \
  --vcf {reference vcf} \
  --chrom_col {chrom name column} \
  --pos_col {bp position column} \
  --effAl_col {effect allele column} \
  --otherAl_col {non effect allele column} \
  --beta_col {beta column} \
  --palin_mode forward \
  --hm_sumstats {output harmonised sumstats} \
  --hm_statfile {output outcome code stats file}
```

#### Todo

- ( ) Check that tbi file exists for the VCF
- ( ) Find way to speed up reference VCF query (currently using tabix which takes up 85% of run time). Possibilities:
  - ( ) [Giggle](https://github.com/ryanlayer/giggle) is reported to be faster than tabix, but I don't know if this is only for multiple querys.
  - ( ) Load required lines from reference VCF into memory
