GWAS summary statistics harmoniser
==================================

Work in progress. Todo:
- Finalise output format
- Add option to assume forward strand
- Add option to use odds ratio (OR), rather than betas
- Find way to speed up reference VCF query (currently using tabix which takes up 85% of run time). Possibilities:
  - [Giggle](https://github.com/ryanlayer/giggle) is reported to be faster than tabix, but I don't know if this is only for multiple querys.
  - Load required lines from reference VCF into memory

#### Requirements

- python3
- [HTSlib](http://www.htslib.org/download/) (for tabix)

```
# Install dependencies into isolated environment
conda env create -n sumstat_harmoniser --file environment.yaml

# Activate environment
source activate sumstat_harmoniser
```

#### Usage

Command line arguments can be viewed with `python sumstat_harmoniser.py --help`

```
usage: sumstat_harmoniser.py [-h] --sumstats <file> --vcf <file> --out <file>
                             --log <file> --rsid_col <str> --chrom_col <str>
                             --pos_col <str> --effAl_col <str> --otherAl_col
                             <str> --beta_col <str> [--eaf_col <str>]
                             [--maf_palin_threshold <float>]
                             [--af_vcf_min <float>] [--af_vcf_field <str>]
                             [--in_sep <str>] [--out_sep <str>]
                             [--infer_strand <bool>] [--infer_palin <bool>]

Summary statistc harmoniser

optional arguments:
  -h, --help            show this help message and exit
  --sumstats <file>     GWAS summary statistics file
  --vcf <file>          Reference VCF file. Use # as chromosome wildcard.
  --out <file>          Harmonised output file. Use 'gz' extension to gzip.
  --log <file>          Log file. Use 'gz' extension to gzip.
  --rsid_col <str>      Rsid column
  --chrom_col <str>     Chromosome column
  --pos_col <str>       Position column
  --effAl_col <str>     Effect allele column
  --otherAl_col <str>   Other allele column
  --beta_col <str>      beta column
  --eaf_col <str>       EAF column
  --only_chrom <str>    Only process provided chromosome.
  --maf_palin_threshold <float>
                        Max MAF that will be used to infer palindrome strand
                        (default: 0.42)
  --af_vcf_min <float>  Min freq of alt allele to be used (default: 0.001)
  --af_vcf_field <str>  VCF info field containing allele freq (default:
                        AF_NFE)
  --in_sep <str>        Input file column separator (default: tab)
  --out_sep <str>       Output file column separator (default: tab)
  --infer_strand <bool>
                        Infer and flip reverse strand alleles? [True|False]
                        (default: True)
  --infer_palin <bool>  Infer and flip palindromic alleles? [True|False]
                        (default: True)
```

#### Pseudo code

```
# Pseudo code:

Load summary stat record (ssrec).

Get records from reference VCF (vcfrec) using ssrec chrom:pos and tabix. (There may be multiple lines in the VCF).
  Discard ssrec if there are no matches in VCF -> discard and break
  Remove alt alleles if allele freq is below af_vcf_min threshold
  Remove alt alleles that don't match ssrec alleles on either strand
  Discard ssrec if there are no records after filtering -> discard and break
  Discard ssrec if there are multiple vcfrecs (>1) after filtering -> discard and break
  Discard ssrec if there are multiple alt alleles (>1) after filtering -> discard and break
  There should now be 1 matching vcf record with 1 matching allele

Check whether ssrec alleles and vcfrec alleles are palindromic
  If palindromic
    If infer_palin==True and infer_strand==True:
      If MAF of either ssrec or vcfrec alleles is > maf_palin_infer_threshold -> discard
      If ssrec and vcfrec EAFs are not concordant (one is >0.5 and other <0.5)
        Flip beta and break
      Else:
        Don't flip beta and break
    Else: -> discard and break

Check whether on the opposite strand: ssrec alleles == reverse_complement(vcfrec alleles)
  If opposite strand
    If infer_strand==True
      If ssrec effect_al == reverse_complement(vcfrec effect_al):
        Don't flip beta and break
      Else:
        Flip beta and break
    Else:
      Assuming forward strand therefore alleles are ambiguous -> discard and break

Check that alleles are on the same strand (they should be by now).
  If same stand:
    If ssrec effect_al == vcfrec effect_al:
      Don't flip beta and break
    Else:
      Flip beta and break

Else error (shouldn't ever get here):
  Raise error
```
