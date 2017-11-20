GWAS summary statistics harmoniser
==================================

Work in progress. Todo:
- Log files
- Finalise output format
- Automatically check test data output for inconsistencies

#### Requirements

- python2

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
                             <str> --eaf_col <str> --beta_col <str>
                             [--maf_palin_threshold <float>]
                             [--af_vcf_field <str>] [--af_vcf_min <float>]
                             [--in_sep <str>] [--out_sep <str>]
                             [--infer_strand <bool>] [--infer_palin <bool>]

Summary statistc harmoniser

optional arguments:
  -h, --help            show this help message and exit
  --sumstats <file>     GWAS summary statistics file
  --vcf <file>          Reference VCF file. Use # as chromosome wildcard.
  --out <file>          Harmonised output file
  --log <file>          Log file
  --rsid_col <str>      Rsid column
  --chrom_col <str>     Chromosome column
  --pos_col <str>       Position column
  --effAl_col <str>     Effect allele column
  --otherAl_col <str>   Other allele column
  --eaf_col <str>       EAF column
  --beta_col <str>      beta column
  --maf_palin_threshold <float>
                        Max MAF that will be used to infer palindrome strand
                        (default: 0.42)
  --af_vcf_field <str>  VCF info field containing allele freq (default:
                        AF_NFE)
  --af_vcf_min <float>  Min freq of alt allele to be included (default: 0.001)
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

Get record from reference VCF (vcfrec) using ssrec chrom:pos and tabix.
  If no record -> discard (code: ??) and break
  Assert that only 1 record was returned (should be the case for gnomad vcfs at least)
  Discard alt alleles if allele freq is below af_vcf_min threshold
  If no alt alleles remaining -> discard (code: ??) and break

As vcfrec could be multialleleic, compare ssrec alleles to vcfrec alleles.
Search for ssrec allele match with each vcfrec ref-alt allele pair and their reverse complements.
Remove alt allele and corresponding INFO entries if not a match.
  If no matches -> discard ssrec as ambiguous (code: ??) and break
  If >1 match -> discard ssrec as ambiguous (code: ??) and break
  If 1 match -> proceed

Check whether ssrec alleles and vcfrec alleles are palindromic
  If palindromic
    If infer_palin==True and infer_strand==True:
      If MAF of either ssrec or vcfrec alleles is > maf_palin_infer_threshold -> discard (error code: ??)
      If ssrec and vcfrec EAFs are not concordant (one is >0.5 and other <0.5)
        Flip beta (code: ??) and break
      Else:
        Don't flip beta (code: ??) and break
    Else: -> discard (code: ??) and break

Check whether on the opposite strand: ssrec alleles == reverse_complement(vcfrec alleles)
  If opposite strand
    If infer_strand==True
      If ssrec effect_al == reverse_complement(vcfrec effect_al):
        Don't flip beta (code: ??) and break
      Else:
        Flip beta (code: ??) and break
    Else:
      Assuming forward strand therefore alleles are ambiguous -> discard (code: ??) and break

Check that alleles are on the same strand (they should be by now).
  If same stand:
    If ssrec effect_al == vcfrec effect_al:
      Don't flip beta (code: ??) and break
    Else:
      Flip beta (code: ??) and break

Else error (shouldn't ever get here):
  Raise error
```
