GWAS summary statistics harmoniser
==================================

Work in progress.

##### Pseudo code

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

##### Set up environment
```
# Install dependencies into isolated environment
conda env create -n sumstat_harmoniser --file environment.yaml

# Activate environment
source activate sumstat_harmoniser
```
