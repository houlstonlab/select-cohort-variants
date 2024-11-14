[![Runs successfully](https://github.com/houlstonlab/select-cohort-variants/actions/workflows/runs-successfully.yml/badge.svg)](https://github.com/houlstonlab/select-cohort-variants/actions/workflows/runs-successfully.yml)

### Introduction

This workflow filters and counts variants meeting certain criteria in an annotated VCF. Qualifying
variants are then tabulated by gene and information about thier frequency in the dataset is 
exported.

### Usage

The typical command looks like the following. `--vcf` and `--cases` are required inputs. 
Different versions of the workflow can be called using `-r` and output directed to `--output_dir`

```bash
nextflow run houlstonlab/test-gene-burden \
    -r main \
    --output_dir results/ \
    --vcf input/*.variants.vcf.gz{,.tbi} \
    --cases input/*.cases.txt
```

### Inputs & Parameters

- `vcf`  : a VCF file with annotated variants and genotypes
- `cases`: a list of the cases to consider
- `genome`: genome version (default is 'hg38')
- `style` : chromosome names style 'UCSC' or 'NCBI'
- Genotype filtering parameters
  - `GQ`: genotype quality (default is 30)
  - `DP`: call depth (default is 10)
  - `VAF`: ratio of reads counts between alt and ref (default is 0.4)
  - `MAF`: minor allele frequency (default is 0.005)
  - `HWE`: hardy weinberg equilibrium test p-value cut off (default is 1e-5)
  - `ExcHet`: excess heterozygousity test p-value cut off (default is 0.9)
- Variant filtering parameters
  - `gnomADe_AF`: allele frequency GnomAD exomes cut off (default is 0.005)
  - `MAX_AF`: max population allele frequency GnomAD exomes cut off (default is 0.005)
  - `DS`: SpliceAI predicted scores cut off (default is 0.8)
    
### Output

- `filtered/`  : filtered variants
- `variants/` : the frequency, annotation and genotypes of qualifying variants
- `aggregate/`: aggregated variant info by gene
- The final ouput in `summary/` is a tsv file with `gene`, `het`, `hom`, and `ch`
