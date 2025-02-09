[![Runs successfully](https://github.com/houlstonlab/select-cohort-variants/actions/workflows/runs-successfully.yml/badge.svg)](https://github.com/houlstonlab/select-cohort-variants/actions/workflows/runs-successfully.yml)

### Introduction

This workflow filters and counts variants meeting certain criteria in an 
annotated VCF. Qualifying variants are then tabulated by gene and information 
about thier frequency in the dataset is exported.

### Usage

The typical command looks like the following. `--cohorts` are required inputs. 
Different versions of the workflow can be called using `-r` and output directed 
to `--output_dir`

```bash
nextflow run houlstonlab/select-cohort-variants \
    -r main \
    --output_dir results/ \
    --cohorts input/cohorts_input.csv
```

### Inputs & Parameters

Other paramters include:
- `genome`    : genome version (default is'hg38') 
- `style`     : chromosome names style 'UCSC' or 'NCBI'
- `coding`    : restrict to coding regions. Default `'true'`
- `categories`: selection categories. One or more of `'Pathogenic,Damaging,Splicing,High,PTV,Stop'`
- `GQ`        : minimum genotyping quality. Default `> 10`
- `DP`        : minimum read depth. Default `> 5`
- `VAF`       : variant allele frequency. Default `> 0.2`
- `MAF`       : maximum allele frequence. Default `> 0.5` 
- `PAF`       : maximum pathogenic variants group allele frequence. Default `> 0.5`
- `HWE`       : cutoff for HWE test. Default `> 1e-5`
- `ExcHet`    : cutoff for excess heterosygosity. Default `> 0.5`
- `gnomADe_AF`: maximum allele frequency. Default `> 0.5`
- `MAX_AF`    : maximum group allele frequence. Default `> 0.5`
- `DS`        : minimum delta score for spliceAI
- `CADD`      : minimum CADD score. Defalt `> 5`

### Output

- `coordiantes/`: coding gene coordinates
- `pheno/`    : split vcf file by chromosome
- `filtered/`   : filtered vcfs
- `combined/`   : a combined vcf
- `variants/`   : extracted variant frequence
- `aggregate/`  : aggregated variant frequency by gene
- `reports/`    : summary reports
