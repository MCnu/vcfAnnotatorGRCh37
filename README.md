## Purpose of this repo
This repo lets users annotate the variants in a VCF file with GRCh37 hgvs info

Primarily, each variant should be annotated the following information:

- Depth of sequence coverage at the site of variation
- Number of reads supporting the variant
- Percentage of reads supporting the variant versus those supporting reference reads
- The gene of the variant (if applicable)
- Type of variation (currently supports: substitution, deletion, insertion, deletion-insertion)
- Effect of variation (missense, silent, intergenic, etc.)
- The frequency of each alt allele (by region, as provided by hgvs, if available)

## How this code works

The code can be broken into the following steps:

- This python script takes a VCF file with reads aligned to GRCh37
- Variant call position is translated into hgvs notation
- Each variant allele is queried using the VEP hgvs API for GRCh37
- Variants are annotated with info from pyvcf and hgvs query
- Each variant and its annotations are output as a row in a tsv file

## How to use the python script

- Create a conda environment using the vcfAnnotator.yml in the condaEnv folder
- Activate the environment
- Navigate to the directory containing vcfAnnotatorGRCh37.py (scripts)
- Execute the following within your CLI, substituting in your own file paths:
```
python vcfAnnotatorGRCh37.py -i {inputVCFPath} -o {outputTSVPath}
```

## Contents of output tsv

Each variant will have a row with the following columns:

| Column Name | Description |
|-:|:-|
| CHROM | Chromosome |
| GRCH37_POS | Position on chromosome from GRCh37 gneome assembly |
| REF | Reference Allele |
| ALT | List of Alternative Alleles |
| TOTCOV | Total read coverage of position |
| VARREADS | List of read count supporting each alt |
| PERCVARREADS | Percentage of total reads supporing each alt |
| VCF_AAF | ALT allele frequency in vcf samples (from pyvcf) |
| VCF_FILTER | List of potential fitlers in vcf |
| VCF_VARTYPE | Variant type (from pyvcf) |
| VCF_VARSUBTYPE | Variant subtype (from pyvcf) |
| HGVS_GENEID | ID of gene containing variant (from hgvs) |
| HGVS_GENESYMBOL | Name or symbol of gene contianing variant (from hgvs) |
| HGVS_VARTYPE | Type of variant query used for hgvs notation |
| HGVS_EFFECT | Most Severe Consequence reported (from hgvs) |
| HGVS_AAF | Alternative allele frequency (from hgvs) |
| HGVS_COLOCVAR | rsID of other known variants at same position (from hgvs) |
