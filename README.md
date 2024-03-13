# Variant effect predictor benchmarking pipeline

This GitHub repository contains the pipeline used in "Benchmarking computational variant effect predictors by their ability to infer human traits" (*under review*).

## Install

Before installing the pipeline, make sure your R version is at least R 4.0.

Next, clone the repository.

```
git clone https://github.com/DanielTabet/VEP_benchmarking.git
```

**Please note**, this repository does NOT contain raw data from either of the UK Biobank or *All of Us* cohorts, these are available upon application. To implement this pipeline as-is, cohort sequencing data is required in pVCF format (e.g, UKB field ID 23157). Alternatively, cohort data can be processed as in [sample_variants_filtered.csv](sample_variants_filtered.csv) and [sample_phenotypes_filtered.csv](sample_phenotypes_filtered.csv) for participant genotypes and phenotypes, respectively. We provide a sample set of scores from all 24 computational variant effect predictors assessed in this analysis in [common/VEP_scores_vFinal_sample.csv](common/VEP_scores_vFinal_sample.csv).

## Confiugre the running parameters

The [config.yaml](config.yaml) file lists all the running parameters that you can congifure.

Here we document all parameters supported by the pipeline and their functions.

| Parameter | Description | Default Value |
| --- | --- | --- |
| gene_list| A file listing all genes included in the study.<br>*See the sample file for format.* | [common/genes.csv](common/genes.csv) |
| phenotype_list | A file with all phenotypes (traits) included in the study. <br>*See the sample file for format.* | [common/phenotypeDescriptions.csv](common/phenotypeDescriptions.csv) |
| rscript_path | The path to the RScript executable. | Rscript |
| input_var_dir | The path to the input variants. | input |
| occurance_cutoff | The occurance cutoff used in the study.<br>*See the Methods section of the manuscript for detail.* | 10 |
| bootstrap_iterations | The number of iterations for the bootstrap resampling process.<br>*See the Methods section of the manuscript for detail.* | 1000 |
| all_variants_path | The filename pattern of the input variant files. | ukb23148_c%s_b%s_v1_filtered_mut.csv |
| unique_variants_path | The filename pattern of the input unique variants and computational predictor scores. | ukb23148_c%s_b%s_v1_all_weights.csv |
| variant_blocks_path | The filename of the pVCF file blocks. UK Biobank exome data is split across many pVCF files. This document itemises the content of these files. Read more about this here: https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=837 | [common/pvcf_blocks.txt](common/pvcf_blocks.txt) |
| withdraw_eids_path | A file with all the participants who have withdrawn their participation from the UK Biobank study.<br>*See the sample file for format. Replace "<withdrawn_eidx>" in the sample file with actual EIDs. We are unable to provide real EIDs due to UK Biobank's data sharing restrictions.* | [common/withdraws.csv](common/withdraws.csv) |
| variant_predictors | A list of variant effect predictors used in the study.<br>*Format: [predictor name]: [predictor column name in [sample_variants_filtered.csv](sample_variants_filtered.csv)]* | <br>VARITY: VARITY_R<br>AlphaMissense: alphaMissense_score<br>... |
| plot_individual_correlations | Whetherindividually plot correlations as scatterplots | FALSE |
| plot_individual_correlations | Whetherindividually plot correlations as scatterplots | FALSE |

## Run the pipeline

Execute.

```
Rscript main.R <configuration-file> <log-dir>
```

There are two required arguments:
1. **configuration-file:** the path to configuration file (e.g. config.yaml)
2. **log-dir:** the directory to the log files (e.g. logs/)

## Output

Outputs should be stored in the [output](output) folder. All pairwise predictor comparisons are included in pvals.csv in the root directory.

## Cohort data processing

Scripts used to process UK Biobank and All of Us cohort data are included in [data_processing](data_processing).

## Contact us

If you have any feedback, suggestions or questions, please reach out via [email](mailto:daniel.tabet@mail.utoronto.ca) or open an issue on github.
