# Neoantigen Pipeline Helper Scripts


- Docker
```bash
  bsub -n 1 -Is -G compute/ -g /evelyn/default -q general-interactive -M 16G -R 'rusage[mem=16G]' -a 'docker('evelyns2000/neoang_scripts')' /bin/bash
```
  
## Get FDA Thresholds

- List of parameters
```
python3 /opt/scripts/get_FDA_thresholds.py --help
usage: get_FDA_thresholds.py [-h] [-WB WB] [-f FIN_RESULTS] [--n_dna N_DNA] [--t_dna T_DNA] [--t_rna T_RNA] [--una_n_dna UNA_N_DNA] [--una_t_dna UNA_T_DNA] [--una_t_rna UNA_T_RNA] [--somalier SOMALIER] [--contam_n CONTAM_N] [--contam_t CONTAM_T]

Process some integers.

optional arguments:
  -h, --help            show this help message and exit
  -WB WB                the path to the gcp_immuno folder of the trial you wish to tun script on, defined as WORKING_BASE in envs.txt
  -f FIN_RESULTS, --fin_results FIN_RESULTS
                        Name of the final results folder in gcp immuno
  --n_dna N_DNA         file path for aligned normal dna FDA report table
  --t_dna T_DNA         file path for aligned tumor dna FDA report table
  --t_rna T_RNA         file path for aligned tumor rna FDA report table
  --una_n_dna UNA_N_DNA
                        file path for unaligned normal dna FDA report table
  --una_t_dna UNA_T_DNA
                        file path for unaligned tumor dna FDA report table
  --una_t_rna UNA_T_RNA
                        file path for unaligned tumor rna FDA report table
  --somalier SOMALIER   file path for Somalier results for sample tumor/normal sample relatedness
  --contam_n CONTAM_N   file path for VerifyBamID results for contamination the normal sample
  --contam_t CONTAM_T   file path for VerifyBamID results for contamination the tumor sample
  ```
### Example Commands
```bash
python3 /opt/scripts/get_FDA_thresholds.py -WB /storage1/fs1/mgriffit/Active/JLF_MCDB/cases/JLF-100-042/gcp_immuno -f final_results_v1
```

```bash
python3 /opt/scripts/get_FDA_thresholds.py --n_dna aligned_normal_dna_table2.csv --t_dna aligned_tumor_dna_table2.csv --t_rna aligned_tumor_rna_table3.csv --una_n_dna unaligned_normal_dna_table1.csv --una_t_dna unaligned_tumor_dna_table1.csv --una_t_rna unaligned_tumor_rna_table1.csv --somalier concordance.somalier.pairs.tsv --contam_n normal.VerifyBamId.selfSM --contam_t tumor.VerifyBamId.selfSM
```

## Get Neoanitgen QC

```
python3 /opt/scripts/get_neoantigen_qc.py --help
usage: get_neoantigen_qc.py [-h] [-WB WB] [-f FIN_RESULTS] [--n_dna N_DNA] [--t_dna T_DNA] [--t_rna T_RNA] [--concordance CONCORDANCE] [--contam_n CONTAM_N] [--contam_t CONTAM_T] [--rna_metrics RNA_METRICS] [--strand_check STRAND_CHECK] --yaml YAML
                            [--fin_variants FIN_VARIANTS]

Get the stats for the basic data QC review in the neoantigen final report.

optional arguments:
  -h, --help            show this help message and exit
  -WB WB                the path to the gcp_immuno folder of the trial you wish to tun script on, defined as WORKING_BASE in envs.txt
  -f FIN_RESULTS, --fin_results FIN_RESULTS
                        Name of the final results folder in gcp immuno
  --n_dna N_DNA         file path for aligned normal dna FDA report table
  --t_dna T_DNA         file path for aligned tumor dna FDA report table
  --t_rna T_RNA         file path for aligned tumor rna FDA report table
  --concordance CONCORDANCE
                        file path for Somalier results for sample tumor/normal sample relatedness
  --contam_n CONTAM_N   file path for VerifyBamID results for contamination the normal sample
  --contam_t CONTAM_T   file path for VerifyBamID results for contamination the tumor sample
  --rna_metrics RNA_METRICS
  --strand_check STRAND_CHECK
  --yaml YAML
  --fin_variants FIN_VARIANTS

```

### Example Commands
```bash
python3 /opt/scripts/get_neoantigen_qc.py -WB /path/to/mcdb048/gcp_immuno -f final_results_v1 --yaml /path/to/gcp_immuno/final_results_v1/workflow_artifacts/mcdb048_immuno_cloud-WDL.yaml
```

```bash
python3 /opt/scripts/get_neoantigen_qc.py --n_dna normal_dna_aligned_metrics.txt --t_dna tumor_dna_aligned_metrics.txt --t_rna tumor_rna_aligned_metrics.txt --concordance concordance.somalier.pairs.tsv --contam_n normal.VerifyBamId.selfSM --contam_t tumor.VerifyBamId.selfSM --rna_metrics rna_metrics.txt --strand_check trimmed_read_1strandness_check.txt --yaml jlf-100-044_immuno_cloud-WDL.yaml --fin_variants variants.final.annotated.tsv
```

## Generate Review Files

```
python3 /opt/scripts/generate_reviews_files.py --help
usage: generate_reviews_files.py [-h] [-a A] [-c C] [-samp SAMP] [-f FIN_RESULTS]

Create the file needed for the neoantigen manuel review

options:
  -h, --help            show this help message and exit
  -a A                  The path to the ITB Reviewed Candidates
  -c C                  The path to annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv from the generate_protein_fasta script
  -samp SAMP            The name of the sample
  -f FIN_RESULTS, --fin_results FIN_RESULTS
                        Name of the final results folder in gcp immuno
```

```bash
python3 /opt/scripts/generate_reviews_files.py -a /Volumes/gillandersw/Active/Project_0001_Clinical_Trials/CTEP/analysis/TWJF-10146-MO011-0021/itb-review-files/10146-0021.Annotated.Neoantigen_Candidates.xlsx -c /Volumes/gillandersw/Active/Project_0001_Clinical_Trials/CTEP/analysis/TWJF-10146-MO011-0021/generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv -samp 10146-0021
```
## Color Peptides 51mer
```
python3 /opt/scripts/color_peptides51mer.py --help
usage: color_peptides51mer.py [-h] -p P -classI CLASSI -classII CLASSII -o O

Color the 51mer peptide

options:
  -h, --help        show this help message and exit
  -p P              The path to the Peptides 51 mer
  -classI CLASSI    The path to the classI all_epitopes.aggregated.tsv used in pVACseq
  -classII CLASSII  The path to the classII all_epitopes.aggregated.tsv used in pVACseq
  -o O              Output location
```

```bash
python3 scripts/color_peptides51mer.py -p /Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/10146-0021_Peptides_51-mer.xlsx -classI /Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/classI.TWJF-10146-0021-Tumor_Lysate.all_epitopes.aggregated.tsv -classII /Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/classII.TWJF-10146-0021-Tumor_Lysate.all_epitopes.aggregated.tsv -o /Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/test_colored_peptide.html
```
## Bold ClassII

```
python3 /opt/scripts/bold_classII.py --help
usage: bold_classII.py [-h] -p P -classI CLASSI -classII CLASSII -o O

Bold the class II pepetides

options:
  -h, --help        show this help message and exit
  -p P              The path to the Peptides 51 mer
  -classI CLASSI    The path to the classI all_epitopes.aggregated.tsv used in pVACseq
  -classII CLASSII  The path to the classII all_epitopes.aggregated.tsv used in pVACseq
  -o O              Output location for the html doc
```

```bash
python3 /opt/scripts/bold_classII.py -p /Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/10146-0021 Peptides 51-mer.xlsx -classI /Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/classI.TWJF-10146-0021-Tumor_Lysate.all_epitopes.aggregated.tsv -classII /Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/classII.TWJF-10146-0021-Tumor_Lysate.all_epitopes.aggregated.tsv -o /Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/test_colored_peptide.html
```
