# Neoantigen Pipeline Helper Scripts
These scripts assist in setting up files ofr manualling reviewing the results for Neoantigen Vaccine Desing results generate from the [Washington University Immuno Pipeline](https://github.com/wustl-oncology/analysis-wdls).
  
## Creating Case Final Report on compute 1

### Before Immunogenomics Tumor Board Review

A written case final report will be created which includes a Genomics Review Report document. This document includes a section of a basic data QC review and a table summarizing values that pass/fail the FDA quality thresholds.

### Basic data QC

Pull the basic data qc from various files. This script will output a file final_results/qc_file.txt and also print the summary to to screen.

```
mkdir $WORKING_BASE/../manual_review
cd $WORKING_BASE/../manual_review

bsub -Is -q oncology-interactive -G $GROUP -a "docker(griffithlab/neoang_scripts)" /bin/bash
python3 /opt/scripts/get_neoantigen_qc.py -WB $WORKING_BASE -f final_results --yaml $WORKING_BASE/yamls/$CLOUD_YAML
```

### FDA Quality Thresholds

This script will output a file final_results/fda_quality_thresholds_report.tsv and also print the summary to to screen.

```
python3 /opt/scripts/get_FDA_thresholds.py -WB  $WORKING_BASE -f final_results
exit
```

### After Immunogenomics Tumor Board Review

After the Immunogenomics Tumor Board Review, both a .tsv and .xlsx file are downloaded from pVACview whihc contains the canidates marked as Accept, Review, Reject, and Pending. These files should be kept in a fould named itb-review-files.

#### Generate Protein Fasta

```bash
cd $WORKING_BASE
mkdir ../generate_protein_fasta
cd ../generate_protein_fasta
mkdir candidates
mkdir all

zcat $WORKING_BASE/final_results/annotated.expression.vcf.gz | less # Get sample ID Found in the #CHROM header of VCF
export SAMPLE_ID="TWJF-10146-0029-0029_Tumor_Lysate"



bsub -Is -q general-interactive -G $GROUP -a "docker(griffithlab/pvactools:4.0.1)" /bin/bash

pvacseq generate_protein_fasta \
  -p $WORKING_BASE/final_results/pVACseq/phase_vcf/phased.vcf.gz \
  --pass-only --mutant-only -d 150 \
  -s $SAMPLE_ID \
  --aggregate-report-evaluation {Accept,Review} \
  --input-tsv ../itb-review-files/*.tsv  \
  $WORKING_BASE/final_results/annotated.expression.vcf.gz \
  25 \
  $WORKING_BASE/../generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa

pvacseq generate_protein_fasta \
  -p $WORKING_BASE/final_results/pVACseq/phase_vcf/phased.vcf.gz \
  --pass-only --mutant-only -d 150 \
  -s $SAMPLE_ID  \
  $WORKING_BASE/final_results/annotated.expression.vcf.gz \
  25  \
  $WORKING_BASE/../generate_protein_fasta/all/annotated_filtered.vcf-pass-51mer.fa

exit 
```

To generate files needed for manual review, save the pVAC results from the Immunogenomics Tumor Board Review meeting as $SAMPLE.revd.Annotated.Neoantigen_Candidates.xlsx (Note: if the file is not saved under this exact name the below command will need to be modified).

```
cd $WORKING_BASE/../manual_review
bsub -Is -q oncology-interactive -G $GROUP -a "docker(griffithlab/neoang_scripts)" /bin/bash

export SAMPLE="TWJF-10146-0029"

python3 /opt/scripts/setup_review.py -WB $WORKING_BASE -a ../itb-review-files/*.xlsx -c $WORKING_BASE/../generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv -samp $SAMPLE  -classI $WORKING_BASE/final_results/pVACseq/mhc_i/*.all_epitopes.aggregated.tsv -classII $WORKING_BASE/final_results/pVACseq/mhc_ii/*.all_epitopes.aggregated.tsv 
```

## Creating Case Final Report locally

### Before Immunogenomics Tumor Board Review

A written case final report will be created which includes a Genomics Review Report document. This document includes a section of a basic data QC review and a table summarizing values that pass/fail the FDA quality thresholds.

### Basic data QC

Pull the basic data qc from various files. This script will output a file final_results/qc_file.txt and also print the summary to to screen.

```
mkdir $WORKING_BASE/../manual_review
cd $WORKING_BASE/../manual_review

docker pull griffithlab/neoang_scripts:latest
docker run -it -v $HOME/:$HOME/ -v $HOME/.config/gcloud:/root/.config/gcloud --env $WORKING_BASE griffithlab/neoang_scripts:latest /bin/bash

cd $WORKING_BASE

python3 /opt/scripts/get_neoantigen_qc.py -WB $WORKING_BASE -f final_results --yaml $WORKING_BASE/yamls/$CLOUD_YAML
```

### FDA Quality Thresholds

This script will output a file final_results/fda_quality_thresholds_report.tsv and also print the summary to to screen.

```
python3 /opt/scripts/get_FDA_thresholds.py -WB  $WORKING_BASE -f final_results
exit
```

### After Immunogenomics Tumor Board Review

After the Immunogenomics Tumor Board Review, both a .tsv and .xlsx file are downloaded from pVACview whihc contains the canidates marked as Accept, Review, Reject, and Pending. These files should be kept in a fould named itb-review-files.

#### Generate Protein Fasta

```bash
cd $WORKING_BASE
mkdir ../generate_protein_fasta
cd ../generate_protein_fasta
mkdir candidates
mkdir all

zcat $WORKING_BASE/final_results/annotated.expression.vcf.gz | less # Get sample ID Found in the #CHROM header of VCF
export SAMPLE_ID="TWJF-10146-0029-0029_Tumor_Lysate"

docker pull griffithlab/pvactools:4.0.5
docker run -it -v $HOME/:$HOME/ --env $WORKING_BASE  --env SAMPLE_ID griffithlab/pvactools:4.0.5 /bin/bash

cd $WORKING_BASE

pvacseq generate_protein_fasta \
  -p $WORKING_BASE/final_results/pVACseq/phase_vcf/phased.vcf.gz \
  --pass-only --mutant-only -d 150 \
  -s $SAMPLE_ID \
  --aggregate-report-evaluation {Accept,Review} \
  --input-tsv ../itb-review-files/*.tsv  \
  $WORKING_BASE/final_results/annotated.expression.vcf.gz \
  25 \
  $WORKING_BASE/../generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa

pvacseq generate_protein_fasta \
  -p $WORKING_BASE/final_results/pVACseq/phase_vcf/phased.vcf.gz \
  --pass-only --mutant-only -d 150 \
  -s $SAMPLE_ID  \
  $WORKING_BASE/final_results/annotated.expression.vcf.gz \
  25  \
  $WORKING_BASE/../generate_protein_fasta/all/annotated_filtered.vcf-pass-51mer.fa

exit 
```

To generate files needed for manual review, save the pVAC results from the Immunogenomics Tumor Board Review meeting as $SAMPLE.revd.Annotated.Neoantigen_Candidates.xlsx (Note: if the file is not saved under this exact name the below command will need to be modified).

```
cd $WORKING_BASE/../manual_review

docker pull griffithlab/neoang_scripts:latest
docker run -it -v $HOME/:$HOME/ -v $HOME/.config/gcloud:/root/.config/gcloud --env $WORKING_BASE griffithlab/neoang_scripts:latest /bin/bash

export SAMPLE="TWJF-10146-0029"

python3 /opt/scripts/setup_review.py -WB $WORKING_BASE -a ../itb-review-files/*.xlsx -c $WORKING_BASE/../generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv -samp $SAMPLE  -classI $WORKING_BASE/final_results/pVACseq/mhc_i/*.all_epitopes.aggregated.tsv -classII $WORKING_BASE/final_results/pVACseq/mhc_ii/*.all_epitopes.aggregated.tsv 
```
Open colored_peptides51mer.html and copy the table into an excel spreadsheet. The formatting should remain. Utilizing the Annotated.Neoantigen_Candidates and colored Peptides_51-mer for manual review.

# Description of Scripts includes

## Get Basic QC

```bash
python3  /opt/scripts/get_neoantigen_qc.py --help
usage: get_neoantigen_qc.py [-h] [-WB WB] [-f FIN_RESULTS] [--n_dna N_DNA] [--t_dna T_DNA] [--t_rna T_RNA] [--concordance CONCORDANCE] [--contam_n CONTAM_N] [--contam_t CONTAM_T]
                            [--rna_metrics RNA_METRICS] [--strand_check STRAND_CHECK] --yaml YAML [--fin_variants FIN_VARIANTS]

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

## GET FDA metrics
python3  /opt/scripts/get_FDA_thresholds.py --help
usage: get_FDA_thresholds.py [-h] [-WB WB] [-f FIN_RESULTS] [--n_dna N_DNA] [--t_dna T_DNA] [--t_rna T_RNA] [--una_n_dna UNA_N_DNA] [--una_t_dna UNA_T_DNA] [--una_t_rna UNA_T_RNA]
                             [--somalier SOMALIER] [--contam_n CONTAM_N] [--contam_t CONTAM_T]

Get FDA qc stats from various files and determine if they pass or fail.

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
  --somalier SOMALIER   file path for Somalier results for sample tumor/normal sample relatedness (concordance.somalier.pairs.tsv)
  --contam_n CONTAM_N   file path for VerifyBamID results for contamination the normal sample
  --contam_t CONTAM_T   file path for VerifyBamID results for contamination the tumor dna sample

## Setup Review

The set up review script runs two other scripts: generate_reviews_files.py and color_peptides51mer.py. The first sets up the Annotated.Neoantige.Canidates spreadsheet and the Peptides 51mer spreadsheet. The second script colors the Peptides 51mer sequences and outputs an html table whihc can be copied into a Microsoft spreadsheet.

```bash
python3  /opt/scripts/setup_review.py --help
usage: setup_review.py [-h] [-WB WB] [-samp SAMP] [-a A] [-c C] -classI CLASSI -classII CLASSII

Sets up manuel review files

optional arguments:
  -h, --help        show this help message and exit
  -WB WB            the path to the gcp_immuno folder of the trial you wish to tun script on, defined as WORKING_BASE in envs.txt
  -samp SAMP        Name of the sample
  -a A              Path to ITB Reviewed Candidates
  -c C              Path to annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv
  -classI CLASSI    Path to classI all_epitopes.aggregated.tsv
  -classII CLASSII  Path to classII all_epitopes.aggregated.tsv

```

## Generate Review Files

```bash
python3  /opt/scripts/generate_reviews_files.py --help
usage: generate_reviews_files.py [-h] [-a A] [-c C] [-samp SAMP] [-WB WB] [-f FIN_RESULTS]

Create the file needed for the neoantigen manuel review

optional arguments:
  -h, --help            show this help message and exit
  -a A                  The path to the ITB Reviewed Candidates
  -c C                  The path to annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv from the generate_protein_fasta script
  -samp SAMP            The name of the sample
  -WB WB                the path to the gcp_immuno folder of the trial you wish to tun script on, defined as WORKING_BASE in envs.txt
  -f FIN_RESULTS, --fin_results FIN_RESULTS
                        Name of the final results folder in gcp immuno
```

## Color Peptides 51mer

```
python3  /opt/scripts/generate_reviews_files.py --help
usage: generate_reviews_files.py [-h] [-a A] [-c C] [-samp SAMP] [-WB WB] [-f FIN_RESULTS]

Create the file needed for the neoantigen manuel review

optional arguments:
  -h, --help            show this help message and exit
  -a A                  The path to the ITB Reviewed Candidates
  -c C                  The path to annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv from the generate_protein_fasta script
  -samp SAMP            The name of the sample
  -WB WB                the path to the gcp_immuno folder of the trial you wish to tun script on, defined as WORKING_BASE in envs.txt
  -f FIN_RESULTS, --fin_results FIN_RESULTS
                        Name of the final results folder in gcp immuno
root@02d92932f0f8:/# python3  /opt/scripts/color_peptides51mer.py --help
usage: color_peptides51mer.py [-h] -p P -classI CLASSI -classII CLASSII [-WB WB] [-samp SAMP]

Color the 51mer peptide

optional arguments:
  -h, --help        show this help message and exit
  -p P              The path to the Peptides 51 mer
  -classI CLASSI    The path to the classI all_epitopes.aggregated.tsv used in pVACseq
  -classII CLASSII  The path to the classII all_epitopes.aggregated.tsv used in pVACseq
  -WB WB            the path to the gcp_immuno folder of the trial you wish to tun script on, defined as WORKING_BASE in envs.txt
  -samp SAMP        Name of the sample
```

## Bold Class II

Bold Class II is not utilized in the current workflow of setting up the manual review. However, it is included as an example of how adding stylization (in this case bold) to certain characters within individual cells of spreadsheets can be accomplished using BeautifulSoup. If you wanted to only bold certain characters (or do any stylzing such as coloring), you can insert a style tag directly into the HTML. However, if you wanted to do formatting inside formatting such as in these review files where there needs to some characters which are red, bold, underlined, or any combination if thise mentioned, BeuaitfulSoup cannot accomplish this. 

```bash
python3 /opt/scripts/bold_classII.py --help
usage: bold_classII.py [-h] -p P -classI CLASSI -classII CLASSII -o O

Bold the class II pepetides

optional arguments:
  -h, --help        show this help message and exit
  -p P              The path to the Peptides 51 mer
  -classI CLASSI    The path to the classI all_epitopes.aggregated.tsv used in pVACseq
  -classII CLASSII  The path to the classII all_epitopes.aggregated.tsv used in pVACseq
  -o O              Output location

```
