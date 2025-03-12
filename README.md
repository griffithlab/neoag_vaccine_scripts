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

bsub -Is -q oncology-interactive -G $GROUP -a "docker(griffithlab/neoang_scripts:version7)" /bin/bash
python3 /opt/scripts/get_neoantigen_qc.py -WB $WORKING_BASE -f final_results --yaml $WORKING_BASE/yamls/$CLOUD_YAML
```

### FDA Quality Thresholds

This script will output a file final_results/fda_quality_thresholds_report.tsv and also print the summary to to screen.

```
python3 /opt/scripts/get_FDA_thresholds.py -WB  $WORKING_BASE -f final_results
```

### HLA Comparison
This script will output a file manual_review/hla_comparison.tsv and also print the summary to to screen.

```
python3 /opt/scripts/hla_comparison.py -WB $WORKING_BASE
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

#generate a protein fasta file using the final annotated/evaluated neoantigen candidates TSV as input
#this will filter down to only those candidates under consideration and use the top transcript

# check the file to find Tumor sample ID in the #CHROM header of VCF

zcat $WORKING_BASE/final_results/annotated.expression.vcf.gz | less
export TUMOR_ID="100-049-BG004667"

bsub -Is -q general-interactive -G $GROUP -a "docker(griffithlab/pvactools:4.0.1)" /bin/bash

pvacseq generate_protein_fasta \
  -p $WORKING_BASE/final_results/pVACseq/phase_vcf/phased.vcf.gz \
  --pass-only --mutant-only -d 150 \
  -s $TUMOR_ID \
  --aggregate-report-evaluation {Accept,Review} \
  --input-tsv ../itb-review-files/*.tsv  \
  $WORKING_BASE/final_results/annotated.expression.vcf.gz \
  25 \
  $WORKING_BASE/../generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa

pvacseq generate_protein_fasta \
  -p $WORKING_BASE/final_results/pVACseq/phase_vcf/phased.vcf.gz \
  --pass-only --mutant-only -d 150 \
  -s $TUMOR_ID  \
  $WORKING_BASE/final_results/annotated.expression.vcf.gz \
  25  \
  $WORKING_BASE/../generate_protein_fasta/all/annotated_filtered.vcf-pass-51mer.fa

exit 
```

To generate files needed for manual review, save the pVAC results from the Immunogenomics Tumor Board Review meeting as $SAMPLE.revd.Annotated.Neoantigen_Candidates.xlsx (Note: if the file is not saved under this exact name the below command will need to be modified).

```
export PATIENT_ID=TWJF-5120-28

bsub -Is -q oncology-interactive -G $GROUP -a "docker(griffithlab/neoang_scripts:version7)" /bin/bash

cd $WORKING_BASE
mkdir ../manual_review

python3 /opt/scripts/generate_reviews_files.py -a ../itb-review-files/*.tsv -c ../generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv -variants final_results/variants.final.annotated.tsv -classI final_results/pVACseq/mhc_i/*.all_epitopes.aggregated.tsv -classII final_results/pVACseq/mhc_ii/*.all_epitopes.aggregated.tsv -samp $PATIENT_ID -o ../manual_review/

# Note: You can change the classI and classI IC50/percentile cutoff for coloring
python3 /opt/scripts/color_peptides51mer.py -p ../manual_review/*Peptides_51-mer.xlsx -probPos C -samp $PATIENT_ID -o ../manual_review/
```

## Creating Case Final Report locally

### Before Immunogenomics Tumor Board Review

A written case final report will be created which includes a Genomics Review Report document. This document includes a section of a basic data QC review and a table summarizing values that pass/fail the FDA quality thresholds.

### Basic data QC

Pull the basic data qc from various files. This script will output a file final_results/qc_file.txt and also print the summary to to screen.

```
cd $WORKING_BASE

docker run -it --env HOME --env WORKING_BASE -v $HOME/:$HOME/ -v $HOME/.config/gcloud:/root/.config/gcloud griffithlab/neoang_scripts:version7 /bin/bash

cd $WORKING_BASE
mkdir ../manual_review
cd ../manual_review

python3 /opt/scripts/get_neoantigen_qc.py -WB $WORKING_BASE -f final_results --yaml $HOME/yamls/${GCS_CASE_NAME}_immuno_cloud-WDL.yaml
python3 /opt/scripts/get_FDA_thresholds.py -WB  $WORKING_BASE -f final_results
python3 /opt/scripts/hla_comparison.py -WB $WORKING_BASE
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

#generate a protein fasta file using the final annotated/evaluated neoantigen candidates TSV as input
#this will filter down to only those candidates under consideration and use the top transcript

# check the file to find Tumor sample ID in the #CHROM header of VCF

gzcat $WORKING_BASE/final_results/annotated.expression.vcf.gz | less
export TUMOR_SAMPLE_ID="100-049-BG004667"

docker pull griffithlab/pvactools:4.0.5
docker run -it -v $HOME/:$HOME/ --env $WORKING_BASE  --env SAMPLE_ID griffithlab/pvactools:4.0.5 /bin/bash

cd $WORKING_BASE

pvacseq generate_protein_fasta \
  -p $WORKING_BASE/final_results/pVACseq/phase_vcf/phased.vcf.gz \
  --pass-only --mutant-only -d 150 \
  -s $TUMOR_SAMPLE_ID \
  --aggregate-report-evaluation {Accept,Review} \
  --input-tsv ../itb-review-files/*.tsv  \
  $WORKING_BASE/final_results/annotated.expression.vcf.gz \
  25 \
  $WORKING_BASE/../generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa

pvacseq generate_protein_fasta \
  -p $WORKING_BASE/final_results/pVACseq/phase_vcf/phased.vcf.gz \
  --pass-only --mutant-only -d 150 \
  -s $TUMOR_SAMPLE_ID  \
  $WORKING_BASE/final_results/annotated.expression.vcf.gz \
  25  \
  $WORKING_BASE/../generate_protein_fasta/all/annotated_filtered.vcf-pass-51mer.fa

exit 
```

To generate files needed for manual review, save the pVAC results from the Immunogenomics Tumor Board Review meeting as $SAMPLE.revd.Annotated.Neoantigen_Candidates.xlsx (Note: if the file is not saved under this exact name the below command will need to be modified).

```
docker pull griffithlab/neoang_scripts
docker run -it --env WORKING_BASE --env PATIENT_ID -v $HOME/:$HOME/ -v $HOME/.config/gcloud:/root/.config/gcloud griffithlab/neoang_scripts:version7 /bin/bash

cd $WORKING_BASE
mkdir manual_review

python3 /opt/scripts/generate_reviews_files.py -a itb-review-files/*.tsv -c generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv -classI final_results/pVACseq/mhc_i/*.all_epitopes.aggregated.tsv -classII final_results/pVACseq/mhc_ii/*.all_epitopes.aggregated.tsv -samp $PATIENT_ID -o manual_review/

python3 /opt/scripts/color_peptides51mer.py -p manual_review/*Peptides_51-mer.xlsx -samp $PATIENT_ID -o manual_review/
```
Open colored_peptides51mer.html and copy the table into an excel spreadsheet. The formatting should remain. Utilizing the Annotated.Neoantigen_Candidates and colored Peptides_51-mer for manual review.

# Description of Scripts

## Get Basic QC

```
python3  /opt/scripts/get_neoantigen_qc.py --help
usage: get_neoantigen_qc.py [-h] [-WB WB] [-f FIN_RESULTS] [--n_dna N_DNA] [--t_dna T_DNA] [--t_rna T_RNA]
                            [--concordance CONCORDANCE] [--contam_n CONTAM_N] [--contam_t CONTAM_T]
                            [--rna_metrics RNA_METRICS] [--strand_check STRAND_CHECK] --yaml YAML
                            [--fin_variants FIN_VARIANTS]

Get the stats for the basic data QC review in the neoantigen final report.

optional arguments:
  -h, --help            show this help message and exit
  -WB WB                The path to the gcp_immuno folder of the trial you wish to run the script on, defined as
                        WORKING_BASE in envs.txt
  -f FIN_RESULTS, --fin_results FIN_RESULTS
                        Name of the final results folder in gcp immuno
  --n_dna N_DNA         File path for aligned normal DNA FDA report table
  --t_dna T_DNA         File path for aligned tumor DNA FDA report table
  --t_rna T_RNA         File path for aligned tumor RNA FDA report table
  --concordance CONCORDANCE
                        File path for Somalier results for sample tumor/normal sample relatedness
  --contam_n CONTAM_N   File path for VerifyBamID results for contamination of the normal sample
  --contam_t CONTAM_T   File path for VerifyBamID results for contamination of the tumor sample
  --rna_metrics RNA_METRICS
                        File path for RNA metrics
  --strand_check STRAND_CHECK
                        File path for strandness check
  --yaml YAML           File path for the pipeline YAML file
  --fin_variants FIN_VARIANTS
                        File path for the final variants file
```

## GET FDA metrics

```
python3  /opt/scripts/get_FDA_thresholds.py --help
usage: get_FDA_thresholds.py [-h] [-WB WB] [-f FIN_RESULTS] [--n_dna N_DNA] [--t_dna T_DNA] [--t_rna T_RNA]
                             [--una_n_dna UNA_N_DNA] [--una_t_dna UNA_T_DNA] [--una_t_rna UNA_T_RNA]
                             [--somalier SOMALIER] [--contam_n CONTAM_N] [--contam_t CONTAM_T]

Get FDA qc stats from various files and determine if they pass or fail.

optional arguments:
  -h, --help            show this help message and exit
  -WB WB                the path to the gcp_immuno folder of the trial you wish to tun script on, defined as
                        WORKING_BASE in envs.txt
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
                        (concordance.somalier.pairs.tsv)
  --contam_n CONTAM_N   file path for VerifyBamID results for contamination the normal sample
  --contam_t CONTAM_T   file path for VerifyBamID results for contamination the tumor dna sample
```

## HLA Comparison
```
python3  /opt/scripts/hla_comparison.py --help
usage: hla_comparison.py [-h] [-WB WB] [-f FIN_RESULTS] [--optitype_n OPTITYPE_N] [--optitype_t OPTITYPE_T]
                         [--phlat_n PHLAT_N] [--phlat_t PHLAT_T] [--clinical CLINICAL] [--o O]

Compare HLA alleles called by phlat, opitype, and clincal data if available.

optional arguments:
  -h, --help            show this help message and exit
  -WB WB                The path to the gcp_immuno folder of the trial you wish to run the script on, defined as
                        WORKING_BASE in envs.txt
  -f FIN_RESULTS, --fin_results FIN_RESULTS
                        Name of the final results folder in gcp immuno
  --optitype_n OPTITYPE_N
                        File path for optitype normal calls
  --optitype_t OPTITYPE_T
                        File path for optitype tumor calls
  --phlat_n PHLAT_N     File path for phlat normal calls
  --phlat_t PHLAT_T     File path for phlat tumor calls
  --clinical CLINICAL   File path for the clinical_calls.txt
  --o O                 Output folder
```

## Generate Review Files

```
python3  /opt/scripts/generate_reviews_files.py --help
usage: generate_reviews_files.py [-h] -a A -c C [-variants VARIANTS] -classI CLASSI -classII CLASSII -samp SAMP
                                 [-o O] [-f FIN_RESULTS]

Create the file needed for the neoantigen manuel review

optional arguments:
  -h, --help            show this help message and exit
  -a A                  The path to the ITB Reviewed Candidates tsv file
  -c C                  The path to candidates annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv from
                        the generate_protein_fasta script
  -variants VARIANTS    The path to the variants.final.annotated.tsv file generated by the pipeline
  -classI CLASSI        The path to the classI all_epitopes.aggregated.tsv used in pVACseq
  -classII CLASSII      The path to the classII all_epitopes.aggregated.tsv used in pVACseq
  -samp SAMP            The name of the sample
  -o O                  the path to output folder
  -f FIN_RESULTS, --fin_results FIN_RESULTS
                        Name of the final results folder in gcp immuno
```

## Color Peptides 51mer

```
python3  /opt/scripts/color_peptides51mer.py --help
usage: color_peptides51mer.py [-h] -p P -samp SAMP [-cIIC50 CIIC50] [-cIpercent CIPERCENT] [-cIIIC50 CIIIC50]
                              [-cIIpercent CIIPERCENT] [-probPos [PROBPOS [PROBPOS ...]]] [-o O]

Color the 51mer peptide

optional arguments:
  -h, --help            show this help message and exit
  -p P                  The path to the Peptides 51 mer
  -samp SAMP            The name of the sample
  -cIIC50 CIIC50        Maximum classI IC50 score to annotate
  -cIpercent CIPERCENT  Maximum classI percentile to annotate
  -cIIIC50 CIIIC50      Maximum classII IC50 score to annotate
  -cIIpercent CIIPERCENT
                        Maximum classII percentile to annotate
  -probPos [PROBPOS [PROBPOS ...]]
                        problematic position to make large
  -o O                  the path to output folder
```



