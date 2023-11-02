# Neoantigen Pipeline Helper Scripts


- Docker
```bash
  bsub -n 1 -Is -G compute/ -g /evelyn/default -q general-interactive -M 16G -R 'rusage[mem=16G]' -a 'docker('evelyns2000/neoang_scripts')' /bin/bash
```
  
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

#### Generate Protein Fasta

```bash
cd $WORKING_BASE
mkdir ../generate_protein_fasta
cd ../generate_protein_fasta
mkdir candidates
mkdir all

zcat $WORKING_BASE/final_results/annotated.expression.vcf.gz | less # Get sample ID Found in the #CHROM header of VCF
export SAMPLE_ID="TWJF-10146-0029-0029_Tumor_Lysate"
export ITB_REVIEW_FILE=10146-0029.Annotated.Neoantigen_Candidates.Revd.tsv


bsub -Is -q general-interactive -G $GROUP -a "docker(griffithlab/pvactools:4.0.1)" /bin/bash

pvacseq generate_protein_fasta \
  -p $WORKING_BASE/final_results/pVACseq/phase_vcf/phased.vcf.gz \
  --pass-only --mutant-only -d 150 \
  -s $SAMPLE_ID \
  --aggregate-report-evaluation {Accept,Review} \
  --input-tsv ../itb-review-files/$ITB_REVIEW_FILE  \
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

## Creating Case Final Report on locally

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
