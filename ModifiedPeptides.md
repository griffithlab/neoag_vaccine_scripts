## Test solubility modifications to peptides

### Preamble
Starting with a list of peptides and proposed amino acid modifications on the N- or C-terminal ends (or automatically selected modifications), aimed at improving the solubility of synthesize long peptides, test for creation of strong binding peptides containing modified amino acids. Summarize these findings and avoid use of modified peptides that lead to predicted strong binding peptide containing these synthetic modifications.

### Local dependencies
The following assumes you have gcloud installed and have authenticated to use the google cloud project below

Set up Google cloud configurations and make sure the right one is activated:
```bash 
export GCS_PROJECT=jlf-rcrf
export GCS_VM_NAME=mg-test-peptide-mods 

#list possible configs that are set up
gcloud config configurations list

#activate the rcrf config
gcloud config configurations activate rcrf

#login if needed (only needs to be done once)
gcloud auth login 

#view active config/login (should show the correct project "jlf-rcrf", zone, and email address)
gcloud config list

```

Configure these configurations to use a specific zone. Once the config is setup and you have logged into at least once the config files should look like this:

`cat ~/.config/gcloud/configurations/config_rcrf`

```
[compute]
region = us-central1
zone = us-central1-c
[core]
account = <email address associated with rcrf account>
disable_usage_reporting = True
project = jlf-rcrf
```

### Launching a Google VM to perform the predictions
Launch GCP instance set up for ad hoc analyses (including docker)

```bash
gcloud compute instances create $GCS_VM_NAME --service-account=cromwell-server@$GCS_PROJECT.iam.gserviceaccount.com --source-machine-image=jlf-adhoc-v1 --network=cloud-workflows --subnet=cloud-workflows-default --boot-disk-size=250GB --boot-disk-type=pd-ssd --machine-type=e2-standard-8
```

### Log into the GCP instance and check status

```bash
gcloud compute ssh $GCS_VM_NAME 

#confirm start up scripts have completed. use <ctrl> <c> to exit
journalctl -u google-startup-scripts -f

#check for expected disk space
df -h 

```

### Configure Docker to work for current user

```bash
sudo usermod -a -G docker $USER
sudo reboot

```

Logout and login to get this change to take effect and test the docker install
```bash
exit

gcloud compute ssh $GCS_VM_NAME 

docker run hello-world

```


## Generating the peptide_table.tsv file

The input is a csv file which has the names and base sequences. The name does not have to be unique.
For example, the gene name to indentify each sequence. The sequences should not contian any modifications. 

```bash
Base Sequence name,Base sequence
CUL9,RMLDYYEEISAGDEGEFRQS
CUL9,RVRMLDYYEEISAGDEGEFRQSN
CUL9,RVRMLDYYEEISAGDEGEFR
EXOC4,SVIRTLSTIDDVEDRENEKGR
EXOC4,ISVIRTLSTIDDVEDRENEKGR
EXOC4,LISVIRTLSTIDDVEDRENEKGR
VPS13B,GLRQGLFRLGISLLGAIAGIVD
VPS13B,GEGLRQGLFRLGISLLGAIAG
VPS13B,SLGEGLRQGLFRLGISLLGAI
DYNC1H1,KRFHATISFDTDTGLKQALET
DYNC1H1,GKRFHATISFDTDTGLKQALET
DYNC1H1,KRFHATISFDTDTGLKQAL
MYO9A,FDWIVFRINHALLNSKVLEHNTK
MYO9A,FDWIVFRINHALLNSKVL
MYO9A,SALFDWIVFRINHALLNSKVLEHN
EPG5,KELPLYLWQPSTSEIAVIRDW
```

```
export HOME=/Users/evelynschmidt/jlf/JLF-100-048/ModifiedPeptides
export HLA_ALLELES=HLA-A*24:02,HLA-A*29:02,HLA-B*14:02,HLA-B*14:02,HLA-C*02:02,HLA-C*08:02
export SAMPLE_NAME=jlf-100-048
```



Using the Docker and executing this command will produce the peptide_table.tsv used by pVAC bind. 
The -n argument is the maximum number of modified peptides and the -m argument is a path to the csv file.

```
docker pull griffithlab/neoang_scripts:latest

docker run -it -v $HOME/:$HOME/ -v $HOME/.config/gcloud:/root/.config/gcloud griffithlab/neoang_scripts:latest /bin/bash

python3 scripts/modify_peptides.py -n 3 --100-043/ModifiedPeptides/sequences.csv  -samp $SAMPLE_NAM -HLA $HLA_ALLELES -WD $HOME
```

For example, if you speficify -n 1 then the modified sequences produced will me:

```bash
CUL9.1.n-term-K	KRMLDYYEEISAGDEGEFRQS	K|RMLDYYEEISAGDEGEFRQS
CUL9.1.c-term-K	RMLDYYEEISAGDEGEFRQSK	RMLDYYEEISAGDEGEFRQS|K
CUL9.1.n-term-R	RRMLDYYEEISAGDEGEFRQS	R|RMLDYYEEISAGDEGEFRQS
CUL9.1.c-term-R	RMLDYYEEISAGDEGEFRQSR	RMLDYYEEISAGDEGEFRQS|R
.
.
.
EPG5.n-term-K	KKELPLYLWQPSTSEIAVIRDW	K|KELPLYLWQPSTSEIAVIRDW
EPG5.c-term-K	KELPLYLWQPSTSEIAVIRDWK	KELPLYLWQPSTSEIAVIRDW|K
EPG5.n-term-R	RKELPLYLWQPSTSEIAVIRDW	R|KELPLYLWQPSTSEIAVIRDW
EPG5.c-term-R	KELPLYLWQPSTSEIAVIRDWR	KELPLYLWQPSTSEIAVIRDW|R
```

### Enter a pVACtools docker environment to run pVACbind on the sub-peptide sequences containing modified AAs

```bash
docker pull griffithlab/pvactools:4.0.5
docker run -it -v $HOME/:$HOME/ --user $(id -u):$(id -g) --env HOME --env SAMPLE_NAME --env HLA_ALLELES griffithlab/pvactools:4.0.5 /bin/bash
cd $HOME

for LENGTH in 8 9 10 11
do 
   #process n-term fasta for this length
   echo "Running pVACbind for length: $LENGTH (n-term sequences)"
   export LENGTH_FASTA=$HOME/n-term/pvacbind_inputs/${LENGTH}-mer-test.fa
   export LENGTH_RESULT_DIR=$HOME/n-term/pvacbind_results/${LENGTH}-mer-test
   pvacbind run $LENGTH_FASTA $SAMPLE_NAME $HLA_ALLELES all_class_i $LENGTH_RESULT_DIR -e1 $LENGTH --n-threads 8 --iedb-install-directory /opt/iedb/ 1>$LENGTH_RESULT_DIR/stdout.txt 2>$LENGTH_RESULT_DIR/stderr.txt

   #process c-term fasta for this length
   echo "Running pVACbind for length: $LENGTH (c-term sequences)"
   export LENGTH_FASTA=$HOME/c-term/pvacbind_inputs/${LENGTH}-mer-test.fa
   export LENGTH_RESULT_DIR=$HOME/c-term/pvacbind_results/${LENGTH}-mer-test
   pvacbind run $LENGTH_FASTA $SAMPLE_NAME $HLA_ALLELES all_class_i $LENGTH_RESULT_DIR -e1 $LENGTH --n-threads 8 --iedb-install-directory /opt/iedb/ 1>$LENGTH_RESULT_DIR/stdout.txt 2>$LENGTH_RESULT_DIR/stderr.txt
done

```


To check for successful completion of all jobs you can check the stdout logs that have been saved. There should be 8 successful jobs total, 4 lengths for n-term modified peptides and 4 lengths for c-term.

```bash
grep "Pipeline finished" */pvacbind_results/*/stdout.txt | wc -l
```

### Combine all the pVACbind results into a single file
Create a combined TSV file by concatenating all the individual "all_epitopes.tsv" files and avoiding redundant headers. Store this file locally (or in a cloud bucket) so that it can be accessed after the VM is destroyed.

```bash
#get the header line
grep -h "^Mutation" --color=never */pvacbind_results/*/MHC_Class_I/${SAMPLE_NAME}.all_epitopes.tsv | sort | uniq > header.tsv

#combine the results from all prediction runs and add the header on
cat */pvacbind_results/*/MHC_Class_I/${SAMPLE_NAME}.all_epitopes.tsv | grep -v "^Mutation" | cat header.tsv - > ${SAMPLE_NAME}.all_epitopes.all_modifications.tsv

```

### Evaluate the proposed modified peptide sequences
The goal of this analysis is to test whether any strong binding peptides are created that include the modified amino acid sequences included to improve solubility. For example, one could require that no such peptides exist where the median binding affinity is < 500nm OR median binding score percentile is < 1%.

For each candidate modified peptide sequence, summarize the number of such potentially problematic peptides. 

```bash

#pull out all the rows that correspond to strong binders according to default criteria (<500nm affinity OR <1 percentile score)
cut -f 1,2,4,5,8 ${SAMPLE_NAME}.all_epitopes.all_modifications.tsv | perl -ne 'chomp; @l=split("\t",$_); $median_affinity=$l[3]; $median_percentile=$l[4]; if ($median_affinity < 500 || $median_percentile < 1){print "$_\n"}' > ${SAMPLE_NAME}.all_epitopes.all_modifications.problematic.tsv

#summarize number of problematic results of each unique candidate proposed peptide
cat ${SAMPLE_NAME}.all_epitopes.all_modifications.problematic.tsv | grep -v "^Mutation" | cut -f 1 | sort | uniq -c | sed 's/^[ ]*//' | tr " " "\t" | awk 'BEGIN {FS="\t"; OFS="\t"} {print $2, $1}' > ${SAMPLE_NAME}.problematic.summary.tsv

#create a list of all unique peptide names for modified peptides to be summarized
cut -f 1 ${SAMPLE_NAME}.all_epitopes.all_modifications.tsv | grep -v "^Mutation" | sort | uniq > peptide_name_list.tsv

#create an output table with a count of problematic binders for all peptides (include 0 if that is the case)
join -t $'\t' -a 1 -a 2 -e'0' -o '0,2.2' peptide_name_list.tsv ${SAMPLE_NAME}.problematic.summary.tsv > ${SAMPLE_NAME}.problematic.summary.complete.tsv

```

### Retrieve final result files to local system

Files to be kept:

- ${SAMPLE_NAME}.all_epitopes.all_modifications.tsv
- ${SAMPLE_NAME}.all_epitopes.all_modifications.problematic.tsv
- ${SAMPLE_NAME}.problematic.summary.complete.tsv

```bash
#leave the GCP VM
exit

export SAMPLE_NAME="jlf-100-026"

mkdir ${SAMPLE_NAME}_modified_peptide_results
cd ${SAMPLE_NAME}_modified_peptide_results

gcloud compute scp $USER@$GCS_VM_NAME:${SAMPLE_NAME}.all_epitopes.all_modifications.tsv ${SAMPLE_NAME}.all_epitopes.all_modifications.tsv

gcloud compute scp $USER@$GCS_VM_NAME:${SAMPLE_NAME}.all_epitopes.all_modifications.problematic.tsv ${SAMPLE_NAME}.all_epitopes.all_modifications.problematic.tsv

gcloud compute scp $USER@$GCS_VM_NAME:${SAMPLE_NAME}.problematic.summary.complete.tsv ${SAMPLE_NAME}.problematic.summary.complete.tsv


```

### Once the analysis is done and results retrieved, destroy the Google VM on GCP to avoid wasting resources

```bash

gcloud compute instances delete $GCS_VM_NAME

```

### Final report generation and interpretation
Use the information in `${SAMPLE_NAME}.all_epitopes.all_modifications.tsv` and `${SAMPLE_NAME}.problematic.summary.complete.tsv` to produce summary spreadsheets.
