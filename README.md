# Neoantigen Pipeline Helper Scripts


- Docker
```bash
  bsub -n 1 -Is -G compute/ -g /evelyn/default -q general-interactive -M 16G -R 'rusage[mem=16G]' -a 'docker('evelyns2000/neoang_scripts')' /bin/bash
```
  
## Get FDA Thresholds

- List of parameters
```bash
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
