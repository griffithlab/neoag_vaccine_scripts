import argparse
import os
import csv
import sys

# Notes
"""
Tested for JLF trials -- wdl sometimes in different place
Tested for CTEP -- wdl/fianl variants work too
Tested for Leidos -- wdl/fianl variants work too

Use:
provide Trial name, then sample name folder, wdl name

ex: $ python3 get_neoantigen_qc.py CTEP TWJF-10146-MO011-0024 10146-MO011-0024_immuno_cloud-WDL.yaml

"""

argv = sys.argv[1:]


# Locations for Various Trilas
CTEP = "/storage1/fs1/gillandersw/Active/Project_0001_Clinical_Trials/CTEP/analysis/"
JLF = "/storage1/fs1/mgriffit/Active/JLF_MCDB/cases/"
LEIDOS = "/storage1/fs1/gillandersw/Active/Project_0001_Clinical_Trials/pancreas_leidos/analysis/"
WARD = "/storage1/fs1/jward2/Active/Project_HRPO_202104143_Clinical_Trial/SCLC/"
COMPASSIONATEUSE = "/storage1/fs1/gillandersw/Active/Project_0001_Clinical_Trials/compassionate_use/analysis/"


# Sample name 
sample_name = argv[1]

# Root path
if argv[0] == "CTEP":
    root_path = CTEP
elif argv[0] == "JLF":
    root_path = JLF
elif argv[0] == "Leidos":
    root_path = LEIDOS
elif argv[0] == "Ward":
    root_path = WARD
else:
    root_path = COMPASSIONATEUSE

wdl_file = argv[2]




# ---- PARSE ARGUMENTS -------------------------------------------------------
# Parses command line arguments
# Enables user help
def parse_arguments():
    # Parse command line arugments
    parser = argparse.ArgumentParser(description='Get the stats for the basic data QC review in the neoantigen final report.')

    parser.add_argument('-WB',
                        help='the path to the gcp_immuno folder of the trial you wish to tun script on, defined as WORKING_BASE in envs.txt')

    # The name of the final results folder 
    parser.add_argument('-f', "--fin_results", help="Name of the final results folder in gcp immuno")

    # Allows user to specify direct path to every file in case of werid naming convention
    parser.add_argument("--n_dna", help="file path for aligned normal dna FDA report table")
    parser.add_argument("--t_dna", help="file path for aligned tumor dna FDA report table")
    parser.add_argument("--t_rna", help="file path for aligned tumor rna FDA report table")
    parser.add_argument("--somalier", help="file path for Somalier results for sample tumor/normal sample relatedness")
    parser.add_argument("--contam_n", help="file path for VerifyBamID results for contamination the normal sample")
    parser.add_argument("--contam_t", help="file path for VerifyBamID results for contamination the tumor sample")
    parser.add_argument("--rna_metrics", help="")
    parser.add_argument("--strand_check", help="")
    parser.add_argument("--yaml", help="")
    parser.add_argument("--fin_variants", help="")



    return(parser.parse_args())

# ---- TOTAL NUMBER OF UNIQELY MAPPING READS and DUPLICATION -----------------
# Summarize total number of uniquely mapping reads generated for Tumor/Normal 
# exome and Tumor RNA-seq (i.e. not counting duplicates that arise from PCR 
# amplification of the same source DNA fragment).
# Summarize duplication rates for tumor/normal DNA samples
def get_read_pairs(normal_dna, tumor_dna, tumor_rna):

    read_pairs_report_string = ""

    file_paths = [normal_dna, tumor_dna, tumor_rna]

    for file_path in file_paths:
        # Open the file in read mode
        with open(file_path, 'r') as file:

            # Search for the string and extract the number
            for line in file:
                if 'Unique Mapped Reads' in line:
                    # Split the line by tab character
                    parts = line.split('\t')
                    if len(parts) > 1:
                        number = parts[1].strip()
                        print(os.path.basename(os.path.normpath(file_path)), "Unique Map Reads:", number)
                        read_pairs_report_string += os.path.basename(os.path.normpath(file_path)) + " Unique Map Reads: " + number + "\n"
    
    for file_path in file_paths:
        # Open the file in read mode
        with open(file_path, 'r') as file:
            for line in file:    
                if 'Mapped Read Duplication' in line and os.path.basename(os.path.normpath(file_path)) != 'tumor_rna_aligned_metrics.txt':
                    # Split the line by tab character
                    parts = line.split('\t')
                    if len(parts) > 1:
                        number = parts[2].strip()
                        print(os.path.basename(os.path.normpath(file_path)), "Mapped Read Duplication Rate:", number)
                        read_pairs_report_string += os.path.basename(os.path.normpath(file_path)) + " Mapped Read Duplication Rate: " + number + "\n"
        
    return read_pairs_report_string


# ---- SAMPLE RELATEDNESS ----------------------------------------------------
# Check Somalier results for sample tumor/normal sample relatedness
def get_relatedness(somalier):
    # Open the TSV file
    with open(somalier, 'r') as file:
        # Create a TSV reader object
        reader = csv.DictReader(file, delimiter='\t')
        
        # Get the first entry in the specified column
        for row in reader:
            first_entry = row['relatedness']
            break
        
        # Print the first entry
        print('relatedness:', first_entry)

        return("relatedness: " + first_entry  + "\n")



# ---- CONTAMINATION ---------------------------------------------------------
# Check VerifyBamID results for contamination of both tumor and normal samples
def get_contaimination(contamination_normal, contamination_tumor):
    contamination_str = ""

    file_paths = [contamination_normal, contamination_tumor]

    for file_path in file_paths:
        with open(file_path, 'r') as file:
            # Create a TSV reader object
            reader = csv.DictReader(file, delimiter='\t')
            
            # Get the first entry in the specified column
            for row in reader:
                first_entry = row['FREEMIX']
                break
            
            # Print the first entry
            print(os.path.basename(os.path.normpath(file_path)) + ' contaimination: ' + first_entry)
            contamination_str += (os.path.basename(os.path.normpath(file_path)) + "contaimination: " + first_entry + "\n")

    return contamination_str


# ---- PERCENT READS ALIGNING TO TRANSCRIPTS ---------------------------------
# Check RNA-seq metrics for % reads aligning to transcripts
def get_rna_alignment(rna_metrics):
    # Open the file
    with open(rna_metrics, 'r') as file:

        line = file.readlines()

        total_lines = len(line)
        for i in range(total_lines - 2):

            current_line = line[i].strip()
            next_line = line[i + 1].strip()

            if current_line.startswith("PF_BASES"):
                keys = current_line.split("\t")
                values = next_line.split("\t")
                break

    for i in range(len(keys)):
        if keys[i] == "PCT_CODING_BASES":
            pct_coding_bases = values[i]
        if keys[i] == "PCT_UTR_BASES":
            pct_utr_bases = values[i]

    rna_reads_mapped = float(pct_coding_bases) + float(pct_utr_bases)

    print("The proportion of RNA reads mapping to cDNA sequence is ", 
        rna_reads_mapped, 
        " (coding (",
        pct_coding_bases,
        ") + UTR (",
        pct_utr_bases, ")")

    return("The proportion of RNA reads mapping to cDNA sequence is " +
        str(rna_reads_mapped) +
        " (coding (" +
        pct_coding_bases +
        ") + UTR (" +
        pct_utr_bases + ")" + "\n")


# ---- STRAND CHECK ----------------------------------------------------------
# Check that the correct RNA strand setting was used in the pipeline YAML file
def check_strand(strandness_check, yaml_file):

    strand_str = ""

    # seek to the end of the file, and move backwards to find a newline
    # file has to be opened in binary mode, otherwise, it will be impossible to seek from the end
    with open(strandness_check, 'rb') as f:
        try:  # catch OSError in case of a one line file 
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        last_line = f.readline().decode()

    # remove extra newline
    print("trimmed_read_1strandness_check.txt: " + last_line.strip())
    strand_str += "trimmed_read_1strandness_check.txt: " + last_line.strip() + "\n"

    with open(yaml_file, 'r') as file:
        for line in file:
            if "strand" in line:
                print("YAML file: " + line.strip())
                strand_str += "YAML file: " + line.strip() + "\n"
    
    return strand_str

# ---- TOTAL VARIANTS---------------------------------------------------------
# Summarize total variants called and total neoantigen variants called 
# (i.e. the subset of variants that could lead to neoantigens).  
# Briefly, total variants is how many somatic mutations were detected 
# in the tumor. Total neoantigen variants is the subset of these that 
# lead to possible neoantigen candidates (i.e. those that cause protein 
# coding changes in known genes). Note that the number of neoantigen 
# candidates selected by the Immunotherapy Tumor Board will be a small 
# subset of this number.
def get_variant_count(final_variants):
    with open(final_variants, 'r') as fp:
        lines = len(fp.readlines())
        print('Total Number of somatic variants called:', lines-1)
        return('Total Number of somatic variants called: ' + str(lines-1) + "\n")




def main():
    # create a text file to store results
    qc_file = open(root_path + sample_name + '/gcp_immuno/final_results/qc_file.txt', 'w')

    args = parse_arguments()

    if args.fin_results:
        final_result = '/' + args.fin_results
    else:
        final_result = "/final_results"


    # Paths for various files needed
    if args.n_dna:
        normal_dna = args.n_dna
    else:
        normal_dna = args.WB + final_result +  '/qc/fda_metrics/aligned_normal_dna/table_metrics/normal_dna_aligned_metrics.txt'
    if args.t_dna:
        tumor_dna = args.t_dna
    else:
        tumor_dna = args.WB + final_result + '/qc/fda_metrics/aligned_tumor_dna/table_metrics/tumor_dna_aligned_metrics.txt'
    if args.t_rna:
        tumor_rna = args.t_rna
    else:
        tumor_rna = args.WB + final_result + '/qc/fda_metrics/aligned_tumor_rna/table_metrics/tumor_rna_aligned_metrics.txt'
    if args.somalier:
        somalier = args.somalier
    else:
        somalier = args.WB + final_result +  '/qc/concordance/concordance.somalier.pairs.tsv'
    if args.contam_n:
        contamination_normal = args.contam_n
    else:
        contamination_normal = args.WB + final_result + '/qc/normal_dna/normal.VerifyBamId.selfSM'
    if args.contam_t:
        contamination_tumor = args.contam_t
    else:
        contamination_tumor = args.WB + final_result + '/qc/tumor_dna/tumor.VerifyBamId.selfSM'
    if args.rna_metrics:
        rna_metrics = args.rna_metrics
    else:
        rna_metrics = args.WB + final_result + '/qc/tumor_rna/rna_metrics.txt'
    if args.strand_check:
        strandness_check = args.strand_check
    else:
        strandness_check = args.WB + final_result + '/qc/tumor_rna/trimmed_read_1strandness_check.txt'
    if args.yaml:
        yaml_file = args.yaml
    else:
        yaml_file = args.WB + final_result + 'workflow_artifacts/*.yaml' 
    if args.final_variants:
        final_variants = args.final_variants
    else:
        final_variants = args.WB + final_result + '/variants.final.annotated.tsv'




    
    qc_file.write(get_read_pairs(normal_dna, tumor_dna, tumor_rna))
    qc_file.write(get_relatedness(somalier))
    qc_file.write(get_contaimination(contamination_normal, contamination_tumor))
    qc_file.write(get_rna_alignment(rna_metrics))
    qc_file.write(check_strand(strandness_check,  yaml_file))
    qc_file.write(get_variant_count(final_variants))



if __name__ == "__main__":
    main()