import argparse
import os
import csv

# Notes
"""
A script to generate the Basic data QC for Genomic Review Reports
This script works with various trials and includes optional argument flags in case of weird naming conventions.

Author: Evelyn Schmidt
Date: August 2023
"""

# ---- PARSE ARGUMENTS -------------------------------------------------------
def parse_arguments():
    """Parses command line arguments and enables user help."""
    parser = argparse.ArgumentParser(description='Get the stats for the basic data QC review in the neoantigen final report.')

    parser.add_argument('-WB', help='The path to the gcp_immuno folder of the trial you wish to run the script on, defined as WORKING_BASE in envs.txt')
    parser.add_argument('-f', "--fin_results", help="Name of the final results folder in gcp immuno")
    parser.add_argument("--n_dna", help="File path for aligned normal DNA FDA report table")
    parser.add_argument("--t_dna", help="File path for aligned tumor DNA FDA report table")
    parser.add_argument("--t_rna", help="File path for aligned tumor RNA FDA report table")
    parser.add_argument("--concordance", help="File path for Somalier results for sample tumor/normal sample relatedness")
    parser.add_argument("--contam_n", help="File path for VerifyBamID results for contamination of the normal sample")
    parser.add_argument("--contam_t", help="File path for VerifyBamID results for contamination of the tumor sample")
    parser.add_argument("--rna_metrics", help="File path for RNA metrics")
    parser.add_argument("--strand_check", help="File path for strandness check")
    parser.add_argument("--yaml", help="File path for the pipeline YAML file", required=True)
    parser.add_argument("--fin_variants", help="File path for the final variants file")

    return parser.parse_args()


# ---- EVALUATION FUNCTIONS --------------------------------------------------
def evaluate_unique_mapped_reads(count):
    """Evaluates the unique mapped reads count."""
    if count < 40000000:
        return "poor"
    elif 40000000 <= count < 50000000:
        return "acceptable"
    elif 50000000 <= count < 100000000:
        return "good"
    else:
        return "excellent"

def evaluate_duplication_rate(rate):

    print(rate)
    """Evaluates the duplication rate."""
    if rate > 75:
        return "very poor"
    elif 50 < rate <= 75:
        return "poor"
    elif 30 < rate <= 50:
        return "acceptable"
    elif 20 < rate <= 30:
        return "good"
    else:
        return "excellent"

def evaluate_relatedness(relatedness):
    """Evaluates the relatedness."""
    if relatedness > 0.975:
        return "excellent"
    elif 0.95 < relatedness <= 0.975:
        return "good"
    elif 0.90 < relatedness <= 0.95:
        return "concerning"
    else:
        return "very concerning"

def evaluate_contamination(rate):
    """Evaluates the contamination rate."""
    if rate < 0.025:
        return "good"
    elif 0.025 <= rate < 0.05:
        return "concerning"
    else:
        return "very concerning"

def evaluate_rna_alignment(alignment):
    """Evaluates the RNA alignment rate."""
    if alignment > 0.90:
        return "excellent"
    elif 0.75 < alignment <= 0.90:
        return "good"
    elif 0.50 < alignment <= 0.75:
        return "acceptable"
    else:
        return "sub-optimal"

# ---- TOTAL NUMBER OF UNIQUELY MAPPING READS and DUPLICATION -----------------
def get_read_pairs(normal_dna, tumor_dna, tumor_rna):
    """Summarizes total number of uniquely mapping reads and duplication rates."""
    read_pairs_report_string = ""
    file_paths = [normal_dna, tumor_dna, tumor_rna]

    for file_path in file_paths:
        try:
            with open(file_path, 'r') as file:
                for line in file:
                    if 'Unique Mapped Reads' in line:
                        parts = line.split('\t')
                        if len(parts) > 1:
                            number = int(parts[1].strip())
                            evaluation = evaluate_unique_mapped_reads(number)
                            print(os.path.basename(file_path), "Unique Map Reads:", f'{number:,d}', "(", evaluation, ")")
                            read_pairs_report_string += os.path.basename(os.path.normpath(file_path)) + " Unique Map Reads: " + f'{number:,d}' + " (" + evaluation + ")\n"
        except FileNotFoundError:
            print(f"File {file_path} not found.")
    
    for file_path in file_paths:
        try:
            with open(file_path, 'r') as file:
                for line in file:
                   if 'Mapped Read Duplication' in line and os.path.basename(os.path.normpath(file_path)) != 'tumor_rna_aligned_metrics.txt':
                        parts = line.split('\t')
                        if len(parts) > 1:
                            number = parts[2].strip()
                            evaluation = evaluate_duplication_rate(float(number.replace(' (%)','')))
                            print(os.path.basename(file_path), "Mapped Read Duplication Rate:", number, "(", evaluation, ")")
                            read_pairs_report_string += os.path.basename(file_path) + " Mapped Read Duplication Rate: " + number + " (" + evaluation + ")\n"
        except FileNotFoundError:
            print(f"File {file_path} not found.")
    return read_pairs_report_string

# ---- SAMPLE RELATEDNESS ----------------------------------------------------
def get_relatedness(concordance):
    """Checks Somalier results for sample tumor/normal sample relatedness."""
    try:
        with open(concordance, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                relatedness = float(row['relatedness'])
                evaluation = evaluate_relatedness(relatedness)
                break
            print('Relatedness:', relatedness, "(", evaluation, ")")
            return "Relatedness: " + str(relatedness) + " (" + evaluation + ")\n"
    except FileNotFoundError:
        print(f"File {concordance} not found.")
        return ""

# ---- CONTAMINATION ---------------------------------------------------------
def get_contamination(contamination_normal, contamination_tumor):
    """Checks VerifyBamID results for contamination of both tumor and normal samples."""
    contamination_str = ""
    file_paths = [contamination_normal, contamination_tumor]

    for file_path in file_paths:
        try:
            with open(file_path, 'r') as file:
                reader = csv.DictReader(file, delimiter='\t')
                for row in reader:
                    contamination_rate = row['FREEMIX']
                    evaluation = evaluate_contamination(float(contamination_rate))
                    print(os.path.basename(file_path), 'Contamination:', contamination_rate, "(", evaluation, ")")
                    contamination_str += os.path.basename(file_path) + " Contamination: " + str(contamination_rate) + " (" + evaluation + ")\n"
        except FileNotFoundError:
            print(f"File {file_path} not found.")
    return contamination_str

# ---- PERCENT READS ALIGNING TO TRANSCRIPTS ---------------------------------
def get_rna_alignment(rna_metrics):
    """Checks RNA-seq metrics for % reads aligning to transcripts."""
    try:
        with open(rna_metrics, 'r') as file:
            lines = file.readlines()
            for i in range(len(lines) - 2):
                current_line = lines[i].strip()
                next_line = lines[i + 1].strip()
                if current_line.startswith("PF_BASES"):
                    keys = current_line.split("\t")
                    values = next_line.split("\t")
                    break
        pct_coding_bases = pct_utr_bases = None
        for i in range(len(keys)):
            if keys[i] == "PCT_CODING_BASES":
                pct_coding_bases = values[i]
            if keys[i] == "PCT_UTR_BASES":
                pct_utr_bases = values[i]

        if pct_coding_bases and pct_utr_bases:
            rna_reads_mapped = float(pct_coding_bases) + float(pct_utr_bases)
            evaluation = evaluate_rna_alignment(rna_reads_mapped)
            print("The proportion of RNA reads mapping to cDNA sequence is",
                  rna_reads_mapped,
                  "(coding (",
                  pct_coding_bases,
                  ") + UTR (",
                  pct_utr_bases, ")", "(", evaluation, ")")
            return ("The proportion of RNA reads mapping to cDNA sequence is " +
                    str(rna_reads_mapped) +
                    " (coding (" +
                    pct_coding_bases +
                    ") + UTR (" +
                    pct_utr_bases + ")) (" + evaluation + ")\n")
    except FileNotFoundError:
        print(f"File {rna_metrics} not found.")
    return ""

# ---- STRAND CHECK ----------------------------------------------------------
def check_strand(strandness_check, yaml_file):
    """Checks the correct RNA strand setting was used in the pipeline YAML file."""
    strand_str = ""
    try:
        with open(strandness_check, 'rb') as f:
            try:
                f.seek(-2, os.SEEK_END)
                while f.read(1) != b'\n':
                    f.seek(-2, os.SEEK_CUR)
            except OSError:
                f.seek(0)
            last_line = f.readline().decode().strip()
        print("trimmed_read_1strandness_check.txt:", last_line)
        strand_str += "trimmed_read_1strandness_check.txt: " + last_line + "\n"

        with open(yaml_file, 'r') as file:
            for line in file:
                if "strand" in line:
                    print("YAML file:", line.strip())
                    strand_str += "YAML file: " + line.strip() + "\n"
    except FileNotFoundError as e:
        print(e)
    return strand_str

# ---- TOTAL VARIANTS---------------------------------------------------------
def get_variant_count(final_variants):
    """Summarizes total variants called and total neoantigen variants called."""
    try:
        with open(final_variants, 'r') as fp:
            lines = len(fp.readlines())
            print('Total Number of somatic variants called:', lines - 1)
            return 'Total Number of somatic variants called: ' + str(lines - 1) + "\n"
    except FileNotFoundError:
        print(f"File {final_variants} not found.")
        return ""

def main():
    args = parse_arguments()

    final_result = f"/{args.fin_results}" if args.fin_results else '/final_results'

    normal_dna = args.n_dna if args.n_dna else f"{args.WB}{final_result}/qc/fda_metrics/aligned_normal_dna/table_metrics/normal_dna_aligned_metrics.txt"
    tumor_dna = args.t_dna if args.t_dna else f"{args.WB}{final_result}/qc/fda_metrics/aligned_tumor_dna/table_metrics/tumor_dna_aligned_metrics.txt"
    tumor_rna = args.t_rna if args.t_rna else f"{args.WB}{final_result}/qc/fda_metrics/aligned_tumor_rna/table_metrics/tumor_rna_aligned_metrics.txt"
    concordance = args.concordance if args.concordance else f"{args.WB}{final_result}/qc/concordance/concordance.somalier.pairs.tsv"
    contamination_normal = args.contam_n if args.contam_n else f"{args.WB}{final_result}/qc/normal_dna/normal.VerifyBamId.selfSM"
    contamination_tumor = args.contam_t if args.contam_t else f"{args.WB}{final_result}/qc/tumor_dna/tumor.VerifyBamId.selfSM"
    rna_metrics = args.rna_metrics if args.rna_metrics else f"{args.WB}{final_result}/qc/tumor_rna/rna_metrics.txt"
    strandness_check = args.strand_check if args.strand_check else f"{args.WB}{final_result}/qc/tumor_rna/trimmed_read_1strandness_check.txt"
  
    # yaml is always in different place
    if args.yaml:
        yaml_file = args.yaml
    
    final_variants = args.fin_variants if args.fin_variants else f"{args.WB}{final_result}/variants.final.annotated.tsv"
        
    # create a text file to store results
    qc_file = open(args.WB + '/../manual_review/qc_file.txt', 'w') if args.WB else pen('qc_file.txt', 'w')
       

    print()
    print()

    qc_file.write(get_read_pairs(normal_dna, tumor_dna, tumor_rna))
    qc_file.write(get_relatedness(concordance))
    qc_file.write(get_contamination(contamination_normal, contamination_tumor))
    qc_file.write(get_rna_alignment(rna_metrics))
    qc_file.write(check_strand(strandness_check,  yaml_file))
    qc_file.write(get_variant_count(final_variants))
    qc_file.write("REMEMBER to visually inspect end bias plot (usually found in qc/tumor_rna/rna_metrics.pdf)")

    print()
    print("REMEMBER to visually inspect end bias plot (usually found in qc/tumor_rna/rna_metrics.pdf)")



if __name__ == "__main__":
    main()
