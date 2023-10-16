import argparse
import subprocess

# Define the command-line arguments
parser = argparse.ArgumentParser(description='Sets up manuel review files')

parser.add_argument('-WB',
                    help='the path to the gcp_immuno folder of the trial you wish to tun script on, defined as WORKING_BASE in envs.txt')

parser.add_argument('-samp', help='Name of the sample')
parser.add_argument('-a', help='Path to ITB Reviewed Candidates')
parser.add_argument('-c', help='Path to annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv')

parser.add_argument('-classI', help='Path to classI all_epitopes.aggregated.tsv', required=True)
parser.add_argument('-classII', help='Path to classII all_epitopes.aggregated.tsv', required=True)



# Parse the command-line arguments
args = parser.parse_args()

# Check if the required arguments are provided
if not args.classI or not args.classII or not args.WB:
    parser.error("Missing required arguments. Please provide -p, -classI, -classII, and -WB.")



command1 = f"python /opt/scripts/generate_reviews_files.py -WB {args.WB} -a {args.a} -c {args.c} -samp {args.samp}"
command2 = f"python /opt/scripts/color_peptides51mer.py -WB {args.WB} -p {args.WB}/../manual_review/{args.samp}_Peptides_51-mer.xlsx -classI {args.classI} -classII {args.classII} -samp {args.samp}"


# Execute the first script
print("Generating Review Files...")

subprocess.run(command1, shell=True)

# Execute the second script
print("Coloring Peptide Sequeces...")

subprocess.run(command2, shell=True)


print("Scripts have been executed successfully.")

