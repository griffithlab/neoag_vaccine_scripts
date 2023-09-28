import argparse
import subprocess

# Define the command-line arguments
parser = argparse.ArgumentParser(description='Sets up manuel review files')

parser.add_argument('-samp', help='Name of the sample')
parser.add_argument('-a', help='Path to ITB Reviewed Candidates')
parser.add_argument('-c', help='Path to annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv')

parser.add_argument('-classI', help='Path to classI all_epitopes.aggregated.tsv', required=True)
parser.add_argument('-classII', help='Path to classII all_epitopes.aggregated.tsv', required=True)
parser.add_argument('-o', help='The name of the html in which the coloring will be printed', required=True)



# Parse the command-line arguments
args = parser.parse_args()

# Check if the required arguments are provided
if not args.classI or not args.classII or not args.o:
    parser.error("Missing required arguments. Please provide -p, -classI, -classII, and -o.")



command1 = f"python /opt/scripts/generate_reviews_files.py -a {args.a} -c {args.c} -samp {args.samp}"
command2 = f"python /opt/scripts/color_peptides51mer.py -p {args.samp}_Peptides_51-mer.xlsx -classI {args.classI} -classII {args.classII} -o {args.o}"


# Execute the first script
print("Generating Review Files...")
subprocess.run(command1, shell=True)

# Execute the second script
print("Coloring Peptide Sequeces...")
subprocess.run(command2, shell=True)

print("Scripts have been executed successfully.")
