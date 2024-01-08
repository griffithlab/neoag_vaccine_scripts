import argparse
from itertools import product, chain
import pandas as pd
import sys
import os
import re
import subprocess



# ---- PARSE ARGUMENTS -------------------------------------------------------
# Parses command line arguments
# Enables user help
def parse_arguments():
    # Parse command line arugments
    parser = argparse.ArgumentParser(description='Create peptide_table.tsv needed for modifying peptide sequences for pVACbind')

    parser.add_argument('-n',
                        help='The maximum number of mofifying peptides to add to the begining or end')
    parser.add_argument('-m',
                        help='A csv file containing the name/indentifer for the sequence which does NOT have to be unique')
    parser.add_argument('-samp',
                        help='sample name')
    parser.add_argument('-HLA',
                        help='a list of the HLA alleles in the format: HLA-A*02:01,HLA-A*24:02,HLA-B*07:02,HLA-B*35:02,HLA-C*04:01,HLA-C*07:02')
    parser.add_argument('-WD',
                        help='The directory in which you would like to run pVACbind')
    
    return(parser.parse_args())

def generate_modifed_peptides(n, name, base_sequence):
    characters = ['K', 'R']
    
    # Generate all possible combinations up to length n
    all_combinations = chain.from_iterable(product(characters, repeat=i) for i in range(1, n))
    possible_modifications = set(''.join(combination) for combination in all_combinations)
    
    peptide_table = []

    for modification in possible_modifications:

        Nterm_dict = {}
        Cterm_dict = {}

        Nterm_sequence_name = name + "." + "n-term" + "-" + modification
        Cterm_sequence_name = name + "." + "c-term" + "-" + modification

        Nterm_sequence = modification + base_sequence
        Cterm_sequence = base_sequence + modification 

        Nterm_parsed_sequence = modification + '|' + base_sequence
        Cterm_parsed_sequence = base_sequence + '|' + modification

        Nterm_dict.update({'sequence_name': Nterm_sequence_name, 'sequence': Nterm_sequence, 'parsed_sequence': Nterm_parsed_sequence})
        Cterm_dict.update({'sequence_name': Cterm_sequence_name, 'sequence': Cterm_sequence, 'parsed_sequence': Cterm_parsed_sequence})

        peptide_table.append(Nterm_dict)
        peptide_table.append(Cterm_dict)

    return(peptide_table)


def assign_unique_numbers(df, column_name):
    counts = df[column_name].value_counts()
    duplicated_values = counts[counts > 1].index

    for value in duplicated_values:
        indices = df.index[df[column_name] == value]
        for i, index in enumerate(indices, start=1):
            df.at[index, column_name] = f"{value}.{i}"

    return df

def main():
    args = parse_arguments()

    max_length = args.n

    peptides = pd.read_csv(args.m, names=["Name", "Sequence"], header=None)
    peptides = peptides[1:]
    peptides = assign_unique_numbers(peptides, "Name")

    max_length = 3
    list = []

    for index, row in peptides.iterrows():
        sequences_list = []
        
        name = row['Name']
        base_sequence = row['Sequence']
        sequences_list  = generate_modifed_peptides(max_length, name, base_sequence)

        list = list + sequences_list

    peptide_table = pd.DataFrame(list) 

    peptide_table.to_csv('peptide_table.tsv', sep="\t", index=False, header=None)


    # Perform checks
    num_entries = len(peptide_table)

    # Check if all entries in the sequence_name column are unique
    are_all_unique = not peptide_table['sequence_name'].duplicated().any()

    if are_all_unique:
        print("All entries in the sequence_name column are unique.")
    else:
        print("There are duplicate entries in the sequence_name column.")
        sys.exit(1)
    
    with open('modified_peptides.fa', 'w') as fasta_file:
        for index, row in peptide_table.iterrows():
            sequence_name = ">" + row['sequence_name']
            sequence = row['parsed_sequence']
            fasta_file.write(f"{sequence_name}\n{sequence}\n")


    # Create dirs for processing the N-term and C-term sequences separately
    os.makedirs("n-term", exist_ok=True)
    os.makedirs("c-term", exist_ok=True)

    with open('modified_peptides.fa', 'r') as input_file, \
            open('n-term/modified_peptides_n-term.fa', 'w') as n_term_file, \
            open('c-term/modified_peptides_c-term.fa', 'w') as c_term_file:
        
        for line in input_file:
            line = line.strip()
            
            # Check for the pattern 'n-term' in the line
            if 'n-term' in line:
                n_term_file.write(line + '\n')
                
                # Ensure there are at least two more lines in the file
                try:
                    next_line = next(input_file).strip()
                    n_term_file.write(next_line + '\n')
                except StopIteration:
                    # Handle the case where there are not enough lines left
                    pass
            
            # Check for the pattern 'c-term' in the line
            if 'c-term' in line:
                c_term_file.write(line + '\n')
                
                # Ensure there are at least two more lines in the file
                try:
                    next_line = next(input_file).strip()
                    c_term_file.write(next_line + '\n')
                except StopIteration:
                    # Handle the case where there are not enough lines left
                    pass

    def get_line_count(file_path):
        with open(file_path, 'r') as file:
            return sum(1 for line in file)

    modified_peptides_count = get_line_count('modified_peptides.fa')
    n_term_count = get_line_count('n-term/modified_peptides_n-term.fa')
    c_term_count = get_line_count('c-term/modified_peptides_c-term.fa')

    if n_term_count+ c_term_count == modified_peptides_count:
        print("Sucessfully split fasta into n-term and c-term")
    else:
        print("N-term and c-term fasta lines to not equal modified fasta lines")


    def  create_subpeptide_fastas_n_term(input_dir, results_dir, infile_path):
        # Define the lengths to iterate over
        lengths = [8, 9, 10, 11]

        # Iterate over each test length
        for LENGTH in lengths:
            # Create input files for each test length
            print(f"Creating input fasta to test peptides of length: {LENGTH}")
            
            # Set the file paths
            LENGTH_FASTA = os.path.join(input_dir , f"{LENGTH}-mer-test.fa")
            LENGTH_RESULT_DIR = os.path.join(results_dir, f"{LENGTH}-mer-test")
            
            # Create the result directory if it doesn't exist
            os.makedirs(LENGTH_RESULT_DIR, exist_ok=True)
            
            with open(infile_path, 'r') as infile, open(LENGTH_FASTA, 'w') as length_fasta:
                for line in infile:
                    line = line.strip()
                    if line.startswith('>'):
                        length_fasta.write(f"{line}\n")
                    elif match := re.match(r'(\w+)\|(\w+)', line):
                        before, after = match.groups()
                        sub = after[:LENGTH - 1]
                        length_fasta.write(f"{before}{sub}\n")
                infile.close()

    def  create_subpeptide_fastas_c_term(input_dir, results_dir, infile_path):
        # Define the lengths to iterate over
        lengths = [8, 9, 10, 11]

        for LENGTH in lengths:
            # Create input files for each test length so that each test sequence will contain at least one modified base
            print(f"Creating input fasta to test peptides of length: {LENGTH}")

            LENGTH_FASTA = os.path.join(input_dir, f"{LENGTH}-mer-test.fa")
            LENGTH_RESULT_DIR = os.path.join(results_dir, f"{LENGTH}-mer-test")

            os.makedirs(LENGTH_RESULT_DIR, exist_ok=True)

            with open(infile_path, 'r') as infile, open(LENGTH_FASTA, 'w') as length_fasta:
                for line in infile:
                    line = line.strip()
                    if line.startswith('>'):
                        length_fasta.write(f"{line}\n")
                    elif match := re.match(r'(\w+)\|(\w+)', line):
                        before, after = match.groups()
                        sub = before[-(LENGTH - 1):]
                        length_fasta.write(f"{sub}{after}\n")
                infile.close()

    # Set up sub-peptide fasta sequences for each target class I prediction length
    working_dir = args.WD

    input_dir = os.path.join(working_dir + "/n-term/pvacbind_inputs")
    results_dir = os.path.join(working_dir + "/n-term/pvacbind_results")
    infile_path = os.path.join(working_dir + "/n-term/modified_peptides_n-term.fa")

    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    print("Creating sub-pepetide fastas for n-term")
    create_subpeptide_fastas_n_term(input_dir, results_dir, infile_path)

    input_dir = os.path.join(working_dir + "/c-term/pvacbind_inputs")
    results_dir = os.path.join(working_dir + "/c-term/pvacbind_results")
    infile_path = os.path.join(working_dir + "/c-term/modified_peptides_c-term.fa")

    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    print("Creating sub-pepetide fastas for c-term")
    create_subpeptide_fastas_c_term(input_dir, results_dir, infile_path)

    # Run pVACtools







if __name__ == "__main__":
    main()
