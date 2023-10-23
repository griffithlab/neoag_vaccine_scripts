import argparse
from itertools import product, chain
import pandas as pd


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
    return(parser.parse_args())

def generate_modifed_peptides(n, name, base_sequence):
    characters = ['K', 'R']
    
    # Generate all possible combinations up to length n
    all_combinations = chain.from_iterable(product(characters, repeat=i) for i in range(1, n + 1))
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

    df = pd.DataFrame(list) 

    df.to_csv('peptide_table.tsv', sep="\t", index=False, header=None)

if __name__ == "__main__":
    main()