import argparse
import csv
import pandas as pd
import sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re

'''
Write a script to create the files for the Case Final Reports
- Sample Peptides 51-mer
- SAMPLE.Annotated.Neoantigen_Candidates.xlsx
'''

# ---- PARSE ARGUMENTS -------------------------------------------------------
# Parses command line arguments
# Enables user help
def parse_arguments():
    # Parse command line arugments
    parser = argparse.ArgumentParser(description='Create the file needed for the neoantigen manuel review')

    parser.add_argument('-a',
                        help='The path to the ITB Reviewed Candidates')
    parser.add_argument('-c',
                        help='The path to annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv from the generate_protein_fasta script')
    parser.add_argument('-samp',
                        help='The name of the sample')
    parser.add_argument('-WB',
                        help='the path to the gcp_immuno folder of the trial you wish to tun script on, defined as WORKING_BASE in envs.txt')

    # The name of the final results folder 
    parser.add_argument('-f', "--fin_results", help="Name of the final results folder in gcp immuno")

    return(parser.parse_args())

# Fucnction to break the pepetides ID on the . to extract gene and AA information
def extract_info(value):
    parts = value.split('.')
    result = '.'.join([parts[2], parts[3], parts[4]])
    return result

# Function to rearrange string so that G518D looks like 518G/D
def rearrange_string(s):
    match = re.match(r'([A-Za-z]+)([\d-]+)([A-Za-z]*)', s)
    if match:
        letters_before = match.group(1)
        numbers = match.group(2)
        letters_after = match.group(3)
                
        return f"{numbers}{letters_before}/{letters_after}"
    else:
        return s
    
# Function to calculate molecular weight---------------------------------------
def calculate_molecular_weight(peptide):
    analyzed_seq = ProteinAnalysis(peptide)
    return analyzed_seq.molecular_weight()

# Function to make id column unique -------------------------------------------
def make_column_unique(df, column_name):
    seen_values = set()
    new_values = []

    for value in df[column_name]:
        if value in seen_values:
            suffix = 1
            while f"{value}.{suffix}" in seen_values:
                suffix += 1
            unique_value = f"{value}.{suffix}"
        else:
            unique_value = value

        seen_values.add(unique_value)
        new_values.append(unique_value)

    df[column_name] = new_values
    return df


def main():

    args = parse_arguments()
    
    reviewed_candidates = pd.read_excel(args.a)
 

    reviewed_candidates.columns = reviewed_candidates.iloc[0]
    reviewed_candidates = reviewed_candidates[1:] # there is a extra row before the col name row
    reviewed_candidates = reviewed_candidates.reset_index(drop=True) # Reset the index of the dataframe
    
    reviewed_candidates = reviewed_candidates[reviewed_candidates.Evaluation != "Pending"]
    reviewed_candidates = reviewed_candidates[reviewed_candidates.Evaluation != "Reject"]

    reviewed_candidates = reviewed_candidates.rename(columns={'Comments':'pVAC Review Comments'})
    reviewed_candidates["Variant Called by CLE Pipeline"] = " "
    reviewed_candidates["IGV Review Comments"] = " "


    # create sorting ID that is gene and transcript to sort in the same order as peptide
    reviewed_candidates['sorting id'] = reviewed_candidates['Gene']  + '.' + reviewed_candidates['Best Transcript']
    # make sure the sorting id column is unique
    reviewed_candidates = make_column_unique(reviewed_candidates, 'sorting id')

    peptides = pd.read_csv(args.c, sep="\t")
    peptides =  peptides.drop(['cterm_7mer_gravy_score', 'cysteine_count', 'n_terminal_asparagine', 'asparagine_proline_bond_count', 
                                 'difficult_n_terminal_residue', 'c_terminal_cysteine', 'c_terminal_proline', 'max_7mer_gravy_score'], axis=1)
    peptides["RESTRICTING HLA ALLELE"] = " "

    peptides["CANDIDATE NEOANTIGEN AMINO ACID SEQUENCE MW (CLIENT)"] = peptides["peptide_sequence"].apply(calculate_molecular_weight)

    peptides = peptides.rename(columns={"id":"ID", "peptide_sequence":"CANDIDATE NEOANTIGEN AMINO ACID SEQUENCE WITH FLANKING RESIDUES"})
    peptides["Comments"] = " "
    peptides["CANDIDATE NEOANTIGEN"] = peptides["ID"].apply(lambda x: '.'.join(x.split('.')[:3]))
    peptides["CANDIDATE NEOANTIGEN"] = args.samp + "." + peptides["CANDIDATE NEOANTIGEN"]

    peptides = peptides[["ID", "CANDIDATE NEOANTIGEN", "CANDIDATE NEOANTIGEN AMINO ACID SEQUENCE WITH FLANKING RESIDUES", 
                           "RESTRICTING HLA ALLELE", "CANDIDATE NEOANTIGEN AMINO ACID SEQUENCE MW (CLIENT)", "Comments"]]
    

    # creating a ID to sort reviewed canidates by the order of the 51mer
    peptides['sorting id'] = peptides['ID'].apply(extract_info)
    # make sure every sorting id is unique
    peptides = make_column_unique(peptides, 'sorting id')

    reviewed_candidates = reviewed_candidates.set_index('sorting id')
    reviewed_candidates = reviewed_candidates.reindex(index=peptides['sorting id'])
    reviewed_candidates = reviewed_candidates.reset_index()

    reviewed_candidates = reviewed_candidates.drop(columns=['sorting id'])
    peptides = peptides.drop(columns=['sorting id'])


    if args.WB:
        Peptide_file_name = args.WB +  '/../manual_review/' + args.samp + "_Peptides_51-mer.xlsx"
    else:
        Peptide_file_name =  args.samp + "_Peptides_51-mer.xlsx"

    peptides.to_excel(Peptide_file_name, index=False)

    if args.WB:
        Neoantigen_Canidates_file_name = args.WB +  '/../manual_review/' + args.samp + ".Annotated.Neoantigen_Candidates.xlsx"
    else:
        Neoantigen_Canidates_file_name =  args.samp + ".Annotated.Neoantigen_Candidates.xlsx"

    reviewed_candidates.to_excel(Neoantigen_Canidates_file_name, index=False)


if __name__ == "__main__":
    main()
