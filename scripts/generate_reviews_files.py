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

Maybe the Sample Genomics Review Report with everything highlighted in yellow


Use:
python3 generate_reviews_files.py -a /Volumes/gillandersw/Active/Project_0001_Clinical_Trials/CTEP/analysis/TWJF-10146-MO011-0021/itb-review-files/10146-0021.Annotated.Neoantigen_Candidates.xlsx -c /Volumes/gillandersw/Active/Project_0001_Clinical_Trials/CTEP/analysis/TWJF-10146-MO011-0021/generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv -samp 10146-0021
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
    
# Function to calculate molecular weight
def calculate_molecular_weight(peptide):
    analyzed_seq = ProteinAnalysis(peptide)
    return analyzed_seq.molecular_weight()

def main():

    # 1. ITB reivew
    # 2. Generate protein Fasta 


    args = parse_arguments()
    
    reviewed_canidates = pd.read_excel(args.a)
 

    reviewed_canidates.columns = reviewed_canidates.iloc[0]
    reviewed_canidates = reviewed_canidates[1:] # there is a extra row before the col name row
    reviewed_canidates = reviewed_canidates.reset_index(drop=True) # Reset the index of the dataframe
    
    reviewed_canidates = reviewed_canidates[reviewed_canidates.Evaluation != "Pending"]
    reviewed_canidates = reviewed_canidates[reviewed_canidates.Evaluation != "Reject"]

    reviewed_canidates = reviewed_canidates.rename(columns={'Comments':'pVAC Review Comments'})
    reviewed_canidates["Variant Called by CLE Pipeline"] = " "
    reviewed_canidates["IGV Review Comments"] = " "


    # create sorting ID that is gene and transcript to sort in the same order as peptide
    reviewed_canidates['sorting id'] = reviewed_canidates['Gene']  + '.' + reviewed_canidates['Best Transcript']


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

    reviewed_canidates = reviewed_canidates.set_index('sorting id')
    reviewed_canidates = reviewed_canidates.reindex(index=peptides['sorting id'])
    reviewed_canidates = reviewed_canidates.reset_index()

    reviewed_canidates = reviewed_canidates.drop(columns=['sorting id'])
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

    reviewed_canidates.to_excel(Neoantigen_Canidates_file_name, index=False)


if __name__ == "__main__":
    main()
