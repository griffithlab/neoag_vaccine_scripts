import argparse
import csv
import pandas as pd
import sys


'''
Read in a excel file with CSBio proposed sequences which have addition amino acids to improve solubility
It is assumed that the additional AA are spearated by a "|" from the rest of the wildtype sequence

Create a files called petide_table.tsv which has three columns (a) unique peptide name, (b) peptide sequence, (c) parsable peptide sequence
(c) is from excel file
(b) just remove "|"
(a) 
    - no gene name, download final 51 mer spreadsheet
        - use canidate neoantigen and restircting HLA allele to join gene name to sheet (ask CSBio to match Canidate Neoantigen ID of SAMPLE.MT.POS.GENE)
    - assemble unique pepetide name
        - gene
        - parse proposed sequence (it would be easier if CSBio gave all possible WT seeunces and we could test any number of additional AA)
            - start counter for gene
            - if you find a | before the fourth position used n-term
            - count how many K/Rs there are
            - 
'''

# ---- PARSE ARGUMENTS -------------------------------------------------------
# Parses command line arguments
# Enables user help
def parse_arguments():
    # Parse command line arugments
    parser = argparse.ArgumentParser(description='Create the file needed for the neoantigen manuel review')

    parser.add_argument('-p',
                        help='The path to the ITB Reviewed Candidates')
    parser.add_argument('-51mer',
                        help='The path to annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv from the generate_protein_fasta script')
   
    # The name of the final results folder 
    parser.add_argument('-f', "--fin_results", help="Name of the final results folder in gcp immuno")

    return(parser.parse_args())

def main():
    proposed = pd.read_excel('/Users/evelynschmidt/jlf/JLF-100-047/proposed.xlsx')
    peptides = pd.read_excel('/Users/evelynschmidt/jlf/JLF-100-047/JLF-100-047_Peptides_51mer.final.xlsx')

    # separate ID and gene name in peptides sheet
    peptides[['Candidate Neoantigen', 'Gene']] = peptides['CANDIDATE NEOANTIGEN'].str.rsplit('.', n=1, expand=True)
    # remove the first period to match format of the proposed sheet
    peptides['Candidate Neoantigen'] = peptides['Candidate Neoantigen'].str.replace('.', '', 1)

    # join peptides with proposed on Canidate Neoantigen persevering only the gene column from peptides
    proposed_full = pd.merge(proposed, peptides[['Candidate Neoantigen', 'Gene']], on='Candidate Neoantigen', how='left')

    # Creating a new dataframe
    print(proposed_full.columns)
    peptide_table = proposed_full[['CSBio Proposed Sequences', 'Gene']].copy()
    # get rid of amy unedited sequences
    peptide_table = peptide_table[peptide_table['CSBio Proposed Sequences'].str.contains('\|')]

    # create a column that has the peptide sequence but with no "|"
    peptide_table['peptide sequence'] =  peptide_table['CSBio Proposed Sequences'].str.replace('|', '')

    # create an ID column
    peptide_table = peptide_table.rename(columns={"CSBio Proposed Sequences":"parsable peptide sequence"})
    peptide_table['ID'] = " "


    # split the parsable peptide sequence and create an ID that is gene and WT sequence
    gene = ""
    inter = 0
    for index, row in peptide_table.iterrows():

        if row['Gene'] == gene:
            inter += 1
        else:
            inter = 1 
            gene = row['Gene']
        
    
        id = row['Gene'] + '.'

        # break the sequence on |
        result_array = row["parsable peptide sequence"].split('|')

        if len(result_array) < 3: # there is only an addition at one end of the sequence

            if len(result_array[0]) > 3: # the addtion is at the end
                id += str(inter)
                id += '.c-term-'
                id += result_array[1]
            else: # the addtion is at the beginning
                id += str(inter)
                id += '.n-term-'
                id += result_array[0]
        else: # there is  an addition at both ends
            id += str(inter)
            id += '.n-term-'
            id += result_array[0]
            id += '.c-term-'
            id += result_array[2]

        peptide_table.at[index, 'ID'] = id
        
    

    final_peptide_table = peptide_table[['ID', 'peptide sequence', 'parsable peptide sequence']].copy()
    print(final_peptide_table)

    final_peptide_table.to_csv('peptide_table.tsv', sep='\t', index=False, header=False)



if __name__ == "__main__":
    main()