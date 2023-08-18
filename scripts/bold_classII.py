import argparse
import csv
import pandas as pd
import sys
import re
from bs4 import BeautifulSoup

# read in class I and class II, pepetides 51mer
# create a key that key that matches the petides 51mer ID
# loop through all rows in peptides 51 mer
# try to bold all class II peptides

'''
Use: 
python3 scripts/bold_classII.py -p /Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/10146-0021 Peptides 51-mer.xlsx -classI /Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/classI.TWJF-10146-0021-Tumor_Lysate.all_epitopes.aggregated.tsv -classII /Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/classII.TWJF-10146-0021-Tumor_Lysate.all_epitopes.aggregated.tsv -o /Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/test_colored_peptide.html
'''
def parse_arguments():
    # Parse command line arugments
    parser = argparse.ArgumentParser(description='Bold the class II pepetides')

    parser.add_argument('-p',
                        help='The path to the Peptides 51 mer', required=True)
    parser.add_argument('-classI',
                        help='The path to the classI all_epitopes.aggregated.tsv used in pVACseq', required=True)
    parser.add_argument('-classII',
                        help='The path to the classII all_epitopes.aggregated.tsv used in pVACseq', required=True)

    parser.add_argument('-o', help="Output location", required=True)

    return(parser.parse_args())

def rearrange_string(s):
    match = re.match(r'([A-Za-z]+)([\d-]+)([A-Za-z]*)', s)
    if match:
        letters_before = match.group(1)
        numbers = match.group(2)
        letters_after = match.group(3)
        
        #return f"{numbers}{letters_before}/{letters_after}"
        # Just use the postion for the key to avoid FS problem
        return f"{numbers}"
    else:
        return s

def insert_around_substring(original_string, target_substring, insert_before, insert_after):

    new_string = original_string
    index = original_string.find(target_substring)
    
    if index != -1:
        new_string = (
            original_string[:index] +
            insert_before +
            target_substring +
            insert_after +
            original_string[index + len(target_substring):]
        )
        
    return new_string


args = parse_arguments()

peptides_51mer = pd.read_excel(args.p)

classI = pd.read_csv(args.classI, sep="\t")

classII = pd.read_csv(args.classII, sep="\t")

# create ID to be shared by all dataframes
peptides_51mer.rename(columns={'ID': 'full ID'}, inplace=True)
peptides_51mer['ID'] = peptides_51mer['full ID']

peptides_51mer['ID'] = peptides_51mer['ID'].apply(lambda x: '.'.join(x.split('.')[1:]))  # Removing before first period, periods will be removed 
peptides_51mer['ID'] = peptides_51mer['ID'].apply(lambda x: '.'.join(x.split('.')[1:]))  # Removing before second period
peptides_51mer['ID'] = peptides_51mer['ID'].apply(lambda x: '.'.join(x.split('.')[:3]) + '.' + '.'.join(x.split('.')[4:]))

for index, row in peptides_51mer.iterrows():
    for i, char in enumerate(row['ID'][::-1]):
        if char.isdigit():
            peptides_51mer.at[index, 'ID'] = row['ID'][:-i]
            break
    else:
        result = row['ID']


# create a key that is gene, transcript, AA change for CLASSI
classII['modified AA Change'] = classII['AA Change'] 
# Apply the function to the 'Value' column
classII['modified AA Change'] = classII['modified AA Change'].apply(rearrange_string)
classII['ID'] = classII['Gene'] + '.' + classII['Best Transcript'] + '.' + classII['modified AA Change'] 


# create a key that is gene, transcript, AA change for CLASSI
classI['modified AA Change'] = classI['AA Change'] 
# Apply the function to the 'Value' column
classI['modified AA Change'] = classI['modified AA Change'].apply(rearrange_string)
classI['ID'] = classI['Gene'] + '.' + classI['Best Transcript'] + '.' + classI['modified AA Change'] 

# Merge the sequences from classI and classII with peptide 51mer
merged_peptide_51mer = pd.merge(peptides_51mer, classII[['ID', 'Best Peptide']], on='ID', how='left')
merged_peptide_51mer.rename(columns = {"Best Peptide":"Best Peptide Class II"}, inplace=True)
merged_peptide_51mer = pd.merge(merged_peptide_51mer, classI[['ID', 'Best Peptide']], on='ID', how='left')
merged_peptide_51mer.rename(columns = {"Best Peptide":"Best Peptide Class I"}, inplace=True)


# convert peptide 51mer to HTML
peptides_51mer_html = peptides_51mer.to_html(index=False) # convert to html

# Creating a BeautifulSoup object and specifying the parser
peptides_51mer_soup = BeautifulSoup(peptides_51mer_html, 'html.parser')


# Loop through all peptide 51mer IDs
for index, row in peptides_51mer.iterrows():

    search_string = row['full ID']

    #classII_sequence 
    classII_peptide = merged_peptide_51mer.loc[merged_peptide_51mer['full ID'] == search_string, 'Best Peptide Class II'].values[0]

    # Find the tag containing the search string
    tag_with_search_string = peptides_51mer_soup.find('td', string=search_string)

    if tag_with_search_string and isinstance(classII_peptide, str):

        # Find the parent <tr> tag of the tag containing the search string
        parent_tr = tag_with_search_string.find_parent('tr')
    
        # Find the next two <td> tags
        next_td_tags = parent_tr.findChildren('td', limit=3)

        sequence = next_td_tags[2].get_text()

        new_string = insert_around_substring(sequence, classII_peptide, '<span style="font-weight:bold;">', '</span>')

        next_td_tags[2].string = new_string

        modified_html = peptides_51mer_soup.prettify(formatter=None)

    else:
        print("Search string not found.")


with open(args.o, "w", encoding = 'utf-8') as file:
    file.write(modified_html)