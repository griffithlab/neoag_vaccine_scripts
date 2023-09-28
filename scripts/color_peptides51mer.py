import numpy as np
import pandas as pd
import re
from bs4 import BeautifulSoup
import argparse


class AminoAcid:

    def __init__(self, nucleotide, bold, color, underline, large, position, open_tag, close_tag):
        self.nucleotide = nucleotide
        self.bold = bold
        self.color = color
        self.underline = underline
        self.large = large
        self.position = position 
        self.open_tag = open_tag
        self.close_tag = close_tag

# ---- PARSE ARGUMENTS -------------------------------------------------------
# Parses command line arguments
# Enables user help
def parse_arguments():
    # Parse command line arugments
    parser = argparse.ArgumentParser(description='Color the 51mer peptide')

    parser.add_argument('-p',
                        help='The path to the Peptides 51 mer', required=True)
    parser.add_argument('-classI',
                        help='The path to the classI all_epitopes.aggregated.tsv used in pVACseq', required=True)
    parser.add_argument('-classII',
                        help='The path to the classII all_epitopes.aggregated.tsv used in pVACseq', required=True)

    parser.add_argument('-o', help="Output location", required=True)

    return(parser.parse_args())

    # Function to rearrange string so that G518D looks like 518G/D
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
        

def annotate_every_nucleotide(sequence, classI_peptide, classII_peptide):

    peptide_sequence = []

    # Make the sequence a list of AminoAcid objects
    for i in range(len(sequence)):
        new_AA = AminoAcid(sequence[i], False, False, False, False, -1, False, False)

        if sequence[i] == 'C':
            new_AA.large = True
        
        peptide_sequence.append(new_AA)

    # CLASS I
    positions = []
    # Get the positions in the peptide_sequence where the classI is located
    for i in range(len(peptide_sequence)):
        for j in range(len(classI_peptide)):
            if peptide_sequence[i].nucleotide == classI_peptide[j]:
                positions.append(i)
                i+=1
            else:
                break

        if len(positions) == len(classI_peptide):
            break
        else:
            positions = []
            
    # set those positions to red
    j = 0
    for i in range(len(peptide_sequence)):
        if j < len(positions) and i == positions[j]:
            peptide_sequence[i].color = True
            j+=1

    # CLASS II
    positions = []
    for i in range(len(peptide_sequence)):
        for j in range(len(classII_peptide)):
            if peptide_sequence[i].nucleotide == classII_peptide[j]:
                positions.append(i)
                i+=1
            else:
                break

        if len(positions) == len(classII_peptide):
            break
        else:
            positions = []

    j = 0
    for i in range(len(peptide_sequence)):
        if j < len(positions) and i == positions[j]:
            peptide_sequence[i].bold = True
            j+=1
    

    return(peptide_sequence)

def set_underline(peptide_sequence, mutant_peptide_pos, row_ID):

    frameshift = False
    classI_position = 0

    if '-' in mutant_peptide_pos:
        positions = mutant_peptide_pos.split("-")
        start_position = int(positions[0])
        end_position = int(positions[1])

        frameshift = True
        

    else:
        mutant_peptide_pos = int(mutant_peptide_pos)

        

    if frameshift:

        continue_underline = False
        
        for i in range(len(peptide_sequence)):

            if peptide_sequence[i].color:
                classI_position += 1
            else:
                 classI_position = 0
                 continue_underline = False

            if classI_position == start_position:
                peptide_sequence[i].underline = True
                continue_underline = True
            elif continue_underline:
                peptide_sequence[i].underline = True
            elif classI_position == end_position:
                peptide_sequence[i].underline = True
                continue_underline = False
            i+=1
    else:
        for i in range(len(peptide_sequence)):

            if peptide_sequence[i].color:
                classI_position += 1
            else:
                 classI_position = 0

            if classI_position == int(mutant_peptide_pos):
                peptide_sequence[i].underline = True
            i+=1
        
def set_span_tags(peptide_sequence):

    currently_bold = False
    currently_red = False
    currently_underlined = False
    currently_large = False
    inside_span = False

    for nucleotide in peptide_sequence:

        if currently_bold != nucleotide.bold or currently_red != nucleotide.color or currently_underlined != nucleotide.underline or currently_large != nucleotide.large:
            
            nucleotide.open_tag = True

            if inside_span:
                nucleotide.close_tag = True # only if its isnide a span tag
            else:
                nucleotide.close_tag = False


            currently_bold = nucleotide.bold 
            currently_red = nucleotide.color
            currently_underlined = nucleotide.underline
            currently_large = nucleotide.large

            inside_span = True
        
    return(peptide_sequence)

def create_stylized_sequence(peptide_sequence):

    new_string = ''

    for nucleotide in peptide_sequence:

        if nucleotide.open_tag or nucleotide.close_tag:
            if nucleotide.close_tag:
                new_string += '</span>'
                

            if nucleotide.open_tag:

                if nucleotide.large: # we are assuming that a cystine is never in the classI and classIi
                    new_string += '<span style="font-size:105%">'
                    new_string += nucleotide.nucleotide

                if nucleotide.bold and nucleotide.color and nucleotide.underline:
                    new_string += '<span style="font-weight:bold;color:#ff0000;text-decoration:underline;">'
                    new_string += nucleotide.nucleotide
                elif nucleotide.bold and not nucleotide.color and not nucleotide.underline:
                    new_string += '<span style="font-weight:bold;">'
                    new_string += nucleotide.nucleotide
                elif not nucleotide.bold and nucleotide.color and not nucleotide.underline:
                    new_string += '<span style="color:#ff0000;">'
                    new_string += nucleotide.nucleotide
                elif not nucleotide.bold and not nucleotide.color and nucleotide.underline:
                    new_string += '<span style="text-decoration:underline;">'
                    new_string += nucleotide.nucleotide
                elif nucleotide.bold and nucleotide.color and not nucleotide.underline:
                    new_string += '<span style="font-weight:bold;color:#ff0000;">'
                    new_string += nucleotide.nucleotide
                elif not nucleotide.bold and nucleotide.color and nucleotide.underline:
                    new_string += '<span style="color:#ff0000;text-decoration:underline;">'
                    new_string += nucleotide.nucleotide
                elif nucleotide.bold and not nucleotide.color and nucleotide.underline:
                    new_string += '<span style="font-weight:bold;text-decoration:underline;">'
                    new_string += nucleotide.nucleotide

            if not nucleotide.large and not nucleotide.bold and not nucleotide.color and not nucleotide.underline:
                new_string += nucleotide.nucleotide
        else:
            new_string += nucleotide.nucleotide
    return(new_string)   

def main():
    args = parse_arguments()

    # read in classI and class II
    #peptides_51mer = pd.read_excel("/Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/10146-0021_Peptides_51-mer.xlsx")
    #classI = pd.read_csv("/Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/classI.TWJF-10146-0021-Tumor_Lysate.all_epitopes.aggregated.tsv", sep="\t")
    #classII = pd.read_csv("/Volumes/mgriffit/Active/griffithlab/gc2596/e.schmidt/neoag_vaccine_scripts/scripts/data_files/classII.TWJF-10146-0021-Tumor_Lysate.all_epitopes.aggregated.tsv", sep="\t")

    peptides_51mer = pd.read_excel(args.p)
    classI = pd.read_csv(args.classI, sep="\t")
    classII = pd.read_csv(args.classII, sep="\t")

    # Create a universal ID by editing the peptide 51mer ID
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

    merged_peptide_51mer = pd.merge(merged_peptide_51mer, classI[['ID', 'Best Peptide', 'Pos']], on='ID', how='left')

    merged_peptide_51mer.rename(columns = {"Best Peptide":"Best Peptide Class I"}, inplace=True)

    # convert peptide 51mer to HTML
    peptides_51mer_html = peptides_51mer.to_html(index=False) # convert to html

    # Creating a BeautifulSoup object and specifying the parser
    peptides_51mer_soup = BeautifulSoup(peptides_51mer_html, 'html.parser')


    for index, row in peptides_51mer.iterrows():

        search_string = row['full ID']

        #classII_sequence 
        classII_peptide = merged_peptide_51mer.loc[merged_peptide_51mer['full ID'] == search_string, 'Best Peptide Class II'].values[0]
        #classI_sequence 
        classI_peptide = merged_peptide_51mer.loc[merged_peptide_51mer['full ID'] == search_string, 'Best Peptide Class I'].values[0]
        
        
        # mutant pepetide position ---
        mutant_peptide_pos = str(merged_peptide_51mer.loc[merged_peptide_51mer['full ID'] == search_string, 'Pos'].values[0])

        # Find the tag containing the search string
        tag_with_search_string = peptides_51mer_soup.find('td', string=search_string)

        if tag_with_search_string and isinstance(classII_peptide, str):

            # Find the parent <tr> tag of the tag containing the search string
            parent_tr = tag_with_search_string.find_parent('tr')    
            # Find the next two <td> tags
            next_td_tags = parent_tr.findChildren('td', limit=3)
            
            sequence = next_td_tags[2].get_text()

            # make sequence the list of objects
            peptide_sequence = annotate_every_nucleotide(sequence, classI_peptide, classII_peptide)

            # actaully lets break class I and classII into two steps and handle the mutated nucleotide in class I function
            # it should be basically like at that position in the class I set 
            
            set_underline(peptide_sequence, mutant_peptide_pos, row['full ID'])

            set_span_tags(peptide_sequence) # pass by reference
            
            new_string = create_stylized_sequence(peptide_sequence)

            next_td_tags[2].string = new_string

            modified_html = peptides_51mer_soup.prettify(formatter=None)

        else:
            print("\nNOT FOUND: ", search_string)
            print("Mutant Peptide Position: ", mutant_peptide_pos)
            print("ClassI: ", classI_peptide)
            print("ClassII: ", classII_peptide, "\n")

    with open(args.o, "w", encoding = 'utf-8') as file:
        file.write(modified_html)



if __name__ == "__main__":
    main()