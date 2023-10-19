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
    parser.add_argument('-WB',
                        help='the path to the gcp_immuno folder of the trial you wish to tun script on, defined as WORKING_BASE in envs.txt')
    parser.add_argument('-samp', help='Name of the sample')


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

                new_string += '<span style="'
                if nucleotide.bold:
                    new_string += 'font-weight:bold;'
                if nucleotide.color:
                    new_string += 'color:#ff0000;'
                if nucleotide.underline:
                    new_string += 'text-decoration:underline;'
                if nucleotide.large:
                     new_string += 'font-size:105%;'
                new_string += '">'
                new_string += nucleotide.nucleotide

            if not nucleotide.large and not nucleotide.bold and not nucleotide.color and not nucleotide.underline:
                new_string += nucleotide.nucleotide
        else:
            new_string += nucleotide.nucleotide
    return(new_string)   

def main():
    args = parse_arguments()

    # read in classI and class II
    peptides_51mer = pd.read_excel(args.p)
    classI = pd.read_csv(args.classI, sep="\t")
    classII = pd.read_csv(args.classII, sep="\t")

    # Create a universal ID by editing the peptide 51mer ID
    peptides_51mer.rename(columns={'ID': 'full ID'}, inplace=True)
    peptides_51mer['51mer ID'] = peptides_51mer['full ID']

    peptides_51mer['51mer ID'] = peptides_51mer['51mer ID'].apply(lambda x: '.'.join(x.split('.')[1:]))  # Removing before first period, periods will be removed 
    
    peptides_51mer['51mer ID'] = peptides_51mer['51mer ID'].apply(lambda x: '.'.join(x.split('.')[1:]))  # Removing before second period
    peptides_51mer['51mer ID'] = peptides_51mer['51mer ID'].apply(lambda x: '.'.join(x.split('.')[:3]) + '.' + '.'.join(x.split('.')[4:]))
    

    for index, row in peptides_51mer.iterrows():
        for i, char in enumerate(row['51mer ID'][::-1]):
            if char.isdigit():
                peptides_51mer.at[index, '51mer ID'] = row['51mer ID'][:-i]
                break
        else:
            result = row['51mer ID']

    # create a dataframe that contains the classI and classII pepetide sequence
    classI.rename(columns = {"Best Peptide":"Best Peptide Class I"}, inplace=True)
    classII.rename(columns = {"Best Peptide":"Best Peptide Class II"}, inplace=True)

    # create a key that is gene, transcript, AA change for ClassI to join to the peptides order form
    classI['modified AA Change'] = classI['AA Change'] 
    classI['modified AA Change'] = classI['modified AA Change'].apply(rearrange_string)
    classI['51mer ID'] = classI['Gene'] + '.' + classI['Best Transcript'] + '.' + classI['modified AA Change'] 

    class_sequences = pd.merge(classI[['ID', 'Best Peptide Class I', '51mer ID', 'Pos']], classII[['ID', 'Best Peptide Class II']], on='ID', how='left')

    # Create a dataframe that has the classI and classII sequence
    merged_peptide_51mer = pd.merge(peptides_51mer, class_sequences, on='51mer ID', how='left')

    # convert peptide 51mer to HTML
    peptides_51mer_html = peptides_51mer.to_html(index=False) # convert to html

    # Creating a BeautifulSoup object and specifying the parser
    peptides_51mer_soup = BeautifulSoup(peptides_51mer_html, 'html.parser')


    for index, row in peptides_51mer.iterrows():

        search_string = row['51mer ID']

        #classII_sequence 
        classII_peptide = merged_peptide_51mer.loc[merged_peptide_51mer['51mer ID'] == search_string, 'Best Peptide Class II'].values[0]
        #classI_sequence 
        classI_peptide = merged_peptide_51mer.loc[merged_peptide_51mer['51mer ID'] == search_string, 'Best Peptide Class I'].values[0]
        
        
        # mutant pepetide position ---
        mutant_peptide_pos = str(merged_peptide_51mer.loc[merged_peptide_51mer['51mer ID'] == search_string, 'Pos'].values[0])

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
            
            set_underline(peptide_sequence, mutant_peptide_pos, row['51mer ID'])

            set_span_tags(peptide_sequence) # pass by reference
            
            new_string = create_stylized_sequence(peptide_sequence)

            next_td_tags[2].string = new_string

            # Remove the tag_with_search_string from the BeautifulSoup tree
            tag_with_search_string.decompose()

            modified_html = peptides_51mer_soup.prettify(formatter=None)

        else:
            print("\nNOT FOUND: ", search_string)
            print("Mutant Peptide Position: ", mutant_peptide_pos)
            print("ClassI: ", classI_peptide)
            print("ClassII: ", classII_peptide, "\n")

    if args.WB:
        html_file_name = args.WB +  '/../manual_review/' + args.samp + ".Colored_Peptides.html" 
    else:
        html_file_name  =  args.samp + ".Colored_Peptides.html"

    with open(html_file_name, "w", encoding = 'utf-8') as file:
        file.write(modified_html)


if __name__ == "__main__":
    main()
