import argparse
import csv
import pandas as pd
import os.path

# python3 scripts/hla_comparison.py -WB "/Users/evelynschmidt/Bioinformatics_tools/neoag_vaccine_scripts/test_data/gcp_immuno"

# ---- PARSE ARGUMENTS -------------------------------------------------------
def parse_arguments():
    """Parses command line arguments and enables user help."""
    parser = argparse.ArgumentParser(description='Compare HLA alleles called by phlat, opitype, and clincal data if available.')

    parser.add_argument('-WB', help='The path to the gcp_immuno folder of the trial you wish to run the script on, defined as WORKING_BASE in envs.txt')
    parser.add_argument('-f', "--fin_results", default="final_results", help="Name of the final results folder in gcp immuno")
    parser.add_argument("--optitype_n", help="File path for optitype normal calls")
    parser.add_argument("--optitype_t", help="File path for optitype tumor calls")
    parser.add_argument("--phlat_n", help="File path for phlat normal calls")
    parser.add_argument("--phlat_t", help="File path for phlat tumor calls")
    parser.add_argument("--clinical", help="File path for the clinical_calls.txt")
    parser.add_argument("--o", help="Output folder")

    return parser.parse_args()


# A function that processes the phlat calls
def process_phlat(input_file, sample_type):
    df = pd.read_csv(input_file, sep='\t')
    
    # Initialize a dictionary to store sorted alleles by locus
    allele_columns = {}

    # Collect and sort alleles for each locus
    for _, row in df.iterrows():
        locus = row['Locus']
        
        # Define a safe function to extract the numeric part
        def get_numeric_part(allele):
            try:
                # Attempt to extract the numeric part after the first ':'
                return int(allele.split(':')[1])
            except (IndexError, ValueError):
                # Return a large number to keep non-numeric alleles at the end
                return float('inf')
        
        # Sort alleles using the safe function
        alleles = sorted([row['Allele1'], row['Allele2']], key=get_numeric_part)
        
        # If the locus is not yet in the dictionary, add it with the sorted alleles
        if locus not in allele_columns:
            allele_columns[locus] = alleles
        else:
            # Extend the list if there are more alleles for this locus
            allele_columns[locus].extend(alleles)
    
    # Create a new DataFrame with a single row
    new_df = pd.DataFrame()

    # Populate the new DataFrame columns with alleles, naming them dynamically
    for locus, alleles in allele_columns.items():
        # Remove 'HLA_' prefix if it's a class II locus with "D" in it
        if "D" in locus:
            locus_parts = locus.split('_') # No prefix
            locus_name = locus_parts[1]
        else:
            locus_name = locus
            
        for i, allele in enumerate(alleles, start=1):
            column_name = f"{locus_name}-{i}"
            # Extract only the first two components of the allele (e.g., "A*01:01")
            new = ':'.join(allele.split(':')[:2])
            new_df.at[0, column_name] = new
    
    new_df['Tool'] = ["phlat"]
    new_df['sample_type'] = [sample_type]
    return new_df


                
# A function that process the opitype calls
def process_optitype(input_file, sample_type):
    df = pd.read_csv(input_file, sep='\t')
    
    # Remove the last two columns
    df = df.drop(df.columns[-2:], axis=1)
    df = df.drop(df.columns[:1], axis=1)
    
    # Rename columns
    df.columns = ['HLA_A-1', 'HLA_A-2', 'HLA_B-1', 'HLA_B-2', 'HLA_C-1', 'HLA_C-2']
    
    # Sort alleles within each HLA type category for each row
    for i, row in df.iterrows():
        # Sort A alleles
        a_alleles = sorted([row['HLA_A-1'], row['HLA_A-2']], key=lambda x: int(x.split(':')[1]))
        # Sort B alleles
        b_alleles = sorted([row['HLA_B-1'], row['HLA_B-2']], key=lambda x: int(x.split(':')[1]))
        # Sort C alleles
        c_alleles = sorted([row['HLA_C-1'], row['HLA_C-2']], key=lambda x: int(x.split(':')[1]))
        
        # Assign back to the DataFrame
        df.loc[i, ['HLA_A-1', 'HLA_A-2']] = a_alleles
        df.loc[i, ['HLA_B-1', 'HLA_B-2']] = b_alleles
        df.loc[i, ['HLA_C-1', 'HLA_C-2']] = c_alleles
    
    df['Tool'] = ["optitype"]
    df['sample_type'] = [sample_type]
    
    return(df)

# A function that process the clincal calls
def process_clinical(input_file):
    # Load the data and transpose to get alleles in a single row
    df = pd.read_csv(input_file, sep=',')
    alleles = df.columns.tolist()
    
    # Dictionary to store sorted alleles by locus
    allele_dict = {}

    # Group alleles by their locus
    for allele in alleles:
        # Split by '*' to get the locus
        locus = allele.split('*')[0]
        locus = locus.replace("-", "_")
        
        # Add allele to the list for this locus
        if locus not in allele_dict:
            allele_dict[locus] = [allele]
        else:
            allele_dict[locus].append(allele)

    # Sort alleles within each locus and create the output DataFrame
    output_df = pd.DataFrame()
    
    for locus, alleles in allele_dict.items():
        # Sort alleles by the numeric part, if available
        def get_numeric_part(allele):
            try:
                # Attempt to extract the numeric part after the first ':'
                return int(allele.split(':')[1])
            except (IndexError, ValueError):
                # Return a large number to keep non-numeric alleles at the end
                return float('inf')
        
        sorted_alleles = sorted(alleles, key=get_numeric_part)
        
        # Create columns for each allele within this locus
        for i, allele in enumerate(sorted_alleles, start=1):
            column_name = f"{locus}-{i}"
            new = allele.split('-')
            output_df.at[0, column_name] = new[1] if len(new) > 1 else new[0]
            
    output_df['Tool'] = ["clinical"]
    
    return output_df

        
def main():
    args = parse_arguments()
    
    if args.WB:
        hla_typing = f"{args.WB}/{args.fin_results}/hla_typing"  
        
        optitype_normal = process_optitype(f"{hla_typing}/optitype_normal_result.tsv", "normal")
        optitype_tumor = process_optitype(f"{hla_typing}/optitype_tumor_result.tsv", "tumor")
        phlat_normal = process_phlat(f"{hla_typing}/phlat_normal_HLA.sum", "normal")
        phlat_tumor = process_phlat(f"{hla_typing}/phlat_tumor_HLA.sum", "tumor")
        
        if os.path.isfile(f"{hla_typing}/clinical_calls.txt"):
            clincal = process_clinical(f"{hla_typing}/clinical_calls.txt")
            df = pd.concat([optitype_normal, optitype_tumor, phlat_normal, phlat_tumor, clincal])
        else:
            df = pd.concat([optitype_normal, optitype_tumor, phlat_normal, phlat_tumor])
    else:
        list_df = []
        if os.path.isfile(args.optitype_n):
            optitype_normal = process_optitype(args.optitype_n, "normal")
            list_df.append(optitype_normal)
        if os.path.isfile(args.optitype_t):
            optitype_tumor = process_optitype(args.optitype_t, "tumor")
            list_df.append(optitype_tumor)
        if os.path.isfile(args.phlat_n):
            phlat_normal = process_phlat(args.phlat_n, "normal")
            list_df.append(phlat_normal)
        if os.path.isfile(args.phlat_t):
            phlat_tumor = process_phlat(args.phlat_t, "tumor")
            list_df.append(phlat_tumor)
        if os.path.isfile(args.phlat_t):
            clinical = process_clinical(args.clinical)
            list_df.append(clinical)
        
        df = pd.concat(list_df)
            
    first_column = df.pop('sample_type') 
    df.insert(0, 'sample_type', first_column) 
    print(df)
    
    first_column = df.pop('Tool') 
    df.insert(0, 'Tool', first_column) 
    print(df)
    
    if args.WB:
        df.to_csv(f"{args.WB}/../manual_review/hla_comparison.tsv", sep='\t', index=False)
    if args.o:
        df.to_csv(f"{args.o}/hla_comparison.tsv", sep='\t', index=False)
    else:
        df.to_csv("hla_comparison.tsv", sep='\t', index=False)
    

if __name__ == "__main__":
    main()