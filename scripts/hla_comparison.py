import argparse
import csv
import pandas as pd
import os.path
import re

# python3 scripts/hla_comparison.py -WB "/Users/evelynschmidt/Bioinformatics_tools/neoag_vaccine_scripts/test_data/gcp_immuno"

# ---- PARSE ARGUMENTS -------------------------------------------------------
def parse_arguments():
    """Parses command line arguments and enables user help."""
    parser = argparse.ArgumentParser(description='Compare HLA alleles called by hlahd, phlat, opitype, and clincal data if available.')

    parser.add_argument('-WB', help='The path to the gcp_immuno folder of the trial you wish to run the script on, defined as WORKING_BASE in envs.txt')
    parser.add_argument('-f', "--fin_results", default="final_results", help="Name of the final results folder in gcp immuno")
    parser.add_argument("--optitype_n", help="File path for optitype normal calls")
    parser.add_argument("--optitype_t", help="File path for optitype tumor calls")
    parser.add_argument("--phlat_n", help="File path for phlat normal calls")
    parser.add_argument("--phlat_t", help="File path for phlat tumor calls")
    parser.add_argument("--hlahd_n", help="File path for hlahd normal calls")
    parser.add_argument("--hlahd_t", help="File path for hla tumor calls")
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
        
        # Sort alleles using the safe function
        alleles = list([row['Allele1'], row['Allele2']])
        
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
def process_clinical(input_file, sample_type):
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
    output_df['sample_type'] = [sample_type]
    
    return output_df

def sort_hla_pairs(df):
    # Make a copy to avoid modifying the original
    sorted_df = df.copy()
    
    # Define the allele pairs to be sorted
    pairs = [
        ('HLA_A-1', 'HLA_A-2'),
        ('HLA_B-1', 'HLA_B-2'),
        ('HLA_C-1', 'HLA_C-2'),
        ('DQA1-1', 'DQA1-2'),
        ('DQB1-1', 'DQB1-2'),
        ('DRB1-1', 'DRB1-2')
    ]
    
    # For each row, sort each pair of columns alphabetically
    def sort_row(row):
        row = row.copy()
        for c1, c2 in pairs:
            if c1 not in row.index or c2 not in row.index:
                continue

            v1, v2 = row[c1], row[c2]
            if pd.isna(v1) or pd.isna(v2):
                continue

            # compare as strings
            if str(v1) > str(v2):
                row[c1], row[c2] = v2, v1
        return row

    return df.apply(sort_row, axis=1)


def _two_field(allele: str) -> str:
    """
    Normalize to 2-field allele resolution:
      HLA-A*02:01:01 -> A*02:01
      HLA-DRB1*11:01:01 -> DRB1*11:01
      DRB3*02:96 -> DRB3*02:96 (already 2-field-ish)
    """
    allele = allele.strip()
    allele = allele.replace("HLA-", "").replace("HLA_", "")

    # If it looks like LOCUS*xx:yy(:...) keep only first two fields after '*'
    if "*" in allele:
        loc, rest = allele.split("*", 1)
        parts = rest.split(":")
        if len(parts) >= 2:
            return f"{loc}*{parts[0]}:{parts[1]}"
        else:
            return f"{loc}*{rest}"
    return allele


def process_hlahd(input_file, sample_type):
    """
    Parse HLA-HD *_DNA_final.result.txt format:
      locus <TAB> allele1 <TAB> allele2 <TAB> ...
    No header. Missing calls can be '-' or 'Not typed'.
    Some loci may have >2 allele tokens (e.g., E).
    """
    # loci you actually want to harmonize with your downstream columns
    # (add/remove here as needed)
    loci_keep = {"A", "B", "C", "DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1"}

    rows = []
    with open(input_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if "Couldn't read result file" in line:
                # keep going; file may still have usable rows above
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                continue

            locus = parts[0].strip()
            if locus not in loci_keep:
                continue

            allele_tokens = [p.strip() for p in parts[1:] if p.strip()]

            # filter missing tokens
            allele_tokens = [
                a for a in allele_tokens
                if a not in {"-", "Not typed", "Not typed\tNot typed"} and a.lower() != "not typed"
            ]

            # normalize + de-dup preserving order
            alleles_norm = []
            for a in allele_tokens:
                norm = _two_field(a)
                if norm not in alleles_norm:
                    alleles_norm.append(norm)

            rows.append((locus, alleles_norm))

    out = pd.DataFrame()

    for locus, alleles in rows:
        # map locus to your column naming scheme:
        # A/B/C -> HLA_A/HLA_B/HLA_C ; class II loci keep as-is
        if locus in {"A", "B", "C"}:
            base = f"HLA_{locus}"
        else:
            base = locus  # e.g. DQA1, DRB1, ...

        # take up to first 2 alleles (your rest of script expects pairs)
        for i, a in enumerate(alleles[:2], start=1):
            out.at[0, f"{base}-{i}"] = a

    out["Tool"] = ["hlahd"]
    out["sample_type"] = [sample_type]
    return out


def find_discrepancies(df):
    """
    For each (sample_type, locus), collapses each tool's allele slots into
    a sorted tuple (order-independent), then flags loci where tools disagree.
    Returns one row per (sample_type, locus, Tool) for loci with disagreement,
    so it's easy to see exactly which tool called what.
    """
    allele_cols = [c for c in df.columns if re.match(r'.+-\d+$', c)]
    tool_order = list(df['Tool'].unique())  # preserve the table's row order (optitype, phlat, hlahd, clinical, ...)

    long = (df.melt(id_vars=['Tool', 'sample_type'], value_vars=allele_cols,
                     var_name='column', value_name='allele')
              .dropna(subset=['allele']))
    long['locus'] = long['column'].str.rsplit('-', n=1).str[0]

    calls = (long.groupby(['sample_type', 'locus', 'Tool'])['allele']
                  .apply(lambda alleles: tuple(sorted(set(alleles))))
                  .reset_index(name='alleles'))
    calls['Tool'] = pd.Categorical(calls['Tool'], categories=tool_order, ordered=True)

    n_distinct_calls = calls.groupby(['sample_type', 'locus'])['alleles'].transform('nunique')
    return calls[n_distinct_calls > 1].sort_values(['sample_type', 'locus', 'Tool']).reset_index(drop=True)


def print_discrepancies(discrepancy_df):
    print("\n" + "=" * 80)
    print("HLA ALLELE DISCREPANCIES")
    print("=" * 80)

    if discrepancy_df.empty:
        print("No discrepancies found between tools.")
        return

    for (sample_type, locus), group in discrepancy_df.groupby(['sample_type', 'locus']):
        print(f"\n[{sample_type}] {locus}")
        for tool, alleles in zip(group['Tool'], group['alleles']):
            print(f"  {tool}: {', '.join(alleles)}")


def get_consensus_alleles(consensus_file):
    """
    Parses consensus_calls.txt and returns ClassI and ClassII alleles.
    ClassI starts with 'HLA-', ClassII starts with 'D'.
    """
    class_i_alleles = []
    class_ii_alleles = []
    
    try:
        with open(consensus_file, 'r') as file:
            content = file.read().strip()
            alleles = content.split(',')
            
            for allele in alleles:
                allele = allele.strip()
                if allele.startswith('HLA-'):
                    class_i_alleles.append(allele)
                elif allele.startswith('D'):
                    class_ii_alleles.append(allele)
    except FileNotFoundError:
        print(f"[WARN] Consensus file not found: {consensus_file}")
    
    return class_i_alleles, class_ii_alleles


def main():
    args = parse_arguments()
    
    if args.WB:
        hla_typing = f"{args.WB}/{args.fin_results}/hla_typing"  
        
        optitype_normal = process_optitype(f"{hla_typing}/optitype_normal_result.tsv", "normal")
        optitype_tumor = process_optitype(f"{hla_typing}/optitype_tumor_result.tsv", "tumor")
        phlat_normal = process_phlat(f"{hla_typing}/phlat_normal_HLA.sum", "normal")
        phlat_tumor = process_phlat(f"{hla_typing}/phlat_tumor_HLA.sum", "tumor")
        list_df = [optitype_normal, optitype_tumor, phlat_normal, phlat_tumor]

        # ---- HLA-HD (optional in older cases) ----
        hlahd_n_path = f"{hla_typing}/hlahd_normal_DNA_final.result.txt"
        hlahd_t_path = f"{hla_typing}/hlahd_tumor_DNA_final.result.txt"

        if os.path.isfile(hlahd_n_path):
            list_df.append(process_hlahd(hlahd_n_path, "normal"))
        else:
            print(f"[WARN] Missing HLA-HD normal file, skipping: {hlahd_n_path}")

        if os.path.isfile(hlahd_t_path):
            list_df.append(process_hlahd(hlahd_t_path, "tumor"))
        else:
            print(f"[WARN] Missing HLA-HD tumor file, skipping: {hlahd_t_path}")

        # ---- clinical (optional) ----
        clinical_path = f"{hla_typing}/clinical_calls.txt"
        if os.path.isfile(clinical_path):
            list_df.append(process_clinical(clinical_path, "normal"))
            
        df = pd.concat(list_df, ignore_index=True)

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
        if os.path.isfile(args.hlahd_n):
            hlahd_normal = process_hlahd(args.hlahd_n, "normal")
            list_df.append(hlahd_normal)
        if os.path.isfile(args.hlahd_t):
            hlahd_tumor = process_hlahd(args.hlahd_t, "tumor")
            list_df.append(hlahd_tumor)
        if os.path.isfile(args.clinical):
            clinical = process_clinical(args.clinical, "normal")
            list_df.append(clinical)
        
        df = pd.concat(list_df)
            
    first_column = df.pop('sample_type') 
    df.insert(0, 'sample_type', first_column) 
    
    first_column = df.pop('Tool') 
    df.insert(0, 'Tool', first_column) 
    
    sorted_df = sort_hla_pairs(df)
    print(sorted_df)

    print_discrepancies(find_discrepancies(sorted_df))

    # ---- Print consensus alleles if available ----
    if args.WB:
        consensus_path = f"{args.WB}/{args.fin_results}/hla_typing/consensus_calls.txt"
        class_i_alleles, class_ii_alleles = get_consensus_alleles(consensus_path)
        
        print("\n" + "=" * 80)
        print("pVAC CONSENSUS ALLELES")
        print("=" * 80)
        
        if class_i_alleles:
            class_i_str = "ClassI consensus alleles used by pVAC predictions: " + ", ".join(class_i_alleles)
            print(class_i_str)
        
        if class_ii_alleles:
            class_ii_str = "ClassII consensus alleles used by pVAC predictions: " + ", ".join(class_ii_alleles)
            print(class_ii_str)
    
    if args.WB:
        sorted_df.to_csv(f"{args.WB}/../manual_review/hla_comparison.tsv", sep='\t', index=False)
    if args.o:
        sorted_df.to_csv(f"{args.o}/hla_comparison.tsv", sep='\t', index=False)
    else:
        sorted_df.to_csv("hla_comparison.tsv", sep='\t', index=False)
    

if __name__ == "__main__":
    main()
