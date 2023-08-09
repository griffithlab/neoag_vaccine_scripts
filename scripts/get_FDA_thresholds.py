import argparse
import csv
import pandas as pd
import sys



# Notes
"""
Generates a tsv file containing the FDA quality thresholds for the 
Gemonics Review Report. This script works with various trials and 
includes optional argument flags in case of werid naming conventions.

ex
python3 get_FDA_thresholds.py JLF mcdb046

python3 get_FDA_thresholds.py -WB /storage1/fs1/mgriffit/Active/JLF_MCDB/cases/JLF-100-042/gcp_immuno -f final_results_v1
python3 get_FDA_thresholds.py --n_dna aligned_normal_dna_table2.csv --t_dna aligned_tumor_dna_table2.csv --t_rna aligned_tumor_rna_table3.csv --una_n_dna unaligned_normal_dna_table1.csv --una_t_dna unaligned_tumor_dna_table1.csv --una_t_rna unaligned_tumor_rna_table1.csv --somalier concordance.somalier.pairs.tsv --contam_n normal.VerifyBamId.selfSM --contam_t tumor.VerifyBamId.selfSM

"""


# ---- PARSE ARGUMENTS -------------------------------------------------------
# Parses command line arguments
# Enables user help
# Future impovements: require the user to enter either the WB OR the list of files
def parse_arguments():
    # Parse command line arugments
    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument('-WB',
                        help='the path to the gcp_immuno folder of the trial you wish to tun script on, defined as WORKING_BASE in envs.txt')

    # The name of the final results folder 
    parser.add_argument('-f', "--fin_results", help="Name of the final results folder in gcp immuno")

    # Allows user to specify direct path to every file in case of werid naming convention
    parser.add_argument("--n_dna", help="file path for aligned normal dna FDA report table")
    parser.add_argument("--t_dna", help="file path for aligned tumor dna FDA report table")
    parser.add_argument("--t_rna", help="file path for aligned tumor rna FDA report table")
    parser.add_argument("--una_n_dna", help="file path for unaligned normal dna FDA report table")
    parser.add_argument("--una_t_dna", help="file path for unaligned tumor dna FDA report table")
    parser.add_argument("--una_t_rna", help="file path for unaligned tumor rna FDA report table")
    parser.add_argument("--somalier", help="file path for Somalier results for sample tumor/normal sample relatedness")
    parser.add_argument("--contam_n", help="file path for VerifyBamID results for contamination the normal sample")
    parser.add_argument("--contam_t", help="file path for VerifyBamID results for contamination the tumor sample")

    return(parser.parse_args())


# ---- RESHAPE QUALITY THRESHOLDS --------------------------------------------
# Reads file with thresholds and sets up a dataframe for recording case specfic values
def reshape_quality_thresholds(fda_quality_thresholds):
        # Split 'Column1' into two columns based on "("
    fda_quality_thresholds[['Criteria_1', 'File']] = fda_quality_thresholds.Criteria.str.split('(', expand=True)
    # Remove parentheses from the 'Column1' values
    fda_quality_thresholds['File'] = fda_quality_thresholds['File'].str.replace(r'\(|\)', '')
    # Split ro
    fda_quality_thresholds = fda_quality_thresholds.drop('File', axis=1) \
            .join(fda_quality_thresholds['File'] \
            .str \
            .split('/', expand=True) \
            .stack() \
            .reset_index(level=1, drop=True).rename('File')) \
            .reset_index(drop=True)

    # Add "DNA" after "tumor" if "DNA" is not already present
    mask = ~fda_quality_thresholds['File'].str.contains('DNA|RNA')
    fda_quality_thresholds.loc[mask, 'File'] = fda_quality_thresholds.loc[mask, 'File'].str.replace('tumor', 'tumor DNA')

    # Replace ''Criteria' values with ''Criteria_1' values
    fda_quality_thresholds['Criteria'] = fda_quality_thresholds['Criteria_1']


    fda_quality_thresholds['File'] = fda_quality_thresholds['File'].str.rstrip(')') # Remove closing parenthese 
    fda_quality_thresholds['Criteria'] = fda_quality_thresholds['Criteria'].str.strip() # Remove extra space at the end of Criteria
    fda_quality_thresholds['File'] = fda_quality_thresholds['File'].str.strip() # Remove extra space at the end of File

    # Delete 'Criteria_1
    qc = fda_quality_thresholds.drop('Criteria_1', axis=1)

    return(qc)



# ---- GET VALUES ------------------------------------------------------------
# Gets values for various files and places them in the dataframe
def get_values(qc, normal_dna, tumor_dna, tumor_rna, unalgined_normal_dna, unalgined_tumor_dna, unalgined_tumor_rna, somalier, contamination_normal, contamination_tumor):
    normal_temp = qc[qc.File == 'normal DNA']
    normal_temp = pd.merge(normal_temp , normal_dna, on="Criteria", how='left')

    tumor_temp = qc[qc.File == 'tumor DNA']
    tumor_temp = pd.merge(tumor_temp , tumor_dna, on="Criteria", how='left')

    rna_temp = qc[qc.File == 'tumor RNA']
    rna_temp = pd.merge(rna_temp, tumor_rna, on="Criteria", how='left')

    qc = pd.concat([normal_temp, tumor_temp, rna_temp], ignore_index =True)

    qc['Pass'] = ''
    qc = qc.astype({'Value':'float'})


    # get total numer of reads
    # form Col[Ssample Name] == Total Number of Reads, get col[jlf-100-042-normal-exome]
    # THIS DOES NOT WORK ALL THE TIME
    # CHANGE TO JUST GRAB FROM THE FIRST COLUMN
    total_number_reads_normal = unalgined_normal_dna.loc[unalgined_normal_dna["Sample Name"] == "Total Number of Reads"].values[0]
    qc.loc[(qc["Criteria"] == "TOTAL_READS") & (qc["File"] == "normal DNA"), "Value"] = total_number_reads_normal[1]

    total_number_reads_tumor = unalgined_tumor_dna.loc[unalgined_tumor_dna["Sample Name"] == "Total Number of Reads"].values[0]
    qc.loc[(qc["Criteria"] == "TOTAL_READS") & (qc["File"] == "tumor DNA"), "Value"] = total_number_reads_tumor[1]

    total_number_reads_tumor_rna = unalgined_tumor_rna.loc[unalgined_tumor_rna["Sample Name"] == "Total Number of Reads"].values[0]
    qc.loc[(qc["Criteria"] == "TOTAL_READS") & (qc["File"] == "tumor RNA"), "Value"] = total_number_reads_tumor_rna[1]

    # Open the TSV file
    with open(somalier, 'r') as file:
        # Create a TSV reader object
        reader = csv.DictReader(file, delimiter='\t')
        
        # Get the first entry in the specified column
        for row in reader:
            first_entry = row['relatedness']
            break

        # maybe get rid of the normal DNA row
        qc.loc[(qc["Criteria"] == "Genotype Concordance"), "Value"] = float(first_entry)


    file_paths = [contamination_normal, contamination_tumor]

    for file_path in file_paths:
        with open(file_path, 'r') as file:
            # Create a TSV reader object
            reader = csv.DictReader(file, delimiter='\t')
            
            # Get the first entry in the specified column
            for row in reader:
                first_entry = row['FREEMIX']
                break

            if file_path == contamination_normal:
                qc.loc[(qc["Criteria"] == "Contamination Estimate") & (qc["File"] == "normal DNA"), "Value"] = float(first_entry)
            else:
                qc.loc[(qc["Criteria"] == "Contamination Estimate") & (qc["File"] == "tumor DNA"), "Value"] = float(first_entry)
    
    return(qc)


# ---- EVALUATE THRESHOLDS ---------------------------------------------------
# Evaluates the value of Criteria and adds PASS or FAIL status to dataframe
# Values of FDA thresholds are currently hard coded 
# -- future iterations could reolsve this
# convert to precentage for visualization -- loses some percision
def evaluate_thresholds(qc):
    
    for index, row in qc.iterrows():

        if (row["File"] == 'normal DNA') and (row["Criteria"] == 'TOTAL_READS'):
            if float(row["Value"])/100000000.00 > 1:
                qc.at[index, "Pass"] = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL"

        if (row["File"] == "tumor DNA") and (row["Criteria"] == "TOTAL_READS"):
            
            if float(row["Value"])/250000000.00 > 1:
                qc.at[index, "Pass"]  = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL"

        if (row["File"] == "tumor RNA") and (row["Criteria"] == "TOTAL_READS"):
            
            if float(row["Value"])/200000000.00 > 1:
                qc.at[index, "Pass"]  = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL"

        if row["Criteria"] == "PCT_USABLE_BASES_ON_TARGET":
            if float(row["Value"])/0.2 > 1:
                qc.at[index, "Pass"] = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL"
            
            
            qc.at[index, "Value"] = "{:.2%}".format(row["Value"])

        if row["Criteria"] == "PCT_EXC_OFF_TARGET":
            if row["Value"]/0.60 < 1:
                qc.at[index, "Pass"] = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL"

        # fix to be decimal and not percentage
        if row["Criteria"] == "PERCENT_DUPLICATION":
            if row["Value"]/40 < 1:
                qc.at[index, "Pass"] = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL"

            qc.at[index, "Value"] = qc.at[index, "Value"]/100


        if (row["File"] == 'normal DNA') and (row["Criteria"] == 'MEAN_TARGET_COVERAGE'):
            if row["Value"]/100.00 > 1:
                qc.at[index, "Pass"] = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL"

        if (row["File"] == "tumor DNA") and (row["Criteria"] == "MEAN_TARGET_COVERAGE"):
            
            if row["Value"]/250.00 > 1:
                qc.at[index, "Pass"] = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL"
        
        if row["Criteria"] == "PCT_TARGET_BASES_20X":
            if row["Value"]/0.95 > 1:
                qc.at[index, "Pass"] = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL"
            

        if row["Criteria"] == "PCT_READS_ALIGNED_IN_PAIRS":
            if row["Value"]/0.95 > 1:
                qc.at[index, "Pass"] = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL"
            
        
        if row["Criteria"] == "MEAN_INSERT_SIZE":
            if row["Value"] >= 125 or row["Value"] <= 300:
                qc.at[index, "Pass"] = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL"

        if row["Criteria"] == "PF_MISMATCH_RATE_1":
            if row["Value"]/0.0075 < 1:
                qc.at[index, "Pass"] = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL"
            

        if row["Criteria"] == "PF_MISMATCH_RATE_2":
            if row["Value"]/0.01 < 1:
                qc.at[index, "Pass"] = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL" 


        # fix to be decimal and not percentage
        if row["Criteria"] == "PCT_PF_READS_ALIGNED":
            if row["File"] != "tumor RNA":
                if row["Value"]/95 > 1:
                    qc.at[index, "Pass"] = "PASS"
                else:
                    qc.at[index, "Pass"] = "FAIL" 
            else:
                if row["Value"]/85 > 1:
                    qc.at[index, "Pass"] = "PASS"
                else:
                    qc.at[index, "Pass"] = "FAIL" 

            qc.at[index, "Value"] = qc.at[index, "Value"]/100
                
        if row["Criteria"] == "Genotype Concordance":
            if row["Value"]/0.95 > 1:
                qc.at[index, "Pass"] = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL" 
        
        if row["Criteria"] == "Contamination Estimate":
            if row["Value"]/0.075 < 1:
                qc.at[index, "Pass"] = "PASS"
            else:
                qc.at[index, "Pass"] = "FAIL" 
            

    return qc


def main():

    args = parse_arguments()


    #sys.exit("Not a known trial\n")

    if args.fin_results:
        final_result = '/' + args.fin_results
    else:
        final_result = "/final_results"



    # Paths for various files needed
    if args.n_dna:
        normal_dna = pd.read_csv(args.n_dna, names=["Criteria", "Value"])
    else:
        normal_dna = pd.read_csv(args.WB + final_result + '/qc/fda_metrics/aligned_normal_dna/aligned_normal_dna_table2.csv', names=["Criteria", "Value"])

    # Reshaping some dataframe for join
    # maybe this could be a function since it is done three timmes
    normal_dna.loc[(normal_dna["Criteria"] == "Total Mapped Reads (%)"), "Criteria"] = "PCT_PF_READS_ALIGNED" # rename column so join works
    normal_dna['Value'] = normal_dna['Value'].str.rstrip(' (%)') # Remove  (%)


    if args.t_dna:
        tumor_dna = pd.read_csv(args.t_dna, names=["Criteria", "Value"])
    else:
        tumor_dna = pd.read_csv(args.WB + final_result + '/qc/fda_metrics/aligned_tumor_dna/aligned_tumor_dna_table2.csv', names=["Criteria", "Value"])

    tumor_dna.loc[(tumor_dna["Criteria"] == "Total Mapped Reads (%)"), "Criteria"] = "PCT_PF_READS_ALIGNED" # rename column so join works
    tumor_dna['Value'] = tumor_dna['Value'].str.rstrip(' (%)') # Remove  (%)

    if args.t_rna:
        tumor_rna = pd.read_csv(args.t_rna, names=["Criteria", "Value"])
    else:
        tumor_rna = pd.read_csv(args.WB + final_result + '/qc/fda_metrics/aligned_tumor_rna/aligned_tumor_rna_table3.csv', names=["Criteria", "Value"])

    tumor_rna.loc[(tumor_rna["Criteria"] == "Total mapped reads (%)"), "Criteria"] = "PCT_PF_READS_ALIGNED" # rename column so join works
    tumor_rna['Value'] = tumor_rna['Value'].str.rstrip(' (%)') # Remove  (%)

    if args.una_n_dna:
        unalgined_normal_dna = pd.read_csv(args.una_n_dna)
    else:
        unalgined_normal_dna = pd.read_csv(args.WB + final_result + '/qc/fda_metrics/unaligned_normal_dna/unaligned_normal_dna_table1.csv')

    if args.una_t_dna:
        unalgined_tumor_dna = pd.read_csv(args.una_t_dna)
    else:
        unalgined_tumor_dna = pd.read_csv(args.WB + final_result + '/qc/fda_metrics/unaligned_tumor_dna/unaligned_tumor_dna_table1.csv')

    if args.una_t_rna:
        unalgined_tumor_rna = pd.read_csv(args.una_t_rna)
    else:
        unalgined_tumor_rna = pd.read_csv(args.WB + final_result + '/qc/fda_metrics/unaligned_tumor_rna/unaligned_tumor_rna_table1.csv')

    if args.somalier:
        somalier = args.somalier
    else:
        somalier = args.WB + final_result + '/qc/concordance/concordance.somalier.pairs.tsv'

    if args.contam_n:
        contamination_normal = args.contam_n
    else:
        contamination_normal = args.WB + final_result + '/qc/normal_dna/normal.VerifyBamId.selfSM'

    if args.contam_t:
        contamination_tumor = args.contam_t
    else:
        contamination_tumor = args.WB + final_result + '/qc/tumor_dna/tumor.VerifyBamId.selfSM'




    fda_quality_thresholds = pd.read_csv("/opt/scripts/fda_quality_thresholds.csv")

    qc = reshape_quality_thresholds(fda_quality_thresholds)

    qc = get_values(qc, normal_dna, tumor_dna, tumor_rna, unalgined_normal_dna, unalgined_tumor_dna, unalgined_tumor_rna, somalier, contamination_normal, contamination_tumor)

    qc = evaluate_thresholds(qc)

    qc = qc.sort_values('Criteria', ignore_index=True)

    if args.WB:
        qc.to_csv(args.WB + final_result + '/qc/fda_quality_thresholds_report.tsv', sep="\t", index=False)
    else:
        qc.to_csv('fda_quality_thresholds_report.tsv', sep="\t", index=False)


    print(qc)



if __name__ == "__main__":
    main()

