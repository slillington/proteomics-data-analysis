"""
This Python file takes the output TSVs from MSGF+ bottom-up proteomic data analysis
and produces a final output similar to Sam Purvine's pipeline at EMSL. Works for
Python 3 only.
"""
import pandas as pd 
import numpy as np 
import os
import errno



def tsv_filter(input_tsv, fdr_tag="DECOY"):
    """
    Input: 
        input_tsv - MSGF+ output after running MzidToTsvConverter. Note that
        for the false discovery rate to be accurate, MSGF+ has to be run with showDecoy = 1
        otherwise the FDR will just be 0.
        fdr_tag - substring for denoting decoy sequences in the protein library (the library sequences
        in reverse, which you know are incorrect). Default for MSGF+ is '%%%'.
    Output: 
        filtered - A pandas dataframe with the filtered data
        fdr - the false discovery rate aka what percentage of hits
        you can expect to be wrong on average 

    This function takes the output from running MSGF+ in .tsv format
    and filters out unneeded columns, but also filters out any peptides
    with Q-values > 0.01 and E-values > 0.01.
    """
    q_cutoff = 0.01
    e_value_cutoff = 0.01
    #First check if input_tsv file exists
    if not os.path.isfile(input_tsv):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), input_tsv)
    elif not '.tsv' in input_tsv:
        raise ValueError('%s is not .tsv and will not be read properly.' %input_tsv)

    #Load input into dataframe
    input_data = pd.read_csv(input_tsv,sep='\t')
    #print(input_data)

    #Filter out rows with Q-value > q_cutoff
    filtered = input_data.loc[input_data['QValue'] < q_cutoff].copy()
    fdr = float(filtered.loc[filtered['Protein'].str.contains(fdr_tag)].shape[0] / filtered.shape[0])

    #Note for multiple filters we would do
    #filtered = input_data.loc[(input_data['QValue'] < q_cutoff) & (input_data['Peptide'].str.contains('+142.03'))]
    #To filter by Qvalue by also filter for a specific PTM if we wanted to do that

    #Filter out columns
    filtered = filtered[['Protein', 'Peptide']]

    return filtered, fdr


def pivot_fxns(filtered_df):
    """
    Input: Filtered Pandas dataframe with proteins and peptides after QC filtering
    Output: Pandas dataframe in the pivot format
    """

    #First count the unique peptides that map to each protein
    piv= filtered_df.pivot_table(values='Peptide', index=['Protein'], aggfunc={'Peptide':[pd.Series.nunique, pd.Series.count]})
    #print(piv)

    return piv





def main():
    my_file = r"C:\Users\steph\Box\OmalleyLabfiles\Data\Proteomics\G1_cellulosome_proteomics_2021\RCG\rcg0525_UniProtlib.tsv"
    output_path = r"C:\Users\steph\Box\OmalleyLabfiles\Data\Proteomics\G1_cellulosome_proteomics_2021\RCG\\"
    new_lib, fdr = tsv_filter(my_file,fdr_tag="DECOY")

    print("Filtered library of shape %s with FDR = %.5f" %(str(new_lib.shape),fdr))
    #print(new_lib)

    piv_table = pivot_fxns(new_lib)

    #Write output to excel file
    filename = os.path.split(my_file)[-1]
    #print(os.path.split(my_file))

    new_filename = output_path + filename[:-4] + "_processed.xlsx"
    piv_table.to_excel(new_filename)



if __name__ == "__main__":
    main()