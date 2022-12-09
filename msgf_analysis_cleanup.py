"""
This Python file takes the output TSVs from MSGF+ bottom-up proteomic data analysis
and produces a final output similar to Sam Purvine's pipeline at EMSL. Works for
Python 3 only.
"""
import pandas as pd 
import numpy as np 
import os
import errno

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment


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
    filtered2 = filtered.loc[~filtered['Protein'].str.contains(fdr_tag)]
    fdr = float(filtered.loc[filtered['Protein'].str.contains(fdr_tag)].shape[0] / filtered.shape[0])

    #Note for multiple filters we would do
    #filtered = input_data.loc[(input_data['QValue'] < q_cutoff) & (input_data['Peptide'].str.contains('+142.03'))]
    #To filter by Qvalue by also filter for a specific PTM if we wanted to do that

    #Filter out columns
    filtered2 = filtered2[['Protein', 'Peptide']]

    return filtered2, fdr


def pivot_fxns(filtered_df):
    """
    Input: Filtered Pandas dataframe with proteins and peptides after QC filtering
    Output: Pandas dataframe in the pivot format
    """

    #First count the unique peptides that map to each protein
    piv= filtered_df.pivot_table(values='Peptide', index=['Protein'], aggfunc={'Peptide':[pd.Series.nunique, pd.Series.count]})
    piv.reset_index(inplace=True)
    #print(piv)

    return piv

def compute_coverage(dataframe, seq_db, outputHTML=False):
    """
     Computes the sequence coverage of all proteins with peptides that map to them.

     Input: tsv_file, type str. Path to the filtered tsv file containing the peptides
     identified by MS.

     Output: A dictionary with keys as protein names and values as computed coverage

    """
    seq_db = SeqIO.index(seq_db, "fasta") #type dict_keyiterator
    ids = seq_db.keys()
    coverage_dict = {}

    if outputHTML:
        html_message = """<html>
        <head></head>
        <title>Sequence coverage map</title>
        <body>\n"""

    for id in ids:
        iseq = seq_db[id]
        #print(iseq.seq)

        #Create boolean array from which to compute coverage
        cov_array = np.zeros(len(iseq.seq))

        #Slice dataframe to only contain peptides that map to given id
        sliced = dataframe.loc[dataframe['Protein'] == id]

        #Loop through Peptide column of sliced and search seq for substring
        for i,pep in sliced.iterrows():
            peptide = pep['Peptide'].replace('.','')
            #print(peptide)

            start = iseq.seq.find(peptide)
            end = start + len(peptide)

            if not start == -1:
                #Change values in cov_array
                cov_array[start:end] = 1

        #If html, add sequence info to html file and color sequence
        #according to coverage
        if outputHTML and np.sum(cov_array) > 1:
            html_message += "<br>>%s<br>" %id

            #Loop through string/coverage array
            for aa,cov in zip(iseq.seq,cov_array):
                if cov:
                    html_message += """<span style="color:#D4D100">%s</span>\n""" %aa
                else:
                    html_message += """<span style="color:#000000">%s</span>\n""" %aa

            #Add last piece to HTML file
            html_message += """</body>\n</html>"""




        coverage = np.sum(cov_array) / len(cov_array)
        coverage_dict[id] = coverage
        #print("Coverage for %s is %.2f" %(id,coverage))

    #Remove entries with 0 coverage
    coverage_dict2 = {k:v for k,v  in coverage_dict.items() if v != 0.}
    #coverage_dictF = dict((v,k) for k,v in coverage_dict2.items())

    if outputHTML:
        return coverage_dict2, html_message
    else:
        return coverage_dict2


def main():
    cb_path = r"C:\Users\steph\Box\OmalleyLabfiles\Data\Proteomics\G1_cellulosome_proteomics_2021\Cellobiose\July2022_data"
    fp_path = r"C:\Users\steph\Box\OmalleyLabfiles\Data\Proteomics\G1_cellulosome_proteomics_2021\FilterPaper\July2022_data"
    rcg_path = r"C:\Users\steph\Box\OmalleyLabfiles\Data\Proteomics\G1_cellulosome_proteomics_2021\RCG\July2022_data"
    output_path = r"C:\Users\steph\Box\OmalleyLabfiles\Data\Proteomics\G1_cellulosome_proteomics_2021\COV_results"

    #Iterate through file paths and run .tsv files
    for fpath in [rcg_path]:
        for tsvFile in os.listdir(fpath):
            if tsvFile.endswith('.tsv'):

                new_lib, fdr = tsv_filter(os.path.join(fpath,tsvFile),fdr_tag="XXX")
                seq_db = r"C:\Users\steph\Box\OmalleyLabfiles\Data\Proteomics\Neosp1_Catalog_doc_scaf_cont.fasta"
                print("Filtered library of shape %s with FDR = %.5f" %(str(new_lib.shape),fdr))
                #print(new_lib)

                
                #print(new_lib.shape)
                piv_table = pivot_fxns(new_lib)
                coverage_dict, coverage_html = compute_coverage(new_lib,seq_db,outputHTML=True)

                if coverage_html:
                    output_fname = os.path.split(tsvFile)[-1][:-4] + "_coverage.html"
                    with open(os.path.join(fpath,output_fname),'w') as f:
                        f.write(coverage_html)

                #print(len([k for k in coverage_dict.keys()]))
                """
                print(piv_table.columns)
                #Add coverages to piv_table

                piv_table['Coverage'] = piv_table['Protein'].map(coverage_dict)
                #piv_table['Coverage'] = piv_table.iloc[:,0].map(coverage_dict)

                #Write output to excel file
                
                filename = os.path.split(tsvFile)[-1]
                #print(os.path.split(my_file))

                new_filename = os.path.join(output_path,filename[:-4] + "_test.xlsx")
                print("Writing to file %s" %new_filename)
                piv_table.to_excel(new_filename)
                """


if __name__ == "__main__":
    main()