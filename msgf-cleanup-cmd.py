"""
This Python file takes the output TSVs from MSGF+ bottom-up proteomic data analysis
and produces a final output similar to Sam Purvine's pipeline at EMSL. Works for
Python 3 only.
"""
import pandas as pd 
import numpy as np 
import os
import sys
import errno
import argparse

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment

Usage = """Takes the tsv output from MzidToTsvConverter after MSGF+ proteomics
    data analysis runs and produces a list of protein IDs with counts (unique
    and non-unique) of spectra that mapped to each protein ID.
    
    msgf-cleanup-cmd.py output_xls [tsv-file-path-list]
     Include file path in brackets even if only converting one file.
    
    **As of July 2022, I have been running MzidToTsvConverter with the following specs:
    MzidToTsvConverter -mzid:<mzid-file-path> -u -sd -MaxE 0.001
    
    """

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

def compute_coverage(tsv_file, seq_db):
    """
     Computes the sequence coverage of all proteins with peptides that map to them.

     Input: tsv_file, type str. Path to the filtered tsv file containing the peptides
     identified by MS.

     Output: A list of protein IDs with % sequence coverage

    """
    df = pd.read_csv(tsv_file,sep='\t')
    seq_db = SeqIO.index(seq_db, "fasta") #type dict_keyiterator
    #seq_db[<id>] will return the Seq object to pull the sequence

    print(next(seq_db)['id'])


def pivot_fxns(filtered_df):
    """
    Input: Filtered Pandas dataframe with proteins and peptides after QC filtering
    Output: Pandas dataframe in the pivot format
    """

    #First count the unique peptides that map to each protein
    piv= filtered_df.pivot_table(values='Peptide', index=['Protein'], aggfunc={'Peptide':[pd.Series.nunique, pd.Series.count]})
    #print(piv)

    return piv



def convert_files(filelist,output_path):
    for f in filelist:
        new_lib, fdr = tsv_filter(f,fdr_tag="DECOY")

        print("Filtered library of shape %s with FDR = %.5f" %(str(new_lib.shape),fdr))
        #print(new_lib)

        piv_table = pivot_fxns(new_lib)

        #Write output to excel file
        filename = os.path.split(f)[-1]
        #print(os.path.split(my_file))

        #Check if output_path ends with '\'
        if not output_path[-1] == '\\':
            output_path += '\\'

        new_filename = output_path + filename[:-4] + "_processed.xlsx"
        #print(new_filename)
        piv_table.to_excel(new_filename)



parser = argparse.ArgumentParser(description=Usage)
parser.add_argument('-p', '--path_list',nargs='+',
                metavar='path_list',
                default=[],
                help="path_list is a list of .tsv files to be converted e.g. -p G1rcg.tsv G1fp.tsv G1cb.tsv")

parser.add_argument('-op', '--output_path',
                metavar='output_path', default=os.getcwd(), type=str,
                help="output_path is the path to which the processed file will be written")

args = parser.parse_args()
tsv_list = args.path_list
outpath = args.output_path
print(outpath)
#Check that files are file
for t in tsv_list:
    if not os.path.isfile(t):
        print("%s cannot be found in current directory: %s"
            %(t,os.getcwd()))
        sys.exit()

if not os.path.isdir(outpath):
    print("Output directory %s is invalid" %outpath)
    sys.exit()

#Run convert_files
#convert_files(tsv_list,outpath)
compute_coverage(tsv_list[0],outpath)






