#Goal: Make a FASTA file with the Top 25 hits by protein sample from G1 cellulosome proteomics data
#1. Load data in from excel file
#Sort by each sample's column and get the list of the 25 Locus numbers
#Get the corresponding protein sequence from the FASTA database
from Bio import SeqIO
from file_parse import *

locus_path = r"C:\Users\steph\Box\OmalleyLabfiles\Data\Proteomics\G1_acetonePrec_cellulosome\ID_007523_2D1C82EB.fasta"
#mycocosm_path = r"C:\Users\steph\Box\OmalleyLabfiles\Data\Proteomics\Haitjema2017reanalysis\Ncaliforniae\Neosp1_GeneCatalog_proteins_20170918.aa"
locus_files = file_parse(locus_path,"EMSL_G1_fasta_")

first_pass_path = r"C:\Users\steph\Box\OmalleyLabfiles\Data\Proteomics\G1_acetonePrec_cellulosome\EMSL_49765_OMall_Iso_FirstHitsResults"

#Load data from first pass Excel file

locus data = SeqIO.parse(locus_path, 'fasta')

#Can't think of a faster way than a linear search
for locus_entry in locus_data:
    if locus_entry.seq 