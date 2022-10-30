'''merge_annotations.py - For Mycocosm Interpro annotations, collapses multiple domain annotations for same protein ID
from multiple lines and cells to one.
example: 
Protein ID  Interpro ID     Interpro desc.
146         IP00452         GTPase
146         IP02732         Dockerin

will be collapsed to:
Protein ID  Interpro ID         Interpro desc.
146         IP00452,IP02732     GTPase,Dockerin
'''
import pandas as pd
#This is easily accomplished with pandas
def merge_annotations(iprfile,idheader,outputfile):
    inputpd = pd.ExcelFile(iprfile)
    df = pd.read_excel(inputpd,0)
    grouped = df.groupby(idheader).agg(lambda x: x.tolist())
    grouped.to_excel(outputfile)