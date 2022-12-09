# proteomics-data-analysis
Takes raw MS/MS peptide-protein mapping output from MSGF+ and performs analysis routines

msgf_analysis_cleanup.py -
 Takes the output TSVs from MSGF+ bottom-up proteomic data analysis
 and produces a final output similar to Sam Purvine's pipeline at EMSL.
 Example output is an excel file with columns: Protein ID, peptide count, unique peptide count, sequence coverage.

 The functions filter the MSGF+ output to a specified false discovery rate.
 **As of July 2022, I have been running MzidToTsvConverter with the following specs:
    MzidToTsvConverter -mzid:<mzid-file-path> -u -sd -MaxE 0.001

 So MzidToTsvConverter filters the spectra matches with an E-value cutoff. See the MzidToTsvConverter
 documentation for more details: https://github.com/PNNL-Comp-Mass-Spec/Mzid-To-Tsv-Converter
