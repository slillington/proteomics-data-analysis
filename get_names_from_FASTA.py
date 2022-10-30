import os
#Reads in a fasta file and outputs a tab-separated list of the protein ID and the protein name
def file_parse(file_name,output_prefix):
    #Parse FASTA file into new files with <5000 protein sequences per file
    master = open(file_name,'r')
    master_text = master.read()
    block_list = master_text.split('>')
    master.close()
    file_list = []
    i = 0
    j = 1
    k = 1
    filename = output_prefix+str(k)+'.fasta'
    part_j = open(filename,'a+')
    while i < len(block_list)-1:
        if not i % 5000 == 0:
            string = ">"+block_list[i]+"\n"
            part_j.write(string)
            j = j+1
        else:
            file_list.append(filename)
            part_j.close()
            k = k+1
            filename = output_prefix + str(k) +'.fasta'
            part_j = open(filename,'a+')
            j = 1
        i = i+1
    part_j.close()
    return file_list
        
            
def sequences_in_file(filename):
    with open(filename,'r') as seqfile:
        txt = seqfile.read()
        ls = txt.split('>')
    return len(ls)

def get_name_from_FASTA(fasta_file,output_file, source):
    #Takes a proteome FASTA file as input and produces a csv output with protein IDs, names, and sequences
    #Mycocosm proteome FASTA files are different from Uniprot
    if source == 'Uniprot':
        with open(fasta_file, 'r') as seqFile:
            txt = seqFile.read().strip()
            ls = txt.split('>')
            with open(output_file, 'a+') as writeFile:
                for fasta_block in ls:
                    if len(fasta_block) > 1:
                        block_list = fasta_block.split('|')
                        #Need to split up block list to isolate entry ID, protein name, and sequence
                        #0 = 'tr' or 'jgi', 1 = EntryID, 2 = 'Entryname_9FUNG ProteinName OS= GN= SV=1\nM...sequence'
                        #prot_ids.append(block_list[1])
                        
                        #Get name and sequence from item 2
                        seq_raw = block_list[2]
                        seq_start = seq_raw.find('\n')
                        #prot_seqs.append(seq_raw[seq_start:])

                        name_start = seq_raw.find('_9FUNG')+7
                        name_end = seq_raw.find('OS=',name_start,seq_start)-1
                        prot_seq = seq_raw[name_start:name_end]
                        line = block_list[1]+","+seq_raw[name_start:name_end]+"\n"
                        writeFile.write(line)           
                    else:
                        print('Empty sequence')
            writeFile.close()
        #Close file to save memory once we don't need it anymore        
        seqFile.close()
    elif source == 'Mycocosm':
        with open(fasta_file, 'r') as seqFile:
            txt = seqFile.read().strip()
            ls = txt.split('>')
            with open(output_file, 'a+') as writeFile:
                for fasta_block in ls:
                    if len(fasta_block) > 1:
                        block_list = fasta_block.split('|')
                        #Need to split up block list to isolate entry ID, protein name, and sequence
                        #0 = 'tr' or 'jgi', 1 = species_ID, 2 = 'Mycocosm protein ID' 3 = string of numbers \nSEQUENCE
                        #prot_ids.append(block_list[1])
                        
                        #Get name and sequence from item 2
                        seq_raw = block_list[3]
                        seq_start = seq_raw.find('\n')
                        #prot_seqs.append(seq_raw[seq_start:])

                        #name_start = seq_raw.find('_9FUNG')+7
                        seq_end = seq_raw.find('*')-1
                        prot_seq = seq_raw[seq_start:seq_end]
                        prot_seq_stripped = prot_seq.replace('\n','')
                        line = block_list[2]+","+prot_seq_stripped+"\n"
                        writeFile.write(line)           
                    else:
                        print('Empty sequence')
            writeFile.close()
        #Close file to save memory once we don't need it anymore

def main():
    input_file = os.getcwd() + r"\Neosp1_GeneCatalog_proteins_20170918.aa.fasta"
    file_list = file_parse(input_file,'G1_Mycocosm_proteome_filtered_Mycocosm')
        
    for f in file_list:
        get_name_from_FASTA(f,"G1_proteome_filtered_Mycocosm.csv",'Mycocosm')
    #print(len(ids),len(sequences),len(cyscounts),len(names))




if __name__ == "__main__":
    main()
