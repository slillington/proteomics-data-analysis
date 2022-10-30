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