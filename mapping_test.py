import sys,os,re

############### FUNCTIONS TO :

## 1/ Check, 

def check(input_file):
    '''
    Check if the input file exists and is in .sam format.
    Verifies that each alignment line (non-header) contains at least the 11 "mandatory" fields
    '''
    if not os.path.exists(input_file): #check existence
        print("No file found")
        sys.exit(1)
    if not input_file.endswith(".sam"): #check extension
        print("File must be in .sam format")
        sys.exit(1)
    print("Exists and format check OK")

    with open(input_file, "r") as file: #open file in read mode
        for line_index, line in enumerate(file,start=1):
            if line.startswith("@"): #check header
                if line.startswith("@SQ"):
                    if "SN:" not in line or "LN:" not in line: #check for mandatory fields
                        print(f"Error at line {line_index}: @SQ must contain SN: and LN:.")
                        sys.exit(1)
                continue

            #Process alignment lines
            columns = line.strip().split("\t")

            #Check for at least 11 mandatory fields
            if len(columns) < 11:
                print(f"Format error at line {line_index}: less than 11 fields.")
                sys.exit(1)

            #Ger the first 11 mandatory fields
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = columns[:11]

            #Check for empty mandatory fields
            if any(field == "" for field in columns[:11]):
                print(f"Format error at line {line_index}: one of the 11 mandatory fields is empty.")
                sys.exit(1)

            #Check if each field is in correct Regexp/Range
            #QNAME
            if not re.fullmatch(r'[!-?A-~]{1,254}', qname):
                print(f"Format error at line {line_index}: QNAME incorrect.")
                sys.exit(1)
            
            #FLAG: first check if of type integer then if in the right range
            try:
                int(flag)
            except ValueError:
                print(f"Format error at line {line_index}: FLAG must be an interger.")
                sys.exit(1)
            
            if not (0 <= int(flag) <= 2**16 - 1):
                print(f"Format error at line {line_index}: FLAG incorrect.")
                sys.exit(1)

            #RNAME
            if not re.fullmatch(r'[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*', rname):
                print(f"Format error at line {line_index}: RNAME incorrect.")
                sys.exit(1)

            #POS: first check if of type integer then if in the right range
            try:
                int(pos)
            except ValueError:
                print(f"Format error at line {line_index}: POS must be an interger.")
                sys.exit(1)

            if not (0 <= int(pos) <= 2**31 - 1):
                print(f"Format error at line {line_index}: POS incorrect.")
                sys.exit(1)

            #MAPQ: first check if of type integer then if in the right range
            try:
                int(mapq)
            except ValueError:
                print(f"Format error at line {line_index}: MAPQ must be an interger.")
                sys.exit(1)

            if not (0 <= int(mapq) <= 2**8 - 1):
                print(f"Format error at line {line_index}: MAPQ incorrect.")
                sys.exit(1)

            #CIGAR
            if not re.fullmatch(r'\*|([0-9]+[MIDNSHPX=])+', cigar):
                print(f"Format error at line {line_index}: CIGAR incorrect.")
                sys.exit(1)

            #RNEXT
            if not re.fullmatch(r'(?:\*|=|[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)', rnext):
                print(f"Format error at line {line_index}: RNEXT incorrect.")
                sys.exit(1)

            #PNEXT: first check if of type integer then if in the right range
            try:
                int(pnext)
            except ValueError:
                print(f"Format error at line {line_index}: PNEXT must be an interger.")
                sys.exit(1)

            if not (0 <= int(pnext) <= 2**31 - 1):
                print(f"Format error at line {line_index}: pNEXT incorrect.")
                sys.exit(1)

            #TLEN: first check if of type integer then if in the right range
            try:
                int(tlen)
            except ValueError:
                print(f"Format error at line {line_index}: TLEN must be an interger.")
                sys.exit(1)

            if not (-2**31 + 1 <= int(tlen) <= 2**31 - 1):
                print(f"Format error at line {line_index}: TLEN incorrect.")
                sys.exit(1)

            #SEQ
            if not re.fullmatch(r'\*|[A-Za-z=.]+', seq):
                print(f"Format error at line {line_index}: SEQ incorrect.")
                sys.exit(1)   

            #QUAL
            if not re.fullmatch(r'[!-~]+', qual):
                print(f"Format error at line {line_index}: QUAL incorrect.")
                sys.exit(1)     



            # #Check for mandatory numeric fields
            # for value, name in [(flag, "FLAG"), (pos, "POS"), (mapq, "MAPQ"), (pnext, "PNEXT"), (tlen, "TLEN")]:
            #     if type(value) != int:
            #         print(f"Format error at line {line_index}: {name} must be an interger.")
            #         sys.exit(1)

            # #Check for mandatory string fields
            # for value, name in [(qname,"QNAME"),(rname, "RNAME"),(cigar,"CIGAR"),(rnext,"RNEXT"),(seq,"SEQ"),(qual,"QUAL")]:
            #     if type(value) != str:
            #         print(f"Format error at line {line_index}: {name} must be of type string.")
            #         sys.exit(1)

            # #for value, name in [(rname, "RNAME"), (seq, "SEQ")]:

    print("Format check OK")

## 2/ Read, 

def reader(input_file):




## 3/ Store,

## 4/ Analyse 

#### Convert the flag into binary ####
def flagBinary(flag) :

    flagB = bin(int(flag)) # Transform the integer into a binary.
    flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagB = list(flagB) 
    if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
        add = 12 - len(flagB) # We compute the difference between the maximal flag size (12) and the length of the binary flag.
        for t in range(add):
            flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.
    return flagB


#### Analyze the unmapped reads (not paired) ####
def unmapped(sam_line):
    
    unmapped_count = 0
    with open ("only_unmapped.fasta", "a+") as unmapped_fasta, open("summary_unmapped.txt", "w") as summary_file:
        for line in sam_line:
            col_line = line.split("\t")
            flag = flagBinary(col_line[1])

            if int(flag[-3]) == 1:
                unmapped_count += 1
                unmapped_fasta.write(toStringOutput(line))

        summary_file.write("Total unmapped reads: " + str(unmapped_count) + "\n") 
        return unmapped_count

#### Analyze the partially mapped reads ####
def partiallyMapped(sam_line):
    
    partially_mapped_count = 0

    with open ("only_partially_mapped.fasta", "a+") as partially_mapped_fasta, open("summary_partially_mapped.txt", "w") as summary_file:
        for line in sam_line:
            col_line = line.split("\t")
            flag = flagBinary(col_line[1]) # We compute the same 

            if int(flag[-2]) == 1: 
                if col_line[5] != "100M":
                    partially_mapped_count += 1
                    partially_mapped_fasta.write(toStringOutput(line))

        summary_file.write("Total partially mapped reads: " + str(partially_mapped_count) + "\n") 
        return partially_mapped_count


### Analyse the CIGAR = regular expression that summarise each read alignment ###
def readCigar(cigar): 
   
    ext = re.findall('\w',cigar) # split cigar 
    key=[] 
    value=[]    
    val=""

    for i in range(0,len(ext)): # For each numeric values or alpha numeric
        if (ext[i] == 'M' or ext[i] == 'I' or ext[i] == 'D' or ext[i] == 'S' or ext[i] == 'H' or ext[i] == "N" or ext[i] == 'P'or ext[i] == 'X'or ext[i] == '=') :
            key.append(ext[i])
            value.append(val)
            val = ""
        else :
            val = "" + val + ext[i]  # Else concatenate in order of arrival
    
    dico = {}
    n = 0
    for k in key:   # Dictionnary contruction in range size lists              
        if k not in dico.keys():    # for each key, insert int value
            dico[k] = int(value[n])   # if key not exist, create and add value
            n += 1
        else:
            dico[k] += int(value[n])  # inf key exist add value
            n += 1
    return dico

### Analyse the CIGAR = regular expression that summarise each read alignment ###
def percentMutation(dico):
        
    totalValue = 0 # Total number of mutations
    for v in dico :
        totalValue += dico[v]

    mutList = ['M','I','D','S','H','N','P','X','=']
    res = ""
    for mut in mutList : # Calculated percent of mutation if mut present in the dictionnary, else, percent of mut = 0
        if mut in dico.keys() :
            res += (str(round((dico[mut] * 100) / totalValue, 2)) + ";")
        else :
            res += ("0.00" + ";")
    return res

def globalPercentCigar():
    """
      Global representation of cigar distribution.
    """
    
    with open ("outpuTable_cigar.txt","r") as outpuTable, open("Final_Cigar_table.txt", "w") as FinalCigar:
        nbReads, M, I, D, S, H, N, P, X, Egal = [0 for n in range(10)]

        for line in outpuTable :
            mutValues = line.split(";")
            nbReads += 2
            M += float(mutValues[2])+float(mutValues[12])
            I += float(mutValues[3])+float(mutValues[13])
            D += float(mutValues[4])+float(mutValues[14])
            S += float(mutValues[5])+float(mutValues[15])
            H += float(mutValues[6])+float(mutValues[16])
            N += float(mutValues[7])+float(mutValues[17])
            P += float(mutValues[8])+float(mutValues[18])
            X += float(mutValues[9])+float(mutValues[19])
            Egal += float(mutValues[10])+float(mutValues[20])

        FinalCigar.write("Global cigar mutation observed :"+"\n"
                        +"Alignlent Match : "+str(round(M/nbReads,2))+"\n"
                        +"Insertion : "+str(round(I/nbReads,2))+"\n"
                        +"Deletion : "+str(round(D/nbReads,2))+"\n"
                        +"Skipped region : "+str(round(S/nbReads,2))+"\n"
                        +"Soft Clipping : "+str(round(H/nbReads,2))+"\n"
                        +"Hard Clipping : "+str(round(N/nbReads,2))+"\n"
                        +"Padding : "+str(round(P/nbReads,2))+"\n"
                        +"Sequence Match : "+str(round(Egal/nbReads,2))+"\n"
                        +"Sequence Mismatch : "+str(round(X/nbReads,2))+"\n")


 
#### Summarise the results ####

def Summary(fileName):

#### Main function ####

def main():
    if len(sys.argv) < 2:
        print("Your entry should be: python3 mapping_test.py file.sam")
        sys.exit(1)

    input_file = sys.argv[1]
    check(input_file)


############### LAUNCH THE SCRIPT ###############

if __name__ == "__main__":
    main()
