############### IMPORT MODULES ###############

import sys,os,re

############### FUNCTIONS TO : ###############

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
                flag_int = int(flag)
            except ValueError:
                print(f"Format error at line {line_index}: FLAG must be an interger.")
                sys.exit(1)
            
            if not (0 <= flag_int <= 2**16 - 1):
                print(f"Format error at line {line_index}: FLAG incorrect.")
                sys.exit(1)

            #RNAME
            if not re.fullmatch(r'\*|[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*', rname):
                print(f"Format error at line {line_index}: RNAME incorrect.")
                sys.exit(1)

            #POS: first check if of type integer then if in the right range
            try:
                pos_int = int(pos)
            except ValueError:
                print(f"Format error at line {line_index}: POS must be an interger.")
                sys.exit(1)

            if not (0 <= pos_int <= 2**31 - 1):
                print(f"Format error at line {line_index}: POS incorrect.")
                sys.exit(1)

            #MAPQ: first check if of type integer then if in the right range
            try:
                mapq_int = int(mapq)
            except ValueError:
                print(f"Format error at line {line_index}: MAPQ must be an interger.")
                sys.exit(1)

            if not (0 <= mapq_int <= 2**8 - 1):
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
                pnext_int = int(pnext)
            except ValueError:
                print(f"Format error at line {line_index}: PNEXT must be an interger.")
                sys.exit(1)

            if not (0 <= pnext_int <= 2**31 - 1):
                print(f"Format error at line {line_index}: pNEXT incorrect.")
                sys.exit(1)

            #TLEN: first check if of type integer then if in the right range
            try:
                tlen_int = int(tlen)
            except ValueError:
                print(f"Format error at line {line_index}: TLEN must be an interger.")
                sys.exit(1)

            if not (-2**31 + 1 <= tlen_int <= 2**31 - 1):
                print(f"Format error at line {line_index}: TLEN incorrect.")
                sys.exit(1)

            #SEQ
            if not re.fullmatch(r'\*|[A-Za-z=.]+', seq):
                print(f"Format error at line {line_index}: SEQ incorrect.")
                sys.exit(1)   

            #QUAL
            if not re.fullmatch(r'\*|[!-~]+', qual):
                print(f"Format error at line {line_index}: QUAL incorrect.")
                sys.exit(1)     
    return True

## 2/ Read, and 3/ Store,

def sam_reader(input_file):
    '''extract useful information and store it in a dictionary'''
    reads_extract = {}
    chromosome_extract = {}

    with open(input_file, "r") as file:
        for line in file:
            if line.startswith("@"): #extract name and length of chromosomes from header, add names as keys and lengths as values
                if line.startswith("@SQ"):
                    sn_ln = parse_header(line)
                    chromosome_extract[sn_ln[0]] = sn_ln[1]
                    reads_extract[sn_ln[0]] = []
                continue
            columns = line.strip().split("\t")
            
            #extract useful fields
            qname = columns[0]
            flag = int(columns[1])
            pos = int(columns[3])
            mapq = int(columns[4])
            cigar = columns[5]

            #extract chromosome from RNAME
            chromosome = columns[2]

            #verify if chromosome key exists in dictionary, if not create it in case not described in @SQ
            if chromosome not in reads_extract:
                reads_extract[chromosome] = []

            reads_extract[chromosome].append((qname, flag, pos, mapq, cigar))

        return reads_extract, chromosome_extract

### Fuction to parse header lines to get the name and length of chromosomes from header lines ###
def parse_header(header_line):
    columns = header_line.strip().split("\t")
    sn_ln = [None, 0]
    for field in columns:
        if field.startswith("SN:"):
            sn_ln[0] = field.split(":")[1]
        elif field.startswith("LN:"):
            sn_ln[1] = int(field.split(":")[1])
    return sn_ln

## 4/ Analyse 

#### 1) Combien de reads sont mappés ? ####

#Convert the flag into binary
def flagBinary(flag) :

    flagB = bin(int(flag)) # Transform the integer into a binary.
    flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagB = list(flagB) 
    if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
        add = 12 - len(flagB) # We compute the difference between the maximal flag size (12) and the length of the binary flag.
        for t in range(add):
            flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.
    return flagB

#compter le nombre de reads en fonction du flag (colonne #2) -> look at bit x4 of flag binary
def readMapped(reads_extract):
    count_mapped = {'mapped': 0, 'unmapped': 0, 'total': 0}
    for chromosome in reads_extract:
        for read in reads_extract[chromosome]:
            flag = read[1]
            flagB = flagBinary(flag)
            count_mapped['total'] += 1
            if int(flagB[-3]) == 0: # check if the read is mapped, if the flag has the bit 4 it means it is unmapped
                count_mapped['mapped'] += 1
            else:
                count_mapped['unmapped'] += 1
    print(count_mapped)
    return count_mapped

#### 2) Comment les reads (et paires de reads) sont-ils mappés ?  # compter le nombre de reads pour chaque flag ####

#function to extract and count flags in a dictionary
def readFlag(reads_extract):
    count_flag = {}
    for chromosome in reads_extract:
        for read in reads_extract[chromosome]:
            if read[1] in count_flag:
                count_flag[read[1]] += 1
            else:
                count_flag[read[1]] = count_flag.get(read[1],0) + 1
    print(count_flag)
    return count_flag

#### 3) Où les reads sont-ils mappés ? L'alignement est-il homogène le long de la séquence de référence ? # compter le nombre de reads par chromosome ####

def readCHROM(reads_extract):
    count_chrom= {key: [0, 0] for key in reads_extract.keys()} #create a counting dico with same keys as reads_extract 
    for chromosome in reads_extract:
        for read in reads_extract[chromosome]:
            flagB = flagBinary(read[1])
            if int(flagB[-3]) == 0: # check if the read is mapped, if the flag has the bit 4 it means it is unmapped
                count_chrom[chromosome][0] += 1 #count_chrom[chromosome][0] is the number of mapped reads
            else:
                count_chrom[chromosome][1] +=1 #count_chrom[chromosome][1] is the number of unmapped reads
    print(count_chrom)
    return count_chrom

# on cherche la position de chaque read sur le chromosome en utilisant la position de départ (POS) et la longueur du CIGAR

def positionsReads(reads_extract):
    positions = {}
    for chromosome in reads_extract:
        for read in reads_extract[chromosome]:
            
            #calculate the position of the read on the reference
            pos, end = read[2], 0        
            cigar = read[4]
            length = lengthRefCigar(cigar)
            end = pos + length - 1

            if chromosome not in positions:
                positions[chromosome] = []
                positions[chromosome].append((pos, end))
            else:
                positions[chromosome].append((pos, end))
    print(positions)
    return positions

            



### Analyse the CIGAR = calculate length consumed on reference ###
def lengthRefCigar(cigar): 
   
    ext = re.findall(r"(\d+)([MIDNSHPX=])",cigar) # split cigar on alpha numeric and numeric, ex. ext = ['1','0','0','M','5','S']
    length = 0

    for nb, op in ext: # for each pair number-operation

        if op in ("M", "D", "N", "X", "="): # these operations consume reference
            length += int(nb)

    return length


#### 4) Avec quelle qualité les reads sont-ils mappés ? # compter le nombre de reads pour chaque valeur de qualité ou par tranche de valeurs (score de mapping) ####

def readMAPQ(reads_extract,threshold):
    above = f"MAPQ above {threshold}"
    below = f"MAPQ below {threshold}"
    count_mapq = {above : 0, below : 0}
    for chromosome in reads_extract:
        for read in reads_extract[chromosome]:
            if read[3] >= threshold :
                count_mapq[above] += 1
            else:
                count_mapq[below] += 1
    print(count_mapq)
    return count_mapq


 
#### Summarise the results ####

def Summary(fileName, count_flag, count_chrom, count_mapped, count_mapq):
    '''create a text file to summarize the results'''
    with open(fileName, "w") as fileSummary: #open file in write mode
        fileSummary.write("============ Summary of SAM file ============\n\n")

        #1) output mapped reads, unmapped reads and total reads
        fileSummary.write("Reads mapping information:\n")
        for key, value in count_mapped.items():
            fileSummary.write(f"{key} reads: {value}\n") #write mapped, unmapped and total reads
        fileSummary.write("\n=============================================\n")

        #2) output reads per flag
        fileSummary.write("Reads per flag:\n")
        for flag, count in count_flag.items():
            fileSummary.write(f"{flag} : {count} reads\n") #write each flag and its count
        fileSummary.write("\n=============================================\n")
        
        #3) output reads per chromosome
        fileSummary.write("Reads per chromosome:\n")
        for chrom, counts in count_chrom.items():
            fileSummary.write(f"{chrom} : {counts[0]} mapped reads, {counts[1]} unmapped reads\n") #write each chromosome and its mapped/unmapped counts
        fileSummary.write("\n=============================================\n")

        #4) output reads per MAPQ
        fileSummary.write("Summary of reads repartition per mapping quality:\n")
        for key, value in count_mapq.items():
            fileSummary.write(f"{key} : {value} reads\n") #write MAPQ summary
        fileSummary.write("\n=============================================\n")

        print(f"Summary written in {fileName}")

#### Main function ####

def main():
    if len(sys.argv) < 2:
        print("Your entry should be: python3 mapping_test.py file.sam")
        sys.exit(1)

    input_file = sys.argv[1]
    if check(input_file):
        print("Format check OK")

        reads_extract = sam_reader(input_file)[0]
        count_mapped = readMapped(reads_extract)
        count_flag = readFlag(reads_extract)
        count_chrom = readCHROM(reads_extract)
        count_mapq = readMAPQ(reads_extract,36)

        Summary("summary.txt", count_flag, count_chrom, count_mapped, count_mapq)

        os.system("cat summary.txt")
            

############### LAUNCH THE SCRIPT ###############

if __name__ == "__main__":
    main()
