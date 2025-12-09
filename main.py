############### IMPORT MODULES ###############

import sys,os,re, matplotlib.pyplot as plt, matplotlib.cm as cm, matplotlib.colors as colors

############### TOOL FUNCTIONS ###############

def isFullyMapped(flag, cigar):
    '''determine if a read is mapped based on flag and cigar'''
    
    #case unmapped
    if flag & 4 == 1: #unmapped
        return False
    
    if cigar == "*" or cigar == None: #unmapped
        return False 
    
    ops = re.findall(r"[MIDNSHPX=]", cigar) # extract operations from CIGAR string

    #case partially mapped
    if "S" in ops or "H" in ops:
        return False    
    
    #case mapped
    if all(op in ["M", "D", "N", "X", "="] for op in ops):
        return True


def lengthRefCigar(cigar): 
    '''calculate the length consumed on the reference sequence based on CIGAR'''
   
    ext = re.findall(r"(\d+)([MIDNSHPX=])",cigar) # split cigar on alpha numeric and numeric, ex. ext = ['100','M','5','S']
    length = 0

    for nb, op in ext: # for each pair number-operation

        if op in ("M", "D", "N", "X", "="): # these operations consume reference
            length += int(nb)

    return length

def nbIndel(cigar):
    '''calculate the number of indel in read based on CIGAR'''

    ext = re.findall(r"(\d+)([MIDNSHPX=])",cigar) # split cigar on alpha numeric and numeric, ex. ext = ['100','M','5','S']
    indel_count = 0

    for nb, op in ext:

        if op in ("I", "D"):
            indel_count += int(nb)

    return indel_count

############### SAM FILE CHECK FUNCTION ###############

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
    print("File exists and format check OK")

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


############### READ FILTERING AND EXTRACTION FUNCTIONS ###############

def sam_reader(input_file, header_parsed, filterMAPQ, fullyMappedOnly):
    '''extract useful information and store it in a dictionary'''
    reads_extract = {chrom: [] for chrom in header_parsed.keys()}

    with open(input_file, "r") as file:
        for line in file:
            if line.startswith("@"): #skip header lines
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

            if fullyMappedOnly:
                #filter fully mapped reads only
                if not isFullyMapped(flag, cigar):
                    continue
            
            if filterMAPQ is not None and mapq < filterMAPQ:
                continue

            if chromosome not in reads_extract:
                reads_extract[chromosome] = []            

            reads_extract[chromosome].append((qname, flag, pos, mapq, cigar))

    return reads_extract


def parse_header(input_file):
    '''parse the header of the SAM file to get the length of each reference sequence {reference_name: length}'''
    length_ref = {}
    with open(input_file, "r") as file:
        for line in file:
            if not line.startswith("@"):
                break
            if line.startswith("@SQ"):
                columns = line.strip().split("\t")
                for field in columns:
                    if field.startswith("SN:"):
                        sn = field.split(":")[1]
                    elif field.startswith("LN:"):
                        ln = int(field.split(":")[1])
                length_ref[sn] = ln
    return length_ref


################ STATISTICS FUNCTIONS ###############

def readMapped(reads_extract):
    '''count the number of mapped, unmapped and total reads {'mapped': x, 'unmapped': y, 'total': z}'''
    count_mapped = {'mapped': 0, 'unmapped': 0, 'total': 0}
    for chromosome in reads_extract:
        for read in reads_extract[chromosome]:
            flag = read[1]
            count_mapped['total'] += 1
            if flag & 4 == 0: # check if the read is mapped, if the flag has the bit 4 it means it is unmapped
                count_mapped['mapped'] += 1
            else:
                count_mapped['unmapped'] += 1
    return count_mapped


def readFlag(reads_extract):
    '''count the number of reads per flag {flag: count}'''
    readsByName = {}
    paired_orientation = {} #count percentage of properly paired reads and of properly oriented pairs

    for chromosome in reads_extract:
        for read in reads_extract[chromosome]:

            qname = read[0]
            flag = read[1]

            if qname not in readsByName: #we fill the dictionnary with flags for each qname (to see properly paired reads)
                readsByName[qname] = []
            readsByName[qname].append(flag)

        total, properly_paired, properly_oriented = 0, 0, 0

        for qname, flags in readsByName.items():
            total += 1

            if len(flags) != 2: #skip unpaired reads
                continue

            f1, f2 = flags

            #check if properly paired (bit 0x40 et 0x80 in either f1 or f2 and 0x2 in both) 
            if ((f1 & 0x80 and f2 & 0x40) or (f2 & 0x80 and f1 & 0x40)) and f1 & 0x2 and f2 & 0x2:
                properly_paired += 1

            #check if properly oriented (RF or FR but FF or RR are misoriented)
            #bit 0x10 indicates R
            if f1 & 0x10 and not f2 & 0x10:
                properly_oriented += 1
                    
            elif f2 & 0x10 and not f1 & 0x10:
                properly_oriented += 1
            
            perc_properly_paired = round(100*properly_paired / total, 2)
            perc_properly_oriented = round(100*properly_oriented / total, 2)
            paired_orientation[chromosome] = [perc_properly_paired, perc_properly_oriented]
        
    return paired_orientation


def readCHROM(reads_extract):
    '''count the number of mapped and unmapped reads per chromosome {chromosome: [mapped, unmapped]}'''
    count_chrom= {key: [0, 0] for key in reads_extract.keys()} #create a counting dico with same keys as reads_extract 
    for chromosome in reads_extract:
        for read in reads_extract[chromosome]:
            if read[1] & 4 == 0: # check if the read is mapped, if the flag has the bit 4 it means it is unmapped
                count_chrom[chromosome][0] += 1 #count_chrom[chromosome][0] is the number of mapped reads
            else:
                count_chrom[chromosome][1] +=1 #count_chrom[chromosome][1] is the number of unmapped reads
    return count_chrom


def readMAPQ(reads_extract, MAPQ_threshold):
    '''count the number of reads per MAPQ {MAPQ above threshold, MAPQ below threshold}'''

    count_mapq = {chromosome: [0, 0] for chromosome in reads_extract.keys()} #dict {chromosome: [nb mapped, nb unmapped]}

    for chromosome in reads_extract:
        for read in reads_extract[chromosome]:
            mapq = read[3]

            if mapq >= MAPQ_threshold :
                count_mapq[chromosome][0] += 1
            else:
                count_mapq[chromosome][1] += 1

    return count_mapq


def statAlignment(reads_extract, short_size, large_size):
    '''return basic statistics on alignment length'''
    lengths = {key: [] for key in reads_extract.keys()}
    stats = {key: [] for key in reads_extract.keys()}

    for chromosome in reads_extract:
        for read in reads_extract[chromosome]:
            cigar = read[4]
            length = lengthRefCigar(cigar)
            lengths[chromosome].append(length)

    for chromosome in lengths:
        under50, over250, mean, total, min_len, max_len = 0, 0, 0, len(lengths[chromosome]), min(lengths[chromosome]), max(lengths[chromosome])
        
        for length in lengths[chromosome]:
            mean += length

            if length <= short_size:
                under50 += 1

            if length >= large_size:
                over250 += 1
        
        mean = round(mean / total, 3)
        stats[chromosome] = (under50, over250, mean, total, min_len, max_len)

    return stats

def statIndel(reads_extract):
    '''calculate percentage of reads w/ at least one indel'''
    indel_dict = {key: 0 for key in reads_extract.keys()}

    for chromosome in reads_extract:
        total, atLeastOne = 0, 0
        for read in reads_extract[chromosome]:
            
            total += 1
            cigar = read[4]
            nb_indel = nbIndel(cigar)

            if nb_indel >= 1:
                atLeastOne += 1
        
        indel_dict[chromosome] += round(atLeastOne / total, 2)
    
    return indel_dict

def Summary(fileName, paired_orientation, count_chrom, count_mapped, count_mapq, stat_alignment, short_size, long_size, MAPQ_threshold, stat_indel):
    '''create a text file to summarize the results'''
    with open(fileName, "w") as fileSummary: #open file in write mode
        fileSummary.write("============ Summary of SAM file ============\n\n")

        #table header
        fileSummary.write(f"CHR_NAME\tTOT_READS\tMAP\tUMAP\tMAPQ-\tMAPQ+\tPAIR%\tRF%\t")
        fileSummary.write(f"<{short_size}BP%\tINT%\t>{long_size}BP%\tMEANL\tMINL\tMAXL\tINDEL%\n")

        #CHROM_NAME
        for chromosome in count_chrom.keys():
            fileSummary.write(f"{chromosome}\t")

            #TOTAL_READS, MAPPED_READS, UNMAPPED_READS
            total_reads = count_chrom[chromosome][0] + count_chrom[chromosome][1]
            mapped_reads = count_chrom[chromosome][0]
            unmapped_reads = count_chrom[chromosome][1]
            fileSummary.write(f"{total_reads}\t{mapped_reads}\t{unmapped_reads}\t")

            #MAPQ<={MAPQ_threshold}, MAPQ>{MAPQ_threshold}
            mapq_below = count_mapq[chromosome][1]
            mapq_above = count_mapq[chromosome][0]
            fileSummary.write(f"{mapq_below}\t{mapq_above}\t")

            #%_PROP_PAIRED, %_PROP_ORIENTED
            perc_properly_paired = paired_orientation[chromosome][0]
            perc_properly_oriented = paired_orientation[chromosome][1]
            fileSummary.write(f"{perc_properly_paired}\t{perc_properly_oriented}\t")

            #<{short_size}BP%, [{short_size};{long_size}BP]%, >{long_size}BP%, MEAN_LENGTH, MIN_LENGTH, MAX_LENGTH
            under = stat_alignment[chromosome][0]
            over = stat_alignment[chromosome][1]
            mean_len = stat_alignment[chromosome][2]
            total_len = stat_alignment[chromosome][3]
            min_len = stat_alignment[chromosome][4]
            max_len = stat_alignment[chromosome][5]

            fileSummary.write(f"{round(under/total_len,2)*100}%\t")
            fileSummary.write(f"{100-(round(over/total_len,2)+round(under/total_len,2))*100}%\t")
            fileSummary.write(f"{round(over/total_len,2)*100}%\t")         
            fileSummary.write(f"{mean_len}\t") 
            fileSummary.write(f"{min_len}\t") 
            fileSummary.write(f"{max_len}\t")

            #percentage of reads w/ at least one indel
            fileSummary.write(f"{stat_indel[chromosome]}\n")

        fileSummary.write(f"\nLEGEND:\nCHR_NAME: Chromosome name\nTOT_READS: Total reads\nMAP: Mapped reads\nUMAP: Unmapped reads\n")
        fileSummary.write(f"MAPQ-: Reads with MAPQ less than or equal to {MAPQ_threshold}\nMAPQ+: Reads with MAPQ greater than {MAPQ_threshold}\n")
        fileSummary.write(f"PAIR%: Percentage of properly paired reads\nRF%: Percentage of properly oriented pairs (forward-reverse or reverse-forward)\n")
        fileSummary.write(f"<{short_size}BP%: Percentage of reads with length less than {short_size} bp\n")
        fileSummary.write(f"INT%: Percentage of reads with intermediate length (between {short_size} and {long_size} bp)\n")
        fileSummary.write(f">{long_size}BP%: Percentage of reads with length greater than {long_size} bp\n")
        fileSummary.write(f"MEANL: Mean length of alignments\nMINL: Minimum length of alignments\nMAXL: Maximum length of alignments\n")
        fileSummary.write(f"INDEL%: Percentage of reads with at least one indel\n\n")


################ STATISTICS ON FILTERED DATA ###############

def positionsReads(reads_extract):
    '''calculate the positions of each mapped read on the reference sequence for each chromosome {[(start1, end1),(start2, end2)...]} and MAPQ'''
    positions = {}
    for chromosome in reads_extract:
        
        if not chromosome == "*": #skip chromosome '*'
            
            for read in reads_extract[chromosome]:
                
                if read[1] & 4 == 0: # only mapped reads

                    mapq = read[3]

                    #calculate the position of the read on the reference
                    start, end = read[2], 0        
                    cigar = read[4]
                    length = lengthRefCigar(cigar)
                    end = start + length - 1

                    if chromosome not in positions:
                        positions[chromosome] = []
                        positions[chromosome].append((start, end, mapq))
                    else:
                        positions[chromosome].append((start, end, mapq))

    return positions


def readsPerWindow(positions, header_parsed, window_size):
    '''calculate the niumber of reads per window on each reference'''
    reads_window = {chrom: [] for chrom in positions.keys()}

    for chrom, reads in positions.items():

        if not reads: #case chromosome has no reads
            continue
        
        length_ref = header_parsed[chrom]

        nb_windows = (length_ref // window_size) + 1 # calculate the number of windows needed to cover the reference
        windows_counts = [0.0] * nb_windows # initialize a list to count reads per window

        # for each read we determine the range of windows it is on
        for start, end, mapq in reads:
            first_window = start // window_size
            last_window = end // window_size
            for window in range(first_window, last_window + 1):
                '''we add to each window's count 1 if fully covered or coverage ratio if not entirely covered (for first and last)'''
                if window == first_window:
                    overlap = (first_window + 1) * window_size - start                    
                    windows_counts[window] += overlap / window_size
                elif window == last_window:
                    overlap = end - last_window * window_size                    
                    windows_counts[window] += overlap / window_size
                else:    
                    windows_counts[window] += 1.0
        
        #round up and add to the dictionnary of repartition per chromosome
        windows_counts = [round(c, 3) for c in windows_counts]
        reads_window[chrom] = windows_counts
    
    return reads_window


def meanMAPQPerWindow(positions, header_parsed, window_size, MAPQ_threshold):
    '''calculate the mean MAPQ per window on each reference'''
    mapq_window = {chrom: [] for chrom in positions.keys()}

    for chrom, reads in positions.items():

        if not reads: #case chromosome has no reads
            continue

        length_ref = header_parsed[chrom]
        nb_windows = (length_ref // window_size) + 1 # calculate the number of windows needed to cover the reference
        windows_mapq = [[] for _ in range(nb_windows)] # initialize a list to store MAPQ per window

        # for each read we determine the range of windows it is on
        for start, end, mapq in reads:

            if mapq < MAPQ_threshold: #skip mapq values below threshold
                continue

            first_window = start // window_size
            last_window = end // window_size
            for window in range(first_window, last_window + 1):

                windows_mapq[window].append(mapq)
        
        #we calculate the mean MAPQ per window and we add to the dictionnary of repartition per chromosome
        mean_mapq_counts = []
        for mapqs in windows_mapq:
            if mapqs:
                mean_mapq = sum(mapqs) / len(mapqs)
                mean_mapq_counts.append(round(mean_mapq, 3))
            else:
                mean_mapq_counts.append(0.0)
        
        mapq_window[chrom] = mean_mapq_counts
    
    return mapq_window


################ PLOTTING FUNCTION ###############

def plotReadsPerWindow(reads_window, mapq_window, window_size):
    '''plot the number of reads per window on each chromosome, colored by mean MAPQ'''

    for chrom, counts in reads_window.items():
        if not counts: #case chromosome has no reads
            continue
        
        mapq_values = mapq_window[chrom]

        # Normalize MAPQ values for color mapping
        min_mapq = min(mapq_values)
        max_mapq = max(mapq_values)

        if max_mapq == min_mapq:
            norm_mapq = [0.5 for x in mapq_values]  # all same color if no variation
        else:
            norm_mapq = [(mapq - min_mapq) / (max_mapq - min_mapq) for mapq in mapq_values]

        #create color map
        colormap = cm.get_cmap('RdYlGn')
        colors_mapped = [colormap(norm) for norm in norm_mapq]

        # plot
        fig, ax = plt.subplots(figsize=(10, 5)) 

        ax.bar(range(len(counts)), counts, width=1.0, color = colors_mapped, edgecolor='none') #bar plot with colored bars red to green
        ax.set_xlabel(f'Windows of size {window_size} bp along {chrom}')
        ax.set_ylabel('Number of reads')
        ax.set_title(f'Read distribution along {chrom} (colored by mean MAPQ)')
        ax.set_xticks(ticks=range(0, len(counts), max(1, len(counts)//10)), 
                   labels=[str(i * window_size) for i in range(0, len(counts), max(1, len(counts)//10))])
        
        #colorbar pour MAPQ
        norm = colors.Normalize(vmin=min_mapq, vmax=max_mapq)
        sm = cm.ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax = ax) #colorbar to show MAPQ scale
        cbar.set_label('Mean MAPQ per window')
        
        plt.grid(axis='y')
        plt.tight_layout()
        plt.show()
 

################ MAIN FUNCTION ###############

def main():
    if len(sys.argv) < 2:
        print("Your entry should be: python3 mapping_test.py file.sam")
        sys.exit(1)

    input_file = sys.argv[1]
    
    ## Check input file ##
    check(input_file)

    ## User inputs ##

    # MAPQ threshold (optional) #
    filterInput = input("Enter a MAPQ threshold to filter reads or press Enter to skip (default 0): ")
    if filterInput == "":
        filterMAPQ = None
    else:
        try:
            filterMAPQ = int(filterInput)
            if not (0 <= filterMAPQ <= 60):
                print("MAPQ threshold must be between 0 and 60.")
                sys.exit(1)
        except ValueError:
            print("MAPQ threshold must be an integer.")
            sys.exit(1)
    
    # Fully mapped reads only (mandatory) #
    fullyMappedInput = input("Do you want to consider only fully mapped reads? (yes/NO): ")
    if fullyMappedInput.lower() in ["yes", "y", "oui", "o", "true", "t"]:
        fullyMappedOnly = True
    else:
        fullyMappedOnly = False

    
    ## Function on SAM file ##
    header_parsed = parse_header(input_file)

    # Window size for read distribution (mandatory) #
    window_size_input = input("Enter window size for read distribution (default 1000): ")
    if window_size_input == "":
        window_size = 1000
    else:
        try:
            window_size = int(window_size_input)
            if window_size <= 0:
                window_size = int(input("Window size must be a positive integer."))
            
            max_size = max(header_parsed.values())/2
            if window_size > max_size:
                window_size = int(input(f"Window size must be smaller than {max_size}."))
                
        except ValueError:
            print("Window size must be an integer.")
            sys.exit(1)

    # Sizes for short and small reads

    short_size_input = input("Enter a threshold size for small reads or press ENTER (default 80): ")
    if short_size_input == "":
        short_size = 80
    else:
        try:
            short_size = int(short_size_input)
        except ValueError:
            print("Size must be an integer.")
            sys.exit(1)

    long_size_input = input("Enter a threshold size for long reads or press ENTER (default 200): ")
    if long_size_input == "":
        long_size = 200
    else:
        try:
            long_size = int(long_size_input)
        except ValueError:
            print("Size must be an integer.")
            sys.exit(1)

    
    ## Analysis of filtered data ##
    reads_extract = sam_reader(input_file, header_parsed, filterMAPQ, fullyMappedOnly) #reads with user filtering

    if filterMAPQ is None:
        MAPQ_threshold = 0 #default threshold
    else:
        MAPQ_threshold = filterMAPQ
           
    positions = positionsReads(reads_extract)
    reads_window = readsPerWindow(positions, header_parsed, window_size)
    mapq_window = meanMAPQPerWindow(positions, header_parsed, window_size, MAPQ_threshold)
    stat_alignment = statAlignment(reads_extract, short_size, long_size)
    stat_indel = statIndel(reads_extract)
    count_mapped = readMapped(reads_extract)
    paired_orientation = readFlag(reads_extract)
    count_chrom = readCHROM(reads_extract)
    count_mapq = readMAPQ(reads_extract, MAPQ_threshold)

    Summary("summary.txt", paired_orientation, count_chrom, count_mapped, count_mapq, stat_alignment, short_size, long_size, MAPQ_threshold, stat_indel)
    plotReadsPerWindow(reads_window, mapq_window, window_size)


    os.system("cat summary.txt")
            

############### LAUNCH THE SCRIPT ###############

if __name__ == "__main__":
    main()
