import sys,os,re

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

            # Process alignment lines
            columns = line.strip().split("\t")

            # Check for at least 11 mandatory fields
            if len(columns) < 11:
                print(f"Format error at line {line_index}: less than 11 fields.")
                sys.exit(1)

            # Unpack the first 11 mandatory fields
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = columns[:11]

            # Check for empty mandatory fields
            if any(field == "" for field in columns[:11]):
                print(f"Format error at line {line_index}: one of the 11 mandatory fields is empty.")
                sys.exit(1)

            # Check if each field is in correct Regexp/Range
            # QNAME
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



            # # Check for mandatory numeric fields
            # for value, name in [(flag, "FLAG"), (pos, "POS"), (mapq, "MAPQ"), (pnext, "PNEXT"), (tlen, "TLEN")]:
            #     if type(value) != int:
            #         print(f"Format error at line {line_index}: {name} must be an interger.")
            #         sys.exit(1)

            # # Check for mandatory string fields
            # for value, name in [(qname,"QNAME"),(rname, "RNAME"),(cigar,"CIGAR"),(rnext,"RNEXT"),(seq,"SEQ"),(qual,"QUAL")]:
            #     if type(value) != str:
            #         print(f"Format error at line {line_index}: {name} must be of type string.")
            #         sys.exit(1)

            # # for value, name in [(rname, "RNAME"), (seq, "SEQ")]:

    print("Format check OK")

def main():
    if len(sys.argv) < 2:
        print("Usage : python3 mapping_test.py fichier.sam")
        sys.exit(1)

    input_file = sys.argv[1]
    check(input_file)


if __name__ == "__main__":
    main()
