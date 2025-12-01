import sys,os

def check(input_file):
    '''
    Check if the input file exists and is in .sam format.
    Verifies that each alignment line (non-header) contains at least the 11 "mandatory" fields
    '''
    if not os.path.exists(input_file):
        print("No file found")
        sys.exit(1)
    if not input_file.endswith(".sam"):
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

            # Check for mandatory numeric fields
            for value, name in [(flag, "FLAG"), (pos, "POS"), (mapq, "MAPQ"), (pnext, "PNEXT"), (tlen, "TLEN")]:
                try:
                    int(value)
                except ValueError:
                    print(f"Format error at line {line_index}: {name} must be an interger.")
                    sys.exit(1)

            # # Check for mandatory string fields
            # for value, name in [(rname, "RNAME"), (seq, "SEQ")]:
                


            



    print("Format check OK")

def main():
    if len(sys.argv) < 2:
        print("Usage : python3 mapping_test.py fichier.sam")
        sys.exit(1)

    input_file = sys.argv[1]
    check(input_file)


#if __name__ == "__main__":
#    main()
