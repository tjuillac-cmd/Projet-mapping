import sys,os

def check(input_file):
    '''
    Check if the input file exists and is in .sam format.
    Also verifies that each alignment line (non-header) contains at least 11 columns as required by the SAM format specification.
    '''
    if not os.path.exists(input_file):
        print("Aucun fichier trouvé")
        sys.exit(1)
        print("Aucun fichier trouvé")
        sys.exit(1)
    if not input_file.endswith(".sam"):
        print("Le fichier doit être au format .sam")
        sys.exit(1)
    print("Exists and format check OK")

    with open(input_file, "r") as file: #ouverture du fichier en lecture
        for line_index, line in enumerate(file,start=1):
            if line.startswith("@"):
                if line.startswith("@SQ"):
                    if "SN:" not in line or "LN:" not in line: # Check for mandatory tags
                        print(f"Erreur de format à la ligne {line_index}: @SQ doit contenir SN: et LN:.")
                        sys.exit(1)
                continue

            # Process alignment lines
            columns = line.strip().split("\t")

            # Check for at least 11 mandatory fields
            if len(columns) < 11:
                print(f"Erreur de format à la ligne {line_index}: moins de 11 champs (seulement {len(columns)}).")
                sys.exit(1)

            # Unpack the first 11 mandatory fields
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = columns[:11]

            # Check for empty mandatory fields
            if any(field == "" for field in columns[:11]):
                print(f"Erreur de format à la ligne {line_index}: un des 11 champs obligatoires est vide.")
                sys.exit(1)

            # Check for mandatory numeric fields
            for value, name in [(flag, "FLAG"), (pos, "POS"), (mapq, "MAPQ"), (pnext, "PNEXT"), (tlen, "TLEN")]:
                try:
                    int(value)
                except ValueError:
                    print(f"Erreur de format à la ligne {line_index}: {name} doit être un entier.")
                    sys.exit(1)

            # Check for mandatory string fields
            for value, name in [(rname, "RNAME"), (seq, "SEQ")]:
                


            



    print("Format check OK")

def main():
    if len(sys.argv) < 2:
        print("Usage : python3 mapping_test.py fichier.sam")
        sys.exit(1)

    input_file = sys.argv[1]
    check(input_file)


if __name__ == "__main__":
    main()
