import sys,os

def check(input_file):
    ''' Check if the input file exists and is in .sam format '''
    if not os.path.exists(input_file):
        print("Aucun fichier trouvé")
        sys.exit(1)
    if not input_file.endswith(".sam"):
        print("Le fichier doit être au format .sam")
        sys.exit(1)
    print("Exists and format check OK")

    with open("input_file", "r") as file: #ouverture du fichier en lecture
        for line, line_index in enumerate(file,start=1):
            if line.startswith("@"):
                if line.startswith("@SQ"):
                    if "SN:" not in line or "LN:" not in line: # Check for mandatory tags
                        print(f"Erreur de format à la ligne {line_index}: @SQ doit contenir SN: et LN:.")
                        sys.exit(1)
            else:
                if line.split("\t") < 11: # Check for number of columns
                    print(f"Erreur de format à la ligne {line_index}: moins de 11 colonnes.")
                    sys.exit(1)

