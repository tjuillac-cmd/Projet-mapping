import sys,os
flag=99

# flagB = bin(int(flag)) # Transform the integer into a binary.
# flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
# flagB = list(flagB) 
# if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
#     add = 12 - len(flagB) # We compute the difference between the maximal flag size (12) and the length of the binary flag.
#     for t in range(add):
#         flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.
# print(flagB)




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
                    if "SN:" not in line or "LN:" not in line:
                        print(f"Erreur de format à la ligne {line_index}: @SQ doit contenir SN: et LN:.")
                        sys.exit(1)
            else:
                if line.split("\t") < 11:
                    print(f"Erreur de format à la ligne {line_index}: moins de 11 colonnes.")
                    sys.exit(1)




    if len(sys.argv) != 2:
        print("L'entrée doit être: python3 SamReader_template.py <sam_file>")
        sys.exit(1)
    elif sys.argv[1].endswith(".sam") == False:
        print("Le fichier doit être au format .sam")
        sys.exit(1)
    else:
        print("Check OK")




sam_file = sys.argv[1]
mapping = open(sam_file, "r")
lines=mapping.readlines() #lecture de toutes les lignes du fichier
print(lines[0:5])