import sys
flag=99

# flagB = bin(int(flag)) # Transform the integer into a binary.
# flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
# flagB = list(flagB) 
# if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
#     add = 12 - len(flagB) # We compute the difference between the maximal flag size (12) and the length of the binary flag.
#     for t in range(add):
#         flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.
# print(flagB)




def check():
    if len(sys.argv) != 2:
        print("L'entrée doit être: python3 SamReader_template.py <sam_file>")
        sys.exit(1)
    elif sys.argv[1] == "":
    elif sys.argv[1].endswith(".sam") == False:
        print("Le fichier doit être au format .sam")
        sys.exit(1)
    else:
        print("Check OK")




sam_file = sys.argv[1]
mapping = open(sam_file, "r")
lines=mapping.readlines() #lecture de toutes les lignes du fichier
print(lines[0:5])