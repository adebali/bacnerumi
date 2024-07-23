import sys

# Check if a filename is provided as a command-line argument
if len(sys.argv) != 3:
    print("Usage: python3 bedpe2bed.py input_file input_sam")
    sys.exit(1)

input_file = sys.argv[1]
input_sam = sys.argv[2]

d = {}

with open(input_sam, 'r') as file:
    for line in file:
        if line[0] != "@":
            ll = line.strip().split('\t')
            if int(ll[1]) in [147, 99]:
                d[ll[0]] = "-"
            elif int(ll[1]) in [163, 83]:
                d[ll[0]] = "+"
            else:
                d[ll[0]] = "."
    
with open(input_file, 'r') as file:
    for line in file:
        ll = line.strip().split('\t')
        name = ll[6]
        strand = d[name]
        print(f"{ll[0]}\t{ll[1]}\t{ll[5]}\t{name}\t.\t{strand}")
