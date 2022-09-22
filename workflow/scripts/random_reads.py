import random
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, help="singleton fasta file")
parser.add_argument("-l", "--length", required=True, help="read length file")
parser.add_argument("-k", "--readlength", required=True, nargs="+", type=int, help="read length list")
parser.add_argument("-o", "--output", default=sys.stdout,type=argparse.FileType('w'), help="FASTQ file with removed duplicates")
args = parser.parse_args()
out = args.output
readlengthList = args.readlength
readLengthSet = set()
for l in readlengthList:
    readLengthSet.add(l)
readLengthDict = {}
lengthDistFile = args.length
lengthIn = open(lengthDistFile, 'r')
for line in lengthIn:
    ll = line.rstrip().split('\t')
    readLength = int(ll[0])
    value = int(ll[1])
    if readLength in readlengthList:
        readLengthDict[readLength] = value

seqList = []
headerList = []
singletonCount = 0

for line in open(args.input, 'r'):
    if line.startswith('>'):
        singletonCount += 1
        header = line.rstrip()
    else:
        seq = line.rstrip()
        headerList.append(header[1:])
        seqList.append(seq)

countDict = {}
for readLength in readLengthDict.keys():
    countDict[readLength] = 0

randomNoList = []
print(readLengthDict)
print(singletonCount)
random.seed(10)

while True:
    if len(seq) in readLengthSet:
        no = random.randint(0, singletonCount-1)
        seq = seqList[no]
        countDict[len(seq)] = countDict.get(len(seq), 0) + 1
        randomNoList.append(no)
        print(countDict)
        if countDict[len(seq)] == readLengthDict[readLength]:
            readLengthSet.remove(readLength)
            if len(readLengthSet) == 0:
                break

for no in randomNoList:
    header = headerList[no]
    hl = header.split('_')
    chromosome = hl[0]
    start = hl[1]
    end = hl[2]
    strand = hl[3]
    lineList = [chromosome, start, end, '.', '.', strand]
    out.write('\t'.join(lineList) + '\n')