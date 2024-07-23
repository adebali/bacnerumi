import random
import argparse
import sys
import pybedtools
import pandas as pd
from functools import reduce

tempdir = "./tmp"
pybedtools.helpers.set_tempdir(tempdir)

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, help="singleton fasta (.fa) or (.bed) file")
parser.add_argument("-b", "--bed", required=True, help="bed file to get the read lengths")
parser.add_argument("-g", "--genes", required=True, help="bed file to get the gene intervals")
parser.add_argument("-o", "--output", default=sys.stdout,type=argparse.FileType('w'), help="bed output file")
args = parser.parse_args()
out = args.output

separator = '|'

bedIn = open(args.bed, 'r')
readLengthDict = {}
readLengthSet = set()
readCount = 0
for line in bedIn:
    ll = line.strip().split('\t')
    readLength = int(ll[2]) - int(ll[1])
    readLengthDict[readLength] = readLengthDict.get(readLength, 0) + 1
    readLengthSet.add(readLength)
    readCount += 1


bedIn.close()
masterDict = {}

if args.input.endswith('.fa'):
    fastaIn = open(args.input, 'r')
    for line in fastaIn:
        header = line.strip()[1:]
        seq = fastaIn.readline().strip()
        seqLen = len(seq)
        masterDict[seqLen] = masterDict.get(seqLen, [])
        masterDict[seqLen].append(header)
    fastaIn.close()
elif args.input.endswith('.bed'):
    for line in open(args.input, 'r'):
        hl = line.rstrip().split('\t')
        header = separator.join([hl[0], hl[1], hl[2], hl[5]])
        seqLen = int(hl[2]) - int(hl[1])
        masterDict[seqLen] = masterDict.get(seqLen, [])
        masterDict[seqLen].append(header)
else:
    assert "Unexpected file input"

# if args.input.endswith('.fa'):
#     for line in open(args.input, 'r'):
#         if line.startswith('>'):
#             singletonCount += 1
#             header = line.rstrip()
#         else:
#             seq = line.rstrip()
#             headerList.append(header[1:])
#             seqList.append(seq)
# elif args.input.endswith('.bed'):
#     for line in open(args.input, 'r'):
#         singletonCount += 1
#         hl = line.rstrip().split('\t')
#         header = separator.join([hl[0], hl[1], hl[2], hl[5]])
#         headerList.append(header)


a = pybedtools.BedTool(args.genes)
TS = []
NTS = []

for i in range(0,5):
    selectedHeaders = []
    for readLength in readLengthDict.keys():
        random.seed(int(readCount/readLength)+i)
        headerList = masterDict[readLength]
        selectedHeaders += random.choices(headerList, k=readLengthDict[readLength])
        
    intervals = []
    for header in selectedHeaders:
        hl = header.split(separator)
        chromosome = hl[0]
        start = int(hl[1])
        end = int(hl[2])
        strand = hl[3]
        intervals.append(pybedtools.Interval(chromosome, start, end, ".", ".", strand))

    b = pybedtools.BedTool(intervals)

    TS_intersect = a.intersect(b, S=True, c=True)
    TS.append(TS_intersect.to_dataframe())
    
    NTS_intersect = a.intersect(b, s=True, c=True)
    NTS.append(NTS_intersect.to_dataframe())


key_columns = ["chrom", "start", "end", "name", "score", "strand"]

for i, df in enumerate(TS):
    suffix = f'_TS{i+1}'
    df.columns = [col if col in key_columns else col + suffix for col in df.columns]

for i, df in enumerate(NTS):
    suffix = f'_NTS{i+1}'
    df.columns = [col if col in key_columns else col + suffix for col in df.columns]

TSm = reduce(lambda left, right: pd.merge(left, right, on=["chrom", "start", "end", "name", "score", "strand"], how='inner'), TS)
TSm['TSmed'] = TSm[['thickStart_TS1', 'thickStart_TS2', 'thickStart_TS3', 'thickStart_TS4', 'thickStart_TS5']].median(axis=1) 
TSmd = TSm.drop(['thickStart_TS1', 'thickStart_TS2', 'thickStart_TS3', 'thickStart_TS4', 'thickStart_TS5'], axis=1)

NTSm = reduce(lambda left, right: pd.merge(left, right, on=["chrom", "start", "end", "name", "score", "strand"], how='inner'), NTS)
NTSm['TSmed'] = NTSm[['thickStart_NTS1', 'thickStart_NTS2', 'thickStart_NTS3', 'thickStart_NTS4', 'thickStart_NTS5']].median(axis=1) 
NTSmd = NTSm.drop(['thickStart_NTS1', 'thickStart_NTS2', 'thickStart_NTS3', 'thickStart_NTS4', 'thickStart_NTS5'], axis=1)

TSNTS = pd.merge(TSmd, NTSmd, on=["chrom", "start", "end", "name", "score", "strand"], suffixes=('_TS', '_NTS'))
TSNTS.to_csv(args.output, index=False, header=False, sep='\t')
