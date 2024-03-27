import argparse
import sys
import re

def getReverseComplement(seq):
    def getComplementBase(base):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
        if base in complement.keys():
            return complement[base.upper()]
        return base
    reverse_complement = "".join(getComplementBase(base) for base in reversed(seq))
    return reverse_complement

def fa2theoreticalReads(filein, out, kmerSizeList, dinucleotide):
    separator = '|'
    seq = ''
    header = ''
    d = {}
    for line in filein:
        if line.strip().startswith('>'):
            if seq != '':
                d[header[1:]] = seq
                seq = ''
            header = line[1:].strip()
        else:
            seq += line.strip()
    d[header] = seq

    for header in d.keys():
        chromosome = header.split(' ')[0]
        seq = d[header]

        for k in kmerSizeList:
            fiveprime_length = str(k-5)
            motifString = "^.{1," + fiveprime_length + "}" + dinucleotide + ".{3}$"
            motif = re.compile(motifString)
            for i in range(0,len(seq)-k):
                read = seq[i:i+k]
                if motif.match(read):
                    out.write('>' + chromosome + separator + str(i) + separator + str(i+k) + separator + '+' + '\n' + read.upper() + '\n')
                if motif.match(getReverseComplement(read)):
                    out.write('>' + chromosome + separator + str(i) + separator + str(i+k) + separator + '-' + '\n' + getReverseComplement(read).upper() + '\n')
    out.close()


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", default=sys.stdin, type=argparse.FileType('r'), help="Genome file in FASTA format")
parser.add_argument("-o", "--output", default=sys.stdout,type=argparse.FileType('w'), help="Theoretical reads in fasta format")
parser.add_argument('-k','--kmer', nargs='+', help='<Required> kmer length list', type=int, required=True)
parser.add_argument('-d','--dinucleotide', default='TT', help='Dipyrimidine nucleotides', type=str)
args = parser.parse_args()

fa2theoreticalReads(args.input, args.output, args.kmer, args.dinucleotide)