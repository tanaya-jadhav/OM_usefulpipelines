from Bio import SeqIO
from Bio.Seq import Seq
import re
import sys


def main(infile, outfile, chrom):
    # infile = 'chr1.fa'
    # outfile = 'chr1_modified.fa'
    fastaheader = '>chr'+str(chrom)+'\n'
    DLE = 'CTTAAG'
    BSPQ = 'GCTCTTC'

    # with open(posfile, 'r') as pf:
    #     lines = pf.readlines()
    # poslist = [int(x.strip('\n')) for x in lines]
    # print(poslist)

    for record in SeqIO.parse(infile, "fasta"):
        ogsequence = record.seq


    # # replace all dle
    ogsequence = str(ogsequence).upper()
    # print(ogsequence)
    nodlesequence = ogsequence.replace('CTTAAG', 'CTTCCG')
    # print(nodlesequence)


    # # replace all bspq
    # nodlebspqsequence = nodlesequence.replace(BSPQ, 'GCCCGGC')
    # nodlebspqsequence = nodlebspqsequence.replace(str(Seq(BSPQ).reverse_complement()), 'GCCCGGC')
    # print(nodlesequence.count(DLE))
    # nodlebspqsequence is a sequence with 0 dle or bspq motifs
    probes = ['TGTAATCCCAGCACTTTGGGAGG', 'TGTAATCCCAGCACTTTGGGTGG', 'TGTAATCCCAGCACTTTGGGCGG', 'TGTAATCCCAGCACTTTGGGGGG']
    probe_replacement = {}
    for probe in probes:
        # print(pos)
        probe_replacement[probe] = 'CTTAAG'
        probe_replacement[str(Seq(probe).reverse_complement())] = 'CTTAAG'

    pattern = re.compile("|".join(probe_replacement.keys()))
    finalsequence = pattern.sub(lambda m: probe_replacement[re.escape(m.group(0))], nodlesequence)

    with open(outfile, 'w') as o:
        o.write(fastaheader)
        o.write(finalsequence + '\n')


if __name__ == '__main__':
    infile = sys.argv[1]
    outfile = sys.argv[2]
    chrom = sys.argv[3]
    main(infile, outfile, chrom)