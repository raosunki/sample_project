from Bio import SeqIO
import random
import numpy as np
import sys

SUCCESS = 0

ECHOICE = {'A': 'CGT',
           'C': 'AGT',
           'G': 'ACT',
           'T': 'ACG'}


def gen_read_quals(read_length):
    """ Returns quality values.
    
    >>> gen_read_quals(read_length=10)
    'IIIIIIIIII'
    """
    
    quals = 'I' * read_length
    return quals


def introduce_error(read_seq, sub_erate):
    """ Replaces read bases with other. 

    Depends on uniform error rate.
    >>> from Bio.Seq import Seq
    >>> read_seq = Seq('ACGTAAGTACGTAAGTACGTAAGTACGTAAGT')
    >>> read_seq_m = introduce_error(read_seq, 0.5)
    >>> read_seq_m != read_seq
    True
    >>> read_seq_m = introduce_error(read_seq, 1e-10)
    >>> read_seq_m != read_seq
    False
    """
    
    pos_to_mutate = []
    for idx, val in enumerate(np.random.uniform(0, 1/sub_erate, len(read_seq))):
        if int(val) == SUCCESS:
            pos_to_mutate.append(idx)
    
    read_seq_m = read_seq.tomutable()
    
    for idx in pos_to_mutate:
        read_seq_m[idx] = random.choice(ECHOICE[read_seq_m[idx]])
    
    return read_seq_m    


def gen_reads(contigs, n, read_length, sub_erate):
    """ Generates reads with a given error rate.
    
    >>> from Bio.Seq import Seq
    >>> contigs = dict()
    >>> contigs['1'] = Seq('ACGATGAGTAGAGACAGAGATAGAGAGATAGAACGATGAGTAGAGAG')
    >>> rec = None
    >>> for i in gen_reads(contigs, 1, 15, 0.1): rec = i
    >>> rec[0] == 'READ1'
    True
    >>> rec[1] == '1'
    True
    >>> 0 <= rec[2] <= len(contigs['1'])-15+1
    True
    >>> len(rec[3]) == 15
    True
    >>> len(rec[4]) == 15
    True
    """
    
    contig_names = list(contigs.keys())
    contig_lens = dict()

    for i in contigs:
        contig_lens[i] = len(contigs[i])
    
    nread = 0
    done = False
    
    while not done:
        contig_name = random.choice(contig_names)
        loc = random.randint(0, contig_lens[contig_name] - read_length + 1)
        read_seq = contigs[contig_name][loc:loc+read_length].upper()

        if 'N' in read_seq:
            continue

        nread += 1
        done = True if nread == n else False

        read_seq_m = introduce_error(read_seq, sub_erate)
        read_qual = gen_read_quals(read_length=read_length)
        
        yield 'READ' + str(nread), contig_name, loc, str(read_seq_m), read_qual


def main(fasta_file):
    contigs = dict()
    with open(fasta_file) as infile:
        for rec in SeqIO.parse(infile, 'fasta'):
            contigs[rec.id] = rec.seq
    
    for rec in gen_reads(contigs, read_length=50, sub_erate=0.01, n=100000):
        header = '%s_%s_%d' % (rec[0], rec[1], rec[2])
        print('@' + header)
        print(rec[3])
        print('+' +  header)
        print(rec[4])

if __name__ == '__main__':
    main(sys.argv[1])
