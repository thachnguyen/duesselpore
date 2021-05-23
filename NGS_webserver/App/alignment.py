'''
This is the draft of barcodes detection.
- Group the sequences contain a specific barcode (8 bps)
- Aligns them use multiple alignments (1st version use Clustal Omega, the best method I used for multiple sequence alignment) I think it work well for DNA sequence also)
I also wrote Python ClustalOmega interface during my PhD
Reference:
Sievers, Fabian, et al. "Fast, scalable generation of high‐quality protein multiple sequence alignments using Clustal Omega." Molecular systems biology 7.1 (2011): 539.

Requirement packages:
- python3 (not recommend python2)
- biopython (1.78)
- clustalo inside csbplus package, the computational biology toolbox developed in my research group.

''' 

from Bio import SeqIO
from clustalo import ClustalO, align

seq = SeqIO.parse('test.fastq', format='fastq')
#!TODO: parser works for Illumina, check for nanopore seq 

sequence_list = list(seq)
# May has problem when there is not enough RAM

def read_barcodes_list(barcode_file = 'barcodes.csv'):
    
    return

def group_by_barcode(barcodes='CTCAGGTA', seq_list = sequence_list):
    """ Select sequence list contains customize barcode

    Args:
        barcodes (str, optional): [description]. Defaults to 'CTCAGGTA'.
        seq_list ([type], optional): [description]. Defaults to sequence_list.

    Returns:
        store_list: contain barcodes sequence list 
    """    

    store_list = []
    for i, seq1 in enumerate(seq_list):
        if barcodes in seq1.seq:
            store_list.append(seq1)

    return store_list

def multiple_alignment(seq_list = store_list):
    '''Use Clustal omega for multiple sequences aligment. 
    '''
    mult_seq = [str(seq.seq) for seq in sequence_list]
    clustalo = ClustalO()
    result = clustalo.align(mult_seq)
    result.to_fasta('sample.fasta')
    #Need to improve, it take ~ 5 minutes 3630 entries (~325 bp per entry) on my PC
    return

def align_by_barcode():
    return


test_list = group_by_barcode()
# multiple_alignment()
