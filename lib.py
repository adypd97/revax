#!/usr/bin/env python3

class NotFASTAFormatError(Exception):
	pass 
class NotAValidKMerLengthError(Exception):
	pass

def is_FASTA(f):
    if f.readline()[0] == ">":
        for n in ''.join(f.read().split("\n")):
            if n not in ['A', 'T', 'G', 'C']:
                return False 
    else:
        return False
    return True

def hamming(seq1, seq2):
    err = sum([0 if s1 == s2 else 1 for s1, s2 in zip(seq1, seq2)])
    return err

def gc_content(seq):
    return (sum([1 for n in seq if n == 'G' or n == 'C'])) / len(seq) * 100

def naive_pattern_match(dna, seq, error=2):
    '''NOTE: SLOW ALGORITHM'''
    import sys
    len_dna = len(dna)
    len_seq = len(seq)
    if len_dna < len_seq:
        print("Genome is smaller than hexamer length")
        sys.exit(1)
    matches = []
    for i in range(0, len_dna - len_seq):
        if hamming(dna[i:i+len_seq], seq) <= error:
            matches.append(i)
    return matches

def Not(n):
    if n == 'A':
        return 'T'
    elif n == 'T':
        return 'A'
    elif n == 'G':
        return 'C'
    else:
        return 'G'

def rev_comp(seq):
    comp = [Not(n) for n in seq]
    return ''.join(comp[::-1])

def _get_seq(filename):
    import sys
    '''
    import argparse
    parser = argparse.ArgumentParser(description="Enter length k for a k-mer and n number of kmers")
    parser.add_argument("-f", "--file", help="Enter a fasta file", required=True)
    args = vars(parser.parse_args())
    '''
    try:
        try:
            with open(filename) as f:
                # get the sequence into 'sq'
                if not is_FASTA(f):
                    raise NotFASTAFormatError
                f.seek(0)
                sq = ''.join(f.read().strip().split("\n")[1:])
        except Exception as e:
            print('Following error occured : ' + str(e))
            sys.exit(1)
    except FileNotFoundError:
        print('No file named "genome.fna" found')
        sys.exit(1)

    return sq


def translate(seq):
    seq = seq.replace("Ψ", "U")
    # Protein translation table
    # AUG is the start codon
    # UAA, UGA and UAG represent stop codons
    tt = { 'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
	'ACA':'U', 'ACC':'U', 'ACG':'U', 'ACU':'U',
	'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',                 
	'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
	'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
	'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
	'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
	'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
	'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
	'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
	'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W',
	}

    assert len(seq) % 3 == 0
    pep_seq = [tt[seq[i:i+3]] for i in range(0, len(seq)-2, 3)]
    return ''.join(pep_seq)


virus_seq = _get_seq('wuhan-sars-cov2-seq.fasta') 

if __name__ == "__main__":

    virus_sig = "AUG UUU GUU UUU CUU GUU UUA UUG CCA CUA GUC UCU AGU CAG UGU GUU".replace(" ", "").replace("U", "T")
    vaccine_sig = "AUG UUC GUG UUC CUG GUG CUG CUG CCU CUG GUG UCC AGC CAG UGU GUG".replace(" ", ""). replace("U", "T")
    print(f'Start of Signal Protein mRNA sequence {naive_pattern_match(virus_seq, virus_sig)}')

