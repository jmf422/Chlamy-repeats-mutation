#!/usr/bin/env python

# check kmers found by tandem repeats finder in the chloroplast assembly to kmers found in the Illumina reads by k-Seek.
# check_cp_kmers.py

def rev_comp(kmer):
    '''Returns reverse complement of k-mer.'''
    result = ''
    nucleotides = ['A', 'C', 'T', 'G']
    for char in kmer[::-1]:
        result += nucleotides[(nucleotides.index(char)+2)%len(nucleotides)]
    return result

def get_rotations(kmer):
    '''Returns all rotated versions of k-mer.'''
    return [(kmer*2)[i:len(kmer)+i] for i in range(0, len(kmer))]

def get_reverse_complements(rotations):
    '''Returns reverse complements of a k-mer's rotations.'''
    return [rev_comp(kmer) for kmer in rotations]

def is_same_kmer(rotations, reverse_complements, query):
    '''Returns True if query is either a rotated version of the target
    k-mer or the reverse complement of one of these rotations.'''
    if (query in rotations) or (query in reverse_complements):
        return True
    else:
        return False


#check all rotations of the cp kmers found with tandem repeat finder against the established chlamy kmers found with kseek

# read in all the chlamy kmers

allkmers=list()
with open('chlamy_kmers.txt', 'r') as f:
	for line in f:
		allkmers.append(line.strip())

    
# read in cp kmers
cpkmers=list()

with open('cp.kmers.txt', 'r') as f:
	for line in f:
		cpkmers.append(line.strip())

# print the kmers that are in both datasets (considering rotations and reverse complements)
for i in cpkmers:
	rotations = get_rotations(i)
	rev_comps = get_reverse_complements(rotations)
	for j in allkmers:
		bool = is_same_kmer(rotations, rev_comps, j)
		if bool == True:
			print (i, j)

			


