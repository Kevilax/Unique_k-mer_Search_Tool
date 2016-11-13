#!/usr/bin/env python
"""
Using the dictionary .DAT file of unique k-mers, we compare each k-mer with
one another to identify the LOWEST hamming distance. Which means that the
k-mer is closely related to another k-mer. This is put into a new dictionary
and saved as a _UKL.DAT file (UNIQUE KMER LIBRARY) that now contains the
lowest hamming distance of each k-mer.

main(input_name, k)
"""
import pickle
import os
import timeit as ti

def scoring_check(input_name, k):
    """Checkpoint to score uniqueness or not"""
    initial_dir = os.getcwd()
    dir = '_LIB//' + input_name
    os.chdir(dir)
    input_name = input_name + ".DAT"
    input_file = open(input_name, 'rb')
    input_dict = pickle.load(input_file)
    input_file.close()
    continue_bool = False
    n = len(input_dict)
    c = pow(n, 2)
    print('{} Unique {}-kmers'.format(n, k))
    print('{} comparisons required to score each {}-mer\'s uniquenss'.format(c, k))
    loop = True
    while loop:
        continue_query = raw_input('DO YOU WISH TO SCORE UNIQUENESS (Y/N): ')
        if continue_query in ['n', 'no']:
            loop = False
            continue_bool = False
        elif continue_query in ['y', 'yes']:
            continue_bool = True
            loop = False

    os.chdir(initial_dir)
    return continue_bool


def hamming_distance(kmer1, kmer2):
    """
    Given two strings, or k-mers, we calculate the number of differences between
    the two k-mers. This returns the Hamming Distance, which we use as the basis
    for our scoring of the 'uniqueness' of each unique k-mer.
    """
    distance = 0
    for base1, base2 in zip(kmer1, kmer2):
        if base1 != base2:
            distance += 1
    return distance


def unique_kmer_library_scored(input_name, k):
    """
    Using the dictionary .DAT file of unique k-mers, we compare each k-mer with
    one another to identify the LOWEST hamming distance. Which means that the
    k-mer is closely related to another k-mer. This is put into a new dictionary
    and saved as a _UKL.DAT file (UNIQUE KMER LIBRARY) that now contains the
    lowest hamming distance of each k-mer.
    """
    initial_dir = os.getcwd()
    dir = '_LIB//' + input_name
    os.chdir(dir)
    new_unique_kmer_library = {} #New Scored Unique kmer Library
    #Input file of unique kmers using pickle
    input_name = input_name + ".DAT"
    input_file = open(input_name, 'rb')
    input_dict = pickle.load(input_file)
    input_file.close()

    #print len(input_dict)
    #Comparison of a kmer with all other kmers
    for kmer_a in input_dict:
        for kmer_b in input_dict:
            if (kmer_a != kmer_b): #To avoid hamming distance of 0
                ham_dist = k #Best Ham Score for uniqueness
                #If the two kmers are from the same reference don't calc ham_dist
                #Else calculate hamming distance
                if input_dict[kmer_a] != input_dict[kmer_b]:
                    ham_dist = hamming_distance(kmer_a, kmer_b)
                #If the kmer entry already in new dictionary, compare the current
                #hamming distance with this new one. Replace if new distance is
                #smaller than the one in the dictionary.
                if kmer_a in new_unique_kmer_library:
                    cur_ham_dist = new_unique_kmer_library[kmer_a][1]
                    if ham_dist < cur_ham_dist:
                        new_unique_kmer_library[kmer_a] = (input_dict[kmer_a], \
                                                           ham_dist)
                #If kmer not in new dictionary, add it.
                elif kmer_a not in new_unique_kmer_library:
                    new_unique_kmer_library[kmer_a] = (input_dict[kmer_a], \
                                                       ham_dist)

    #Output into a new _UKL.DAT file.close via pickle dictionary dump
    output_name = input_name[: input_name.rfind('.')] + "_UKL.DAT"
    output_file = open(output_name, 'wb')
    pickle.dump(new_unique_kmer_library, output_file)
    output_file.close()
    os.chdir(initial_dir)


def main(input_name, k):
    """
    Simple way of calling the function needed in the inputs script.
    """
    bool = scoring_check(input_name, k)
    if bool:
        startTime = ti.default_timer()
        print ("Processing uniqueness for unique {}-mer library . . .".format(k)),
        unique_kmer_library_scored(input_name, k)
        print ("FINISHED")
        stopTime = ti.default_timer()
        print('\t\t\t\t\tRUN TIME: {:.2f} s'.format((stopTime - \
                                                     startTime)))
'''
start = ti.default_timer()
main('jj', 12)
stop = ti.default_timer()
print stop - start
'''

#Docstrings
__author__ = "Kevin Cheng"
__credits__ = ["Kevin Cheng", "Aaron Jex"]
__version__ = "python.1.0"
__email__ = "kkcheng.au@gmail.com"
