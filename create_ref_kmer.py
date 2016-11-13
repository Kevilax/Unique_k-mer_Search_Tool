#!/usr/bin/env python
"""
This Python script is coded to take a list of reference files. Each reference
file has the sequence extracted using BioPython, then the sequence is divided
into k-mers. K-mers are put into a dictionary with the value equal to the file
name. If a k-mer appears twice, the value is set to 0 to identify a non-unique
k-mer. These k-mer dictionaries are saved to the hard disk using the pickle
package for later use, where we consolidate the k-mers into a single unique
library. These k-mer dictionaries can also be used for later use, in conjunction
with other k-mer dictionaries.

main(file_list, output_name, k)
"""
import pickle
import os
import timeit as ti
from Bio import SeqIO

#Global Variables for Script
global_kmer_file_list = []
global_kmer_dictionary = {}
unique_kmer_library = {}

'''
#Temporary Testers
#file_name_list = ['fasta_short_1.fna', 'fasta_short_2.fna', \
#                  'fasta_short_3.fna', 'fasta_short_4.fna']
file_name_list = ['Ebola.fasta', 'Pyrobaculum.fasta', 'Rabies.fasta', \
                  'Simian.fasta', 'Zika.fasta']
'''

def kmer_dir_check(k):
    """Folder for all exisiting k-mers.DAT files"""
    if not os.path.exists('_KMERS'):
        os.makedirs('_KMERS')
    kmer_dir = '_KMERS//' + str(k) + '-MERS'
    if not os.path.exists(kmer_dir):
        os.makedirs(kmer_dir)


def library_dir_check(dir, k):
    """Unique k-mer library directory"""
    if not os.path.exists('_LIB'):
        os.makedirs('_LIB')
    dir = '_LIB//' + dir
    if not os.path.exists(dir):
        os.makedirs(dir)
    

def set_kmer_dictionary(file_name, k):
    """
    Function creates a hash table of k-mers with their counts. This hash table
    for each reference file is saved to the disk for later and future uses.
    """
    global global_kmer_dictionary
    #Use SeqIO to parse sequence from FASTA, convert to a hash table of k-mers
    #where the hash table contains the file name and reference ID for a unique
    #k-mer within the reference file itself, or a 0 for non unique k-mers
    sequence = SeqIO.parse(open(file_name), 'fasta')
    for line in sequence:
        #ref_id = str(line.id)
        line = str(line.seq)
        n = len(line)
        for i in range(0, n - k + 1):
            kmer = line[i: i + k]
            if len(kmer) == k:
                if kmer in global_kmer_dictionary:
                    global_kmer_dictionary[kmer] = 0
                elif kmer not in global_kmer_dictionary:
                    global_kmer_dictionary[kmer] = file_name

    #Saving kmer hash table object onto the hard disk using pickle package
    initial_dir = os.getcwd()
    dir = '_KMERS//' + str(k) + '-MERS//'
    os.chdir(dir)
    output_name = (file_name[: file_name.rfind('.')] + \
                   '_' + str(k) + '-MER.DAT').upper()
    output_file = open(output_name, 'wb')
    pickle.dump(global_kmer_dictionary, output_file)
    output_file.close()
    global_kmer_dictionary = {}
    os.chdir(initial_dir) #Change back to initial directory

def set_kmer_dictionary_all(file_list, k):
    """
    This function takes each of the reference file and checks if the K-MER.DAT
    file exists. If it exists, we skip the creation of k-mer hash table step.
    Else we will still create the k-mer hash table and store this information.
    """
    initial_dir = os.getcwd()
    dir = '_KMERS//' + str(k) + '-MERS//'
    os.chdir(dir)
    i = 1
    n = len(file_list)
    for file in file_list:
        kmer_dat_file = (file[: file.rfind('.')] + \
                        '_' + str(k) + '-MER.DAT').upper()
        if os.path.isfile(kmer_dat_file) == True:
            print ("Processing {} of {} FASTA to {}-mers . . . . .".\
                   format(str(i), str(n), str(k))),
            print ("FILE EXISTS")
        else:
            print ("Processing {} of {} FASTA to {}-mers . . . . .".\
                   format(str(i), str(n), str(k))),
            os.chdir(initial_dir)
            set_kmer_dictionary(file, k)
            os.chdir(dir)
            print ("FINISHED")
        i += 1
    os.chdir(initial_dir)
            

def set_kmer_file_list(file_list, k):
    """
    Given a list of reference files, the appropriate _K-MER.DAT file is named.
    """
    kmer_file_list = []
    for reference in file_list:
        kmer_file_list.append((reference[: reference.rfind('.')] + \
                              '_' + str(k) + '-MER.DAT').upper())
    global global_kmer_file_list
    global_kmer_file_list = kmer_file_list


def set_unique_library(file_list, output_name, k):
    """
    This function takes our FILE_NAME_K_KMER.DAT files created with the function
    set_kmer_dictionary(file_name, k) for each of our reference files. In this
    function these hash table information are consolidated into a single hash
    table of k-mers. Unique k-mers will contain values of their source file,
    whereas non-unique k-mers are labelled with 0, and these k-mers are deleted
    when we prepare to save the hash tables of unique k-mers to the disk using
    the pickle package.
    """
    #print ("Processing unqiue {}-mers library . . . . .".format(str(k)))
    #Initialize k-mer file names in the kmer_file_list
    global global_kmer_file_list
    kmer_file_list = global_kmer_file_list

    #For each of the k-mer .DAT files, the hash table information is extracted
    #and compared with the unique_kmer_library. If it is unique, information is
    #retained from the .DAT file. Else, the k-mer is marked as non-unique with
    #the value 0.
    global unique_kmer_library
    initial_dir = os.getcwd()
    dir = '_KMERS//' + str(k) + '-MERS//'
    os.chdir(dir)
    i = 1
    n = len(kmer_file_list)
    for input_name in kmer_file_list:
        print ("Consolidating {} of {} k-mer data files".format(str(i), \
                str(n)) + ". . ."),
        input_file = open(input_name, 'rb')
        input_dict = pickle.load(input_file)
        for kmer in input_dict:
            if kmer in unique_kmer_library:
                unique_kmer_library[kmer] = 0
            elif kmer not in unique_kmer_library:
                unique_kmer_library[kmer] = input_dict[kmer]
        print ("FINISHED")
        input_file.close()
        i += 1
    os.chdir(initial_dir)

    #Remove non-unique k-mers
    for kmer in unique_kmer_library.keys():
        if unique_kmer_library[kmer] == 0:
            del unique_kmer_library[kmer]
        #else:
        #    print unique_kmer_library[kmer], kmer

    #Write unique kmer library hash table to the disk with pickle
    dir = '_LIB//' + output_name
    library_dir_check(output_name, k)
    os.chdir(dir)
    output_name = output_name + ".DAT"
    output_file = open(output_name, 'wb')
    pickle.dump(unique_kmer_library, output_file)
    output_file.close()
    output_name = output_name[: output_name.rfind('.')]
    output_list = output_name + ".txt"
    output_list_file = open(output_list, 'w')
    output_list_file.write("{} REFERENCE {}-MER FILES IN {}.DAT\n".\
                               format(len(file_list), str(k), output_name))
    output_list_file.write("@{}\n".format(str(k)))
    file_list = sorted(file_list)
    for file in file_list:
        output_list_file.write('>' + file + '\n')
    output_list_file.close()
    unique_kmer_library = {} #Reset unique_kmer_library
    os.chdir(initial_dir)
    print ("Processing unique {}-mer library . . . . . FINISHED".format(str(k)))


def reset_globals():
    """
    Precautionary function to make sure that the global variables are reset
    for the next time this Python script is called in the same run of the
    program.
    """
    global global_kmer_file_list
    global_kmer_file_list = []
    global global_kmer_dictionary
    global_kmer_dictionary = {}
    global unique_kmer_library
    unique_kmer_library = {}

def main(file_list, output_name, k):
    """
    This would be the public function for the complete execution of creating a
    unique k-mer library.
    """
    #startTime = ti.default_timer()
    #Sets up k-mer directories
    kmer_dir_check(k)
    #Sets up the .DAT file names as a list
    set_kmer_file_list(file_list, k)
    #Creates and saves the .DAT files for k-mers
    set_kmer_dictionary_all(file_list, k)
    #Creates the unique k-mer library
    set_unique_library(file_list, output_name, k)
    reset_globals()
    #stopTime = ti.default_timer()
    #print('\t\t\t\t\tRUN TIME: {:.2f} s'.format((stopTime - \
    #                                             startTime)))
'''
#Testing
file_name_list = ['Ebola.fasta', 'Rabies.fasta']
start = ti.default_timer()    
main(file_name_list, 'JJ', 6)
stop = ti.default_timer()
print stop - start
'''

#Docstrings
__author__ = "Kevin Cheng"
__credits__ = ["Kevin Cheng", "Aaron Jex"]
__version__ = "python.1.0"
__email__ = "kkcheng.au@gmail.com"
