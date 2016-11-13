#!/usr/bin/env python
"""
This Python script creates a dictionary of k-mers and saves it to the hard drive.
Only the k-mers in the FASTQ file from all the reads that pass the quality checks
such as the mean and the quality window, will be added to the dictionary of kmers
The dictionary format will be k-mer as the keys followed by a tuple with the
lowest phred quality score and the average score. e.g. {kmer: (min, avg), ... }

main(file_name, k, mean_threshold, quality_threshold, quality_window)
"""
import os
import pickle
import timeit as ti
from Bio import SeqIO

count = 0 #Number of k-mers in the FASTQ file from all the reads
count_passed = 0 #Number of k-mers that pass the quality checks
total_chunks = 0

def folder_check(file_name, k):
    """Checks and/or creates the _QUERY folder to store query k-mer data"""
    if not os.path.exists('_QUERY'):
        os.makedirs('_QUERY')
    query_dir = '_QUERY//' + file_name[: file_name.rfind('.')] + '_' + \
                str(k) + '-MER'
    if os.path.exists(query_dir):
        return True
    elif not os.path.exists(query_dir):
        os.makedirs(query_dir)
        return False


def kmer_dictionary(file_name, k, mean_threshold, \
                    quality_threshold, quality_window):
    """
    Produces a k-mer dictionary and saves as a pickle. As {k-mer: [values]}.
    Only k-mers that pass the kmer_quality_check() will be stored.
    """
    
    #Initialize the k-mer dictionary with the k-mer as keys, quality as [value]
    kmer_dict = {}
    sequence = SeqIO.parse(open(file_name), 'fastq') #Read FASTQ
    global count
    global count_passed
    count = 0
    count_passed = 0
    #Iterate every read in the FASTQ
    chunk_i = 0; chunk = 0
    initial_dir = os.getcwd()
    dir = '_QUERY//' + file_name[: file_name.rfind('.')] + '_' + \
          str(k) + '-MER'
    os.chdir(dir)
    
    for read in sequence:
        n = len(str(read.seq))
        #Find every k-mer and its quality scores
        i = 0
        while i < (n - k + 1):
        #for i in range(0, n - k + 1):
            kmer = str(read.seq)[i: i + k]; count += 1
            quality = read.letter_annotations["phred_quality"][i: i + k]
            #If kmer quality does not meet thresholds, move on to next k-mer
            #at the position where it fails the quality check.
            check = kmer_quality_check(quality, mean_threshold, \
                                       quality_threshold, quality_window)
            #print n, quality, (check),
            if check[0] == False:
                i += check[1]
                # print i
            #If passes the checks then append to dictionary
            else:
                if kmer in kmer_dict:
                    kmer_dict[kmer].append((min(quality), check[2]))
                elif kmer not in kmer_dict:
                    kmer_dict[kmer] = [(min(quality), check[2])]
                count_passed += 1                
                # print i
            i += 1
        chunk_i += 1
        if chunk_i == 10000:
            chunk_i = 0
            chunk += 1
            output_name = (file_name[: file_name.rfind('.')] + \
                           '_' + str(k) + '-MER_' + str(chunk) + '.DAT').upper()
            print output_name
            output_file = open(output_name, 'wb')
            pickle.dump(kmer_dict, output_file)
            kmer_dict = {}
            output_file.close()
    #Last chunk
    chunk += 1
    output_name = (file_name[: file_name.rfind('.')] + \
                   '_' + str(k) + '-MER_'+ str(chunk) + '.DAT').upper()
    print output_name
    output_file = open(output_name, 'wb')
    pickle.dump(kmer_dict, output_file)
    kmer_dict = {}
    output_file.close()
    global total_chunks
    total_chunks = chunk
    sequence.close()
    os.chdir(initial_dir)
   
    

def kmer_quality_check(quality, mean_threshold, \
                       quality_threshold, quality_window):
    """
    This checks the quality of the k-mer, to ensure that the quality passes
    the thresholds. The mean is very straight forward formula to calculate
    the average of the quality scores. Whereas the quality window check utilizes
    an algorithm, where we only check each base once to find a consequtive seq
    of low quality k-mers, rather than creating a k-mer of a k-mer and checking
    every possible window again and again. This way we avoid a polynomial Big-Oh,
    but instead, in theory, obtain a linear computation time. 
    """
    
    #Mean quality must pass threshold
    mean = float(sum(quality))/len(quality) #Mean of quality
    bool = (True, 0, mean) #Initialization will pass quality check == True
    if mean < quality_threshold:
        bool = (False, 0, mean)
        return bool

    #Quality Window must pass threshold in the given quality k-mer
    i = 0
    while i < (len(quality) - quality_window):
    #for i in range(0, len(quality) - quality_window + 1):
        #Condition at first quality doesn't meet threshold
        if quality[i] < quality_threshold:
            #Create a kmer_window to iterate, with legnth m
            kmer_window = quality[i: i + quality_window]; m = len(kmer_window)
            #Iterate from 1, as we have already determined 0-index has failed
            for j in range(1, m):
                #At each index, if quality fails
                if kmer_window[j] < quality_threshold:
                    #And we are at the end of the window, FAIL the quality check
                    if j == (m - 1):
                        bool = (False, i + j, mean)
                        return bool
                    #If not at the end of the window, keep iterating as the next
                    #quality may or may not pass threshold. Cannot fail yet.
                    else:
                        continue
                #Else where quality passes in this window, update the index to
                #where we ended the place that passes threshold.
                else:
                    #Update i to i + the j-th position we last observed a passing
                    #quality. This leads to moving to the next i + j + 1 position
                    #at the end of the loop so we do not check the same index 2x.
                    i += j
                    break
        i += 1

    #After the two checks, if passes then should return True. Just in case bool
    #was set to False at any point but wasn't returned, then returns False.
    return bool
'''
def main(file_name, k, mean_threshold, quality_threshold, quality_window):
    """
    This would be the public function for the complete execution of creating a
    dictionary of k-mers and their quality scores in the FASTQ file.
    """
    global total_chunks
    global count; global count_passed
    
    bool = folder_check(file_name, k)
    if bool == False:
        print ("Processing {} into {}-mers chunks. . . . .\n".format(file_name, k)),
        kmer_dictionary(file_name, k, mean_threshold, \
                        quality_threshold, quality_window)
        print ("Processed {}-mers in {} chunks".format(k, total_chunks)),
        print (" . . . . . FINISHED")

        initial_dir = os.getcwd()
        dir = '_QUERY//' + file_name[: file_name.rfind('.')] + '_' + \
              str(k) + '-MER'
        os.chdir(dir)

        prop_file = open("Properties.txt", "w")
        prop_file.write(str(k) + '\n')
        prop_file.write(str(count) + '\n')
        prop_file.write(str(count_passed) + '\n')
        prop_file.write(str(mean_threshold) + '\n')
        prop_file.write(str(quality_threshold) + '\n')
        prop_file.write(str(quality_window))
        prop_file.close()

        os.chdir(initial_dir)

        return(count, count_passed, mean_threshold, quality_threshold, quality_window)
    elif bool == True:
        initial_dir = os.getcwd()
        dir = '_QUERY//' + file_name[: file_name.rfind('.')] + '_' + \
              str(k) + '-MER'
        os.chdir(dir)
        
        prop_file = open("Properties.txt", "r")
        lines = prop_file.readlines()
        print lines
        count = lines[1]
        count_passed = lines[2]
        mean_threshold_2 = lines[3]
        quality_threshold_2 = lines[4]
        quality_window_2 = lines[5]
        prop_file.close()

        os.chdir(initial_dir)

        loop = True
        while loop:
            print("Existing query {}-mer file exists".format(k))
            print("Quality Threshold: {}".format(mean_threshold_2))
            input = raw_input("USE EXISTING QUERY DATA (Y/N): ")

            if input.lower() in ['n', 'no']:
                print ("Processing {} into {}-mers chunks. . . . .\n".format(file_name, k)),
                kmer_dictionary(file_name, int(k), int(mean_threshold), \
                                int(quality_threshold), int(quality_window))
                print ("Processed {}-mers in {} chunks".format(k, total_chunks)),
                print (" . . . . . FINISHED")

                initial_dir = os.getcwd()
                dir = '_QUERY//' + file_name[: file_name.rfind('.')] + '_' + \
                      str(k) + '-MER'
                os.chdir(dir)

                prop_file = open("Properties.txt", "w")
                prop_file.write(str(k) + '\n')
                prop_file.write(str(count) + '\n')
                prop_file.write(str(count_passed) + '\n')
                prop_file.write(str(mean_threshold) + '\n')
                prop_file.write(str(quality_threshold) + '\n')
                prop_file.write(str(quality_window))
                prop_file.close()

                os.chdir(initial_dir)

                return(count, count_passed, mean_threshold, quality_threshold, \
                       quality_window)
            
            elif input.lower() in ['y', 'yes']:
                return(int(count[:-1]), int(count_passed[:-1]), int(mean_threshold_2[:-1]),\
                       int(quality_threshold_2[:-1]), int(quality_window_2))
'''            
'''
print main("Query.fastq", 13, 20, 20, 2)
'''        

def main(file_name, k, mean_threshold, quality_threshold, quality_window):
    """
    This would be the public function for the complete execution of creating a
    dictionary of k-mers and their quality scores in the FASTQ file.
    """
    print ("Processing {} into {}-mers chunks. . . . .\n".format(file_name, k)),
    folder_check(file_name, k)
    kmer_dictionary(file_name, k, mean_threshold, \
                    quality_threshold, quality_window)
    global total_chunks
    print ("Processed {}-mers in {} chunks".format(k, total_chunks)),
    print (" . . . . . FINISHED")
    global count; global count_passed
    return(count, count_passed)


'''
start = ti.default_timer()
print main("SRR1660062.fastq", 6, 30, 30, 2)
stop = ti.default_timer()
print stop - start
'''
'''
#Testing
start = ti.default_timer()
print main("Query.fastq", 12, 20, 20, 2)
stop = ti.default_timer()
print stop - start
'''

#Docstrings
__author__ = "Kevin Cheng"
__credits__ = ["Kevin Cheng", "Aaron Jex"]
__version__ = "python.1.0"
__email__ = "kkcheng.au@gmail.com"
