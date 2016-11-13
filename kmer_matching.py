#!/usr/bin/env python
"""
Using the unique_kmer_library, we process the number of kmers in the query file
that matches into the unique_kmer_library. This information first put into the
dictionary, which is then further processed into two text files for the user.
A match log for all the k-mer matches, and just a match text file with the order
of top matches and a scoring value.

match_output(file_list, unique_kmer_lib, query_name, k, count, mean_threshold)
match_output_scored(file_list, unique_kmer_lib, query_name, k, count, \
                    mean_threshold)
"""
from operator import itemgetter as ig
import math
import os
import pickle
import timeit as ti

counts = {}

#file_list = ['Ebola.fasta', 'Pyrobaculum.fasta', 'Rabies.fasta', \
#                  'Simian.fasta', 'Zika.fasta']
#file_list = ['fasta_short_1.fna', 'fasta_short_2.fna', \
#                  'fasta_short_3.fna', 'fasta_short_4.fna']
    
def Pb(N,p,k):
    '''N = # unique k-mers, p = probability of that k-mer, k = # matched'''
    p = 1.0 / pow(4, p)
    q=1-p
    lt1=[q]*(N-k)
    gt1=list(map(lambda x: p*(N-k+x)/x, range(1,k+1)))
    Pb=1.0
    while (len(lt1) + len(gt1)) > 0:
        if Pb>1:
            if len(lt1)>0:
                Pb*=lt1.pop()
            else:
                if len(gt1)>0:
                    Pb*=gt1.pop()
        else:
            if len(gt1)>0:
                Pb*=gt1.pop()
            else:
                if len(lt1)>0:
                    Pb*=lt1.pop()
    return Pb


def result_file(project_name):
    '''Check and create project files'''
    if not os.path.exists('_RESULTS'):
        os.makedirs('_RESULTS')
    dir = '_RESULTS//' + project_name
    if not os.path.exists(dir):
        os.makedirs(dir)

def score_library_present(unique_kmer_lib):
    """Finds if the UKL scored library exists, return True or False"""
    initial_dir = os.getcwd()
    dir = '_LIB//' + unique_kmer_lib
    try:
        os.chdir(dir)
    except WindowsError:
        os.chdir(initial_dir)
        return False
    ukl = unique_kmer_lib.upper() + '_UKL.DAT'
    if os.path.exists(ukl):
        os.chdir(initial_dir)
        return True
    else:
        os.chdir(initial_dir)
        return False


def ref_kmer_count(file_list, unique_kmer_lib):
    """
    Counts the unique k-mers for each of the reference files, in the unique lib.
    """
    global counts
    for file in file_list:
        counts[file.lower()] = 0
    initial_dir = os.getcwd()
    dir = '_LIB//' + unique_kmer_lib
    os.chdir(dir)
    unique_kmer_lib = unique_kmer_lib + ".DAT"
    library_file = open(unique_kmer_lib, 'rb')
    unique_dict = pickle.load(library_file)
    library_file.close()
    for value in unique_dict.values():
        counts[value.lower()] += 1
    os.chdir(initial_dir)


def matching_noscore(file_list, unique_kmer_lib, query_name, k):
    """
    Using the unique_kmer_library file and the DAT file of query k-mers with the
    min and average phred scores. We look for the perfect matches of k-mers.
    Matches are stored in a dictionary with the source file as the key,
    e.g. Ebola.fastq. The values will be and array of arrays containing the
    min, avg phred scores followed by the kmer itself.
    > {'Ebola.fastq': [[min_phred, avg_phred, kmer] . . . ] . . . }
    """    
    #Initialize dictionary of k-mer matches
    match_dictionary = {}; formatted_dict = {}
    for file in file_list:
        match_dictionary[file.lower()] = {}
        formatted_dict[file.lower()] = []

    initial_dir = os.getcwd()

    dir = '_LIB//' + unique_kmer_lib
    os.chdir(dir)
    #Load Unique K-Mer Library
    unique_kmer_lib = unique_kmer_lib + ".DAT"
    library_file = open(unique_kmer_lib, 'rb')
    unique_dict = pickle.load(library_file)
    #print match_dictionary
    library_file.close()
    os.chdir(initial_dir)

    
    dir = '_QUERY//' + query_name[: query_name.rfind('.')] + \
          '_' + str(k) + '-MER'
    os.chdir(dir)
    all_chunks = os.listdir(os.getcwd())
    len_chunks = len(all_chunks)
    j = 0
    for query_name in all_chunks:
        j += 1
        print('Processing {} of {}, {}-mer chunks . . .'.format(j, len_chunks, k)),
        query_file = open(query_name, 'rb')
        query_dict = pickle.load(query_file) #Query Dictionary with (min, avg)
        query_file.close()
        #os.chdir(initial_dir)

        #Check if query k-mer in unique library
        for query_kmer in query_dict:
            if query_kmer in unique_dict:
                source = unique_dict[query_kmer].lower()
                for i in range(0, len(query_dict[query_kmer])):
                    min_phred = query_dict[query_kmer][i][0]
                    avg_phred = query_dict[query_kmer][i][1]
                    if query_kmer in match_dictionary[source]:
                        cur_min = match_dictionary[source][query_kmer][0]
                        if cur_min < min_phred:
                            match_dictionary[source][query_kmer] = (min_phred, \
                                                                    avg_phred)
                    elif query_kmer not in match_dictionary[source]:
                        match_dictionary[source][query_kmer] = (min_phred, \
                                                                avg_phred)
        print('FINISHED')

    for key in match_dictionary:
        for kmer in match_dictionary[key]:
            format = [match_dictionary[key][kmer][0], \
                      match_dictionary[key][kmer][1], \
                      kmer]
            formatted_dict[key].append(format)
    
    os.chdir(initial_dir)
    #print formatted_dict
    return formatted_dict
    #return match_dictionary


def match_output(file_list, unique_kmer_lib, query_name, \
                 k, count, mean_threshold):
    """
    Consolidates the information of matches, and outputs the information
    in a log file that is viewable in a text reader application.
    NOTE: These matches are ordered by the proportion of matches
    """
    #Dictionary of matches
    match_dictionary = matching_noscore(file_list, unique_kmer_lib, \
                                        query_name, k)
    
    #Initialization Statistics
    global counts
    ref_kmer_count(file_list, unique_kmer_lib)
    C = counts
    L = len(match_dictionary.keys())
    P = math.pow(10, (float(-mean_threshold)/10)) * 100
    M = {}
    for file in file_list:
        M[file] = 0

    initial_dir = os.getcwd()
    dir = '_RESULTS//' + unique_kmer_lib
    os.chdir(dir)
    
    #Match Log
    log_name = ("_MATCH_LOG_" + unique_kmer_lib + ".LOG").upper()
    match_log = open(log_name, 'w')
    #Match Log Writing
    match_log.write("Source: {}\n".format(query_name))
    match_log.write("[K-MER SIZE: {} | ".format(k) + \
                    "SOURCE K-MER COUNT: {} | ".format(count) + \
                    "MINIMUM QUALITY: {} ({}%)]".format(mean_threshold, P))
    match_log.write("\n\n")
    match_log.write("ANALYZED AGAINST {} UNIQUE {}-MER ".format(L, k) + \
                    "REFERENCE SEQUENCES\n")
    #Loop for writing all the matches
    for key in match_dictionary.keys():
        #Sort in decending order
        match_dictionary[key] = sorted(match_dictionary[key], \
                                       key = ig(1), reverse = True)
        match_log.write("\n")
        match_log.write(">{}\n".format(key.upper()))
        #match_log.write("{}\t{}\t{}\n".format("MIN", "AVG", "KMER"))
        for val in match_dictionary[key]:
            M[key] += 1
            match_log.write("{}\t{:.2f}\t{}\n".format(val[0], val[1], val[2]))
        try:
            match_log.write("^{}/{} matches ".format(M[key], C[key]) + \
                            "({:.3E})\n".format(float(M[key]) / C[key]))
        except ZeroDivisionError:
            match_log.write("^{}/{} matches ".format(M[key], C[key]) + \
                            "({:.3E})\n".format(float(0)))
    match_log.close()
    print('{} . . . SAVED IN _RESULTS/{}'.format(log_name, unique_kmer_lib))

    #Match Text
    match = ("_MATCH_" + unique_kmer_lib + ".txt").upper()
    match_txt = open(match, 'w')
    match_txt.write("Source: {}\n".format(query_name))
    match_txt.write("[K-MER SIZE: {} | ".format(k) + \
                    "SOURCE K-MER COUNT: {} | ".format(count) + \
                    "MINIMUM QUALITY: {} ({}%)]".format(mean_threshold, P))
    match_txt.write("\n\n")
    match_txt.write("ANALYZED AGAINST {} UNIQUE {}-MER ".format(L, k) + \
                    "REFERENCE SEQUENCES\n")
    match_txt.write("[ORDER OF TOP MATCHES]\n\n")
    match_txt.write("Order\tReference\tMatches\tMatches".format() + \
                    "/k-mers\tProportion\tP-Value\n".format())
    #Working out the order of scoring
    Score = []
    for file in M:
        try:
            Score.append((file, float(M[file])/C[file]))
        except ZeroDivisionError:
            Score.append((file, float(0)))
    Score = sorted(Score, key = ig(1), reverse = True)
    #print Score
    #Writing order of top matches
    i = 1
    for item in Score:
        source = item[0]
        pbval = Pb(M[source], k, C[source])
        match_txt.write("{:02d}.\t{}\t".format(i, source))
        match_txt.write("{0}\t{0}/{1}\t".format(M[source], C[source]))
        match_txt.write("{:.3E}\t{:.3E}\n".format(item[1], pbval))
        i += 1
    match_txt.close()
    os.chdir(initial_dir)
    print('{} . . . SAVED IN _RESULTS/{}'.format(match, unique_kmer_lib))

'''
file_list = ['ebola.fasta', 'pyrobaculum.fasta', 'rabies.fasta', 'simian.fasta', 'zika.fasta']
unique_kmer_lib = 'EPRSZ_TEST_6'
result_file(unique_kmer_lib)
query_name = 'SRR1660062.fastq'
k = 6
count = 300033
mean_threshold = 30
match_output(file_list, unique_kmer_lib, query_name, k, count, mean_threshold)
'''
    
###SCORING ALGORITHMS BELOW###
def matching_score(file_list, unique_kmer_lib, query_name, k):
    """
    Using the scored unique kmer file and the DAT file of query k-mers with the
    hamming score, min and average phred scores. We look for the perfect matches
    of k-mers. Matches are stored in a dictionary with the source file as the key
    e.g. Ebola.fastq. The values will be and array of arrays containing the
    hamming distance, min, a score, followed by the kmer itself.
    > {'Ebola.fastq': [[ham_dist, min_phred, avg_phred, score, kmer] ... ] ... }
    """    
    #Initialize dictionary of k-mer matches
    match_dictionary = {}; formatted_dict = {}
    for file in file_list:
        match_dictionary[file.lower()] = {}
        formatted_dict[file.lower()] = []

    #Load query kmer .DAT file
    initial_dir = os.getcwd()

    #Load Unique K-Mer Library
    dir = '_LIB//' + unique_kmer_lib
    os.chdir(dir)
    unique_kmer_lib = unique_kmer_lib + "_UKL.DAT"
    library_file = open(unique_kmer_lib, 'rb')
    unique_dict = pickle.load(library_file)
    library_file.close()
    os.chdir(initial_dir)
    
    dir = '_QUERY//' + query_name[: query_name.rfind('.')] + \
          '_' + str(k) + '-MER'
    print dir
    os.chdir(dir)
    all_chunks = os.listdir(os.getcwd())
    len_chunks = len(all_chunks)
    j = 0
    for query_name in all_chunks:
        j += 1
        print('Processing {} of {}, {}-mer chunks . . .'.format(j, len_chunks, k)),
        query_file = open(query_name, 'rb')
        query_dict = pickle.load(query_file) #Query Dictionary with (min, avg) qual.
        query_file.close()   

        #Check if query k-mer in unique library
        for query_kmer in query_dict:
            if query_kmer in unique_dict:
                source = unique_dict[query_kmer][0].lower()
                #print source
                for i in range(0, len(query_dict[query_kmer])):
                    ham_dist = unique_dict[query_kmer][1]
                    min_phred = query_dict[query_kmer][i][0]
                    avg_phred = query_dict[query_kmer][i][1]
                    err_rate = math.pow(10, (float(-min_phred)/10)) * 100
                    score = 0
                    try:
                        score = (float(ham_dist)/k) / err_rate
                    except ZeroDivisionError:
                        score = 0
                    if query_kmer in match_dictionary[source]:
                        cur_score = match_dictionary[source][query_kmer][3]
                        if cur_score < score:
                            match_dictionary[source][query_kmer] = [ham_dist, \
                                                                    min_phred, \
                                                                    avg_phred, \
                                                                    score]
                    elif query_kmer not in match_dictionary[source]:
                        match_dictionary[source][query_kmer] = [ham_dist, \
                                                                min_phred, \
                                                                avg_phred, \
                                                                score]
        print('FINISHED')
    

    for key in match_dictionary:
        for kmer in match_dictionary[key]:
            format = [match_dictionary[key][kmer][0], \
                      match_dictionary[key][kmer][1], \
                      match_dictionary[key][kmer][2], \
                      match_dictionary[key][kmer][3], \
                      kmer]
            formatted_dict[key].append(format)

    os.chdir(initial_dir)
    return formatted_dict

def match_output_scored(file_list, unique_kmer_lib, query_name,
                        k, count, mean_threshold):
    """
    Consolidates the information of matches, and outputs the information
    in a log file that is viewable in a text reader application. With the
    hamming score, and our scoring algorithm.
    NOTE: That the scoring of each k-mer is based on the uniqueness in this case.
    """
    #Dictionary of matches
    match_dictionary = matching_score(file_list, unique_kmer_lib, \
                                        query_name, k)
    
    #Initialization of Statistics
    global counts
    ref_kmer_count(file_list, unique_kmer_lib)
    C = counts
    L = len(match_dictionary.keys())
    P = math.pow(10, (float(-mean_threshold)/10)) * 100
    M = {}; S = {}
    for file in file_list:
        M[file] = 0; S[file] = 0

    initial_dir = os.getcwd()
    dir = '_RESULTS//' + unique_kmer_lib
    os.chdir(dir)
    
    #Match Log
    log_name = ("_MATCH_LOG_" + unique_kmer_lib + ".LOG").upper()
    match_log = open(log_name, 'w')
    #Match Log Writing
    match_log.write("Source: {}\n".format(query_name))
    match_log.write("[K-MER SIZE: {} | ".format(k) + \
                    "SOURCE K-MER COUNT: {} | ".format(count) + \
                    "MINIMUM QUALITY: {} ({}%)]".format(mean_threshold, P))
    match_log.write("\n\n")
    match_log.write("ANALYZED AGAINST {} UNIQUE {}-MER ".format(L, k) + \
                    "REFERENCE SEQUENCES\n")
    #Loop for writing all the matches
    for key in match_dictionary.keys():
        score = 0
        #Sort in decending order
        match_dictionary[key] = sorted(match_dictionary[key], \
                                       key = ig(3), reverse = True)
        match_log.write("\n")
        match_log.write(">{}\n".format(key.upper()))
        #match_log.write("{}\t{}\t{}\n".format("MIN", "AVG", "KMER"))
        for val in match_dictionary[key]:
            M[key] += 1
            S[key] += val[3]
            match_log.write("{:.6f}\t{}\t{:.2f}\t".format(val[3], val[4], val[2]) + \
                            "{}\t{}\n".format(val[1], val[0]))
        match_log.write("[Score: {:.2f} || Matches: ".format(S[key]))
        try:
            match_log.write("{}/{} ".format(M[key], C[key]) + \
                            "({:.3E})]\n".format(float(M[key]) / C[key]))
        except ZeroDivisionError:
            match_log.write("0/0 ({:.3E})]\n".format(float(0)))
    match_log.close()
    print('{} . . . SAVED IN _RESULTS/{}'.format(log_name, unique_kmer_lib))

    #Match Text
    match = ("_MATCH_SCORED_" + unique_kmer_lib + ".txt").upper()
    match_txt = open(match, 'w')
    match_txt.write("Source: {}\n".format(query_name))
    match_txt.write("[K-MER SIZE: {} | ".format(k) + \
                    "SOURCE K-MER COUNT: {} | ".format(count) + \
                    "MINIMUM QUALITY: {} ({}%)]".format(mean_threshold, P))
    match_txt.write("\n\n")
    match_txt.write("ANALYZED AGAINST {} UNIQUE {}-MER ".format(L, k) + \
                    "REFERENCE SEQUENCES\n")
    match_txt.write("[ORDER OF TOP SCORING MATCHES]\n\n")
    match_txt.write("Order\tReference\tScore\tMatches\tMatches".format() + \
                    "/k-mers\tProportions\tP-Value\n".format())
    #Working out the order of scoring
    Score = []
    for file in S:
        Score.append((file, S[file]))
    Score = sorted(Score, key = ig(1), reverse = True)
    #Writing order of top matches
    i = 1
    for item in Score:
        source = item[0]; score_sum = item[1]
        match_txt.write("{:02d}.\t{}\t{:.0f}\t".format(i, source, score_sum))
        match_txt.write("{0}\t{0}/{1}\t".format(M[source], C[source]))
        try:
            match_txt.write("{:.3E}\t{:.3E}\n".format((float(M[source]) / C[source]), Pb(M[source], k, C[source])))
        except ZeroDivisionError:
            match_txt.write("{:.3E}\t{:.3E}\n".format(float(0), float(0)))
        i += 1
    match_txt.close()
    os.chdir(initial_dir)
    print('{} . . . SAVED IN _RESULTS/{}'.format(match, unique_kmer_lib))

def main(file_list_source, unique_kmer_lib, query_file, k, count, q):
    '''main method'''
    result_file(unique_kmer_lib)
    file_list = []
    for file in file_list_source:
        file_list.append(file.lower())
    
    bool = False
    if score_library_present(unique_kmer_lib):
        loop = True
        while loop:
            print('\nUniqueness Scores Present')
            input = raw_input('USE UNIQUENESS SCORE LIBRARY (Y/N): ')
            if input in ['y', 'yes']:
                bool = True
                loop = False
            elif input in ['n', 'no']:
                bool = False
                loop = False
            else:
                continue

    start = ti.default_timer()
    if bool == False:
        match_output(file_list, unique_kmer_lib, query_file, k, count, q)
    else:
        match_output_scored(file_list, unique_kmer_lib, query_file, k, count, q)
    stop = ti.default_timer()
    #print('\t\t\t\t\tRUN TIME: {:.2f} s'.format((stop - start)))

'''   
file_list = ['ebola.fasta', 'pyrobaculum.fasta', 'rabies.fasta', 'simian.fasta', 'zika.fasta']
unique_kmer_lib = 'EPRSZ_TEST_1'
result_file(unique_kmer_lib)
query_name = 'SRR1660062.fastq'
k = 7
count = 300033
mean_threshold = 30
match_output_scored(file_list, unique_kmer_lib, query_name, k, count, mean_threshold)
'''
'''    
#Testing                            
#matching_noscore(file_list, 'TESTING1', 'Query.fastq', 6)
#print matching_score(file_list, 'TESTING1', 'Query.fastq', 6)
match_output(file_list, 'TEST12', 'Query.fastq', 12, 729, 30)
#match_output_scored(file_list, 'TEST12', 'Query.fastq', 12, 729, 30)
'''

#Docstrings
__author__ = "Kevin Cheng"
__credits__ = ["Kevin Cheng", "Aaron Jex"]
__version__ = "python.1.0"
__email__ = "kkcheng.au@gmail.com"
