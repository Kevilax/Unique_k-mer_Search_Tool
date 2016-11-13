#!/usr/bin/env python
"""
This Python script ties every component used in the creation of unique k-mer
libraries, creation of k-mers, then the matching of unique k-mers to the query
k-mers.
"""
import os
import sys
import timeit as ti
import create_ref_kmer as crk
import unique_kmer_lib as crk_score
import create_query_kmer as cqk
import kmer_matching as km

#Global Instance Variables
help = "'addRef' or '+'     - add multiple files to reference library.\n" + \
       "'match' or '='      - query k-mer matching to unique k-mers.\n" + \
       "'refLib' or 'L'     - displays current files in reference library.\n" + \
       "'createLib' or 'C'  - creates and save reference k-mer library.\n" + \
       "'openLib' or 'O'    - open reference library list.\n" + \
       "\nRefer to README.txt for further help and example."
invalid = ['\'', '/', ':', '*', '?', '"', '<', '>', '|', '.']

reference_file_set = set()

current_lib = ''
count_passed = 0
global_k = 0
query_file = ''

#Creation of the library
def create_method():
    """Create the unique kmer library"""
    global global_k
    global count_passed
    global query_file
    global current_lib
    global reference_file_set
    loop = True
    #Reprint current reference library
    n = len(reference_file_set)
    if n == 1 or n == 0:
        print(str(n) + " REFERENCE LIBRARY")
        if n == 0:
            print("REQUIRES AT LEAST 1 REFERENCE SEQUENCE TO CREATE"),
            print("A UNIQUE K-MER LIBRARY.")
            loop = False
    else:
        print(str(n) + " REFERENCE LIBRARIES")
    for reference in reference_file_set:
        print reference
    #Query to create k-mer library
    while loop:
        query = raw_input('Create unique k-mer library (Y/N): ')                
        if query.lower() in ['n', 'no']:
            loop = False
            return False
        elif query.lower() in ['y', 'yes']:
            #Unique K-mer Library Name
            project_name = raw_input('\nUNIQUE K-MER LIBRARY NAME: ')
            try:                    
                for char in invalid:
                    if char in project_name:
                        raise ValueError
            except ValueError:
                print('LIBRARY NAME MUST NOT CONTAIN ' + \
                      '\\, /, ., :, *, ?, ", <, >, |')
                loop = False
                break
            #K-MER SIZE
            kmer_size = raw_input('K-MER SIZE: ')
            try:
                kmer_size = int(kmer_size)
                if kmer_size < 1:
                    raise ValueError

                dir = '_LIB//' + project_name
                if os.path.exists(dir):
                    loop2 = True
                    while loop2:
                        ow = \
                           raw_input('OVERWRITE EXISTING LIBRARY (Y/N): ')
                        if ow.lower() in ['n', 'no']:
                            loop = False
                            loop2 = False
                        elif ow.lower() in ['y', 'yes']:
                            loop2 = False
                    if loop == False:
                        break
                    global_k = kmer_size
                    current_lib = project_name
            except ValueError:
                print('K-MER SIZE MUST BE POSITIVE INTEGER')
                loop = False
                break
            #Creatation of k-mer library
            print
            startTime = ti.default_timer()
            crk.main(reference_file_set, project_name, kmer_size)
            stopTime = ti.default_timer()
            print('\t\t\t\t\tRUN TIME: {:.2f} s'.format((stopTime - \
                                                         startTime)))
            #Scored Unique k-mer library
            crk_score.main(project_name, kmer_size)
            
            
            current_lib = project_name
            global_k = kmer_size
            
            loop = False
    return True

###MAIN###
#START
print "UNIQUE K-MER MATCHING TOOL v.Python.1.1b\n"
print "'help' or '?'       - displays available commands."
print help
running = True
#Loop to keep program running
while running:
    #User raw inputs as a string
    input = raw_input('\n> ')
    try:
        input_array = input.split(); command = input_array[0]
    except IndexError:
        continue

    ###INPUT COMMANDS###
    #Help
    if command.lower() in ['help', '?']:
        print help
        
    #Exit
    elif command.lower() in ['exit', 'close', 'q']:
        loop = True
        while loop:
            exit_input = raw_input('\tEnd Program (Y/N): ')
            if exit_input.lower() in ['n', 'no']:
                loop = False
            elif exit_input.lower() in ['y', 'yes', '']:
                print '\tPROGRAM TERMINATED'
                sys.exit()
                
    #Add to reference_file_set
    elif command.lower() in ['addref', '+', 'add']:
        #If in the argument format
        if len(input_array) > 1:
            for arg in input_array[1:]:
                arg = arg.lower()
                if os.path.exists(arg):
                    if (('.fasta' in arg.lower()) or ('.fna' in arg.lower())):
                        reference_file_set.add(arg)
                        print("{} . . . FILE ADDED".format(arg))
                    else:
                        print("{} . . . FILE TYPE ERROR".format(arg))
                else:
                    print("{} . . . FILE NOT FOUND".format(arg))
        #If no arguments
        else:
            loop = True
            while loop:
                arg = raw_input('Reference File Name: ')
                arg = arg.lower()
                if arg == '' or arg.lower() == 'exit':
                    loop = False
                elif os.path.exists(arg):
                    if (('.fasta' in arg.lower()) or ('.fna' in arg.lower())):
                        reference_file_set.add(arg)
                        print("{} . . . FILE ADDED\n".format(arg))
                    else:
                        print("{} . . . FILE TYPE ERROR\n".format(arg))
                else:
                    print("{} . . . FILE NOT FOUND\n".format(arg))

    #Reprint current reference library
    elif command.lower() in ['reflib', 'l', 'ref', 'lib', 'r']:
        n = len(reference_file_set)
        if n == 1 or n == 0:
            print(str(n) + " REFERENCE LIBRARY")
        else:
            print(str(n) + " REFERENCE LIBRARIES")
        for reference in reference_file_set:
            print reference

    #Create Unique K-Mer Library with current Reference Library List
    elif command.lower() in ['createlib', 'c', 'savelib', 's', 'v']:
        create_method()

    #Open an existing library
    elif command.lower() in ['openlib', 'o', 'loadlib', 'load', '^']:
        #If in the argument format
        if len(input_array) > 1:
            arg = "_LIB//" + input_array[1] + "//" + input_array[1] + '.txt'
            if os.path.exists(arg):
                reference_file_set = set()
                current_lib = input_array[1]
                input_file = open(arg, 'r')
                for line in input_file:
                    if line[0] == '@':
                        global_k = int(line[1:])
                    elif line[0] == '>':
                        reference_file_set.add(line[1: -1])
                input_file.close()
                print
                n = len(reference_file_set)
                if n == 1 or n == 0:
                    print(str(n) + " REFERENCE LIBRARY")
                else:
                    print(str(n) + " REFERENCE LIBRARIES")
                for reference in reference_file_set:
                    print reference
        #If no arguments
        else:
            loop = True
            while loop:
                arg = raw_input('Unique Library Name: ')
                path = "_LIB//" + arg + "//" + arg + ".txt"
                if arg == '' or arg.lower() == 'exit':
                    loop = False
                elif os.path.exists(path):
                    reference_file_set = set()
                    current_lib = arg
                    input_file = open(path, 'r')
                    for line in input_file:
                        if line[0] == '@':
                            global_k = int(line[1:])
                        elif line[0] == '>':
                            reference_file_set.add(line[1: -1])
                    input_file.close()
                    loop = False
                    print
                    n = len(reference_file_set)
                    if n == 1 or n == 0:
                        print(str(n) + " REFERENCE LIBRARY")
                    else:
                        print(str(n) + " REFERENCE LIBRARIES")
                    for reference in reference_file_set:
                        print reference
                else:
                    print("{} . . . FILE NOT FOUND\n".format(arg))

    #Matching Algorithm
    elif command.lower() in ['match', '=']:
        q_threshold = 0
        q_window = 0
        q_window_t = 0
        #Use current library?
        if global_k != 0 or len(current_lib) > 0:
            loop_1 = True
            while loop_1:
                step_2 = False
                #Reprint Reference Library
                print('CURRENT UNIQUE {}-MER LIBRARY: {}'.format(global_k, \
                                                                 current_lib))
                n = len(reference_file_set)
                if n == 1 or n == 0:
                    print(str(n) + " REFERENCE LIBRARY")
                else:
                    print(str(n) + " REFERENCE LIBRARIES")
                for reference in reference_file_set:
                    print reference

                print('\n[1] MATCH WITH CURRENT K-MER LIBRARY')
                print('[2] CHANGE K-MER SIZE AND MATCH')
                print('[3] EDIT K-MER LIBRARY')
                use_input = raw_input('INPUT: ')
                try:
                    use_input = int(use_input)
                except ValueError:
                    print
                
                if use_input == 2:
                    use_input = int(use_input)
                    loop_1 = False
                    print
                    create_method()
                    print
                    step_2 = True
                elif use_input == 1:
                    use_input = int(use_input)
                    loop_1 = False
                    step_2 = True
                elif use_input == '':
                    loop_1 = False

                step_3 = False
                if step_2 == True:
                    while step_2:
                        query_file = raw_input('QUERY FILE: ')
                        if query_file == '' or query_file.lower() == 'exit':
                            step_2 = False
                        elif os.path.exists(query_file):
                            step_2 = False
                            step_3 = True
                        else:
                             print("{} . . . FILE NOT FOUND\n".format(query_file))

                step_4 = False
                if step_3 == True:
                    while step_3:
                        if ('.fastq' in query_file.lower()):
                            step_3 = False
                            step_4 = True
                        else:
                            step_3 = False
                            print("{} . . . FILE TYPE ERROR".format(query_file))

                step_5 = False
                if step_4 == True:
                    while step_4:
                        q_threshold = raw_input('QUALITY THRESHOLD [0-40]: ')
                        try:
                            q_threshold = int(q_threshold)
                            if q_threshold < 0 or q_threshold > 40:
                                print('QUALITY MUST BE BETWEEN 0 AND 40\n')
                            elif q_threshold == '':
                                step_4 = False
                            elif q_threshold >= 0 or q_threshold <= 40:
                                q_threshold = int(q_threshold)
                                step_4 = False
                                step_5 = True
                        except ValueError:
                            print('QUALITY MUST BE AN INTEGER [0-40]\n')

                step_6 = False
                if step_5 == True:
                    while step_5:
                        q_window = raw_input('QUALITY WINDOW [0-' + \
                                             '{}]: '.format(global_k))
                        try:
                            q_window = int(q_window)
                            if q_window < 0 or q_window > global_k:
                                print('QUALITY WINDOW MUST BE BETWEEN 0 AND k\n')
                            elif q_window == '':
                                step_5 = False
                            elif q_window >= 0 or q_window <= global_k:
                                q_window = int(q_window)
                                step_5 = False
                                step_6 = True
                        except ValueError:
                            print('QUALITY WINDOW MUST BE AN INTEGER [0-' + \
                                  '{}]\n'.format(global_k))

                step_7 = False
                if step_6 == True:
                    while step_6:
                        q_window_t = raw_input('QUALITY WINDOW THRESHOLD' + \
                                               '[0-40]: ')
                        try:
                            q_window_t = int(q_window_t)
                            if q_window_t < 0 or q_window_t > 40:
                                print('QUALITY MUST BE BETWEEN 0 AND 40\n')
                            elif q_window_t == '':
                                step_6 = False
                            elif q_window_t >= 0 or q_window_t <= 40:
                                q_window_t = int(q_window_t)
                                step_6 = False
                                step_7 = True
                        except ValueError:
                            print('QUALITY MUST BE AN INTEGER [0-40]\n')

                if step_7 == True:
                    start = ti.default_timer()
                    count_passed = cqk.main(query_file, global_k, q_threshold, \
                                   q_window_t, q_window)
                    '''
                    q_t = count_passed[2]
                    q_w_t = count_passed[3]
                    q_w = count_passed[4]
                    c_p = count_passed[1]
                    '''
                    #print count_passed
                    km.main(reference_file_set, current_lib, query_file, \
                            global_k, count_passed, q_threshold)
                    stop = ti.default_timer()
                    print('\t\t\t\t\tRUN TIME: {:.2f} s'.format((stop - start)))

                
        #If empty library
        elif len(reference_file_set) == 0:
            loop_1 = False
            print('NO REFERENCE LIBRARIES')
            print('USE \'+\' or \'addRef\' to add reference sequences')
        #If no preset unique k-mer library
        else:
            bool = create_method()
            step_2 = True
            if bool == False:
                step_2 = False
            else:
                print
            step_3 = False
            if step_2 == True:
                while step_2:
                    query_file = raw_input('QUERY FILE: ')
                    if query_file == '' or query_file.lower() == 'exit':
                        step_2 = False
                    elif os.path.exists(query_file):
                        step_2 = False
                        step_3 = True
                    else:
                         print("{} . . . FILE NOT FOUND\n".format(query_file))

            step_4 = False
            if step_3 == True:
                while step_3:
                    if ('.fastq' in query_file.lower()):
                        step_3 = False
                        step_4 = True
                    else:
                        step_3 = False
                        print("{} . . . FILE TYPE ERROR".format(query_file))

            step_5 = False
            if step_4 == True:
                while step_4:
                    q_threshold = raw_input('QUALITY THRESHOLD [0-40]: ')
                    try:
                        q_threshold = int(q_threshold)
                        if q_threshold < 0 or q_threshold > 40:
                            print('QUALITY MUST BE BETWEEN 0 AND 40\n')
                        elif q_threshold == '':
                            step_4 = False
                        elif q_threshold >= 0 or q_threshold <= 40:
                            q_threshold = int(q_threshold)
                            step_4 = False
                            step_5 = True
                    except ValueError:
                        print('QUALITY MUST BE AN INTEGER [0-40]\n')

            step_6 = False
            if step_5 == True:
                while step_5:
                    q_window = raw_input('QUALITY WINDOW [0-' + \
                                         '{}]: '.format(global_k))
                    try:
                        q_window = int(q_window)
                        if q_window < 0 or q_window > global_k:
                            print('QUALITY WINDOW MUST BE BETWEEN 0 AND k\n')
                        elif q_window == '':
                            step_5 = False
                        elif q_window >= 0 or q_window <= global_k:
                            q_window = int(q_window)
                            step_5 = False
                            step_6 = True
                    except ValueError:
                        print('QUALITY WINDOW MUST BE AN INTEGER [0-' + \
                              '{}]\n'.format(global_k))

            step_7 = False
            if step_6 == True:
                while step_6:
                    q_window_t = raw_input('QUALITY WINDOW THRESHOLD' + \
                                           '[0-40]: ')
                    try:
                        q_window_t = int(q_window_t)
                        if q_window_t < 0 or q_window_t > 40:
                            print('QUALITY MUST BE BETWEEN 0 AND 40\n')
                        elif q_window_t == '':
                            step_6 = False
                        elif q_window_t >= 0 or q_window_t <= 40:
                            q_window_t = int(q_window_t)
                            step_6 = False
                            step_7 = True
                    except ValueError:
                        print('QUALITY MUST BE AN INTEGER [0-40]\n')

            if step_7 == True:
                start = ti.default_timer()
                count_passed = cqk.main(query_file, global_k, q_threshold, \
                               q_window_t, q_window)
                #print count_passed
                km.main(reference_file_set, current_lib, query_file, \
                        global_k, count_passed, q_threshold)
                stop = ti.default_timer()
                print('\t\t\t\t\tRUN TIME: {:.2f} s'.format((stop - start)))

    #Removing Reference Libraries
    elif command.lower() in ['-', 'remove', 'rm']:
        loop = True
        while loop:
            n = len(reference_file_set)
            list_ref = []
            if n == 1 or n == 0:
                print(str(n) + " REFERENCE LIBRARY")
            else:
                print(str(n) + " REFERENCE LIBRARIES")
            i = 0
            for reference in reference_file_set:
                list_ref.append(reference)
                print('[{:02}] {}'.format(i, list_ref[i]))
                i += 1
                
            arg = raw_input('REMOVE WHICH REFERENCE LIBRARY: ')
            if arg == '' or arg.lower() == 'exit':
                loop = False
            else:
                try:
                    arg = int(arg)
                    popped = list_ref.pop(arg)
                    print('{} . . . FILE REMOVED'.format(popped))
                    loop = False
                    reference_file_set = set(list_ref)
                    print
                    create_method()
                except ValueError:
                    print('VALUE MUST BE AN INTEGER\n')
                except IndexError:
                    print('INVALID SELECTION\n')

#Docstrings
__author__ = "Kevin Cheng"
__credits__ = ["Kevin Cheng", "Aaron Jex"]
__version__ = "python.1.0"
__email__ = "kkcheng.au@gmail.com"
