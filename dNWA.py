"""
Compute the optimal global alignment of two strings with NWA for all sequences in the input file.
"""

# import modules
import itertools
import argparse
import math
import logging
import re
import os
import time
import sys

# global variable definitions
number_of_alignments = 0
number_of_lines_in_infile = 0
number_of_lines_in_interfile = 0
number_of_lines_in_outfile = 0
number_of_zero_matches = 0
try:
    timestamp = time.strftime("%Y%m%d-%H%M%S")
except ValueError as v:
    print("Time format error:",v)
    sys.exit(1)
    
# GAP values
PENALTY_GAP_OPENING = -2
PENALTY_GAP_EXTENSION = -1

# Scoring Matrix
WEIGHT_MATCH_F_GROUP = 8
WEIGHT_MATCH_T_GROUP = 4
WEIGHT_MATCH_H_GROUP = 2
WEIGHT_MATCH_X_GROUP = 1
WEIGHT_MATCH_NO_FID = 0
WEIGHT_MISMATCH_FIDS = -5
WEIGHT_MISMATCH_NOFID_FID = -1

# build pointer matrix for an alignment
ptr_matrix = None

# parsing arguments from command line
parser = argparse.ArgumentParser(
    description='aligns all input sequences pairwise with a modified Needleman-Wunsch algorithm.',
    epilog='Generates two output files per default with provided filename at the end.'
    'Output format: >ClusID of seq1,accession of seq1,>ClusID of seq2,accession of seq2,'
    'Alignmentscore, length of alignment'
    )
parser.add_argument('filename', help='please provide an input file.')
parser.add_argument('-v', '--verbose',
                    help='writes a second output file that is more human-readable.',
                    action='store_true')
parser.add_argument('-l', '--log', help='set loglevel',action='store', const='INFO', nargs='?',
                    choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'])
parser.add_argument('-t', '--temp', help='keep temporary file with reformatted input',
                    action='store_true')
parser.add_argument('-a','--noselfalignment', 
                    help='Deactivate the alignment of every sequence with itself',
                    action='store_true')
parser.add_argument('-g','--gapextension',
                    help='Differentiate between gap openings and gap extensions',
                    action='store_true')
parser.add_argument('-s','--score',help='provide a score after keyword -s', action='store',
                    nargs=9, type=int)
args = parser.parse_args()

# activate logging if --log is set. Set log level if provided, otherwise set INFO as default
if args.log:
    numeric_level = args.log
    splitext_list = os.path.splitext(args.filename)
    file_extension = splitext_list[1]
    filename_wo_extension = os.path.basename(args.filename).replace(file_extension,'')
    logging.basicConfig(filename=timestamp+ '_basic_log_'+ filename_wo_extension+ '.log',
                        encoding='utf-8',
                        filemode='w',
                        format='%(levelname)s:%(message)s (%(asctime)s)',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=numeric_level)

# define functions
def convertInput():
    """
    Transform the given input into a format, that NWA can handle, and store it in a txt-file."""

    # open input file and close it in the end
    with open(args.filename, "r", encoding="utf-8") as infile:

        # Initialize variables
        clusID = None
        cluster_count = 0
        numberOfLinesRead = 0
        totalAccessions = 0
        counted_accessions = 0
        global number_of_lines_in_infile
        global number_of_lines_in_interfile
        filename = timestamp+'_'+'formatted_input_'+os.path.basename(args.filename)

        # Read the file line by line
        line = infile.readline()
        while line:
            numberOfLinesRead = numberOfLinesRead +1

            # Skip the iteration if the line is empty, but continue with the next iteration
            if line == '':
                continue

            # if line contains cluster infos, store them
            if line.startswith(">"):
                if clusID is not None:
                    if args.log:
                        logging.info('The number of read lines, that should be accessions, '
                                     'in cluster %s is %s.', str(clusID),str(counted_accessions))

                # set counter variables
                cluster_count = cluster_count+1
                counted_accessions = 0

                # If the line starts with ">", store the clusID in the clusID variable
                # find all numbers in line and store them into a list of strings
                temp = re.findall(r'\d+', line)
                # map strings in list temp to int
                numbersInLine = list(map(int, temp))
                if len(numbersInLine) < 2:
                    print('An error occured in convertInput(). The following description line is '
                          'malformed: ',line)
                    if args.log:
                        logging.critical('An error occured in convertInput(). The following '
                                         'description line is malformed: %s.',line)
                    sys.exit(1)
                # store clusID (first number)
                clusID = numbersInLine[0]
                # store number of sequences in the cluster
                totalAccessions = numbersInLine[1]
                if args.log:
                    logging.info('The number of accessions in cluster %s is indicated in the file '
                                 'as %s.',str(clusID),str(totalAccessions))

            # else split the line into parts using the clusIDentifier as separator
            # and write a transformed line into an intermediate file
            else:

                #assertion, that cluster infos are provided
                if clusID is None:
                    if args.log:
                        logging.error('Cluster clusID is empty.')
                    print('A cluster clusID is empty, please check input file.')

                # count variable for checking number of accessions
                counted_accessions = counted_accessions+1

                # stores the information of an accession in the list parts
                # example format of parts: ['#accession', '(0,0,START)', '(0,1,FID)', '(1,1,END)']
                parts = line.split(" ")

                # Check whether the input format has been complied with
                if len(parts) < 3:
                    s = ''.join(parts)
                    print('An error occured in convertInput(). The following sequence '
                          'representation is malformed: ',s)
                    if args.log:
                        logging.critical('An error occured in convertInput(). The following '
                                         'sequence representation is malformed: %s.',s)
                    sys.exit(1)
                # extracts unique proteinID of sequence
                first_entry = parts[0] # example format of first_entry: #accession
                # extracts last entry of parts
                last_entry = parts[-1] # example format of last_entry: (1,1,END)
                # proteinID = first_entry

                # store the length of the protein sequence, the x as found in (x,y,END)
                parts_length = int(last_entry.split(',')[0][1:])

                # Write the output to a file, for every sequence there are 2 lines
                try:
                    with open(filename, 'a', encoding="utf-8") as outfile:

                        # Write the first line with "clusID,name of the sequence"
                        # to store this info with the sequence
                        outfile.write('>' + str(clusID) + ',' + first_entry + '\n')

                        # declare array for the second line with zeros (as int)
                        sequence_list = [0]* parts_length

                        # Write the information of a domain into the list part to store it
                        # example format of entry: (0,10,A.A.A.A)
                        for entry in parts[2:-1]:
                            part = entry.split(',') # example format of part: ['(0', '10', 'A.A.A.A)']
                            # take the last item in the list part
                            # without the last char (closing bracket)
                            motif_number = part[-1][:-1] # example format of motif_number: A.A.A.A
                            # store the first item in the list part (start of the protein domain)
                            # but without the first char (opening bracket)
                            first_number = int(part[0][1:])
                            # store the second item in the list part (end of the protein domain)
                            second_number = int(part[1])

                        # write the domain F-FIDs into the sequence_list
                            for index in range(first_number,second_number):
                                sequence_list[index] = motif_number

                        # write data into a file
                        output_string = ','.join(map(str, sequence_list))
                        outfile.write(output_string + '\n')
                        # counter for logging purpose
                        number_of_lines_in_interfile = number_of_lines_in_interfile+1

                # handle exception
                except ValueError as v:
                    if args.log:
                        logging.critical('Input file format error. ValueError in convertInput(): %s.',v)
                    print("Input file format error. ValueError in convertInput():",v)
                    sys.exit(1)
                except IOError as e:
                    if args.log:
                        logging.critical('IOError in convertInput(): %s.',e)
                    print("IOError in convertInput():",e)
                    sys.exit(1)

            # Read the next line
            line = infile.readline()

    #store number global for checks
    number_of_lines_in_infile = numberOfLinesRead

    # write infos to log if logging activated
    if args.log:
        logging.info('The number of read lines, that should be accessions, in cluster %s is %s.',
                     str(clusID),str(counted_accessions))
        logging.info('Input transformed')
        logging.info('#cluster= %s',str(cluster_count))
        logging.info('#lines read = %s',str(numberOfLinesRead))
        if totalAccessions != counted_accessions:
            logging.warning('In cluster %s the number of accessions read does '
                            'not match the number of accessions specified in the file',str(clusID))
        if numberOfLinesRead <2:
            logging.error('Less than 2 lines read')

    # if file contains less than 2 lines it is not a useful input for the algorithm
    if numberOfLinesRead <2:
        print('The provided file contains less than 2 lines. Please check input file.')
        if args.log:
            logging.critical('The provided file contains less than 2 lines. Please check input '
                             'file.')
        sys.exit(1)

def readInput():
    """Read the file into a list and return that list."""
    filename = timestamp+'_'+'formatted_input_'+os.path.basename(args.filename)
    try:
        with open(filename, "r", encoding="utf-8") as f:

            # read file into a list (with newlines); every line of the read file is one list entry
            lines_from_file = f.readlines()

            # strip newline from entries and store entries in a list
            lines = [lines_from_file[x].strip('\n') for x in range(len(lines_from_file))]

            # read all transformed sequences into a 2D list
            # so that NWA can align all of them in all combinations
            list_with_all_sequences = [lines[i:i+2] for i in range(0, len(lines), 2)]

            return list_with_all_sequences
    except IOError as e:
        print("The intermediate file cannot be opened by readInput():",e)
        if args.log:
            logging.critical('Error in function readInput().')
        sys.exit(1)

def scoring(i,j):
    """Scoring function, assign the values for match calculation and return the value.
    
    scoring is called in scorematrix() for every position in the matrix.
    arguments:
    i -- value_list2[i-1]
    j -- value_list1[j-1]
    i and j are the counting variables of matrix and i-1 and j-1 the corresponding
    values in the sequence strings
    
    return value: the calculated score
    """

    # variables
    score = 0

    # scoring logic

    # 0 indicates "no labeled proteindomain"
    if i == '0' and j == '0':
        score = WEIGHT_MATCH_NO_FID
    # match
    elif i == j:
        score = WEIGHT_MATCH_F_GROUP
    # one sequence with NOFID, but not both
    elif i == '0' or j == '0':
        score = WEIGHT_MISMATCH_NOFID_FID
    # Partial matches of the F-IDs
    else:
        # split label to compare the different groups
        if len(i) == 0 or len(j) == 0:
            print('Empty F-ID found. Please check input file.')
            if args.log:
                logging.critical('Empty F-ID found. Please check input file.')
            sys.exit(1)
        list_with_splitted_groups_i = i.split('.')
        list_with_splitted_groups_j = j.split('.')
        if len(list_with_splitted_groups_i[0]) == 0 or len(list_with_splitted_groups_j[0]) == 0:
            print('Empty F-ID found. Please check input file.')
            if args.log:
                logging.critical('Empty F-ID found. Please check input file.')
            sys.exit(1)
        if len(list_with_splitted_groups_i)<4:
            list_with_splitted_groups_i.append(None)
        if len(list_with_splitted_groups_j)<4:
            list_with_splitted_groups_j.append(None)

        # compare the F-ID at position 0,1,2 and 3
        if list_with_splitted_groups_i[0] != list_with_splitted_groups_j[0]:
            score = WEIGHT_MISMATCH_FIDS
        elif list_with_splitted_groups_i[1] != list_with_splitted_groups_j[1]:
            score = WEIGHT_MATCH_X_GROUP
        elif list_with_splitted_groups_i[1] is None and list_with_splitted_groups_j[1] is None:
            score = WEIGHT_MATCH_X_GROUP
        elif list_with_splitted_groups_i[2] != list_with_splitted_groups_j[2]:
            score = WEIGHT_MATCH_H_GROUP
        elif list_with_splitted_groups_i[2] is None and list_with_splitted_groups_j[2] is None:
            score = WEIGHT_MATCH_H_GROUP
        elif list_with_splitted_groups_i[3] != list_with_splitted_groups_j[3]:
            score = WEIGHT_MATCH_T_GROUP
        elif list_with_splitted_groups_i[3] is None and list_with_splitted_groups_j[3] is None:
            score = WEIGHT_MATCH_T_GROUP
        else:
            print("Something went wrong with scoring(), no score could be assigned.")
    return score

def build_matrix(sequence_length1,sequence_length2):
    """
    Build a matrix of the correct length for the scorematrix function.

    arguments:
    sequence_length1 -- length of the first sequence
    sequence_length2 -- length of the second sequence
    
    Dimensions of the matrix are len(seq)+1 for each seq.
    Initiate the matrix with 0. 
    return value: matrix as 2d list
    """
    # compute matrix size
    row = sequence_length2+1
    column = sequence_length1+1

    # initialize matrix with 0
    try:
        matrix = [[0 for column in range(column)] for row in range(row)]
    except Exception as e:
        print('An error occured in build_matrix():', e)
        if args.log:
            logging.critical('Error in build_matrix(): %s', e)
        sys.exit(1)
    if args.log:
        logging.debug('build_matrix executed with a matrix of size %s x %s.',sequence_length1,sequence_length2)
    return matrix

def scorematrix(sequence_length1,sequence_length2,value_list1,value_list2):
    """
    Build and return a scorematrix.

    arguments:
    sequence_length1 -- length of the first sequence
    sequence_length2 -- length of the second sequence
    value_list1 -- first sequence as list
    value_list2 -- second sequence as list
   
    Stores computed scores in the scorematrix with dimensions
    m + 1, n + 1, if m and n are the lengths of the sequences.
    
    return value: matrix
    """

    matrix = build_matrix(sequence_length1,sequence_length2)

    #initialize matrix
    # set amount of rows
    MATRIX_COLUMN_N = sequence_length2+1 # length of the sequence in the column plus init row
    # set amount of columns
    MATRIX_ROW_N = sequence_length1+1 # length of the sequence in the row plus init column
    # i and j have to be 1 at the beginning, because M[0][0]=0
    for i in range(1,MATRIX_COLUMN_N):
        matrix[i][0] = PENALTY_GAP_OPENING * i
    for j in range(1,MATRIX_ROW_N):
        matrix[0][j] = PENALTY_GAP_OPENING * j

    # fill the score matrix
    if not args.gapextension:
        for i in range(1,MATRIX_COLUMN_N):
            for j in range(1,MATRIX_ROW_N):
                # [i-1] and [j-1] because of the start of i,j at 1 but start of the list at 0.
                # Match calculation: F(i-1,j-1) + s(a_i,b_j)
                match = matrix[i-1][j-1] + scoring(value_list2[i-1], value_list1[j-1])
                insert = matrix[i][j-1] + PENALTY_GAP_OPENING
                delete = matrix[i-1][j] + PENALTY_GAP_OPENING
                if match >= insert and match >= delete:
                    matrix[i][j] = match
                    ptr_matrix[i][j] = 1 # stores "diag"
                elif insert >= delete:
                    matrix[i][j] = insert
                    ptr_matrix[i][j] = 2 # stores "left"
                else:
                    matrix[i][j] = delete
                    ptr_matrix[i][j] = 3 # stores "up"

    # fill score matrix (gapextension version)
    if args.gapextension:
        # variables for tracing gap opening and extension
        previous_is_insertion = False
        previous_is_deletion =  [False for i in range(MATRIX_ROW_N)]

        for i in range(1,MATRIX_COLUMN_N):
            for j in range(1,MATRIX_ROW_N):
                match = matrix[i-1][j-1] + scoring(value_list2[i-1], value_list1[j-1])
                insert = matrix[i][j-1] + PENALTY_GAP_OPENING
                delete = matrix[i-1][j] + PENALTY_GAP_OPENING
                if match >= insert and match >= delete:
                    matrix[i][j] = match
                    ptr_matrix[i][j] = 1 # stores "diag" pointer
                    # store the value for gap o/e
                    previous_is_deletion[j] = False
                    previous_is_insertion = False
                elif insert >= delete:
                    # check if insertion is extending a gap or not
                    if previous_is_insertion is False:
                        insert = matrix[i][j-1] + PENALTY_GAP_OPENING
                    else:
                        insert = matrix[i][j-1] + PENALTY_GAP_EXTENSION
                    matrix[i][j] = insert
                    ptr_matrix[i][j] = 2 # stores "left" pointer
                    # store the value for gap o/e
                    previous_is_insertion = True
                    previous_is_deletion[j] = False
                else:
                    # check if deletion is extending a gap or not
                    if previous_is_deletion[j] is False:
                        delete = matrix[i-1][j] + PENALTY_GAP_OPENING
                    else:
                        delete = matrix[i-1][j] + PENALTY_GAP_EXTENSION
                    matrix[i][j] = delete
                    ptr_matrix[i][j] = 3 # stores "up" pointer
                    # store the value for gap o/e
                    previous_is_deletion[j] = True
                    previous_is_insertion = False

    if args.log:
        logging.debug('Function scorematrix() is executed, with args.gapextension = %s.',args.gapextension)

    return matrix

def traceback(name1,name2,value_list1,value_list2, matrix):
    """Go back through the matrix with information about the matches and align the sequences.

    arguments:
    name1 -- Name of sequence 1
    name2 -- Name of sequence 2
    value_list1 -- String of sequence 1
    value_list2 -- String of sequence 2
    matrix -- matrix filled with scores by scorematrix()
    
    write the alignment and additional information in a file
    """

    # initialize variables
    AlignmentA = ''
    AlignmentB = ''
    GAP_CHARACTER = '-'
    list_length1 = len(value_list1)
    list_length2 = len(value_list2)
    i = list_length2
    j = list_length1
    list_index_1 = i-1
    list_index_2 = j-1
    global number_of_lines_in_outfile
    global number_of_zero_matches
    # compute alignment as in Needleman-Wunsch algorithm
    while i>0 or j>0:
        if i>0 and j > 0 and ptr_matrix[i][j] == 1: # match, go diag
            AlignmentA = str(value_list2[list_index_1]) + ',' + AlignmentA
            AlignmentB = str(value_list1[list_index_2]) + ',' + AlignmentB
            if value_list2[list_index_1] == '0' and value_list1[list_index_2] == '0':
                number_of_zero_matches += 1
            list_index_1 = list_index_1-1
            i = i - 1
            list_index_2 = list_index_2-1
            j = j - 1
        elif i>0 and (ptr_matrix[i][j] == 3 or j==0): # deletion or edge of the matrix
            AlignmentA = str(value_list2[list_index_1]) + ',' + AlignmentA
            AlignmentB = GAP_CHARACTER + ',' + AlignmentB
            list_index_1 = list_index_1-1
            i = i - 1
        elif j>0 and (ptr_matrix[i][j] == 2 or i==0): # insertion or edge of the matrix
            AlignmentA = GAP_CHARACTER + ',' + AlignmentA
            AlignmentB = str(value_list1[list_index_2]) + ',' + AlignmentB
            j = j - 1
            list_index_2 = list_index_2-1
        else:
            print('error in traceback')
            if args.log:
                logging.critical('error in traceback')
            sys.exit(1)

    # prepare Strings for writing
    AlignmentA = AlignmentA.strip(',')
    AlignmentB = AlignmentB.strip(',')
    score_alignment = str(matrix[list_length2][list_length1])

    # build lists for verbose output and length computation
    list_of_AlignmentA = AlignmentA.split(',')
    list_of_AlignmentB = AlignmentB.split(',')
    alignmentlength = len(list_of_AlignmentA)

    # normalize the alignment score over length of alignment
    # (divide score_alignment by length of alignment minus number of "0 matches")
    try:
        score_alignment_normalized = int(score_alignment) / (int(alignmentlength-number_of_zero_matches))
    except ZeroDivisionError as e:
        print('An error occurred in traceback(): Sequence without protein domain '
              'found. Please run the program again with reasonable input.')
        if args.log:
            logging.critical('An error occured in traceback(): Sequence without '
                             'protein domain found. Please run the program again with reasonable '
                             'input. Error: %s.',e)
        sys.exit(1)
        

    # build a string with names, score and length of alignment, set filename
    first_line = name1 + ',' + name2 + ',' + score_alignment + ',' + str(alignmentlength) + ',' + str(score_alignment_normalized)
    filename = timestamp+'_'+'alignments_'+os.path.basename(args.filename)

    # write 3 lines in an outputfile (for every alignment)
    with open(filename, 'a', encoding="utf-8") as file_alignments:
        file_alignments.write(first_line + '\n' + AlignmentB + '\n' + AlignmentA + '\n')

    # write log info
    number_of_lines_in_outfile = number_of_lines_in_outfile+3

    # check if option -verbose is set and if so, create the second output file
    if args.verbose:
        writeHumanreadableOutput(name1,name2,score_alignment,alignmentlength, 
                                list_of_AlignmentA,list_of_AlignmentB)
    number_of_zero_matches = 0

def needleman_wunsch(subset):
    """Function with all calls for NWA, except selfalignments."""

    # global variables
    global number_of_alignments
    global ptr_matrix

    #input variables
    name1 = subset[0][0]
    name2 = subset[1][0]
    value_string1 = subset[0][1]
    value_string2 = subset[1][1]
    value_list1 = value_string1.split(',')
    value_list2 = value_string2.split(',')

    sequence_length1 = len(value_list1)
    sequence_length2 = len(value_list2)

    # build scorematrix for this alignment
    ptr_matrix = build_matrix(sequence_length1,sequence_length2)
    matrix = scorematrix(sequence_length1,sequence_length2,value_list1,value_list2)

    # compute traceback for this alignment
    traceback(name1,name2,value_list1,value_list2, matrix)

    # global count for logging
    number_of_alignments = number_of_alignments+1

def needleman_wunschSelf(entry):
    """Function with all selfalingment calls for NWA."""

    # global variables
    global number_of_alignments
    global ptr_matrix

    #input variables
    name1 = entry[0]
    name2 = name1
    value_string = entry[1]
    value_list1 = value_string.split(',')
    value_list2 = value_list1

    sequence_length1 = len(value_list1)
    sequence_length2 = sequence_length1

    # build scorematrix for this alignment
    ptr_matrix = build_matrix(sequence_length1,sequence_length2)
    matrix = scorematrix(sequence_length1,sequence_length2,value_list1,value_list2)

    # compute traceback for this alignment
    traceback(name1,name2,value_list1,value_list2, matrix)

    # global count for logging
    number_of_alignments = number_of_alignments+1

def alignAllPairs(list_with_all_sequences):
    """Calls NWA for every possible pair of sequences."""

    for subset in itertools.combinations(list_with_all_sequences, 2):
        needleman_wunsch(subset)

    # if it is not explicitly stated that selfalignments should not be calculated,
    # calculate selfalignments.
    if not args.noselfalignment:
        if args.log:
            logging.info('Compute self-alignments. '
                         'For an output file without self-alignments, check program options.'
                         )
        for entry in list_with_all_sequences:
            needleman_wunschSelf(entry) # appended to the entries of needleman_wunsch()

def writeHumanreadableOutput(name1,name2,score_alignment,alignmentlength,
list_of_AlignmentA,list_of_AlignmentB):
    """Writes a more human readable output if -v is active and inserts padding into the alignment."""

    # variables
    entry_splitted=[]
    transformedListA=[0]*alignmentlength
    transformedListB=[0]*alignmentlength
    iterate_A = 0
    iterate_B = 0

    # add padding to both Strings
    for entry in list_of_AlignmentB:
        if entry == '0':
            transformedListB[iterate_B] = "0000.0000.0000.0000"
        elif entry == "-":
            transformedListB[iterate_B] = "----.----.----.----"
        else:
            entry_splitted = entry.split('.')
            split_list = [0]*len(entry_splitted)
            j = 0
            for label in entry_splitted:
                split_list[j] = label.ljust(4, '_')
                j = j+1

            a_string = '.'.join(split_list)
            if len(split_list)==3:
                a_string = a_string + '.____'
            elif len(split_list)==2:
                a_string = a_string + '.____.____'
            elif len(split_list)==1:
                a_string = a_string + '.____.____.____'
            transformedListB[iterate_B] = a_string
            entry_splitted=[]
        iterate_B = iterate_B +1

    for entry in list_of_AlignmentA:
        if entry == '0':
            transformedListA[iterate_A] = "0000.0000.0000.0000"
        elif entry == "-":
            transformedListA[iterate_A] = "----.----.----.----"
        else:
            entry_splitted = entry.split('.')
            split_list = [0]*len(entry_splitted)
            j = 0
            for label in entry_splitted:
                if len(label)>4:
                    logging.warning(
                        "Found a label (F-ID) with a group with len >4, code asserts len 4,"
                        "please update code accordingly."
                        )
                split_list[j] = label.ljust(4, '_')
                j = j+1

            a_string = '.'.join(split_list)
            if len(split_list)==3:
                a_string = a_string + '.____'
            elif len(split_list)==2:
                a_string = a_string + '.____.____'
            elif len(split_list)==1:
                a_string = a_string + '.____.____.____'
            transformedListA[iterate_A] = a_string
            entry_splitted=[]
        iterate_A = iterate_A +1

    # write sequences as strings
    hAlignmentA = ','.join(map(str, transformedListA))
    hAlignmentB = ','.join(map(str, transformedListB))

    # normalize the alignment score over length of alignment
    #(divide score_alignment by length of alignment minus number of "0 matches")
    try:
        score_alignment_normalized = int(score_alignment) / (int(len(list_of_AlignmentA)-number_of_zero_matches))
    except ZeroDivisionError as e:
        print('An error occured: Only 0 but no protein domain was found in an alignment.'
        'Please run the program again with reasonable input.')
        if args.log:
            logging.error('An error occured in writeHumanreadableOutput(): %s.',e)
        sys.exit(1)
    # build annotation string with names and score
    first_line = name1 + ',' + name2 + ',' + score_alignment + ',' + str(score_alignment_normalized)

    # write data in an output file
    with open(timestamp+'_'+'file_alignments_verbose_'+ os.path.basename(args.filename), 'a', encoding="utf-8") as file_alignments:
        file_alignments.write(first_line + '\n' + hAlignmentB + '\n' + hAlignmentA + '\n')

    if args.log:
        logging.debug('Additional file with better humanreadable output created.')

def setScore():
    """Check if user provided a list with scores and set the score values accordingly."""

    # Scoring Matrix
    global WEIGHT_MATCH_F_GROUP
    global WEIGHT_MATCH_T_GROUP
    global WEIGHT_MATCH_H_GROUP
    global WEIGHT_MATCH_X_GROUP
    global WEIGHT_MATCH_NO_FID
    global WEIGHT_MISMATCH_FIDS
    global WEIGHT_MISMATCH_NOFID_FID

    # Gap penalties
    global PENALTY_GAP_OPENING
    global PENALTY_GAP_EXTENSION

    # set values
    WEIGHT_MATCH_F_GROUP = args.score[0]
    WEIGHT_MATCH_T_GROUP = args.score[1]
    WEIGHT_MATCH_H_GROUP = args.score[2]
    WEIGHT_MATCH_X_GROUP = args.score[3]
    WEIGHT_MATCH_NO_FID = args.score[4]
    WEIGHT_MISMATCH_FIDS = args.score[5]
    WEIGHT_MISMATCH_NOFID_FID = args.score[6]
    PENALTY_GAP_OPENING = args.score[7]
    PENALTY_GAP_EXTENSION = args.score[8]

    if args.log:
        logging.info('The scores were set to %s by input.', str(args.score))


def checkInput():
    """Check if given file exists, is a file and not empty."""

    # check if given file exists
    file_exists = os.path.exists(args.filename)
    if not file_exists:
        logging.critical('The input file does not exist.')
    assert file_exists, "The provided file does not exist."

    # check if given path points to a file
    path_is_file = os.path.isfile(args.filename)
    if not path_is_file:
        logging.critical('The specified path does not point to a file.')
    assert path_is_file, "The specified path does not point to a file."

    # check if given file is not empty
    if os.stat(args.filename).st_size == 0:
        logging.critical('The file is empty.')
    assert os.stat(args.filename).st_size > 0, "The file is empty."

def deleteTempFiles():
    """Delete intermediate file."""

    if os.path.exists(timestamp+'_'+'formatted_input_'+ os.path.basename(args.filename)):
        os.remove(timestamp+'_'+'formatted_input_'+ os.path.basename(args.filename))
    else:
        print('No file found to delete')
        if args.log:
            logging.error('The intermediate file could not be found.')

def logStats():
    """Logs some basic metrics."""

    logging.info('The number of lines read in the input file is %s.',
                 str(number_of_lines_in_infile))
    logging.info('The number of lines written to the intermediate file is %s'
                 ', the number of lines in the output file is %s.'
                 ,str(number_of_lines_in_interfile), str(number_of_lines_in_outfile)
                 )
    logging.info('Number of alignments processed is %s.'
                 ,str(number_of_alignments))
    binomial = math.comb(number_of_lines_in_interfile,2)
    logging.info('Number of alignments should be %s plus #sequences (= %s, for selfalignments)'
                 ' = %s.'
                 ,str(binomial)
                 ,str(number_of_lines_in_interfile)
                 ,str(binomial+number_of_lines_in_interfile)
                 )
    expected_lines_in_alignment = (number_of_alignments*2)+ number_of_alignments
    logging.info('Number of lines in out file (= %s'
                 ') should match #alignments times two + #alignments (= %s).'
                 ,str(number_of_lines_in_outfile), str(expected_lines_in_alignment)
                 )
    if number_of_lines_in_outfile != expected_lines_in_alignment:
        logging.warning('There are not the expected number of lines in the alignment file.')
    if number_of_alignments != binomial+number_of_lines_in_interfile:
        logging.warning('The algorithm did not calculate the expected number of alignments.')

def main():
    """Main function which calls all functions"""
    try:
        # check if input file is okay
        checkInput()

        # convert the given input in proper input for needleman wunsch algorithm
        convertInput()

        # check if alternative scores are provided and if yes, store them
        if args.score:
            setScore()

        # read input data and align all sequences
        alignAllPairs(readInput())

        # check if option -temp is set and if so, delete file with transformed input
        if not args.temp:
            deleteTempFiles()

        # log info with gapextension value
        if args.log:
            if args.gapextension:
                infostring = (
                    'Different penalties for gap extensions and openings are used. '
                    'Gap extension penalty is '
                    + str(PENALTY_GAP_EXTENSION) + ' and gap opening penalty is '
                    + str(PENALTY_GAP_OPENING) + '.'
                )
                logging.info(infostring)

        # write log with basic stats
        if args.log:
            logStats()
    except Exception as e:
        print('An error occured in main():', e)
        if args.log:
            logging.critical('An error occured in main(): %s.',e)
        sys.exit(1)

# call main function to start program
main()