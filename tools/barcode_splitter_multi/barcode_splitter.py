#!/usr/bin/env python
""" Split one or more fastq files based on barcode sequence.
"""
from __future__ import division
import gzip
import argparse
import os
import re
import sys
from collections import defaultdict
import subprocess

__version__      = "0.14"
__author__       = "Lance Parsons & Robert Leach"
__author_email__ = "lparsons@princeton.edu,rleach@princeton.edu"
__copyright__    = "Copyright 2011, Lance Parsons & Robert leach"
__license__      = ("BSD 2-Clause License "
                    "http://www.opensource.org/licenses/BSD-2-Clause")

UNMATCHED    = 'unmatched'
MATCHED      = 'matched'
MULTIMATCHED = 'multimatched'
ANY          = 'any'       #Used in the N-dimensional dictionary to help
                           #diagnose which barcode set is matching/not matching

def main (argv=None):
    if argv is None:
        argv = sys.argv
    
    parser = argparse.ArgumentParser(description=globals()['__doc__'])
    parser.add_argument('--version', action='version',
                        version='%(prog)s version ' + globals()['__version__'])

    required_group = parser.add_argument_group("Barcodes")
    required_group.add_argument('--bcfile', metavar='FILE',
                                help='REQUIRED: Tab delimited file: '
                                '"Sample_ID <tab> Barcode_Sequence"')
    required_group.add_argument ('--idxread', metavar='READNUM', type=int,
                                 nargs='+', help="REQUIRED: Indicate in "
                                 "which read file(s) to search for the "
                                 "corresponding column of barcode sequences, "
                                 "e.g. if the first column of barcodes is in "
                                 "the second sequence read file and the "
                                 "second column's barcodes are in the third "
                                 "sequence read file, you'd supply "
                                 "`--idxread 2 3`")
    required_group.add_argument('--mismatches', default=0, type=int,
                                help='Number of mismatches allowed in '
                                'barcode matching')
    required_group.add_argument('--barcodes_at_end', action='store_true',
                                default=False, help='Barcodes are at the end '
                                'of the index read (default is at the '
                                'beginning)')

    output_group = parser.add_argument_group("Output Options")
    output_group.add_argument('--prefix', default='',
                              help='Prefix for output files')
    output_group.add_argument('--suffix', default=None,
                              help='Suffix for output files (default based '
                              'on --format)')
    output_group.add_argument('--galaxy', action='store_true', default=False,
                              help='Produce "Galaxy safe" filenames by '
                              'removing underscores (default: %(default)s)')
    output_group.add_argument('-v', '--verbose', action='store_true',
                              default=False, help='verbose output')
    output_group.add_argument('--gzipout', action='store_true', default=False,
                              help='Output files in compressed gzip format '
                              '(default is uncompressed)')

    input_group = parser.add_argument_group("Input format")
    input_group.add_argument('--format', default='fastq',
                             help='Specify format for sequence files (fasta '
                             'or fastq)')
    input_group.add_argument('--gzipin', action='store_true', default=False,
                             help='Assume input files are in gzip format, '
                             'despite file extension (default is auto based '
                             'on input file extension)')

    seqs_group = parser.add_argument_group("Sequence Inputs")
    seqs_group.add_argument('fastq_files', metavar='FILE', type=str, nargs='+',
                            help='A series of 1 or more [optionally zipped] '
                            'fastq files.')

    #Error-check the values provided from the command line
    try:
        options = parser.parse_args()
	if options.fastq_files is None or options.idxread is None:
	    parser.error('Sequence files and at least one number indicating '
                         'the indexed file(s) (--idxread) is required')
        if len(options.fastq_files) < 1:
            parser.error('Must specify at least one sequence file')
        if len(options.fastq_files) < len(options.idxread):
            parser.error('Must specify at least one sequence file')
        if not options.bcfile:
            parser.error('Must specify a barcodes file with "--bcfile" option')
        if ((min(options.idxread) < 1) or
            (max(options.idxread) > len(options.fastq_files))):
            parser.error('Invalid index read number ("--idxread"), must be '
                         'between 1 and %s (the number of supplied sequence '
                         'files)' % len(options.fastq_files))
    except SystemExit: # Prevent exit when calling as function
        return 2

    #Read barcodes files into n-dimensional dictionary
    (barcode_dict,approx_bc_dict) = read_barcodes(options.bcfile,
                                                  len(options.idxread),
                                                  options.mismatches)
    barcode_sizes = getBarcodeSizes(barcode_dict)
    #print(barcode_dict)
    #print(approx_bc_dict)

    #initialize the read counts for each barcode
    total_read_count = 0
    matched_counts   = multi_dimensions(len(options.idxread), int)
    unmatched_counts = multi_dimensions(len(options.idxread), int)
    #TODO Verbose: print barcode_dict

    #Determine if we should use gzip for input
    basename, extension = os.path.splitext(options.fastq_files[0])
    if extension == '.gz':
        options.gzipin = True

    #Set the suffix (before determination of gzip output mode)
    if options.suffix is not None:
        suffix = options.suffix
    elif options.format is not None:
        suffix = '.%s' % options.format
    elif extension != '.gz':
        suffix = '.%s' % extension
    else:
        suffix = 'fastq'

    #Determine if we should use gzip for output
    if (options.gzipout is True or
        (options.suffix is not None and options.suffix == '.gz')):
        options.gzipout = True
        if options.suffix is None or suffix != '.gz':
            suffix = '%s.gz' % suffix

    #Determine whether to use the gzip command line tool or not
    use_gzip_cmd_line = False
    if options.gzipin is True or options.gzipout is True:
        use_gzip_cmd_line = exeExists('gzip')

    #Determine if we need to edit the prefix for galaxy mode
    prefix = options.prefix
    if options.galaxy:
        prefix = options.prefix.replace("_","-")
        
    #Open input filehandles for each read
    #print "Opening the input file handles..."
    inputs = {}
    for i in xrange(0,len(options.fastq_files)):
        if options.gzipin:
            if(use_gzip_cmd_line is True):
                inputs[i] = run_gzip_in(options.fastq_files[i])
            else:
                inputs[i] = gzip.open(options.fastq_files[i], 'rb')
        else:
            inputs[i] = open(options.fastq_files[i], 'rb')

    #Open output filehandles for each barcode set and sequence file
    #print "Opening the output barcode file handles..."
    outputs = multi_dimensions((len(options.idxread) + 1), file)
    openOutfiles(outputs, barcode_dict, prefix, suffix, inputs, options.galaxy,
                 options.gzipout,use_gzip_cmd_line)

    #Open the file handle for the unmatched file
    #print "Opening the output un/multi-matched file handles..."
    unmatchedOutputs    = defaultdict(file)
    multimatchedOutputs = defaultdict(file)
    for i in xrange(0,len(inputs)):
        unmf = '%s%s-read-%s%s' %(prefix, UNMATCHED, i+1, suffix)
        mmf  = '%s%s-read-%s%s' %(prefix, MULTIMATCHED, i+1, suffix)
	if options.gzipout:
            if use_gzip_cmd_line is True:
                unmatchedOutputs[i] = run_gzip_out(unmf)
                multimatchedOutputs[i] = run_gzip_out(mmf)
            else:
                unmatchedOutputs[i] = gzip.open(unmf,'wb')
                multimatchedOutputs[i] = gzip.open(mmf, 'wb')
	else:
	    unmatchedOutputs[i] = open(unmf, 'wb')
	    multimatchedOutputs[i] = open(mmf, 'wb')

    readidxs     = [i - 1 for i in options.idxread]
    prim_readidx = readidxs[0]

    #Debug mode reads the 1st 10000 reads of each file when prefix has "debug"
    debug_mode = False
    if re.match('debug', prefix):
	debug_mode = True

    # For each input line in index read, get index sequence
    for prim_index_read in read_fastq(inputs[prim_readidx]):

        total_read_count += 1

        #print "Processing read number %s\r" % (total_read_count),

        #For debugging any file
        if debug_mode and total_read_count > 10000:
            break

        #Determine the ID format using the first record
        if total_read_count == 1:
            id_format = determine_id_format(prim_index_read['seq_id'][1:])

	#Get 1 from each file & assert that their IDs match
	cur_reads = {}
	for read_index in xrange(0,len(inputs)):
	    if(read_index is not prim_readidx):
		cur_reads[read_index] = read_fastq(inputs[read_index]).next()
                #If the IDs don't match between all the files, issue a warning
                try:
                    assert(match_id(prim_index_read['seq_id'],
                                    cur_reads[read_index]['seq_id'],id_format))
                except AssertionError:
                    sys.stderr.write("WARNING: Index ID mismatch: %s does "
                                     "not match %s\n"
                                     %(prim_index_read['seq_id'],
                                       cur_reads[read_index]['seq_id']))
	    else:
		cur_reads[read_index] = prim_index_read

        #Get the index sequences from the indexed reads
        index_seqs = []
        if options.barcodes_at_end:
            index_seqs = [cur_reads[readidxs[i]]['seq'][-barcode_sizes[i]:]
                          for i in xrange(0,len(readidxs))]
        else:
            index_seqs = [cur_reads[readidxs[i]]['seq'][0:barcode_sizes[i]]
                          for i in xrange(0,len(readidxs))]

        barcode_path = getNDDictVal(approx_bc_dict,index_seqs)

        if(barcode_path is None):
            cur_outputs = unmatchedOutputs
            unmatched_path = getBarcodeMatchPath(approx_bc_dict,index_seqs)
	    incrementNDDictInt(unmatched_counts,unmatched_path)
            if (not UNMATCHED in unmatched_path and
                not MULTIMATCHED in unmatched_path):
                sys.stderr.write('WARNING: Sequences match barcodes on '
                                 'different rows: %s for sequence ID: %s\n'
                                 %(index_seqs, prim_index_read['seq_id']))
        elif(MULTIMATCHED in barcode_path):
            cur_outputs = multimatchedOutputs
            unmatched_path = getBarcodeMatchPath(approx_bc_dict,barcode_path)
	    incrementNDDictInt(unmatched_counts,unmatched_path)
            sys.stderr.write('WARNING: More than one barcode matches %s, '
                             'moving to %s category\n'
                             %(prim_index_read['seq_id'], MULTIMATCHED))
        else:
            cur_outputs = getNDDictVal(outputs,barcode_path)
	    incrementNDDictInt(matched_counts,barcode_path)

        #Write each sequence to the matched or unmatched output handle
        for readnum in xrange(0,len(inputs)):
            cur_outputs[readnum].write(fastq_string(cur_reads[readnum]))

    #Report the final matched/unmatched counts in a table to STDOUT
    printCounts(barcode_dict,matched_counts,unmatched_counts,total_read_count,
                len(options.idxread))

    return 0

def read_barcodes(filename, num_dims, mismatches):
    '''Read barcodes file into dictionary'''

    #Declare a variable-dimension dictionary
    #This was previously a static multi-dimensional dictionary, but it would
    #not have worked...
    l = lambda:defaultdict(l)
    real_bc_dict = l()
    approx_bc_dict = l()
    approx_bc_any_dict = l()

    linenum = 0
    filehandle = open(filename, 'rb')
    for line in filehandle:
        
        linenum += 1
        line = line.strip()
        cols = []

        #If the line has something on it and isn't commented
        if line and line[0] != '#':
            
            cols = re.split('\t+',line)

            #If we got the expected number of columns
            if len(cols) == (num_dims + 1):
                sample_id         = cols.pop(0)
                barcode_sequences = cols
            else:
                sys.stderr.write("Unable to parse line in barcode file: [%s]: "
                                 "[%s].  Must contain %s columns (the number "
                                 "of values supplied to --idxread plus 1 for "
                                 "the sample ID), but found [%s].  Skipping.\n"
                                 %(filename, line, (num_dims + 1), len(cols)))
                continue

            #If the sample ID is not none and the number of defined barcodes is
            #equal to the number of expected barcodes
            defined_barcodes = [e for e in barcode_sequences if e is not None]
            if (sample_id is not None) and len(defined_barcodes) == num_dims:

                #Set the sample ID in the N-Dimensional dictionary
                setNDDictVal(real_bc_dict, defined_barcodes, sample_id)

                #This creates a dictionary containing all possible barcode
                #combos, with MULTIMATCHED values pre-computed
                #print ("Calling setNDApproxDictVal with %s and dict: "
                #       "[%s]") % (defined_barcodes, approx_bc_dict)
                setNDApproxDictVal(approx_bc_dict,
                                   defined_barcodes,
                                   mismatches)

                #Make sure we can get all the barcodes for each dict dimension
                dim = 0
                for barcode in defined_barcodes:
                    dim += 1
                    real_bc_dict[ANY][str(dim)][barcode] = ""

                    #This creates a dictionary containing all possible barcode
                    #combos, with MULTIMATCHED values pre-computed
                    setNDApproxDictVal(approx_bc_any_dict[ANY][str(dim)],
                                       [barcode],
                                       mismatches)

            elif ((sample_id is not None) and
                  (len(defined_barcodes) > num_dims)):
                raise Exception("More barcode indexes were found in the "
                                "barcodes file on line %s: '%s' than were "
                                "expected: [%s] (the number of values "
                                "supplied to --idxread)."
                                %(linenum, line, num_dims))

            else:
                raise Exception("Unable to parse barcode(s) from line %s: "
                                "'%s'. Expected a sample ID followed by [%s] "
                                "tab-delimited barcodes"
                                %(linenum, line, num_dims))

    if len(real_bc_dict) == 0:
	raise Exception("Unable to parse any barcodes from barcode file [%s]."
                        %(filename))

    #print(approx_bc_dict)

    reduceNDDictPaths(approx_bc_dict)

    for dim in real_bc_dict[ANY]:
        reduceNDDictPaths(approx_bc_any_dict[ANY][dim])
        approx_bc_dict[ANY][dim] = approx_bc_any_dict[ANY][dim]

    return(real_bc_dict,approx_bc_dict)

def openOutfiles(outputs,barcode_dict,prefix,suffix,inputs,galaxy,gzip_mode,
                 use_gzip_cmd_line):
    '''Opens all the output files for all barcode matches (unmatched
    barcode output files not opened)'''
    if isinstance(barcode_dict,dict):
	for barcode in getRealBarcodes(barcode_dict.keys()):
            openOutfiles(outputs[barcode], barcode_dict[barcode], prefix,
                         suffix, inputs, galaxy, gzip_mode, use_gzip_cmd_line)
    else:
        sample_id = barcode_dict
        if galaxy:
            #Replace underscore to allow this to work with Galaxy
            sample_id = barcode_dict.replace("_","-")
        for i in xrange(0,len(inputs)):
            fn = '%s%s-read-%s%s' % (prefix, sample_id, i+1, suffix)
            if gzip_mode:
                if use_gzip_cmd_line is True:
                    outputs[i] = run_gzip_out(fn)
                else:
                    outputs[i] = gzip.open(fn, 'wb')
            else:
                outputs[i] = open(fn, 'wb')

def getBarcodeSizes(barcode_dict):
    '''Uses the first key in each level of a dictionary to build a list of key
    sizes (assuming all keys are the same size)'''

    if (isinstance(barcode_dict,dict) and len(barcode_dict.keys()) > 0 and
        isinstance(barcode_dict[barcode_dict.keys()[0]],dict) and
        len(barcode_dict[barcode_dict.keys()[0]].keys()) > 0):

	return([len(barcode_dict.keys()[0])] +
               getBarcodeSizes(barcode_dict[barcode_dict.keys()[0]]))

    elif isinstance(barcode_dict,dict) and len(barcode_dict.keys()) > 0:

	return([len(barcode_dict.keys()[0])])

    else:

        return(None)

def setNDDictVal(in_dict, keys_list, val):
    '''Sets the value of an N-Dimensional dictionary to the supplied value
    using the list of keys'''
    cur_dict = in_dict
    for cur_key in keys_list[0:-1]:
	cur_dict = cur_dict[cur_key]
    cur_dict[keys_list[-1]] = val

def incrementNDDictInt(in_dict, keys_list):
    '''Sets the value of an N-Dimensional dictionary to the supplied value
    using the list of keys'''
    cur_dict = in_dict
    for cur_key in keys_list[0:-1]:
	cur_dict = cur_dict[cur_key]
    cur_dict[keys_list[-1]] += 1

def getNDDictVal(in_dict, keys_list):
    '''Gets the value of an N-Dimensional dictionary using the list of keys'''
    cur_dict = in_dict
    for cur_key in keys_list:
        if cur_key in cur_dict:
            cur_dict = cur_dict[cur_key]
        else:
            return None
    return cur_dict

def multi_dimensions(n, type):
    """Creates an N-dimensional dictionary containing values at the leaves of
    type 'type'. E.g. A 2-dimensional dictionary of strs would have 2 levels of
    keys that hold string values (mydict['key1']['key2'] = 'my string value')
    """

    #n <= 0 is correct. E.g. if called with 2 (i.e. a 2D dictionary - like a 2D
    #array/matrix - where array[2][1] = 5 or dict['key1']['key2'] = 'str' are
    #2-dimensional data structures), the the first type returned is a dict (for
    #holding 'key1') and the second thing returned is a defaultdict for holding
    #'key2' and the last thing returned is the type of the data held at the
    #leaves of the dictionary
    if n <= 0:
        return type()

    return defaultdict(lambda:multi_dimensions(n-1, type))

def getRealBarcodes(barcodes):
    '''Return a list of barcodes that exclude generic "matched", "unmatched",
    and "any" pseudo-barcodes'''
    return [bc for bc in barcodes if (bc is not UNMATCHED and bc is not MATCHED
                                      and bc is not ANY and bc is not
                                      MULTIMATCHED and bc is not None)]

def matchBarcodes(sequence, barcodes, mismatches):
    '''Find closest match(es) in barcodes to specified sequence with max number
    of mismatches'''
    best_distance = mismatches
    results = []
    for barcode in barcodes:
        if mismatches == 0:
            if (sequence == barcode):
                results.append(barcode)
        else:
            distance = hamming_distance(sequence, barcode)
            if (distance <= best_distance):
                best_distance = distance
                results.append(barcode)
    return results
    
def hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))    

''' Supported types of Fastq IDs '''
ILLUMINA = 'illumina' # CASVA 1.8+, match up to space
STRIPONE = 'stripone' # Illumina CASAVA 1.7 and lower (/1 /2) and NCBI SRA
                      #   (/f /r), match all but last character
OTHER    = 'other'    # Other, match exactly

def determine_id_format(seq_id):
    '''Determine if the id is new illumina, old illumina (/1 /2 ...), sanger
    (/f /r), or other'''
    
    id_format = None
    # Illumina CASAVA 1.8+ fastq headers use new format
    read_id_regex = re.compile(r'(?P<instrument>[a-zA-Z0-9_-]+):'
                               '(?P<run_number>[0-9]+):'
                               '(?P<flowcell_id>[a-zA-Z0-9]+):'
                               '(?P<lane>[0-9]+):(?P<tile>[0-9]+):'
                               '(?P<x_pos>[0-9]+):'
                               '(?P<y_pos>[0-9]+) (?P<read>[0-9]+):'
                               '(?P<is_filtered>[YN]):'
                               '(?P<control_number>[0-9]+):'
                               '(?P<index_sequence>[ACGT]+){0,1}')
    # Old illumina and sanger reads use /1 /2 /3 or /f /r
    strip_one_endings = ['/1', '/2', '/3', '/f', '/r']
    
    if read_id_regex.match(seq_id):
        id_format = ILLUMINA
    elif (seq_id[-2:] in strip_one_endings):
        id_format = STRIPONE
    else:
        id_format = OTHER
    return id_format

def strip_read_from_id(seq_id, id_format=None):
    new_id = seq_id
    if not id_format:
        id_format = determine_id_format(seq_id)
    elif id_format == STRIPONE:
        new_id = seq_id[0:-1]
    elif id_format == ILLUMINA:
        new_id = seq_id.split(' ')[0]
    return new_id

def strip_read_from_id_stripone(seq_id):
    return seq_id[0:-1]

def strip_read_from_id_illumina(seq_id):
    return seq_id.split(' ')[0]

def match_id(id1, id2, id_format=OTHER):
    ''' Return true if IDs match using rules for specified format '''
    if id_format == STRIPONE:
        if id1[0:-1] == id2[0:-1]:
            return True
        else:
            return False
    elif id_format == ILLUMINA:
        if (id1.split(' ')[0] == id2.split(' ')[0]):
            return True
        else:
            return False
    elif id1 == id2:
        return True
    else:
        return False

def read_fastq(filehandle):
    ''' Return dictionary with "seq_id", "seq", "qual_id", and "qual" '''
    record_line = 0
    read_number = 0
    fastq_record = dict()
    for line in filehandle:
        record_line += 1
        if record_line == 1:
            fastq_record['seq_id'] = line.strip()
        elif record_line == 2:
            fastq_record['seq'] = line.strip()
        elif record_line == 3:
            fastq_record['qual_id'] = line.strip()
        elif record_line == 4:
            record_line = 0
            fastq_record['qual'] = line.strip()
            read_number += 1
            yield fastq_record

def fastq_string(record):
    return "%s\n%s\n%s\n%s\n" % (record['seq_id'], record['seq'],
                                 record['qual_id'], record['qual'])

def printCounts(barcode_dict,matched_counts,unmatched_counts,total_read_count,
                numBarcodes):
    '''Prints the header of the counts table and calls the 2 recursive
    functions for printing barcode and unmatched barcode count data'''

    #Print the table header
    print "Sample\t",
    for b in xrange(1,numBarcodes + 1):
	print "Barcode%s\t" %(b),
    print "Count\tPercent"

    if total_read_count == 0:
        print "ERROR: Total read count is 0.  No reads were parsed."
        return()

    #Print the matched counts, sorting the lines
    for line in sorted(getMatchedCountsTable(barcode_dict,matched_counts,
                                             total_read_count,[])):
        print(line),

    #Print the unmatched counts (default sorted by string)
    printUnmatchedCounts(unmatched_counts,total_read_count,[])

def getMatchedCountsTable(barcode_dict,matched_counts,total_read_count,
                          barcode_path):
    '''Recursively builds a "barcode path" and prints a line in the barcode
    counts table for every leaf in the barcode dictionary'''

    ret_list = []
    if isinstance(barcode_dict,dict):
	for barcode in sorted(barcode_dict, key=barcode_dict.get):
	    if barcode is not ANY:
		if len(barcode_path) > 0:
		    new_barcode_path = barcode_path[:]
		else:
		    new_barcode_path = []
		new_barcode_path.append(barcode)
		ret_list += getMatchedCountsTable(barcode_dict[barcode],
                                                  matched_counts[barcode],
                                                  total_read_count,
                                                  new_barcode_path)
    else:
	ret_str = barcode_dict + "\t"
	for bc in barcode_path:
	    ret_str += bc + "\t"
        ret_str += "%s\t%.2f%%\n" % (matched_counts,
                                     (matched_counts/total_read_count)*100)
        ret_list += [ret_str]

    return(ret_list)

def printUnmatchedCounts(unmatched_counts,total_read_count,barcode_path):
    '''Recursively builds a "barcode path" (consisting of "MATCHED" and
    "UNMATCHED" values) and prints a line in the barcode counts table for every
    leaf in the unmatched counts dictionary'''
    if isinstance(unmatched_counts,dict):
	for barcode in sorted(unmatched_counts, key=unmatched_counts.get):
	    if len(barcode_path) > 0:
		new_barcode_path = barcode_path[:]
	    else:
		new_barcode_path = []
	    new_barcode_path.append(barcode)
	    printUnmatchedCounts(unmatched_counts[barcode],total_read_count,
                                 new_barcode_path)
    else:
        label = UNMATCHED
        if (UNMATCHED not in barcode_path and MULTIMATCHED in barcode_path):
            label = MULTIMATCHED
        print "%s\t" %(label),
	for bc in barcode_path:
	    print "%s\t" %(bc),
	print "%s\t%.2f%%" %(unmatched_counts,
                             (unmatched_counts/total_read_count)*100)

def getBarcodeMatchPath(in_dict, keys_list):
    '''Generates a path of pseudo-barcodes (consisting of "MATCHED",
    "UNMATCHED", and "MULTIMATCHED" values) to use in the summary output.  The
    intent is to give the user useful feedback about which level(s) of barcodes
    are not matrching instead of listing all sequences as simply "UNMATCHED".
    It also reduces complexity of the unmatched output by not including
    specific barcodes, which do not matter in unmatched/multimatched cases.'''
    dim = 0
    path = []
    for cur_key in keys_list:
        dim += 1
        if cur_key is MATCHED or cur_key is MULTIMATCHED:
            path.append(cur_key)
        elif cur_key in in_dict[ANY][str(dim)]:
            path.append(MATCHED)
        else:
            path.append(UNMATCHED)
    return path

def setNDApproxDictVal(in_dict, keys_list, mismatches):
    '''Creates (or builds upon) an N-Dimensional ("ND") dictionary (in_dict)
    whose keys are mismatch variants (constructed using the keys_list and the
    number of mismatches provided).  The values at each level (initially:
    changed by reduceNDDictPaths) is a 2 member list containing the number of
    mismatches the key has and the dictionary for the next level.  The values
    at the end of the ND dictionary are (initially: changed by
    reduceNDDictPaths) a 2 member list containing the number of mismatches and
    a list of "original key paths".  (Each path is a list of barcodes.)  The
    keys_list is a series of barcodes representing a row from the barcodes file
    (i.e. 1 barcode from each dimension, in order).
    Note: in_dict must be "l()" where "l" = lambda:defaultdict(l).'''
    setNDApproxDictValHelper(in_dict, keys_list, mismatches, [])

def setNDApproxDictValHelper(in_dict, keys_list, mismatches, real_path):
    '''A recursive helper of setNDApproxDictVal which tracks the "real barcode
    path"  to use in the list of real barcode paths at the end of the
    dictionary.  Note that barcodes may have a 1:many relationship, indicated
    by a barcode being repeated in earlier columns of the barcode file.
    Different paths may also intersect given the number of allowed mismatches.
    '''
    l = lambda:defaultdict(l)
    real_path.append(keys_list[0])
    if len(keys_list) == 1:
        #print "Doing last key %s" %(keys_list[0])
        for bcl in getMismatchedBarcodes(keys_list[0],mismatches):
            #print "Processing mismatched barcode record: %s" % (bcl)
            nmm = bcl[0]
            bc  = bcl[1]
            if (bc in in_dict and type(in_dict[bc][1]) is list
                and len(in_dict[bc][1]) > 0
                and not real_path in in_dict[bc][1]):

                #print ("Adding real path %s to the one existing at barcode "
                #       "%s: %s whose type should be a list: %s"
                #       ) % (real_path, bc, in_dict[bc], type(in_dict[bc]))

                #If the number of mismatches is the same for this barcode
                #variant as it is for the real paths that have already been
                #added to the real paths list, append this real path
                if nmm == in_dict[bc][0]:
                    in_dict[bc][1] += [real_path]
                #Else if the number of mismatches is less for this barcode
                #variant as it is for the real paths that have already been
                #added to the real paths list, overwrite the old list
                elif nmm < in_dict[bc][0]:
                    in_dict[bc] = [nmm, [real_path]]
            elif not (bc in in_dict and type(in_dict[bc][1]) is list
                      and len(in_dict[bc][1]) > 0):
                #print ("Setting new real path %s to real path list for "
                #       "barcode %s: %s whose type should be a list: %s"
                #       ) % (real_path, bc, in_dict[bc], type(in_dict[bc]))
                in_dict[bc] = [nmm, [real_path]]
            #Else: skip poorer match
    else:
        #print "Doing first key %s" %(keys_list[0])
        for bcl in getMismatchedBarcodes(keys_list[0],mismatches):
            nmm = bcl[0]
            bc  = bcl[1]
            if bc not in in_dict:
                in_dict[bc] = [nmm, l()]
                setNDApproxDictValHelper(in_dict[bc][1],
                                         keys_list[1:],
                                         mismatches,
                                         real_path[:])
            elif (bc in in_dict and nmm == in_dict[bc][0]):
                setNDApproxDictValHelper(in_dict[bc][1],
                                         keys_list[1:],
                                         mismatches,
                                         real_path[:])
            elif (bc in in_dict and nmm < in_dict[bc][0]):
                #print "Type of in_dict[%s] is: %s Resetting list." % (bc, type(in_dict[bc]))
                in_dict[bc] = [nmm, l()]
                #in_dict[bc][0] = nmm
                #in_dict[bc][1] = l()
                setNDApproxDictValHelper(in_dict[bc][1],
                                         keys_list[1:],
                                         mismatches,
                                         real_path[:])

def reduceNDDictPaths(in_dict):
    '''After an approximate barcode dictionary has been fully created by
    setNDApproxDictVal, this method reduces the lists of barcode paths to
    individual discrete paths or individual generic MATCHED/MULTIMATCHED paths.
    '''
    if isinstance(in_dict,dict) and len(in_dict) > 0:
        arb_bc = in_dict.iterkeys().next()
        #print "Checking arbitrary barcode to see what level of the dictionary we're at: %s" % (arb_bc)
        #print "The keys of this level of the dictionary are: %s" % (in_dict.keys())
    if (isinstance(in_dict,dict) and len(in_dict) > 0 and
        isinstance(in_dict[arb_bc][1],dict)):
        for bc in in_dict:
            #Reduction (skipping the list...)
            in_dict[bc] = in_dict[bc][1]
            reduceNDDictPaths(in_dict[bc])
    #Else if the value is a list of lists
    elif (isinstance(in_dict,dict) and len(in_dict) > 0 and
          isinstance(in_dict[in_dict.keys()[0]][1],list)):
        #For each last barcode in the dictionary
        for bc in in_dict:
            #Reduction, skipping the outer list
            in_dict[bc] = in_dict[bc][1]
            #If the value at the leaf of the dict is a list of lists
            if (len(in_dict[bc]) > 0 and isinstance(in_dict[bc],list) and
                len(in_dict[bc][0]) > 0 and isinstance(in_dict[bc][0],list)):
                #If the number of paths is 1
                if len(in_dict[bc]) == 1:
                    in_dict[bc] = in_dict[bc][0]
                else:
                    inner_len = len(in_dict[bc][0])
                    ref_path = in_dict[bc][0]
                    new_path = []
                    for bci in xrange(0,inner_len):
                        all_matched = True
                        for path in in_dict[bc][1:]:
                            if path[bci] != ref_path[bci]:
                                all_matched = False
                                break
                        if all_matched is True:
                            new_path.append(MATCHED)
                        else:
                            new_path.append(MULTIMATCHED)
                    in_dict[bc] = new_path

def getMismatchedBarcodes(barcode,mismatches):
    '''Given a barcode sequence and a number of allowed mismatches, generate a
    list of barcodes representing all possible matching sequences with the
    allowed number of mismatches.'''
    return(getMismatchedBarcodesHelper(barcode,mismatches,[]))

def getMismatchedBarcodesHelper(barcode,mismatches,skips):
    offbys = [[0, barcode]]
    if mismatches <= 0:
        return(offbys)
    for p in xrange(0,len(barcode)):
        if p in skips:
            continue
        #Keep track of bases that have already been mutated
        newskips = skips + [p]
        if barcode[p] != 'N':
            offbys.append([len(newskips),
                           barcode[:p] + 'N' + barcode[p+1:]])
        if barcode[p] != 'A':
            offbys.append([len(newskips),
                           barcode[:p] + 'A' + barcode[p+1:]])
        if barcode[p] != 'T':
            offbys.append([len(newskips),
                           barcode[:p] + 'T' + barcode[p+1:]])
        if barcode[p] != 'G':
            offbys.append([len(newskips),
                           barcode[:p] + 'G' + barcode[p+1:]])
        if barcode[p] != 'C':
            offbys.append([len(newskips),
                           barcode[:p] + 'C' + barcode[p+1:]])
        #Recurse to add next allowed mismatch
        if mismatches > 1:
            if barcode[p] != 'N':
                offbys += getMismatchedBarcodesHelper(barcode[:p] + 'N' +
                                                      barcode[p+1:],
                                                      mismatches - 1,
                                                      newskips)
            if barcode[p] != 'A':
                offbys += getMismatchedBarcodesHelper(barcode[:p] + 'A' +
                                                      barcode[p+1:],
                                                      mismatches - 1,
                                                      newskips)
            if barcode[p] != 'T':
                offbys += getMismatchedBarcodesHelper(barcode[:p] + 'T' +
                                                      barcode[p+1:],
                                                      mismatches - 1,
                                                      newskips)
            if barcode[p] != 'G':
                offbys += getMismatchedBarcodesHelper(barcode[:p] + 'G' +
                                                      barcode[p+1:],
                                                      mismatches - 1,
                                                      newskips)
            if barcode[p] != 'C':
                offbys += getMismatchedBarcodesHelper(barcode[:p] + 'C' +
                                                      barcode[p+1:],
                                                      mismatches - 1,
                                                      newskips)
    return(offbys)

def exeExists(program):
    '''Determines whether an executable exists (i.e. is in the user\'s path)'''
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True

    return False

def run_gzip_in(infile):
    '''Takes a gzipped file and returns an iterator on the unzipped lines of
    the file using a system call to gzip.'''
    p = subprocess.Popen(["gzip", "-dc", infile], stdout=subprocess.PIPE,
                         bufsize=1)
    return iter(p.stdout.readline, b'')

def run_gzip_out(outfile):
    p = subprocess.Popen(["gzip"], stdin=subprocess.PIPE,
                         stdout=open(outfile,'wb'), bufsize=-1)
    return p.stdin

if __name__ == '__main__':
    sys.exit(main())
