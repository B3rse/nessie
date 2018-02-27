#!/usr/bin/env python


######################################################################################
#
##	NessieOutParser.py
#		Script to format the output from NeSSie related to Palindromes, Mirrors or
#		Triplexes into a tab formatted file.
#
#	Author: Michele Berselli
#		University of Padova
#		berselli.michele@gmail.com
#
##	LICENSE:
#		Copyright (C) 2017  Michele Berselli
#
#		This program is free software: you can redistribute it and/or modify
#		it under the terms of the GNU General Public License as published by
#		the Free Software Foundation.
#
#		This program is distributed in the hope that it will be useful,
#		but WITHOUT ANY WARRANTY; without even the implied warranty of
#		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#		GNU General Public License for more details.
#
#		You should have received a copy of the GNU General Public License
#		along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
######################################################################################


# Import libraries
import sys, argparse, re


# Functions definition
def routine_score(seq_len, align):

	# Variables
	score_increment, score_idls, score_mm = -2, -1, -1
	score = seq_len

	if align:
		aln, last_aln = '', ''
		for i in range(0, len(align), 2):
			aln += align[i]
			aln += align[i + 1]

			if aln == '00':
				if last_aln == '01':
					score += score_increment
					#score_increment -= 1
				else:
					score += score_idls
				#end if
			elif aln == '11':
				if last_aln == '01':
					score += score_increment
					#score_increment -= 1
				else:
					score += score_idls
				#end if
			elif aln == '10':
				if last_aln == '01':
					score += score_increment
					#score_increment -= 1
				else:
					score += score_mm
				#end if
			#end if

			last_aln = aln
			aln = ''
		#end for
	#end if

	return score
#end def routine_score

def print_tab(dict_to_print, of, strand, seq_type, flag, counts_only, indexes_only):

	if flag == 0: #0=indexes, 1=counts, 2=score
		for (key, value) in sorted(dict_to_print.iteritems(), key=lambda(x, y): y['indexes'][0]):
			routine_print_tab(key, value, of, strand, seq_type, counts_only, indexes_only)
		#end for
	elif flag == 1:
		for (key, value) in sorted(dict_to_print.iteritems(), key=lambda(x, y): y['counts'], reverse=True):
			routine_print_tab(key, value, of, strand, seq_type, counts_only, indexes_only)
		#end for
	else:
		for (key, value) in sorted(dict_to_print.iteritems(), key=lambda(x, y): y['score'], reverse=True):
			routine_print_tab(key, value, of, strand, seq_type, counts_only, indexes_only)
		#end for
	#end if
#end def print_tab

def routine_print_tab(key, value, of, strand, seq_type, counts_only, indexes_only):

	#fasta_ID\tmotif\tmotif_type\tstrand\tmotif_length\tscore\tcounts\tindexes\n'

	of.write('{0}\t{1}\t'.format(value['seq_ID'], key))
	if seq_type == 0: #0=mirror, 1=palindrome, 2=triplex
		of.write('m\t{0}\t'.format(strand))
	elif seq_type == 1:
		of.write('p\t{0}\t'.format(strand))
	else:
		of.write('t\t{0}\t'.format(strand))
	#end if
	of.write('{0}\t{1}'.format(value['seq_len'], value['score']))
	if not indexes_only:
		of.write('\t{0}'.format(value['counts']))
	#end if
	if not counts_only:
		of.write('\t')
		of.write(','.join(str(idx) for idx in value['indexes']))
	#end if
	of.write('\n')
#end def routine_print_tab

def main(args):

	## Variables
	flag_score = True if args['score'] else False
	flag_counts = True if args['orderbycounts'] else False
	counts_only = False
	indexes_only = False
	seq_type = 0 #0=mirror, 1=palindrome, 2=triplex
	strand = '+'

	## Init data structures
	dict_hits = {} # {seq: {seq_ID: 'i', seq_len: s, score: n, counts: c, indexes: [i, ...]}, ...}

	## Opening output file
	of = open(args['outputfile'], 'w')

	## Reading input file
	with open(args['inputfile']) as fi:
		first_seq = True

		for line in fi:
			if line.startswith('#'):
				command_line = line.rstrip().split()
				if ('-P' in command_line) or ('--palindrome' in command_line):
					seq_type = 1
				elif ('-T' in command_line) or ('--triplex' in command_line):
					seq_type = 2
				#end if
				if ('-C' in command_line) or ('--complement' in command_line):
					strand = '-'
				#end if
				if ('-c' in command_line) or ('--counts' in command_line):
					flag_counts = True # having only counts hits will be ordered by counts 
					counts_only = True
				elif ('-i' in command_line) or ('--indexes' in command_line):
					if flag_counts:
						print '\nTo order by counts provide output with counts, results will be ordered by indexes!\n'
					#end if
					indexes_only = True
					flag_counts = False # having only indexes hits will be ordered by indexes 
				#end if

				## Decide ordering used to print
				if flag_score:
					flag = 2 #0=indexes, 1=counts, 2=score
				else:
					if flag_counts:
						flag = 1
					else:
						flag = 0
					#end if
				#end if
			elif line.startswith('>') and first_seq:
				seq_ID_space = re.sub(r'[^\w\s]', '', line.rstrip())
				seq_ID = re.sub('\s+', '_', seq_ID_space)
				first_seq = False

				## Writing header
				if counts_only:
					of.write('#fasta_ID\tmotif\tmotif_type\tstrand\tmotif_length\tscore\tcounts\n')
				elif indexes_only:
					of.write('#fasta_ID\tmotif\tmotif_type\tstrand\tmotif_length\tscore\tindexes\n')
				else:
					of.write('#fasta_ID\tmotif\tmotif_type\tstrand\tmotif_length\tscore\tcounts\tindexes\n')
				#end if
			elif line.startswith('>') and not first_seq:
				print_tab(dict_hits, of, strand, seq_type, flag, counts_only, indexes_only) #<----
				seq_ID_space = re.sub(r'[^\w\s]', '', line.rstrip())
				seq_ID = re.sub('\s+', '_', seq_ID_space)
				dict_hits = {}	# re-initializing dict_hits for new sequence
			elif line.startswith('$'):
				seq_len, seq = int(line.split('|')[1]), line.split('|')[2].rstrip()
				if len(line.split('|')) == 4:
					align = line.split('|')[3].rstrip()
				else:
					align = None
				#end if
				dict_hits.setdefault(seq, {'seq_ID': seq_ID, 'seq_len' : seq_len, 'score': routine_score(seq_len, align), 'counts': 0, 'indexes': []})
			elif line.startswith('@counts'):
				counts = int(line.rstrip().split()[1])
				dict_hits[seq]['counts'] += counts
			elif line.startswith('@indexes'):
				indexes = map(lambda x: int(x), line.rstrip().split()[1].split('|')[:-1])
				dict_hits[seq]['indexes'] += indexes
			#end if
		#end for
		print_tab(dict_hits, of, strand, seq_type, flag, counts_only, indexes_only) #<----
		#end if
	#end with

	## Closing output file
	of.close()

#end def main

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Script to format the output from NeSSIe')

	parser.add_argument('-i','--inputfile', help='output file from NeSSIe as input', required=True)
	parser.add_argument('-o','--outputfile', help='file to store formatted output', required=True)
	parser.add_argument('-c','--orderbycounts', help='order the results by counts and not by indexes', action='store_true', required=False)
	parser.add_argument('-s','--score', help='order by score', action='store_true', required=False)

	args = vars(parser.parse_args())

	main(args)
#end if

