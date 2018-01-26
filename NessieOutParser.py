#!/usr/bin/env python


######################################################################################
#
##	NessieOutParser.py
#		Script to format the output from NeSSie related to Palindromes, Mirrors or
#		Triplexes. 
#		The -s flag allows to assign a score to the retrieved motifs, that are ordered
#		from highest to lowest. 
#		Alternatively the results can be ordered by indexes (default), or by counts (-c). 
#		The script allows also to visualize the optimal alignment using the -a flag.
#		Finally -g will generate a gff file to easy the visualization of the results.
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
import sys, argparse


# Functions definition
def routine_score(seq, seq_len, align):

	# Variables
	score_increment, score_idls, score_mm = -2, -1, -1
	score = seq_len

	if align:
		align_string = ''
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
				align_string += 'u'
			elif aln == '11':
				if last_aln == '01':
					score += score_increment
					#score_increment -= 1
				else:
					score += score_idls
				#end if
				align_string += 'l'
			elif aln == '01':
				align_string += 'M'
			else:
				if last_aln == '01':
					score += score_increment
					#score_increment -= 1
				else:
					score += score_mm
				#end if
				align_string += 'm'
			#end if

			last_aln = aln
			aln = ''
		#end for
	else:
		align_string = ''
	#end if

	return score, align_string
#end def routine_score

def routine_print_score_align(of, seq, seq_len, align, palindrome, print_aln):

	score, align_string = routine_score(seq, seq_len, align)

	# Writing score
	of.write('{0}\n'.format(score))

	if print_aln:
		routine_print_align(of, seq, seq_len, align, palindrome, align_string)
	#end if
#end def routine_print_score_align

def routine_print_align(of, seq, seq_len, align, palindrome, align_string):

	# Variables
	equal_len = False
	dict_conversion = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}

	if align_string:
		# Write first sequence
		of.write('@aln: ')
		i = 0
		for l in align_string:
			if l == 'M' or l == 'm' or l == 'u':
				of.write('{0}'.format(seq[i]))
				i += 1
			else:
				of.write('-')
			#end if
		#end for
		of.write('\n')

		# Write align
		of.write('@aln: ')
		for l in align_string:
			if l == 'M':
				of.write('|')
			else:
				of.write(' ')
			#end if
		#end for
		of.write('\n')

		# Write second sequence
		of.write('@aln: ')
		i = 0
		for l in align_string:
			if l == 'M' or l == 'm' or l == 'l':
				of.write('{0}'.format(seq[seq_len - i - 1]))
				i += 1
			else:
				of.write('-')
			#end if
		#end for
		of.write('\n')

	else:
		first_half = seq_len / 2
		second_half = seq_len - first_half
		
		if first_half == second_half:
			equal_len = True
		#end if
		
		# Write first sequence
		of.write('@aln: ')
		for i in range(first_half):
			of.write('{0}'.format(seq[i]))
		#end for
		if not equal_len:
			of.write('-')
		#end if
		of.write('\n')
		
		# Write align
		of.write('@aln: ')
		if not palindrome:
			for i in range(first_half):
				if seq[i] == seq[seq_len - i - 1]:
					of.write('|')
				else:
					of.write(' ')
				#end if
			#end for
		else:
			for i in range(first_half):
				if seq[i] == dict_conversion[seq[seq_len - i - 1]]:
					of.write('|')
				else:
					of.write(' ')
				#end if
			#end for
		#end if

		of.write('\n')

		# Write align
		of.write('@aln: ')
		for i in range(second_half):
			of.write('{0}'.format(seq[seq_len - i - 1]))
		#end for
		of.write('\n')
	#end if
#end def routine_print_align

def routine_print_ordered(key, value, of, counts_only, indexes_only, show_align, palindrome):

	of.write('$|{0}|{1}|'.format(value['seq_len'], key))

	if show_align:
		routine_print_score_align(of, key, value['seq_len'], value['align'], palindrome, True)
	else:
		routine_print_score_align(of, key, value['seq_len'], value['align'], palindrome, False)
	#end if

	if counts_only:
		of.write('@counts: {0}\n'.format(value['counts']))
	elif indexes_only:
		of.write('@indexes: ')
		[of.write('{0}|'.format(i)) for i in value['indexes']]
		of.write('\n')
	else:
		of.write('@counts: {0}\n'.format(value['counts'])) 
		of.write('@indexes: ')
		[of.write('{0}|'.format(i)) for i in value['indexes']]
		of.write('\n')
	#end if
#end def routine_print_ordered


#def routine_print_gff(dict_to_print, of, counts_only, indexes_only, show_align, palindrome, seq_ID, seq_type, strand):

	## Variables
	#if seq_type == 0:
		#tipo = 'mirror'
	#elif seq_type == 1:
		#tipo = 'palindrome'
	#else:
		#tipo = 'triplex'
	##end if

	#dict_idxs_key = {}
	#for (key, value) in dict_to_print.iteritems():
		#for idx in value['indexes']:
			#dict_idxs_key.setdefault(idx, key)
		##end for
	##end for

	#for (idx, key) in sorted(dict_idxs_key.iteritems()):
		#score, _ = routine_score(key, dict_to_print[key]['seq_len'], dict_to_print[key]['align'])
		#of.write('{0}\tNeSSie\t{1}\t{2}\t{3}\t{4}\t{5}\t.\tID={6}\n'.format(
																			#seq_ID, 
																			#tipo,
																			#idx,
																			#idx + dict_to_print[key]['seq_len'],
																			#score,
																			#strand,
																			#key
																				#))
	##end for
##end def routine_print_gff

def routine_print_gff(dict_to_print, of, counts_only, indexes_only, show_align, palindrome, seq_ID, seq_type, strand):

	# Variables
	if seq_type == 0:
		tipo = 'mirror'
	elif seq_type == 1:
		tipo = 'palindrome'
	else:
		tipo = 'triplex'
	#end if

	dict_idxs_key = {}
	max_len = 1
	for (key, value) in dict_to_print.iteritems():
		if value['seq_len'] > max_len:
			max_len = value['seq_len']
		#end if 
		for idx in value['indexes']:
			dict_idxs_key.setdefault(idx, key)
		#end for
	#end for

	# Color scale
	color_dict = {0: 'FF0000', 1: 'FF0000', 2: 'FFEE50', 3: '0080FF', 4: '00FF80', 5: '00FF80', 6: '00FF80'}
	incr = max_len / 5

	for (idx, key) in sorted(dict_idxs_key.iteritems()):
		score, _ = routine_score(key, dict_to_print[key]['seq_len'], dict_to_print[key]['align'])
		of.write('{0}\tNeSSie\t{1}\t{2}\t{3}\t{4}\t{5}\t.\tID={6};color=#{7};score={8}\n'.format(
																			seq_ID, 
																			tipo,
																			idx,
																			idx + dict_to_print[key]['seq_len'],
																			score,
																			strand,
																			key,
																			color_dict[score / incr],
																			score
																				))
	#end for
#end def routine_print_gff

def print_ordered(dict_to_print, flag_counts, of, counts_only, indexes_only, show_align, palindrome, flag_gff, seq_ID, seq_type, strand):

	if not flag_gff:
		if flag_counts:
			for (key, value) in sorted(dict_to_print.iteritems(), key=lambda(x, y): y['counts'], reverse=True):
				routine_print_ordered(key, value, of, counts_only, indexes_only, show_align, palindrome)
			#end for
		else:
			for (key, value) in sorted(dict_to_print.iteritems(), key=lambda(x, y): y['indexes'][0]):
				routine_print_ordered(key, value, of, counts_only, indexes_only, show_align, palindrome)
			#end for
		#end if
	else:
		routine_print_gff(dict_to_print, of, counts_only, indexes_only, show_align, palindrome, seq_ID, seq_type, strand)
	#end if
#end def print_ordered

def main(args):

	## Variables
	flag_counts = True if args['orderbycounts'] else False
	show_align = True if args['alignshow'] else False
	flag_gff = True if args['gff'] else False
	counts_only = False
	indexes_only = False
	palindrome = False
	seq_type = 0 #0=mirror, 1=palindrome, 2=triplex
	strand = '+'

	## Init data structures
	dict_hits = {} # {seq: {seq_len: s, align: 'a', counts: c, intervals: [(i, i + s), ...], indexes: [i, ...]}, ...}
	dict_sorted = {}

	## Opening output file
	of = open(args['outputfile'], 'w')

	## Reading input file
	with open(args['inputfile']) as fi:
		first_seq = True

		for line in fi:
			if line.startswith('#'):
				command_line = line.rstrip().split()
				if ('-P' in command_line) or ('--palindrome' in command_line):
					palindrome = True
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
			elif line.startswith('>') and first_seq:
				seq_ID = line.rstrip().split()[0][1:]
				if not flag_gff:
					of.write('{0}\n'.format(line.rstrip()))
				#end if
				first_seq = False
			elif line.startswith('>') and not first_seq:
				print_ordered(dict_hits, flag_counts, of, counts_only, indexes_only, show_align, palindrome, flag_gff, seq_ID, seq_type, strand)
				seq_ID = line.rstrip().split()[0][1:]
				if not flag_gff:
					of.write('{0}\n'.format(line.rstrip()))
				#end if
				dict_hits = {}	# re-initializing dict_hits for new sequence
			elif line.startswith('$'):
				seq_len, seq = int(line.split('|')[1]), line.split('|')[2].rstrip()
				if len(line.split('|')) == 4:
					align = line.split('|')[3].rstrip()
				else:
					align = None
				#end if
				dict_hits.setdefault(seq, {'seq_len' : seq_len, 'align': align, 'counts': 0, 'indexes': []})
			elif line.startswith('@counts'):
				counts = int(line.rstrip().split()[1])
				dict_hits[seq]['counts'] += counts
			elif line.startswith('@indexes'):
				indexes = map(lambda x: int(x), line.rstrip().split()[1].split('|')[:-1])
				dict_hits[seq]['indexes'] += indexes
			#end if
		#end for
		print_ordered(dict_hits, flag_counts, of, counts_only, indexes_only, show_align, palindrome, flag_gff, seq_ID, seq_type, strand)
		#end if
	#end with

	## Closing output file
	of.close()

	## Order by score
	dict_to_order_score = {}

	#if args['score']:
		#if flag_gff:
			#print '\nNot possible to order a gff file by score, a normal gff have been produced!\n'
		#else:
			#old_of = open(args['outputfile'], 'r')
			#for line in old_of:
				#if line.startswith('#'):
					#command_line = line.rstrip()
				#elif line.startswith('>'):
					#seq_code = line.rstrip()[1:]
					#dict_to_order_score.setdefault(seq_code, {})
				#elif line.startswith('$'):
					#seq_len, seq, score = int(line.split('|')[1]), line.split('|')[2].rstrip(), int(line.split('|')[3].rstrip())
					#dict_to_order_score[seq_code].setdefault(seq, {'seq_len' : seq_len, 'align': align, 'counts': 0, 'indexes': [], 'score': score})
				#elif line.startswith('@counts'):
					#counts = int(line.rstrip().split()[1])
					#dict_to_order_score[seq_code][seq]['counts'] += counts
				#elif line.startswith('@indexes'):
					#indexes = map(lambda x: int(x), line.rstrip().split()[1].split('|')[:-1])
					#dict_to_order_score[seq_code][seq]['indexes'] += indexes
				##end if
			##end for
			#old_of.close()
			
			## Write
			#of = open(args['outputfile'], 'w')
			#for key, dict_key in dict_to_order_score.iteritems():
				#of.write('>{0}\n'.format(key))
				#for key_i, value_i in sorted(dict_key.iteritems(), key=lambda(x, y): y['score'], reverse = True):
					#of.write('$|{0}|{1}|{2}\n'.format(value_i['seq_len'], key_i, value_i['score']))
					#if value_i['counts']:
						#of.write('@counts: {0}\n'.format(value_i['counts']))
					##end if
					#if value_i['indexes']:
						#of.write('@indexes: ')
						#[of.write('{0}|'.format(index)) for index in value_i['indexes']]
						#of.write('\n')
					##end if
				##end for
			##end for
		##end if
	##end if

	if args['score']:
		if flag_gff:
			print '\nNot possible to order a gff file by score, a normal gff have been produced!\n'
		else:
			old_of = open(args['outputfile'], 'r')
			idx_seq = 0
			c = 1
			for line in old_of:
				if line.startswith('#'):
					command_line = line.rstrip()
				elif line.startswith('>'):
					idx_seq += 1
					seq_code = line.rstrip()[1:]
					dict_to_order_score.setdefault((seq_code, idx_seq), {})
				elif line.startswith('$'):
					seq_len, seq, score = int(line.split('|')[1]), line.split('|')[2].rstrip(), int(line.split('|')[3].rstrip())
					dict_to_order_score[(seq_code, idx_seq)].setdefault(seq, {'seq_len' : seq_len, 'align_1': '', 'align_2': '', 'align_3': '', 'counts': 0, 'indexes': [], 'score': score})
				elif line.startswith('@counts'):
					counts = int(line.rstrip().split()[1])
					dict_to_order_score[(seq_code, idx_seq)][seq]['counts'] += counts
				elif line.startswith('@indexes'):
					indexes = map(lambda x: int(x), line.rstrip().split()[1].split('|')[:-1])
					dict_to_order_score[(seq_code, idx_seq)][seq]['indexes'] += indexes
				elif line.startswith('@aln') and c == 1:
					c += 1
					dict_to_order_score[(seq_code, idx_seq)][seq]['align_1'] = line.rstrip()
				elif line.startswith('@aln') and c == 2:
					c += 1
					dict_to_order_score[(seq_code, idx_seq)][seq]['align_2'] = line.rstrip()
				elif line.startswith('@aln') and c == 3:
					c = 1
					dict_to_order_score[(seq_code, idx_seq)][seq]['align_3'] = line.rstrip()
				#end if
			#end for
			old_of.close()
			
			# Write
			of = open(args['outputfile'], 'w')
			for (key, idx_key), dict_key in sorted(dict_to_order_score.iteritems(), key=lambda(x,y): x[1]):
				of.write('>{0}\n'.format(key))
				for key_i, value_i in sorted(dict_key.iteritems(), key=lambda(x, y): y['score'], reverse = True):
					of.write('$|{0}|{1}|{2}\n'.format(value_i['seq_len'], key_i, value_i['score']))
					if value_i['align_1']:
						of.write('{0}\n'.format(value_i['align_1']))
						of.write('{0}\n'.format(value_i['align_2']))
						of.write('{0}\n'.format(value_i['align_3']))
					if value_i['counts']:
						of.write('@counts: {0}\n'.format(value_i['counts']))
					#end if
					if value_i['indexes']:
						of.write('@indexes: ')
						[of.write('{0}|'.format(index)) for index in value_i['indexes']]
						of.write('\n')
					#end if
				#end for
			#end for
		#end if
	#end if

#end def main

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Script to format the output from NeSSIe')

	parser.add_argument('-i','--inputfile', help='output file from NeSSIe as input', required=True)
	parser.add_argument('-o','--outputfile', help='file to store formatted output', required=True)
	parser.add_argument('-c','--orderbycounts', help='order the results by counts and not by indexes', action='store_true', required=False)
	parser.add_argument('-a','--alignshow', help='shows optimal alignment', action='store_true', required=False)
	parser.add_argument('-g','--gff', help='generate GFF file', action='store_true', required=False)
	parser.add_argument('-s','--score', help='order by score', action='store_true', required=False)

	args = vars(parser.parse_args())

	main(args)
#end if

