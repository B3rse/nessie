#!/usr/bin/env python


######################################################################################
#
##	to_wig.py
#		Script to format the output from NeSSie related to complexity and entropy 
#		into wig format file. 
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


import sys, argparse,re
 
 
def main(args): # use as args['name']

	# Variables:
	dict_hits = {}

	## Reading input file
	with open(args['inputfile'], 'r') as fi:
		first_seq = True
		
		for line in fi:
			if line.startswith('#'): #tipo = 0 linguistic, tipo = 1 entropy
				command_line = line.rstrip()
				if ('-E' in command_line) or ('--entropy' in command_line) or ('-L' in command_line) or ('--linguistic' in command_line):
					regex = re.search(r"l\s(\d+)", command_line)
					win_len = int(regex.group(1))

					regex_1 = re.search(r"s\s(\d+)", command_line)
					if regex_1:
						shift = int(regex_1.group(1))
					else:
						shift = 1
					#end if

					if ('-E' in command_line) or ('--entropy' in command_line):
						tipo = 'Shannon entropy'
						tipo_1 = 'Entropy scores'
					else:
						tipo = 'linguistic complexity'
						tipo_1 = 'Complexity scores'
				else:
					print 'Provide an output file generated from complexity or entropy analysis'
					raise SystemExit
				#end if
			elif line.startswith('>'): 
				seq_name = line.rstrip().split()[0][1:]
				dict_hits.setdefault(seq_name, [])
			elif line.startswith('@'):
				start_idx = int(line.rstrip().split('-')[0][1:])
			else:
				idx, score_raw = int(line.split()[0]), line.split()[1]
				if score_raw == '-0':
					score = 0.0
				else:
					score = float(score_raw)
				#dict_hits[seq_name].append((idx + start_idx, score))
				dict_hits[seq_name].append((idx + start_idx + (win_len / 2), score))
			#end if
		#end for
	#end with

	#Print output to file
	with open(args['outputfile'], 'w') as fo:
		for key, value in dict_hits.iteritems():
			fo.write('track type=wiggle_0 name="{0} {1}, windows {2} - shift {3}"'.format(key, tipo_1, win_len, shift)) 
			fo.write(' description="{0} scores" visibility=full color=50,150,255\n'.format(tipo))
			#fo.write(' visibility=full autoScale=off viewLimits=0.0:25.0 color=50,150,255 yLineMark=11.76 yLineOnOff=on priority=10')
			#fo.write('variableStep chrom={0} span={1}\n'.format(key, shift))
			fo.write('variableStep chrom={0}\n'.format(key))
			for index, score in value:
				fo.write('{0} {1}\n'.format(index, score))
			#end for
	#end with
# end def main
 
 
if __name__ == '__main__':
 
    parser = argparse.ArgumentParser(description='Description of your program')
 
    parser.add_argument('-i','--inputfile', help='input file', required=True)
    parser.add_argument('-o','--outputfile', help='output file', required=True)
 
    args = vars(parser.parse_args())
 
    main(args)
 
# end if
