/**************************************************************************************
*
**	MAIN (main.cpp)
*
*	Author: Michele Berselli
*		University of Padova
*		berselli.michele@gmail.com
*
**	LICENSE:
*   	Copyright (C) 2017  Michele Berselli
*
*   	This program is free software: you can redistribute it and/or modify
*   	it under the terms of the GNU General Public License as published by
*   	the Free Software Foundation.
*
*  	 	This program is distributed in the hope that it will be useful,
*   	but WITHOUT ANY WARRANTY; without even the implied warranty of
*   	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   	GNU General Public License for more details.
*
*   	You should have received a copy of the GNU General Public License
*   	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
**************************************************************************************/


// INCLUDE
#include <iostream>
#include <stdint.h>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include "Nessie.h"
#include "FastaUtilities.h"

#define TEST
#undef TEST


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//													COMMAND LINE ARGUMENTS FUNCTIONS										//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
//		print_h
/////////////////////////////////////////////////////////////////////////////////////
void print_h(ostream &pout = std::cout){

	pout << "-------- Documentation -------- " << std::endl;
	pout << std::endl;
	pout << "Basic arguments (one is required)" << std::endl;
	pout << "  -h/--help:  help" << std::endl;
	pout << "  -P/--palindrome:  search for palindromic sequences" << std::endl;
	pout << "  -M/--mirror:  search for mirror sequences" << std::endl;
	pout << "  -N/--motif FILEPATH:  search for motif/motifs, FILEPATH to the file containing motifs to be searched (fasta/multifasta)" << std::endl;
	pout << "  -E/--entropy:  Shannon entropy calculation" << std::endl;
	pout << "  -L/--linguistic:  linguistic complexity calculation" << std::endl;
	pout << "  -T/--triplex:  search for possible triplex forming sequences" << std::endl;
	pout << "  -G/--quadruplex:  search for possible G4-DNA forming sequences" << std::endl;
	pout << "  -A/--allkmer:  search for all kmers in the sequence" << std::endl;

	pout << std::endl;
	pout << "I/O arguments (both are required)" << std::endl;
	pout << "  -I/--input FILEPATH:  input file (fasta/multifasta file)" << std::endl;
	pout << "  -O/--output FILEPATH:  output file" << std::endl;

	pout << std::endl;
	pout << "Additional arguments for all searches" << std::endl;
	pout << "  -b/--begin N:  begin of the interval to search" << std::endl;
	pout << "  -e/--end N:  end of the interval to search" << std::endl;
	pout << "  -c/--counts:  print counts only" << std::endl;
	pout << "  -i/--indexes:  print indexes only" << std::endl;
	pout << "  -C/--complement:  search in the reverse complement of the sequence" << std::endl;

	pout << std::endl;
	pout << "Additional arguments for -P/-M/-A/-L/-T" << std::endl;
	pout << "  -k/--kmin N:  minimum kmer length (required for -P/-M/-A/-T)" << std::endl;
	pout << "  -K/--kmax N:  maximum kmer length (required for -MAX)" << std::endl;

	pout << std::endl;
	pout << "Additional arguments for -P/-M/-T" << std::endl;
	pout << "  -m/--mismatch N:  percentage of permitted mismatches" << std::endl;
	pout << "  -g/--gap N:  percentage of permitted gaps" << std::endl;
	pout << "  -t/--total N:  percentage of permitted gaps and mismatches in total" << std::endl;
	pout << "  -MAX/--MAX:  flag that sets the search for the longest kmers with the selected symmetry in range [kmin..kmax]" << std::endl;

	pout << std::endl;
	pout << "Additional arguments for -E/-L" << std::endl;
	pout << "  -l/--interval N:  interval length" << std::endl;
	pout << "  -s/--shift N:  shift step for the interval" << std::endl;

	pout << std::endl;
	pout << "Additional arguments for -T" << std::endl;
	pout << "  -p/--purine N:  percentage of permitted non purine bases" << std::endl;
	pout << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////
//		print_basic
/////////////////////////////////////////////////////////////////////////////////////
void print_basic(ostream &pout = std::cout){

	pout << std::endl;
	pout << "Usage:" << std::endl;
	pout << "  nessie -h" << std::endl;
	pout << "  nessie -I inputFile -O outputFile -N motifsFile [ADDITIONAL ARGUMENTS]" << std::endl;
	//         0       1  2           3  4            5
	pout << "  nessie -I inputFile -O outputFile {-E | -L | -G} [ADDITIONAL ARGUMENTS]" << std::endl;
	//         0       1  2           3  4            5
	pout << "  nessie -I inputFile -O outputFile {-P | -M | -A | -T} -k N [ADDITIONAL ARGUMENTS]" << std::endl;
	//         0       1  2           3  4            5       6
	pout << std::endl;
	pout << "The shown arguments are required in the correct order, the additional arguments can be specified in any order instead!" << std::endl;
	std::cout << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////
//		startswith
/////////////////////////////////////////////////////////////////////////////////////
bool startswith(const char *a, const char *b){
   if(strncmp(a, b, strlen(b)) == 0) return 1;
   return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
//		init_basic_arg_I_O
/////////////////////////////////////////////////////////////////////////////////////
void init_basic_arg_I_O(char *argv[], std::string &inFile, std::string &outFile){

	if (("-I" == (std::string) argv[1] || "--input" == (std::string) argv[1]) && !startswith(argv[2], "-")){
		inFile = (std::string) argv[2];

		if (("-O" == (std::string) argv[3] || "--output" == (std::string) argv[3]) && !startswith(argv[4], "-")){
			outFile = (std::string) argv[4];
		}
		else{
			throw std::invalid_argument("[-O/--output FILEPATH] missing, call [-h] for documentation");
		}
	}
	else{
		throw std::invalid_argument("[-I/--input FILEPATH] missing, call [-h] for documentation");
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//		parsing_additional_arg_p_m_a_t
/////////////////////////////////////////////////////////////////////////////////////
void parsing_additional_arg_p_m_a_t(int &i, char *argv[], size_t &begin, size_t &end, bool &indexes, bool &counts,
		size_t &kmax, size_t &perc, size_t &perc_gap, bool &complement, bool &MAX, size_t &perc_purine, size_t &perc_gapmm){

	if (("-b" == (std::string) argv[i] || "--begin" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		begin = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}
	else if (("-e" == (std::string) argv[i] || "--end" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		end = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}
	else if ("-c" == (std::string) argv[i] || "--counts" == (std::string) argv[i]){
		indexes = false;
		i += 1;
	}
	else if ("-i" == (std::string) argv[i] || "--indexes" == (std::string) argv[i]){
		counts = false;
		i += 1;
	}
	else if (("-K" == (std::string) argv[i] || "--kmax" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		kmax = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}
	else if (("-m" == (std::string) argv[i] || "--mismatch" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		perc = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}
	else if ("-C" == (std::string) argv[i] || "--complement" == (std::string) argv[i]){
		complement = true;
		i += 1;
	}
	else if (("-g" == (std::string) argv[i] || "--gap" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		perc_gap = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}
	else if ("-MAX" == (std::string) argv[i] || "--MAX" == (std::string) argv[i]){
		MAX = true;
		i += 1;
	}
	else if (("-p" == (std::string) argv[i] || "--purine" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		perc_purine = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}
	else if (("-t" == (std::string) argv[i] || "--total" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		perc_gapmm = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}
	else{
		throw std::invalid_argument("non-recognized additional argument, call [-h] for documentation");
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//		parsing_additional_arg_e_l
/////////////////////////////////////////////////////////////////////////////////////
void parsing_additional_arg_e_l(int &i, char *argv[], size_t &begin, size_t &end, bool &indexes, bool &counts,
		size_t &kmin, size_t &kmax, size_t &interval, size_t &shift, bool &complement){

	if (("-b" == (std::string) argv[i] || "--begin" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		begin = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}
	else if (("-e" == (std::string) argv[i] || "--end" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		end = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}
	else if ("-c" == (std::string) argv[i] || "--counts" == (std::string) argv[i]){
		indexes = false;
		i += 1;
	}
	else if ("-i" == (std::string) argv[i] || "--indexes" == (std::string) argv[i]){
		counts = false;
		i += 1;
	}
	else if (("-k" == (std::string) argv[i] || "--kmin" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		kmin = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}
	else if (("-K" == (std::string) argv[i] || "--kmax" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		kmax = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}
	else if (("-l" == (std::string) argv[i] || "--interval" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		interval = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}

	else if (("-s" == (std::string) argv[i] || "--shift" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		shift = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}
	else if ("-C" == (std::string) argv[i] || "--complement" == (std::string) argv[i]){
		complement = true;
		i += 1;
	}
	else{
		throw std::invalid_argument("non-recognized additional argument, call [-h] for documentation");
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//		parsing_additional_arg_n
/////////////////////////////////////////////////////////////////////////////////////
void parsing_additional_arg_n(int &i, char *argv[], size_t &begin, size_t &end, bool &indexes, bool &counts, bool &complement){

	if (("-b" == (std::string) argv[i] || "--begin" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		begin = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}
	else if (("-e" == (std::string) argv[i] || "--end" == (std::string) argv[i]) && !startswith(argv[i + 1], "-")){
		end = strtoll(argv[i + 1], NULL, 10);
		i += 2;
	}
	else if ("-c" == (std::string) argv[i] || "--counts" == (std::string) argv[i]){
		indexes = false;
		i += 1;
	}
	else if ("-i" == (std::string) argv[i] || "--indexes" == (std::string) argv[i]){
		counts = false;
		i += 1;
	}
	else if ("-C" == (std::string) argv[i] || "--complement" == (std::string) argv[i]){
		complement = true;
		i += 1;
	}
	else{
		throw std::invalid_argument("non-recognized additional argument, call [-h] for documentation");
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//		calling_function
/////////////////////////////////////////////////////////////////////////////////////
void calling_function(std::ofstream &out, Fasta &fasta, int mode,
					  size_t begin, size_t end, bool counts, bool indexes,
					  size_t kmin, size_t kmax,
					  size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t modulo_purine, bool MAX,
					  size_t interval, size_t shift,
					  bool complement, std::vector<Fasta> &fasta_vector){

	// Variables
	const char *fasta_sequence_ptr = fasta.get_sequence().c_str();
	size_t fasta_sequence_len = fasta.get_sequence().length();
	if (!end){ end = fasta_sequence_len - 1; }

	// Check Variables
	if (end > fasta_sequence_len - 1){ throw std::invalid_argument("selected ending index [-e] is larger than sequence ending index"); }
	if (begin > end){ throw std::invalid_argument("selected begin index [-b] is larger than selected ending index [-e] or ending index has not been selected"); }
	if (interval > fasta_sequence_len){ throw std::invalid_argument("selected sliding interval [-l] is longer than sequence length"); }
	if (interval > (end - begin + 1)){ throw std::invalid_argument("selected sliding interval [-l] is longer than selected sequence interval"); }
	if (kmax > fasta_sequence_len){ throw std::invalid_argument("selected kmax [-K] is longer than sequence length"); }
	if (kmax > (end - begin + 1)){ throw std::invalid_argument("selected kmax [-K] is larger than selected selected sequence interval"); }
	if ((1 == mode || 2 == mode || 7 == mode || 8 == mode)){
		if (kmax && (kmax > fasta_sequence_len)){ throw std::invalid_argument("selected kmax [-K] is longer than sequence length"); }
		if (!kmax && (kmin > fasta_sequence_len)){ throw std::invalid_argument("selected kmin [-k] is longer than sequence length"); }
		if (kmax && (kmax > (end - begin + 1))){ throw std::invalid_argument("selected kmax [-K] is larger than selected selected sequence interval"); }
		if (!kmax && (kmin > (end - begin + 1))){ throw std::invalid_argument("selected kmin [-k] is larger than selected selected sequence interval"); }
	}
//	if (kmax && (kmin > kmax)){ throw std::invalid_argument("selected kmin [-k] is larger than selected kmax [-K]"); }
//	if ((interval && !shift) || (shift && !interval)){ throw std::invalid_argument("missing interval [-l] or shift [-s] argument"); }
//	if ((!kmax && interval) && (interval < 20)){ throw std::invalid_argument("default kmax [-K] is larger than selected sliding interval [-l], pls select a smaller kmax"); }
//	if (interval && (kmax > interval)){ throw std::invalid_argument("selected kmax [-K] is larger than selected sliding interval [-l]"); }

	// Other Variables
	std::vector<size_t> tmp_idx_unknown;
	std::vector<Fasta>::iterator IT;
	std::vector<size_t>::iterator it;
	size_t begin_i = begin;	//keep stored the current position on the sequence, starting from the first canonical base
	size_t end_i = end;	//stores the latest possible ending index (later indexes correspond to non-canonical bases)

	// Calculating begin_i, end_i and building tmp_idx_unknown for the interval [begin_i..end_i]
	if (fasta.get_idx_unknown().size()){
		// Retrieve index first canonical base from begin
		it = fasta.get_idx_unknown().begin();
		while(*it == begin_i) {
			++begin_i;
			++it;
		} //std::cout << begin_i << std::endl;

		// Retrieve index last canonical base from end
		it = fasta.get_idx_unknown().end(); --it;
		while(*it == end_i) {
			--end_i;
			--it;
		} //std::cout << end_i << std::endl;

		// Check non-canonical bases indexes relative to selected search interval [begin_i..end_i]
		for (it = fasta.get_idx_unknown().begin(); it != fasta.get_idx_unknown().end();){
			if (*it > begin_i && *it < end_i){
				tmp_idx_unknown.push_back(*it); //std::cout << *it << "|";
			}
			++it;
		} //std::cout << std::endl;
	}

	// Printing fasta sequence id
	out << '>' << fasta.get_id() << std::endl;

	// Checking mode
	switch(mode){
		case 1:	//palindrome
		{
			if (!kmax){ kmax = kmin; }
			for (it = tmp_idx_unknown.begin(); it != tmp_idx_unknown.end();){	// working sub-intervals skipping non-canonical bases
				if ((*it - begin_i) >= kmax){ //std::cout << begin_i << " - " << (*it - 1) << std::endl;
					class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, (*it - 1), false, complement);
					if (!MAX){ fasta_sequence_nessie.print_kmers_palindrome_gap(kmin, kmax, modulo, modulo_gap, modulo_gapmm, 0, 0, out, counts, indexes, begin_i); }
					else if (MAX){ fasta_sequence_nessie.print_max_kmers_palindrome_gap(kmax, kmin, modulo, modulo_gap, modulo_gapmm, 0, 0, out, counts, indexes, begin_i); }
				}
				begin_i = *it + 1;
				++it;
			}
			// working last sub-interval skipping non-canonical bases
			if ((end_i - begin_i + 1) >= kmax){ //std::cout << begin_i << " - " << end_i << std::endl;
				class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, end_i, false, complement);
				if (!MAX){ fasta_sequence_nessie.print_kmers_palindrome_gap(kmin, kmax, modulo, modulo_gap, modulo_gapmm, 0, 0, out, counts, indexes, begin_i); }
				else if (MAX){ fasta_sequence_nessie.print_max_kmers_palindrome_gap(kmax, kmin, modulo, modulo_gap, modulo_gapmm, 0, 0, out, counts, indexes, begin_i); }
			}
			break;
		}
		case 2:	//mirror
		{
			if (!kmax){ kmax = kmin; }
			for (it = tmp_idx_unknown.begin(); it != tmp_idx_unknown.end();){	// working sub-intervals skipping non-canonical bases
				if ((*it - begin_i) >= kmax){ //std::cout << begin_i << " - " << (*it - 1) << std::endl;
					class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, (*it - 1), false, complement);
					if (!MAX){ fasta_sequence_nessie.print_kmers_mirror_gap(kmin, kmax, modulo, modulo_gap, modulo_gapmm, 0, 0, out, counts, indexes, begin_i); }
					else if (MAX){ fasta_sequence_nessie.print_max_kmers_mirror_gap(kmax, kmin, modulo, modulo_gap, modulo_gapmm, 0, 0, out, counts, indexes, begin_i); }
				}
				begin_i = *it + 1;
				++it;
			}
			// working last sub-interval skipping non-canonical bases
			if ((end_i - begin_i + 1) >= kmax){ //std::cout << begin_i << " - " << end_i << std::endl;
				class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, end_i, false, complement);
				if (!MAX){ fasta_sequence_nessie.print_kmers_mirror_gap(kmin, kmax, modulo, modulo_gap, modulo_gapmm, 0, 0, out, counts, indexes, begin_i); }
				else if (MAX){ fasta_sequence_nessie.print_max_kmers_mirror_gap(kmax, kmin, modulo, modulo_gap, modulo_gapmm, 0, 0, out, counts, indexes, begin_i); }
			}
			break;
		}
		case 3:	//motif
		{
			class Nessie fasta_sequence_nessie(1, fasta_sequence_ptr, fasta_sequence_len, begin, end, complement);
			for (IT = fasta_vector.begin(); IT != fasta_vector.end();){
				try {
					Kmer *tmp_motif = fasta_sequence_nessie.check_kmer_char(IT->get_sequence().c_str(), IT->get_sequence().length(), 0, 0, complement);
					out << '!' << IT->get_id() << std::endl;
					tmp_motif->print(out, counts, indexes);
					++IT;
				}
				catch (exception &e){
					++IT;
				}
			}
			break;
		}
		case 4:	//entropy
		{
			if (interval){
				for (it = tmp_idx_unknown.begin(); it != tmp_idx_unknown.end();){	// working sub-intervals skipping non-canonical bases
					if ((*it - begin_i) >= interval){
						class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, (*it - 1), false, complement);
						out << '@' << begin_i << '-' << (*it - 1) << std::endl;
						fasta_sequence_nessie.print_shannon_entropy_sliding(interval, shift, out, 0, 0);
					}
					begin_i = *it + 1;
					++it;
				}
				// working last sub-interval skipping non-canonical bases
				if ((end_i - begin_i + 1) >= interval){
					class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, end_i, false, complement);
					out << '@' << begin_i << '-' << end_i << std::endl;
					fasta_sequence_nessie.print_shannon_entropy_sliding(interval, shift, out, 0, 0);
				}
			}
			else{
				for (it = tmp_idx_unknown.begin(); it != tmp_idx_unknown.end();){	// working sub-intervals skipping non-canonical bases
					if ((*it - 1) > begin_i){
						class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, (*it - 1), false, complement);
						out << '@' << begin_i << '-' << (*it - 1) << ": ";
						fasta_sequence_nessie.print_shannon_entropy_interval(out, 0, 0);
					}
					begin_i = *it + 1;
					++it;
				}
				// working last sub-interval skipping non-canonical bases
				if (end_i > begin_i){
					class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, end_i, false, complement);
					out << '@' << begin_i << '-' << end_i << ": ";
					fasta_sequence_nessie.print_shannon_entropy_interval(out, 0, 0);
				}
			}
			break;
		}
		case 5:	//linguistic
		{
			if (interval){
				for (it = tmp_idx_unknown.begin(); it != tmp_idx_unknown.end();){	// working sub-intervals skipping non-canonical bases
					if ((*it - begin_i) >= interval){
						class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, (*it - 1), false, complement);
						out << '@' << begin_i << '-' << (*it - 1) << std::endl;
						fasta_sequence_nessie.print_linguistic_complexity_sliding(interval, shift, out, 0, 0, kmin, kmax);
					}
					begin_i = *it + 1;
					++it;
				}
				// working last sub-interval skipping non-canonical bases
				if ((end_i - begin_i + 1) >= interval){
					class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, end_i, false, complement);
					out << '@' << begin_i << '-' << end_i << std::endl;
					fasta_sequence_nessie.print_linguistic_complexity_sliding(interval, shift, out, 0, 0, kmin, kmax);
				}
			}
			else{
				if (!kmax){ kmax = (fasta_sequence_len < 20) ? fasta_sequence_len : 20; }
				for (it = tmp_idx_unknown.begin(); it != tmp_idx_unknown.end();){	// working sub-intervals skipping non-canonical bases
					if ((*it - begin_i) >= kmax ){
						class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, (*it - 1), false, complement);
						out << '@' << begin_i << '-' << (*it - 1) << ": ";
						fasta_sequence_nessie.print_linguistic_complexity_interval(out, 0, 0, kmin, kmax);
					}
					begin_i = *it + 1;
					++it;
				}
				// working last sub-interval skipping non-canonical bases
				if ((end_i - begin_i + 1) >= kmax){
					class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, end_i, false, complement);
					out << '@' << begin_i << '-' << end_i << ": ";
					fasta_sequence_nessie.print_linguistic_complexity_interval(out, 0, 0, kmin, kmax);
				}
			}
			break;
		}
		case 6:	//quadruplex
		{
			break;
		}
		case 7:	//allkmer
		{
			if (!kmax){ kmax = kmin; }
			for (it = tmp_idx_unknown.begin(); it != tmp_idx_unknown.end();){	// working sub-intervals skipping non-canonical bases
				if ((*it - begin_i) >= kmax){ //std::cout << begin_i << " - " << (*it - 1) << std::endl;
					class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, (*it - 1), false, complement);
					fasta_sequence_nessie.print_kmers(kmin, kmax, 0, 0, out, counts, indexes, begin_i);
				}
				begin_i = *it + 1;
				++it;
			}
			// working last sub-interval skipping non-canonical bases
			if ((end_i - begin_i + 1) >= kmax){ //std::cout << begin_i << " - " << end_i << std::endl;
				class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, end_i, false, complement);
				fasta_sequence_nessie.print_kmers(kmin, kmax, 0, 0, out, counts, indexes, begin_i);
			}
			break;
		}
		case 8:	//triplex
		{
			if (!kmax){ kmax = kmin; }
			for (it = tmp_idx_unknown.begin(); it != tmp_idx_unknown.end();){	// working sub-intervals skipping non-canonical bases
				if ((*it - begin_i) >= kmax){ //std::cout << begin_i << " - " << (*it - 1) << std::endl;
					class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, (*it - 1), false, complement);
					if (!MAX){ fasta_sequence_nessie.print_kmers_triplex_gap(kmin, kmax, modulo, modulo_gap, modulo_gapmm, modulo_purine, 0, 0, out, counts, indexes, begin_i); }
					else if (MAX){ fasta_sequence_nessie.print_max_kmers_triplex_gap(kmax, kmin, modulo, modulo_gap, modulo_gapmm, modulo_purine, 0, 0, out, counts, indexes, begin_i); }
				}
				begin_i = *it + 1;
				++it;
			}
			// working last sub-interval skipping non-canonical bases
			if ((end_i - begin_i + 1) >= kmax){ //std::cout << begin_i << " - " << end_i << std::endl;
				class Nessie fasta_sequence_nessie(fasta_sequence_ptr, fasta_sequence_len, begin_i, end_i, false, complement);
				if (!MAX){ fasta_sequence_nessie.print_kmers_triplex_gap(kmin, kmax, modulo, modulo_gap, modulo_gapmm, modulo_purine, 0, 0, out, counts, indexes, begin_i); }
				else if (MAX){ fasta_sequence_nessie.print_max_kmers_triplex_gap(kmax, kmin, modulo, modulo_gap, modulo_gapmm, modulo_purine, 0, 0, out, counts, indexes, begin_i); }
			}
			break;
		}
		default:
		{
			throw std::runtime_error("invalid search mode");
			break;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//													MAIN																	//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){

	// Variables
	std::string inFile_path;
	std::string motifsFile_path;
	std::string outFile_path;
	std::ifstream inFile;
	std::ifstream motifsFile;
	std::ofstream outFile;
	std::ofstream logFile;

	int mode = 0; //palindrome = 1, mirror = 2, motif = 3, entropy = 4, linguistic = 5, quadruplex = 6, allkmer = 7, triplex = 8;

	// Varibales for additional arguments
	size_t begin = 0, end = 0;	//additional arguments for all searches
	bool counts = true, indexes = true; bool MAX = false;
	size_t kmin = 0, kmax = 0;	//additional arguments for -P/-M/-L
	size_t perc = 0, perc_gap = 0, perc_gapmm = 0; //additional arguments for -P/-M/-T
	size_t modulo = 0, modulo_gap = 0, modulo_gapmm = 0;
	size_t interval = 0, shift = 0;	//additional arguments for -E/-L
	bool complement = false;	//additional arguments for -N
	size_t perc_purine = 0;
	size_t modulo_purine = 0;	//additional arguments for -T
	MultiFasta motifs;	//class to store the motifs to be searched with -N flag

	// Credits
	std::cout << std::endl;
	std::cout << "NeSSIe - "
			  << "patterns search and sequence analysis on DNA strings using the Nessie library." << std::endl << std::endl
			  << "NeSSIe Copyright (C) 2017  Michele Berselli." << std::endl
			  << "This program comes with ABSOLUTELY NO WARRANTY." << std::endl
			  << "This is free software, and you are welcome to redistribute it under certain conditions;" << std::endl
			  << "see the GNU General Public License for more details." << std::endl;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//													Parsing command line													//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Parsing
	if (1 == argc){
		std::cout << std::endl;
		std::cout << "Check the README or call [-h] for documentation" << std::endl;
		std::cout << std::endl;
		return 0;
	}
	else if (argc < 6){	//error in command line or just -h called
		if ("-h" == (std::string) argv[1] || "--help" == (std::string) argv[1]){
			print_basic();
			print_h();
			return 0;
		}
		else{
			std::cerr << std::endl;
			std::cerr << "at least 3 arguments are required for basic analysis with [-E/-L/-G] flag, call [-h] for documentation" << std::endl;
			print_basic(std::cerr);
			return 1;
		}
	}
	else if (6 == argc){	//minimal command line for -N/-E/-L/-G or wrong command line for -P/-M
		if ("-P" == (std::string) argv[5] || "--palindrome" == (std::string) argv[5] || //-P/-M require the additional parameters -k/--kmin
			"-M" == (std::string) argv[5] || "--mirror" == (std::string) argv[5] ||
			"-A" == (std::string) argv[5] || "--allkmer" == (std::string) argv[5] ||
			"-N" == (std::string) argv[5] || "--motif" == (std::string) argv[5] ||
			"-T" == (std::string) argv[5] || "--triplex" == (std::string) argv[5]){
				std::cerr << std::endl;
				std::cerr << "at least 4 arguments are required for basic analysis with [-P/-M/-A/-N/-T] flag, call [-h] for documentation" << std::endl;
				print_basic(std::cerr);
				return 1;
		}
		else{	//minimal command line for -N/-E/-L/-G
			// Parsing first 5 arguments (I/O)
			try{
				init_basic_arg_I_O(argv, inFile_path, outFile_path);
			}
			catch (exception &e){
				std::cerr << std::endl;
				std::cerr << e.what() << std::endl;
				print_basic(std::cerr);
				return 1;
			}

			// Parsing 6th argument (flag for the search mode)
			if ("-E" == (std::string) argv[5] || "--entropy" == (std::string) argv[5]){
				mode = 4;
			}
			else if ("-L" == (std::string) argv[5] || "--linguistic" == (std::string) argv[5]){
				mode = 5;
			}
			else if ("-G" == (std::string) argv[5] || "--quadruplex" == (std::string) argv[5]){
				mode = 6;
				std::cout << std::endl;
				std::cout << "This option is currently under development, stay tuned!" << std::endl;
				std::cout << std::endl;
				return 1;
			}
			else{
				std::cerr << std::endl;
				std::cerr << "non-recognized argument for the search mode [-E/-L/-G], call [-h] for documentation" << std::endl;
				print_basic(std::cerr);
				return 1;
			}
		}
	}
	else{	//minimal command line for -P/-M/-A/-N or command line with additional arguments
		// Parsing first 5 arguments (I/O)
		try{
			init_basic_arg_I_O(argv, inFile_path, outFile_path);
		}
		catch (exception &e){
			std::cerr << std::endl;
			std::cerr << e.what() << std::endl;
			print_basic(std::cerr);
			return 1;
		}

		// Parsing 6th argument (flag for the search mode)
		if ("-P" == (std::string) argv[5] || "--palindrome" == (std::string) argv[5] ||
			"-M" == (std::string) argv[5] || "--mirror" == (std::string) argv[5] ||
			"-A" == (std::string) argv[5] || "--allkmer" == (std::string) argv[5] ||
			"-T" == (std::string) argv[5] || "--triplex" == (std::string) argv[5]){

			if ("-P" == (std::string) argv[5] || "--palindrome" == (std::string) argv[5]){
				mode = 1;
			}
			else if ("-M" == (std::string) argv[5] || "--mirror" == (std::string) argv[5]){
				mode = 2;
			}
			else if ("-T" == (std::string) argv[5] || "--triplex" == (std::string) argv[5]){
				mode = 8;
			}
			else{	//allkmer
				mode = 7;
			}

			// Parsing 7th, kmin
			if (("-k" == (std::string) argv[6] || "--kmin" == (std::string) argv[6]) && !startswith(argv[7], "-")){
				kmin = strtoll(argv[7], NULL, 10);
			}
			else{
				std::cerr << std::endl;
				std::cerr << "missing [-k/--kmin N] argument, call [-h] for documentation" << std::endl;
				print_basic(std::cerr);
				return 1;
			}

			//Parsing remaining arguments
			int i = 8;
			while (i < argc){
				try{
					parsing_additional_arg_p_m_a_t(i, argv, begin, end, indexes, counts, kmax, perc, perc_gap, complement, MAX, perc_purine, perc_gapmm);
				}
				catch (exception &e){
					std::cerr << std::endl;
					std::cerr << "non-recognized argument for the search mode [-P/-M/-A/-T], call [-h] for documentation" << std::endl;
					print_basic(std::cerr);
					return 1;
				}
			}
		}
		else{
			if ("-N" == (std::string) argv[5] || "--motif" == (std::string) argv[5]){
				mode = 3;

				// Parsing 7th, kmin
				if (!startswith(argv[6], "-")){
					motifsFile_path = (std::string) argv[6];
				}
				else{
					std::cerr << std::endl;
					std::cerr << "missing [-N FILEPATH] argument, call [-h] for documentation" << std::endl;
					print_basic(std::cerr);
					return 1;
				}

				int i = 7;
				while (i < argc){
					try{
						parsing_additional_arg_n(i, argv, begin, end, indexes, counts, complement);
					}
					catch (exception &e){
						std::cerr << std::endl;
						std::cerr << e.what() << std::endl;
						print_basic(std::cerr);
						return 1;
					}
				}
			}
			else if ("-E" == (std::string) argv[5] || "--entropy" == (std::string) argv[5]){
				mode = 4;

				int i = 6;
				while (i < argc){
					try{
						parsing_additional_arg_e_l(i, argv, begin, end, indexes, counts, kmin, kmax, interval, shift, complement);
					}
					catch (exception &e){
						std::cerr << std::endl;
						std::cerr << e.what() << std::endl;
						print_basic(std::cerr);
						return 1;
					}
				}
			}
			else if ("-L" == (std::string) argv[5] || "--linguistic" == (std::string) argv[5]){
				mode = 5;

				int i = 6;
				while (i < argc){
					try{
						parsing_additional_arg_e_l(i, argv, begin, end, indexes, counts, kmin, kmax, interval, shift, complement);
					}
					catch (exception &e){
						std::cerr << std::endl;
						std::cerr << e.what() << std::endl;
						print_basic(std::cerr);
						return 1;
					}
				}
			}
			else if ("-G" == (std::string) argv[5] || "--quadruplex" == (std::string) argv[5]){
				mode = 6;
				std::cout << std::endl;
				std::cout << "This option is currently under development, stay tuned!" << std::endl;
				std::cout << std::endl;
				return 1;
			}
			else{
				std::cerr << std::endl;
				std::cerr << "non-recognized argument for the search mode [-N/-E/-L/-G], call [-h] for documentation" << std::endl;
				print_basic(std::cerr);
				return 1;
			}
		}
	}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//													Check Variables															//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Check
	if (end && (begin > end)){
		std::cerr << std::endl;
		std::cerr << "selected begin index [-b] is larger than selected ending index [-e]" << std::endl;
		print_basic(std::cerr); return 1;
	}
	else if (end && (interval > (end - begin + 1))){
			std::cerr << std::endl;
			std::cerr << "selected sliding interval [-l] is longer than selected sequence interval" << std::endl;
			print_basic(std::cerr); return 1;
	}
	else if (end && (kmax > (end - begin + 1))){
			std::cerr << std::endl;
			std::cerr << "selected kmax [-K] is larger than selected selected sequence interval" << std::endl;
			print_basic(std::cerr); return 1;
	}

	if (kmax && (kmin > kmax)){
		std::cerr << std::endl;
		std::cerr << "selected kmin [-k] is larger than selected kmax [-K]" << std::endl;
		print_basic(std::cerr); return 1;
	}

	if ((interval && !shift) || (shift && !interval)){
		std::cerr << std::endl;
		std::cerr << "missing interval [-l] or shift [-s] argument" << std::endl;
		print_basic(std::cerr); return 1;
	}

	if (interval && (kmax > interval)){
		std::cerr << std::endl;
		std::cerr << "selected kmax [-K] is larger than selected sliding interval [-l]" << std::endl;
		print_basic(std::cerr); return 1;
	}

	if (MAX && !kmax){
		std::cerr << std::endl;
		std::cerr << "kmax [-K] is needed to search patterns with -MAX flag" << std::endl;
		print_basic(std::cerr); return 1;
	}

	if (perc_gapmm){
		if (perc > perc_gapmm){
			std::cerr << std::endl;
			std::cerr << "allowed mismatches percentage [-m] is higher than allowed mismatches and gaps total percentage [-t]" << std::endl;
			print_basic(std::cerr); return 1;
		}
		else if (perc_gap > perc_gapmm) {
			std::cerr << std::endl;
			std::cerr << "allowed gaps percentage [-g] is higher than allowed mismatches and gaps total percentage [-t]" << std::endl;
			print_basic(std::cerr); return 1;
		}
		else if ((perc + perc_gap) > perc_gapmm){
			std::cerr << std::endl;
			std::cerr << "allowed mismatches and gaps percentage [-m, -g] is higher than allowed mismatches and gaps total percentage [-t]" << std::endl;
			print_basic(std::cerr); return 1;
		}
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//													Reading input file and writing to output								//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Setting variables
	if (perc){ modulo = perc; }
	if (perc_gap){ modulo_gap = perc_gap; }
	if (perc_purine){ modulo_purine = perc_purine; }
	if (perc_gapmm){ modulo_gapmm = perc_gapmm; }

	// Open files
	inFile.open(inFile_path.c_str(), ios::in);
	outFile.open(outFile_path.c_str(), ios::out);
	logFile.open("logfile.txt", ios::out);
	if (3 == mode){
		motifsFile.open(motifsFile_path.c_str(), ios::in);
		motifs.get_data(motifsFile);
	}

	// Printing command line
	outFile << "#Command ";
	for (int i = 1; i < argc; ++i){
		outFile << argv[i] << " ";
	}
	outFile << std::endl;

	// Variables
	char c;	//current read char
	bool first = true;
	bool add_id = false;
	size_t curr_base_idx = 0;

	// Variables Fasta
	std::string tmp_id;
	std::string tmp_sequence;
	std::vector<size_t> tmp_idx_unknown;
	std::vector<Fasta>::iterator IT;

	// Reding inFile by char (fasta or multifasta)
	while (inFile.get(c)){
		if (('>' == c) && first){
			add_id = true;
			first = false;
		}
		else if (('>' == c) && !first){
			Fasta fasta_seq(tmp_id, tmp_sequence, tmp_idx_unknown);
			//TODO --> call operations here
			try{
				calling_function(outFile, fasta_seq, mode,
								  begin, end, counts, indexes,
								  kmin, kmax,
								  modulo, modulo_gap, modulo_gapmm, modulo_purine, MAX,
								  interval, shift,
								  complement, motifs.get_sequences_vector());
			}
			catch (exception &e){
				logFile << '>' << tmp_id << std::endl;
				logFile << e.what() << std::endl;
			}

			// Resetting tmp variables for the new fasta entry
			tmp_id.clear();
			tmp_sequence.clear();
			tmp_idx_unknown.clear();
			add_id = true;
			curr_base_idx = 0;
		}
		else if ('\n' == c){
			add_id = false;
		}
		else{
			if (add_id){
				tmp_id += c;
			}
			else{
				switch (c){
					case 'A':
					case 'a':
					case 'C':
					case 'c':
					case 'G':
					case 'g':
					case 'T':
					case 't':
						++curr_base_idx;
						break;
					default:
						tmp_idx_unknown.push_back(curr_base_idx);
						++curr_base_idx;
				}
				tmp_sequence += c;
			}
		}
	}
	// TODO --> call operations here
	Fasta fasta_seq(tmp_id, tmp_sequence, tmp_idx_unknown);
	try{
		calling_function(outFile, fasta_seq, mode,
						  begin, end, counts, indexes,
						  kmin, kmax,
						  modulo, modulo_gap, modulo_gapmm, modulo_purine, MAX,
						  interval, shift,
						  complement, motifs.get_sequences_vector());
	}
	catch (exception &e){
		logFile << '>' << tmp_id << std::endl;
		logFile << e.what() << std::endl;
	}

	std::cout << std::endl;
	// Closing files
	inFile.close();
	outFile.close();
	logFile.close();

	std::cout << "Analysis completed successfully! Check the log file for errors report."<< std::endl;


#ifdef TEST
	// CHECK command line
	std::cout << std::endl;
	std::cout << "---RECAP---" << std::endl;
	std::cout << "inFile_path: " << inFile_path << std::endl;
	std::cout << "outFile_path: " << outFile_path << std::endl;
	std::cout << "motifsFile_path: " << motifsFile_path << std::endl;
	std::cout << "mode: " << mode << std::endl;

	std::cout << std::endl;
	std::cout << "---RECAP_BIS---" << std::endl;
	std::cout << "begin: " << begin << std::endl;
	std::cout << "end: " << end << std::endl;
	std::cout << "kmin: " << kmin << std::endl;
	std::cout << "kmax: " << kmax << std::endl;
	std::cout << "counts: " << counts << std::endl;
	std::cout << "indexes: " << indexes << std::endl;
	std::cout << "modulo: " << modulo << std::endl;
	std::cout << "modulo_gap: " << modulo_gap << std::endl;
	std::cout << "modulo_purine: " << modulo_purine << std::endl;
	std::cout << "interval: " << interval << std::endl;
	std::cout << "shift: " << shift << std::endl;
	std::cout << "complement: " << complement << std::endl;
	std::cout << std::endl;

	std::string prova = "AGGAAAAGGAatgcatgcGGGCTCCTTCTTCTTCAAACTCTTTatgcatgGGGCTCCTTACTTCTTCCTCTTTatgcatgcAGGAAGAAGGA";//"CATAGCTTTATG";//"CATGATCTCAGATCTATATG"; //"AATCGCGTTATT"; //"ATTCGATAATA"; //"actgaaCCTGTCATGACAGGaacATCATGAT";
	Nessie nessie_p(prova.c_str(), prova.length());

#endif TEST

	return 0;
}

