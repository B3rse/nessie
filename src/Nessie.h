/**************************************************************************************
*
**	CLASS (Nessie.h) 'nucleic-acids elements of sequence symmetry identification'
*		Nessie is a class that implements a reasonably compact and efficient library of functions to perform patterns search analysis in a DNA string of length n.
*		The string is eventually stored as a bit set in a structure object called dna_bitset (uint8_t array), being stored as it is or as the reverse complement.
*		The library contains also functions to perform other analysis of sequence features such as the sequence entropy and complexity.
*
*		note: bases different from A, C, G, T, a, c, g and t are not handled throwing an error.
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


// INCLUDE CHECK
#ifndef __IOSTREAM_H_INCLUDED
#define __IOSTREAM_H_INCLUDED
#include <iostream>
#endif /* __IOSTREAM_H_INCLUDED */

#ifndef __STDEXCEPT_H_INCLUDED
#define __STDEXCEPT_H_INCLUDED
#include <stdexcept>
#endif /* __STDEXCEPT_H_INCLUDED */

#ifndef __IOMANIP_H_INCLUDED
#define __IOMANIP_H_INCLUDED
#include <iomanip>
#endif /* __IOMANIP_H_INCLUDED */

#ifndef __STDINT_H_INCLUDED
#define __STDINT_H_INCLUDED
#include <stdint.h>
#endif /*__STDINT_H_INCLUDED */

#ifndef __MATH_H_INCLUDED
#define __MATH_H_INCLUDED
#include <math.h>
#endif /* __MATH_H_INCLUDED */

#ifndef __CSTRING_H_INCLUDED
#define __CSTRING_H_INCLUDED
#include <cstring>
#endif /*__CSTRING_H_INCLUDED */

#ifndef __VECTOR_H_INCLUDED
#define __VECTOR_H_INCLUDED
#include <vector>
#endif /* __VECTOR_H_INCLUDED */

#ifndef __LIST_H_INCLUDED
#define __LIST_H_INCLUDED
#include <list>
#endif /*__LIST_H_INCLUDED */

#ifndef __FUNCTIONS_H_INCLUDED
#define __FUNCTIONS_H_INCLUDED
#include "Functions.h"
#endif /* __FUNCTIONS_H_INCLUDED */

#ifndef __LINKEDLISTKMER_H_INCLUDED
#define __LINKEDLISTKMER_H_INCLUDED
#include "LinkedlistKmer.h"
#endif /* __LINKEDLISTKMER_H_INCLUDED */

#ifndef __HASHTABLE_H_INCLUDED
#define __HASHTABLE_H_INCLUDED
#include "HashTable.h"
#endif /* __HASHTABLE_H_INCLUDED */

#ifndef __BITSCAN_H_INCLUDED
#define __BITSCAN_H_INCLUDED
#include "bitscan/bitscan.h"
#endif /* __BITSCAN_H_INCLUDED */

#ifndef __BITARRAY_H_INCLUDED
#define __BITARRAY_H_INCLUDED
#include "BitArray/bit_array.h"
#endif /* __BITARRAY_H_INCLUDED */


//CONSTANTS DEFINITION
#ifndef DNA_BIT_MASKS
#define DNA_BIT_MASKS

// Definition of a standard bit mask
#define BASE_MASK 0x3	// binary: 11
#define BASE_MASK_SHORT 0x1	// binary: 01
#define BASE_MASK_LONG 0x0F	// binary: 1111

// Definition of a bit encoding for DNA
enum
{
	ENCODING_A = 0x0,	// binary: 00
	ENCODING_C = 0x1,	// binary: 01
	ENCODING_G = 0x2,	// binary: 10
	ENCODING_T = 0x3	// binary: 11
};
#endif /* DNA_BIT_MASKS */


// CLASS
#ifndef NESSIE_H
#define NESSIE_H

/////////////////////////////////////////////////////////////////////////////////////
//
//	STRUCT dna_bitset DEFINITION
//		Structure to store a DNA string of length n into a bit set of length 2*n bit (uint8_t array).
//		Every four DNA bases are stored per 8 bit uint8_t (i.e. 1 byte unsigned integer).
//		The structure also store information on bases statistics for the full DNA string.
//
/////////////////////////////////////////////////////////////////////////////////////
struct dna_bitset{

    uint8_t *data_ptr;	// ptr to the uint8_t array in which to store DNA string as bit set
    size_t data_len;	// length of the DNA string stored (n)
    char *dna_str_ptr;	// ptr to a char array that can store the result of the to_string function
    size_t *array_counts_UP_ptr;	// ptr to size_t array that stores counts of upper case bases in the DNA string
    size_t *array_counts_LOW_ptr;	// ptr to size_t array that stores counts of lower case bases in the DNA string
    double *array_frequencies_UP_ptr;	// ptr to size_t array that stores counts frequencies of upper case bases, p = base count / string length
    double *array_frequencies_LOW_ptr;	// ptr to size_t array that stores counts frequencies of lower case bases, p = base count / string length
};

/////////////////////////////////////////////////////////////////////////////////////
//
//	CLASS Nessie DEFINITION
//		Nessie -- class constructor for encoding DNA and build data structure if needed, non-canonical bases throws an error
//		Nessie -- class constructor used only to build data structure, non-canonical bases are skipped
//		~Nessie -- class destructor
//
//		get_indexes_checked_ptr -- returns a ptr to the bitarray that stores indexes informations on kmers checked
//
//		to_string -- converts the DNA stored as bit in the dna_bitset structure back to a upper case string, returns a ptr to the char array containing the DNA string
//		print_interval_to_string --	prints the DNA stored as bit in the dna_bitset structure as string for an interval
//
//		check_kmer_bit -- returns a ptr to a Kmer object containing informations on a kmer encoded as bit (uint8_t array) in an interval
//		check_kmer_char -- returns a ptr to a Kmer object containing informations on a kmer encoded as a char array in an interval
//		check_kmer_list -- TODO
//		get_kmers_mirror -- returns a ptr to a LinkedlistKmer that stores all the Kmers with mirror symmetry of length [k_min..k_max] in the interval
//		get_kmers_mirror_gap -- returns a ptr to a LinkedlistKmer that stores all the Kmers with mirror symmetry of length [k_min..k_max] in the interval allowing for gaps
//		get_max_kmers_mirror_gap -- returns a ptr to a LinkedlistKmer that stores all the maximum Kmers with mirror symmetry (maximum length max_k and minimum length min_k) in the interval, allows for gaps
//		get_kmers_palindrome -- returns a ptr to a LinkedlistKmer that stores all the Kmers with palindrome symmetry of length [k_min..k_max] in the interval
//		get_kmers_palindrome_gap -- returns a ptr to a LinkedlistKmer that stores all the Kmers with palindrome symmetry of length [k_min..k_max] in the interval allowing for gaps
//		get_max_kmers_palindrome_gap -- returns a ptr to a LinkedlistKmer that stores all the maximum Kmers with palindrome symmetry (maximum length max_k and minimum length min_k) in the interval, allows for gaps
//		print_kmers_mirror --
//		print_kmers_mirror_gap --
//		print_max_kmers_mirror_gap --
//		print_kmers_palindrome --
//		print_kmers_palindrome_gap --
//		print_max_kmers_palindrome_gap --
//		print_kmers -- prints all the Kmers of length [k_min..k_max] in the interval
//
//		shannon_entropy_interval -- returns the Shannon entropy score for an interval
//		print_shannon_entropy_interval --
//		shannon_entropy_sliding -- returns a ptr to a std::vector<double> storing the Shannon entropy scores for a sliding interval
//		print_shannon_entropy_sliding --
//		linguistic_complexity_interval -- returns the Linguistic complexity for an interval
//		print_linguistic_complexity_interval --
//		linguistic_complexity_sliding -- returns a ptr to a std::vector<double> storing the Linguistic complexity for a sliding interval
//		print_linguistic_complexity_sliding --
//
//		get_kmers_triplex_gap -- returns a ptr to a LinkedlistKmer that stores all the Kmers with triplex forming potential of length [k_min..k_max] in the interval, allows for gaps
//		print_kmers_triplex_gap --
//		get_max_kmers_triplex_gap -- returns a ptr to a LinkedlistKmer that stores all the maximum Kmers with with triplex forming potential (maximum length max_k and minimum length min_k) in the interval, allows for gaps
//		print_max_kmers_triplex_gap --
//
/////////////////////////////////////////////////////////////////////////////////////
class Nessie{

private:
	dna_bitset string_bit;	// initializing dna_bitset structure
	dna_bitset *string_bit_ptr;
	sparse_bitarray *array_dimer[16];	// array to store pointers to bit sets (sparse_bitarrays) encoding indexes for dimers in the string
	bitarray *indexes_checked_ptr;	// ptr to a bitarray of length n bit that stores indexes information on kmers checked
									// bits are initialized to 1 and set to 0 when a kmer starting at that index is found
public:
	// Basic functions
	Nessie(const char *dna_str_ptr, size_t dna_str_len, size_t start = 0, size_t end = 0, bool build_structure = false, bool reverse_complement = false);
	Nessie(int, const char *dna_str_ptr, size_t dna_str_len, size_t start = 0, size_t end = 0, bool reverse_complement = false);
	~Nessie();
	void routine_set_array_dimer(sparse_bitarray *array_dimer[], uint8_t *dimer_mask_ptr, bool dimer_set, size_t i, size_t encode);
	bitarray *get_indexes_checked_ptr();

	// Functions for dna_bitset
	char *to_string();
	void print_interval_to_string(size_t start, size_t end, std::ostream &fout = std::cout);

	// Functions for searching a kmer on the Nessie data structure
	Kmer *check_kmer_bit(uint8_t *kmer_bit_ptr, size_t k, size_t start = 0, size_t end = 0);
	void routine_check_kmer_bit(Kmer *kmer_ptr, size_t k, size_t idx, size_t *kmer_dimers_ptr, size_t start, size_t end);
	Kmer *check_kmer_char(const char *kmer_str_ptr, size_t k, size_t start = 0, size_t end = 0, bool reverse_complement = false);
	LinkedlistKmer *check_kmer_list(); //TODO

	// Functions to search kmers with symmetries
	bool check_mirror_symmetry(BIT_ARRAY **monomer_bitarray_ptr, size_t k, size_t max_mm);
	bool check_palindrome_symmetry(BIT_ARRAY **monomer_bitarray_ptr, size_t k, size_t max_mm);
	bool routine_check_mirror_symmetry_interval(BIT_ARRAY **monomer_bitarray_ptr, size_t k, size_t max_mm, size_t start, size_t end, size_t mm_c = 0);
	bool routine_check_palindrome_symmetry_interval(BIT_ARRAY **monomer_bitarray_ptr, size_t k, size_t max_mm, size_t start, size_t end, size_t mm_c = 0);
	bool routine_check_global_alignment(uint8_t *mask_kmer_ptr, std::vector<bool> *align_vector_ptr, size_t k, size_t max_mm, size_t max_gap, size_t max_gapmm, int type);
	bool routine_check_global_alignment_interval(uint8_t *mask_kmer_ptr, std::vector<bool> *align_vector_ptr, size_t k, size_t max_mm, size_t max_gap, size_t max_gapmm, int type, size_t start, size_t end);
	int routine_compare_bases(uint8_t *mask_kmer_ptr, size_t i, size_t j, int m, int mm);
	int routine_compare_bases_complement(uint8_t *mask_kmer_ptr, size_t i, size_t j, int m, int mm);
	HashTable *routine_get_kmers_k_mirror(size_t k, size_t max_mm, size_t start, size_t end);
	HashTable *routine_get_max_kmer_mirror(size_t k_max, size_t k_min, size_t modulo, size_t start, size_t end);
	HashTable *routine_get_kmers_k_gap(size_t k, size_t max_mm, size_t max_gap, size_t max_gapmm, size_t start, size_t end, int type);
	HashTable *routine_get_max_kmer_gap(size_t k_max, size_t k_min, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t start, size_t end, int type);
	HashTable *routine_get_kmers_k_palindrome(size_t k, size_t max_mm, size_t start, size_t end);
	HashTable *routine_get_max_kmer_palindrome(size_t k_max, size_t k_min, size_t modulo, size_t start, size_t end);
	void routine_init_bitarray_and_mask(BIT_ARRAY **monomer_bitarray_ptr, uint8_t *mask_kmer_ptr, size_t k, size_t start);
	void routine_init_mask(uint8_t *mask_kmer_ptr, size_t k, size_t start);
	void routine_shift_bitarray_and_mask(BIT_ARRAY **monomer_bitarray_ptr, uint8_t *mask_kmer_ptr, size_t mask_kmer_len, size_t k, size_t i);
	void routine_shift_mask(uint8_t *mask_kmer_ptr, size_t mask_kmer_len, size_t k, size_t i);
	LinkedlistKmer *get_kmers_mirror(size_t k_min, size_t k_max = 0, size_t modulo = 0, size_t start = 0, size_t end = 0);
	LinkedlistKmer *get_kmers_mirror_gap(size_t k_min, size_t k_max = 0, size_t modulo = 0, size_t modulo_gap = 0, size_t modulo_gapmm = 0, size_t start = 0, size_t end = 0);
	LinkedlistKmer *get_max_kmers_mirror_gap(size_t k_max, size_t k_min = 0, size_t modulo = 0, size_t modulo_gap = 0, size_t modulo_gapmm = 0, size_t start = 0, size_t end = 0);
	LinkedlistKmer *get_kmers_palindrome(size_t k_min, size_t k_max = 0, size_t modulo = 0, size_t start = 0, size_t end = 0);
	LinkedlistKmer *get_kmers_palindrome_gap(size_t k_min, size_t k_max = 0, size_t modulo = 0, size_t modulo_gap = 0, size_t modulo_gapmm = 0, size_t start = 0, size_t end = 0);
	LinkedlistKmer *get_max_kmers_palindrome_gap(size_t k_max, size_t k_min = 0, size_t modulo = 0, size_t modulo_gap = 0, size_t modulo_gapmm = 0, size_t start = 0, size_t end = 0);
	void print_kmers_mirror(size_t k_min, size_t k_max = 0, size_t modulo = 0, size_t start = 0, size_t end = 0, std::ostream &fout = std::cout, bool counts = true, bool indexes = true, size_t start_idx = 0);
	void print_kmers_mirror_gap(size_t k_min, size_t k_max = 0, size_t modulo = 0, size_t modulo_gap = 0, size_t modulo_gapmm = 0, size_t start = 0, size_t end = 0, std::ostream &fout = std::cout, bool counts = true, bool indexes = true, size_t start_idx = 0);
	void print_max_kmers_mirror_gap(size_t k_max, size_t k_min = 0, size_t modulo = 0, size_t modulo_gap = 0, size_t modulo_gapmm = 0, size_t start = 0, size_t end = 0, std::ostream &fout = std::cout, bool counts = true, bool indexes = true, size_t start_idx = 0);
	void print_kmers_palindrome(size_t k_min, size_t k_max = 0, size_t modulo = 0, size_t start = 0, size_t end = 0, std::ostream &fout = std::cout, bool counts = true, bool indexes = true, size_t start_idx = 0);
	void print_kmers_palindrome_gap(size_t k_min, size_t k_max = 0, size_t modulo = 0, size_t modulo_gap = 0, size_t modulo_gapmm = 0, size_t start = 0, size_t end = 0, std::ostream &fout = std::cout, bool counts = true, bool indexes = true, size_t start_idx = 0);
	void print_max_kmers_palindrome_gap(size_t k_max, size_t k_min = 0, size_t modulo = 0, size_t modulo_gap = 0, size_t modulo_gapmm = 0, size_t start = 0, size_t end = 0, std::ostream &fout = std::cout, bool counts = true, bool indexes = true, size_t start_idx = 0);
	HashTable *routine_get_kmers_k(size_t k, size_t start, size_t end);
	void print_kmers(size_t k_min, size_t k_max = 0, size_t start = 0, size_t end = 0, std::ostream &fout = std::cout, bool counts = true, bool indexes = true, size_t start_idx = 0);

	// Shannon entropy
	void routine_init_counts(size_t *array_counts_ptr, size_t start, size_t end);
	void routine_shift_counts(size_t *array_counts_ptr, size_t interval_len, size_t shift, size_t start);
	double routine_shannon_entropy(size_t *array_counts_ptr, size_t sequence_len);
	std::vector<double> *shannon_entropy_sliding(size_t interval_len, size_t shift, size_t start = 0, size_t end = 0);
	void print_shannon_entropy_sliding(size_t interval_len, size_t shift, std::ostream &fout = std::cout, size_t start = 0, size_t end = 0);
	double shannon_entropy_interval(size_t start = 0, size_t end = 0);
	void print_shannon_entropy_interval(std::ostream &fout = std::cout, size_t start = 0, size_t end = 0);

	// Linguistic complexity
	double linguistic_complexity_interval(size_t start = 0, size_t end = 0, size_t k_min = 0, size_t k_max = 0);
	void print_linguistic_complexity_interval(std::ostream &fout = std::cout, size_t start = 0, size_t end = 0, size_t k_min = 0, size_t k_max = 0);
	std::vector<double> *linguistic_complexity_sliding(size_t interval_len = 0, size_t shift = 0, size_t start = 0, size_t end = 0, size_t k_min = 0, size_t k_max = 0);
	void print_linguistic_complexity_sliding(size_t interval_len = 0, size_t shift = 0, std::ostream &fout = std::cout, size_t start = 0, size_t end = 0, size_t k_min = 0, size_t k_max = 0);
	uint64_t routine_convert_to_uint64(uint8_t *mask_kmer_ptr, size_t mask_kmer_len);
	std::list<uint64_t> *routine_get_kmers_k_unique(size_t k, size_t start, size_t end);

	// Quadruplex
	//TODO

	// Triplex
	bool routine_check_triplex_forming(BIT_ARRAY **monomer_bitarray_ptr, size_t k, size_t max_mm, size_t max_purine);
	bool routine_check_triplex_forming_interval(BIT_ARRAY **monomer_bitarray_ptr, size_t k, size_t max_mm, size_t max_purine, size_t start, size_t end);
	bool routine_check_triplex_forming_gap(uint8_t *mask_kmer_ptr, std::vector<bool> *align_vector_ptr, size_t k, size_t max_mm, size_t max_gap, size_t max_gapmm, size_t max_purine);
	bool routine_check_triplex_forming_gap_interval(uint8_t *mask_kmer_ptr, std::vector<bool> *align_vector_ptr, size_t k, size_t max_mm, size_t max_gap, size_t max_gapmm, size_t max_purine, size_t start, size_t end);
	HashTable *routine_get_kmers_k_triplex(size_t k, size_t max_mm, size_t max_purine, size_t start, size_t end);
	HashTable *routine_get_kmers_k_triplex_gap(size_t k, size_t max_mm, size_t max_gap, size_t max_gapmm, size_t max_purine, size_t start, size_t end);
	HashTable *routine_get_max_kmer_triplex(size_t k_max, size_t k_min, size_t modulo, size_t modulo_purine, size_t start, size_t end);
	HashTable *routine_get_max_kmer_triplex_gap(size_t k_max, size_t k_min, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t modulo_purine, size_t start, size_t end);
	LinkedlistKmer *get_kmers_triplex_gap(size_t k_min, size_t k_max = 0, size_t modulo = 0, size_t modulo_gap = 0, size_t modulo_gapmm = 0, size_t modulo_purine = 0, size_t start = 0, size_t end = 0);
	LinkedlistKmer *get_max_kmers_triplex_gap(size_t k_max, size_t k_min = 0, size_t modulo = 0, size_t modulo_gap = 0, size_t modulo_gapmm = 0, size_t modulo_purine = 0, size_t start = 0, size_t end = 0);
	void print_kmers_triplex_gap(size_t k_min, size_t k_max = 0, size_t modulo = 0, size_t modulo_gap = 0, size_t modulo_gapmm = 0, size_t modulo_purine = 0, size_t start = 0, size_t end = 0,
								 std::ostream &fout = std::cout, bool counts = true, bool indexes = true, size_t start_idx = 0);
	void print_max_kmers_triplex_gap(size_t k_max, size_t k_min = 0, size_t modulo = 0, size_t modulo_gap = 0, size_t modulo_gapmm = 0, size_t modulo_purine = 0, size_t start = 0, size_t end = 0,
									 std::ostream &fout = std::cout, bool counts = true, bool indexes = true, size_t start_idx = 0);

};

#endif /* NESSIE_H */

