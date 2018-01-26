/**************************************************************************************
*
**	FUNCTIONS (Nessie.cpp)
*		Implements the functions of the Nessie class.
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



#include "Nessie.h"

/////////////////////////////////////////////////////////////////////////////////////
//
//	Nessie (constructor): converts the DNA string into a bit set (dna_bitset) and initializes the Nessie data structure if needed,
//						  does not accept non-canonical bases
//
//	parameters:
//		dna_str_ptr - a ptr to a char array containing the DNA string
//		dna_str_len - length of the DNA string (i.e. n)
//		start - starting index of the interval to be stored
//		end - ending index of the interval to be stored
//		build_structure - bool value, if true array_dimer is initialized and bits are set / if false bits are not set [false]
//		reverse_complement - bool value, if false stores the normal string / if true stores the string in reverse complement [false]
//
/////////////////////////////////////////////////////////////////////////////////////
Nessie::Nessie(const char *dna_str_ptr, size_t dna_str_len, size_t start, size_t end, bool build_structure, bool reverse_complement){

	if (!end){ end = dna_str_len - 1; }
	if (start > end){ throw std::invalid_argument("Nessie build: starting index is larger than ending index"); }

	// Defining some variables
	uint8_t shift_DNA;
	uint8_t dimer_mask = 0x0;	// bit set mask to keep track of the dimers to initialize the indexes
	uint8_t *dimer_mask_ptr = &dimer_mask;
	string_bit_ptr = &string_bit;
	size_t dna_len = end - start + 1;

	// Bytes necessary to store the DNA string as a bit set
	size_t dna_bytes = (dna_len >> 2) + (0 != (dna_len & ((1 << 2) - 1)));	// (dna_len / 4) + (0 != (dna_len % 4))
	    																	// (0 != dna_len % 4) evaluate to 1 when dna_len is not multiple of 4
	// Initializing Nessie data structure
	for(int i = 0; i < 16; ++i){	// array_dimer initialization
		array_dimer[i] = (build_structure) ? new sparse_bitarray(dna_len) : NULL; // creating the i-bit set (sparse_bitarray) of length n bit if build structure
	}

	//
	if (build_structure) {
		indexes_checked_ptr = new bitarray(dna_len);	// indexes_checked_ptr bit set (bitarray) initialization
		indexes_checked_ptr->set_bit(0, (dna_len - 1));	// set to 1 all bits
	}
	else{
		indexes_checked_ptr = NULL;
	}

	// Initializing string_bit
	string_bit_ptr->data_len = dna_len;
	string_bit_ptr->data_ptr = new uint8_t[dna_bytes];	// string_bit uint8_t array in which to store the DNA string encoded as bit
	std::memset(string_bit_ptr->data_ptr, 0, dna_bytes);	// initializing every bit of the string_bit array to 0
	string_bit_ptr->dna_str_ptr = NULL;	// initializing to NULL the ptr to the char array that stores to_string output

	// Arrays string_bit
	string_bit_ptr->array_counts_UP_ptr = new size_t[4]();
	string_bit_ptr->array_counts_LOW_ptr = new size_t[4]();
	string_bit_ptr->array_frequencies_UP_ptr = new double[4]();
	string_bit_ptr->array_frequencies_LOW_ptr = new double[4]();

	// Storing DNA into string_bit bit set and filling Nessie data structure
	if (reverse_complement){	// storing the DNA string in reverse complement
		bool dimer_set = false;	// becomes true after first iteration, the first dimer is available only while iterating trough the second base

		for (size_t i = 0; i < dna_len; ++i){
			shift_DNA = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
													// 2 * (i % 4)
			switch (dna_str_ptr[(end - i)]){ // index reversed to pick bases from the end of the DNA string
			case 'A':
				++string_bit_ptr->array_counts_UP_ptr[ENCODING_T];	// incrementing base counts
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_T << shift_DNA;	// encoding base as bit
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_T);	// setting bit corresponding to dimer index
				break;
			case 'a':
				++string_bit_ptr->array_counts_LOW_ptr[ENCODING_T];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_T << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_T);
				break;
			case 'C':
				++string_bit_ptr->array_counts_UP_ptr[ENCODING_G];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_G << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_G);
				break;
			case 'c':
				++string_bit_ptr->array_counts_LOW_ptr[ENCODING_G];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_G << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_G);
				break;
			case 'G':
				++string_bit_ptr->array_counts_UP_ptr[ENCODING_C];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_C << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_C);
				break;
			case 'g':
				++string_bit_ptr->array_counts_LOW_ptr[ENCODING_C];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_C << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_C);
				break;
			case 'T':
				++string_bit_ptr->array_counts_UP_ptr[ENCODING_A];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_A << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_A);
				break;
			case 't':
				++string_bit_ptr->array_counts_LOW_ptr[ENCODING_A];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_A << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_A);
				break;
			default:
				throw std::invalid_argument("invalid DNA base");
			}
			dimer_set = (build_structure) ? true : false;
		}
	}
	else{	// storing the DNA string as it is
		bool dimer_set = false;

		for (size_t i = 0; i < dna_len; ++i){
			shift_DNA = (i & ((1 << 2) - 1)) << 1;

			switch (dna_str_ptr[(start + i)]){
			case 'A':
				++string_bit_ptr->array_counts_UP_ptr[ENCODING_A];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_A << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_A);
				break;
			case 'a':
				++string_bit_ptr->array_counts_LOW_ptr[ENCODING_A];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_A << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_A);
				break;
			case 'C':
				++string_bit_ptr->array_counts_UP_ptr[ENCODING_C];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_C << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_C);
				break;
			case 'c':
				++string_bit_ptr->array_counts_LOW_ptr[ENCODING_C];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_C << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_C);
				break;
			case 'G':
				++string_bit_ptr->array_counts_UP_ptr[ENCODING_G];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_G << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_G);
				break;
			case 'g':
				++string_bit_ptr->array_counts_LOW_ptr[ENCODING_G];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_G << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_G);
				break;
			case 'T':
				++string_bit_ptr->array_counts_UP_ptr[ENCODING_T];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_T << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_T);
				break;
			case 't':
				++string_bit_ptr->array_counts_LOW_ptr[ENCODING_T];
				string_bit_ptr->data_ptr[i >> 2] |= ENCODING_T << shift_DNA;
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_T);
				break;
			default:
				throw std::invalid_argument("invalid DNA base");
			}
			dimer_set = (build_structure) ? true : false;
		}
	}

	// Calculating bases frequencies and storing in string_bit structure
	for (size_t i = 0; i < 4; ++i){
		string_bit_ptr->array_frequencies_UP_ptr[i] = (string_bit_ptr->array_counts_UP_ptr[i])
													? (double) string_bit_ptr->array_counts_UP_ptr[i] / dna_len
													: 0;
		string_bit_ptr->array_frequencies_LOW_ptr[i] = (string_bit_ptr->array_counts_LOW_ptr[i])
													? (double) string_bit_ptr->array_counts_LOW_ptr[i] / dna_len
													: 0;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	Nessie (second constructor): build only Nessie data structure,
//								 skips non-canonical bases
//
//	parameters:
//		int - just an integer to discriminate between the two constructors
//		dna_str_ptr - a ptr to a char array containing the DNA string
//		dna_str_len - length of the DNA string (i.e. n)
//		start - starting index of the interval to be stored
//		end - ending index of the interval to be stored
//		dna_len - length of the DNA string (i.e. n)
//		reverse_complement - bool value, if false stores the normal string / if true stores the string in reverse complement [false]
//
/////////////////////////////////////////////////////////////////////////////////////
Nessie::Nessie(int, const char *dna_str_ptr, size_t dna_str_len, size_t start, size_t end, bool reverse_complement){

	if (!end){ end = dna_str_len - 1; }
	if (start > end){ throw std::invalid_argument("Nessie build: starting index is larger than ending index"); }

	// Defining some variables
	uint8_t shift_DNA;
	uint8_t dimer_mask = 0x0;	// bit set mask to keep track of the dimers to initialize the indexes
	uint8_t *dimer_mask_ptr = &dimer_mask;
	string_bit_ptr = &string_bit;
	size_t dna_len = end - start + 1;

	// Initializing Nessie data structure
	for(int i = 0; i < 16; ++i){	// array_dimer initialization
		array_dimer[i] = new sparse_bitarray(dna_len);	// creating the i-bit set (sparse_bitarray) of length n bit
	}

	//
	indexes_checked_ptr = new bitarray(dna_len);	// indexes_checked_ptr bit set (bitarray) initialization
	indexes_checked_ptr->set_bit(0, (dna_len - 1));	// set to 1 all bits

	// Initializing string_bit
	string_bit_ptr->data_len = dna_len;
	string_bit_ptr->data_ptr = NULL;	// string_bit uint8_t array in which to store the DNA string encoded as bit
	string_bit_ptr->dna_str_ptr = NULL;	// initializing to NULL the ptr to the char array that stores to_string output

	// Arrays string_bit
	string_bit_ptr->array_counts_UP_ptr = NULL;
	string_bit_ptr->array_counts_LOW_ptr = NULL;
	string_bit_ptr->array_frequencies_UP_ptr = NULL;
	string_bit_ptr->array_frequencies_LOW_ptr = NULL;

	// Storing DNA into string_bit bit set and filling Nessie data structure
	if (reverse_complement){	// storing the DNA string in reverse complement
		bool dimer_set = false;	// becomes true after first iteration, the first dimer is available only while iterating trough the second base

		for (size_t i = 0; i < dna_len; ++i){
			shift_DNA = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
													// 2 * (i % 4)
			switch (dna_str_ptr[(end - i)]){ // index reversed to pick bases from the end of the DNA string
			case 'A':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_T);	// setting bit corresponding to dimer index
				dimer_set = true;
				break;
			case 'a':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_T);
				dimer_set = true;
				break;
			case 'C':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_G);
				dimer_set = true;
				break;
			case 'c':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_G);
				dimer_set = true;
				break;
			case 'G':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_C);
				dimer_set = true;
				break;
			case 'g':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_C);
				dimer_set = true;
				break;
			case 'T':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_A);
				dimer_set = true;
				break;
			case 't':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_A);
				dimer_set = true;
				break;
			default:
				if (dimer_set){	//skipping non-canonical bases and reinitializing dimer mask and dimer_set
					dimer_mask = 0x0;
					dimer_set = false;
				}
			}
		}
	}
	else{	// storing the DNA string as it is
		bool dimer_set = false;

		for (size_t i = 0; i < dna_len; ++i){
			shift_DNA = (i & ((1 << 2) - 1)) << 1;

			switch (dna_str_ptr[(start + i)]){
			case 'A':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_A);
				dimer_set = true;
				break;
			case 'a':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_A);
				dimer_set = true;
				break;
			case 'C':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_C);
				dimer_set = true;
				break;
			case 'c':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_C);
				dimer_set = true;
				break;
			case 'G':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_G);
				dimer_set = true;
				break;
			case 'g':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_G);
				dimer_set = true;
				break;
			case 'T':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_T);
				dimer_set = true;
				break;
			case 't':
				Nessie::routine_set_array_dimer(array_dimer, dimer_mask_ptr, dimer_set, i, ENCODING_T);
				dimer_set = true;
				break;
			default:
				if (dimer_set){
					dimer_mask = 0x0;
					dimer_set = false;
				}
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	~Nessie (destructor): destructs the Nessie data structure and the dna_bitset structure containing the DNA string encoded as a bit set
//
/////////////////////////////////////////////////////////////////////////////////////
Nessie::~Nessie(){

	// Destructing the string_bit uint8_t array
	if (string_bit_ptr->data_ptr){
		delete[] string_bit_ptr->data_ptr;
		string_bit_ptr->data_ptr = NULL;
		//std::cout << "DELETE BIT ARRAY" << std::endl;
	}

	// Destructing counts and frequencies size_t arrays
	if (string_bit_ptr->array_counts_UP_ptr){
		delete[] string_bit_ptr->array_counts_UP_ptr;
		string_bit_ptr->array_counts_UP_ptr = NULL;
	}

	if (string_bit_ptr->array_counts_LOW_ptr){
		delete[] string_bit_ptr->array_counts_LOW_ptr;
		string_bit_ptr->array_counts_LOW_ptr = NULL;
	}

	if (string_bit_ptr->array_frequencies_UP_ptr){
		delete[] string_bit_ptr->array_frequencies_UP_ptr;
		string_bit_ptr->array_frequencies_UP_ptr = NULL;
	}

	if (string_bit_ptr->array_frequencies_LOW_ptr){
		delete[] string_bit_ptr->array_frequencies_LOW_ptr;
		string_bit_ptr->array_frequencies_LOW_ptr = NULL;
	}
	//std::cout << "DELETE COUNTS ARRAY" << std::endl;

	// Destructing string_bit char array containing the DNA string if existent
	if (string_bit_ptr->dna_str_ptr){
		delete[] string_bit_ptr->dna_str_ptr;
		string_bit_ptr->dna_str_ptr = NULL;
		//std::cout << "DELETE CHAR ARRAY" << std::endl;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_set_array_dimer
//
//	parameters:
//		array_dimer - array that stores dimers indexes information
//		dimer_mask_ptr - ptr to a mask that stores the current dimer
//		dimer_set - bool value that check if the dimer corresponding index has to be set
//		i - index from the loop
//		encode - bit encoding for the base as size_t
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::routine_set_array_dimer(sparse_bitarray *array_dimer[], uint8_t *dimer_mask_ptr, bool dimer_set, size_t i, size_t encode){

	if (dimer_set){ // Initializing dimers
		*dimer_mask_ptr >>= 2;	// 0000xx00 -> 000000xx
		*dimer_mask_ptr |= encode << 2;	// 0000yyxx
		array_dimer[*dimer_mask_ptr]->set_bit(i - 1);	// setting dimer starting index
	}
	else{
		*dimer_mask_ptr |= encode << 2;	// first iteration no dimers discovered yet, filling 0000xx00
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	get_indexes_checked_ptr: returns a ptr to the bitarray that stores indexes informations on kmers checked
//
/////////////////////////////////////////////////////////////////////////////////////
bitarray *Nessie::get_indexes_checked_ptr(){

	return indexes_checked_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	to_string: converts the DNA stored as bit in the dna_bitset structure back to a upper case string,
//			   returns a ptr to the char array containing the DNA string
//
/////////////////////////////////////////////////////////////////////////////////////
char *Nessie::to_string(){

	// Creating a ptr to the new char array that will store the DNA as string
	string_bit_ptr->dna_str_ptr = new char[string_bit_ptr->data_len + 1];

	// Reading back DNA to string
	for (size_t i = 0; i < string_bit_ptr->data_len; ++i){
		uint8_t shift_DNA = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
														// 2 * (i % 4);
		uint8_t base = (string_bit_ptr->data_ptr[i >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING

		switch(base){
		case ENCODING_A:
			string_bit_ptr->dna_str_ptr[i] = 'A';
			break;
		case ENCODING_C:
			string_bit_ptr->dna_str_ptr[i] = 'C';
			break;
		case ENCODING_G:
			string_bit_ptr->dna_str_ptr[i] = 'G';
			break;
		case ENCODING_T:
			string_bit_ptr->dna_str_ptr[i] = 'T';
			break;
		default:
			throw std::runtime_error("invalid DNA base");
		}
	}

	string_bit_ptr->dna_str_ptr[string_bit_ptr->data_len] = '\0';	// set last char to \0

	return string_bit_ptr->dna_str_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_interval_to_string: prints the DNA stored as bit in the dna_bitset structure as string for an interval
//
//	parameters
//		fout - ostream element to be used for printing [cout]
//		start - starting index of the interval to print
//		end - ending index of the interval to print
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::print_interval_to_string(size_t start, size_t end, std::ostream &fout){

	if (start > end){ throw std::invalid_argument("Print interval: starting index is larger than ending index"); }

	for (size_t i = start; i <= end; ++i){
		uint8_t shift_DNA = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
														// 2 * (i % 4);
		uint8_t base = (string_bit_ptr->data_ptr[i >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING

		switch(base){
		case ENCODING_A:
			fout << 'A';
			break;
		case ENCODING_C:
			fout << 'C';
			break;
		case ENCODING_G:
			fout << 'G';
			break;
		case ENCODING_T:
			fout << 'T';
			break;
		default:
			throw std::runtime_error("invalid DNA base");
		}
	}
	fout << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	check_kmer_bit -- returns a ptr to a Kmer object containing informations on a kmer encoded as bit (uint8_t array) in an interval
//
//	parameters:
//		kmer_bit_ptr - a ptr to a kmer encoded as a uint8_t array to search for
//		k - length of the kmer
//		start - starting index of the interval considered [0]
//	 	end - ending index of the interval considered [0]
//
/////////////////////////////////////////////////////////////////////////////////////
Kmer *Nessie::check_kmer_bit(uint8_t *kmer_bit_ptr, size_t k, size_t start, size_t end){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }
	if (start > end){ throw std::invalid_argument("Check kmer bit: starting index is larger than ending index"); }
	if (k > (end - start + 1)){ throw std::invalid_argument("Check kmer bit: motif is longer than sequence interval"); }

	// Defining some variables
	size_t idx = k >> 1;	// dimers index
	size_t kmer_dimers[idx];	// look up array to retrieve subsequent dimers identity (to be used as indexes within the array_dimer)
	size_t *kmer_dimers_ptr = kmer_dimers;
	uint8_t shift_monomer;
	uint8_t shift_dimer;

	// Bytes necessary to store the kmer as a bit set
	size_t dna_bytes = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

	// Defining Kmer object
	Kmer *kmer_ptr = new Kmer(k);

	// Initializing kmer_dimers array
	if (k & 1){	// odd kmer
		for (size_t i = 0; i < idx; ++i){
			shift_dimer = (i & ((1 << 1) - 1)) << 2;
			kmer_dimers_ptr[i] = (kmer_bit_ptr[i >> 1] & (BASE_MASK_LONG << shift_dimer)) >> shift_dimer;
		}
		uint8_t last_dimer_mask = 0x0;

		shift_monomer = ((k - 1) & ((1 << 2) - 1)) << 1;	// 2 * ((k-1) % 4)
		last_dimer_mask |= ((kmer_bit_ptr[(k - 1) >> 2] & (BASE_MASK << shift_monomer)) >> shift_monomer);

		shift_monomer = ((k - 2) & ((1 << 2) - 1)) << 1;
		kmer_dimers_ptr[idx] = ((last_dimer_mask << 2) | ((kmer_bit_ptr[(k - 2) >> 2] & (BASE_MASK << shift_monomer)) >> shift_monomer));
	}
	else{	// even kmer
		for (size_t i = 0; i < idx; ++i){
			shift_dimer = (i & ((1 << 1) - 1)) << 2;
			kmer_dimers_ptr[i] = (kmer_bit_ptr[i >> 1] & (BASE_MASK_LONG << shift_dimer)) >> shift_dimer;
		}
	}

	// Checking windows
	size_t step = 5000000;
	if ((end - start + 1) > step){ //std::cout << "STEP " << std::endl;
		size_t i = start;	// index
		while(true){
			size_t end_i = (i + step + k - 2 > end) ? end : i + step + k - 2;
			if (end - end_i > step + k - 1){
				Nessie::routine_check_kmer_bit(kmer_ptr, k, idx, kmer_dimers_ptr, i, end_i);
			}
			else{
				Nessie::routine_check_kmer_bit(kmer_ptr, k, idx, kmer_dimers_ptr, i , end);
				break;
			}
			i += step;
		}
	}
	else{
		Nessie::routine_check_kmer_bit(kmer_ptr, k, idx, kmer_dimers_ptr, start , end);
	}

	// Assigning remaining informations to Kmer object
	kmer_ptr->counts = kmer_ptr->indexes.size();
	if (kmer_ptr->counts){
		copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, kmer_bit_ptr, dna_bytes);
	}
	else{
		throw std::runtime_error("kmer not found");
	}

	return kmer_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_check_kmer_bit
//
//	parameters:
//		kmer_ptr - a ptr to the Kmer object for the kmer searched
//		k - length of the kmer
//		idx - dimers index
//		kmer_dimers_ptr - ptr to the array containing the decomposition of the kmer in dimers
//		start - starting index of the interval considered
//		end - ending index of the interval considered
//
// 	nota: per velocizzare dovrei riuscire a tirare loop fuori dal while
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::routine_check_kmer_bit(Kmer *kmer_ptr, size_t k, size_t idx, size_t *kmer_dimers_ptr, size_t start, size_t end){

	// Variables
	size_t nBit=EMPTY_ELEM;

	// Reading indexes associated to the first dimer of the kmer and saving into vector kmer_ptr.indexes if a match for the full kmer is found
	//bool first = true;	// allows to keep set the first index for the kmer in indexes_checked_ptr
	if (0 == start){	// searching from the beginning
		if (!(k & 1)){	// even kmer
			array_dimer[kmer_dimers_ptr[0]]->init_scan(bbo::NON_DESTRUCTIVE);
			while(true){
				bool add = true;
				nBit=array_dimer[kmer_dimers_ptr[0]]->next_bit();
				if((nBit==EMPTY_ELEM) || (nBit > (end - k + 1))) break;

				for (size_t i = 1; i < idx; ++i){
					if (!(array_dimer[kmer_dimers_ptr[i]]->is_bit(nBit + (i << 1)))){
						add = false;
						break;
					}
				}

				if (add){
					//if (first){
						kmer_ptr->indexes.push_back(nBit);
						//first = false;
					//}
					//else{
						//kmer_ptr->indexes.push_back(nBit);
						indexes_checked_ptr->erase_bit((nBit-start));
					//}
				}
			}
		}
		else{	// odd kmer
			array_dimer[kmer_dimers_ptr[0]]->init_scan(bbo::NON_DESTRUCTIVE);
			while(true){
				bool add = true;
				nBit=array_dimer[kmer_dimers_ptr[0]]->next_bit();
				if((nBit==EMPTY_ELEM) || (nBit > (end - k + 1))) break;

				for (size_t i = 1; i <= idx; ++i){
					if (i != idx){
						if (!(array_dimer[kmer_dimers_ptr[i]]->is_bit(nBit + (i << 1)))){
							add = false;
							break;
						}
					}
					else{
						if (!(array_dimer[kmer_dimers_ptr[i]]->is_bit(nBit + (i << 1) - 1))){
							add = false;
							break;
						}
					}
				}

				if (add){
					//if (first){
						kmer_ptr->indexes.push_back(nBit);
						//first = false;
					//}
					//else{
						//kmer_ptr->indexes.push_back(nBit);
						indexes_checked_ptr->erase_bit((nBit-start));
					//}
				}
			}
		}
	}
	else{	// searching from an index
		if (!(k & 1)){	// even kmer
			array_dimer[kmer_dimers_ptr[0]]->init_scan_from((start - 1), bbo::NON_DESTRUCTIVE);
			while(true){
				bool add = true;
				nBit=array_dimer[kmer_dimers_ptr[0]]->next_bit();
				if((nBit==EMPTY_ELEM) || (nBit > (end - k + 1))) break;

				for (size_t i = 1; i < idx; ++i){
					if (!(array_dimer[kmer_dimers_ptr[i]]->is_bit(nBit + (i << 1)))){
						add = false;
						break;
					}
				}

				if (add){
					//if (first){
						kmer_ptr->indexes.push_back(nBit);
						//first = false;
					//}
					//else{
						//kmer_ptr->indexes.push_back(nBit);
						indexes_checked_ptr->erase_bit((nBit-start));
					//}
				}
			}
		}
		else{	// odd kmer
			array_dimer[kmer_dimers_ptr[0]]->init_scan_from((start - 1), bbo::NON_DESTRUCTIVE);
			while(true){
				bool add = true;
				nBit=array_dimer[kmer_dimers_ptr[0]]->next_bit();
				if((nBit==EMPTY_ELEM) || (nBit > (end - k + 1))) break;

				for (size_t i = 1; i <= idx; ++i){
					if (i != idx){
						if (!(array_dimer[kmer_dimers_ptr[i]]->is_bit(nBit + (i << 1)))){
							add = false;
							break;
						}
					}
					else{
						if (!(array_dimer[kmer_dimers_ptr[i]]->is_bit(nBit + (i << 1) - 1))){
							add = false;
							break;
						}
					}
				}

				if (add){
					//if (first){
						kmer_ptr->indexes.push_back(nBit);
						//first = false;
					//}
					//else{
						//kmer_ptr->indexes.push_back(nBit);
						indexes_checked_ptr->erase_bit((nBit-start));
					//}
				}
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	check_kmer_char -- returns a ptr to a Kmer object containing informations on a kmer encoded as a char array in an interval
//
//	parameters:
//		kmer_str_ptr - a ptr to the char array containing the kmer to search for
// 		k - length of the kmer
//		start - starting index of the interval considered [0]
//		end - ending index of the interval considered [0]
//		reverse_complement - bool value, if false stores the normal string / if true stores the string in reverse complement [false]
//
//	note: this function basically convert the kmer from string to a bit set (uint8_t array) and call check_kmer_bit on it
//
/////////////////////////////////////////////////////////////////////////////////////
Kmer *Nessie::check_kmer_char(const char *kmer_str_ptr, size_t k, size_t start, size_t end, bool reverse_complement){

	// Bytes necessary to store the kmer as a bit set
	size_t dna_bytes = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	//	(k / 4) + (0 != (k % 4))

	// Initializing uint8_t array in which to store the kmer
	uint8_t kmer_bit[dna_bytes];	// defining kmer_bit array
	uint8_t *kmer_bit_ptr = kmer_bit;	// defining ptr to kmer_bit array
	std::memset(kmer_bit_ptr, 0, dna_bytes);	// set every bit of the uint8_t array to 0

	// Converting kmer into a bit set (uint8_t array)
	if (reverse_complement){	// storing the DNA string in reverse complement
		for (size_t i = 0; i < k; ++i){
			uint8_t shift_kmer = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array

			switch (kmer_str_ptr[(k - i - 1)]){ // index reversed to pick bases from the end of the DNA string
			case 'A':
				kmer_bit_ptr[i >> 2] |= ENCODING_T << shift_kmer;
				break;
			case 'a':
				kmer_bit_ptr[i >> 2] |= ENCODING_T << shift_kmer;
				break;
			case 'C':
				kmer_bit_ptr[i >> 2] |= ENCODING_G << shift_kmer;
				break;
			case 'c':
				kmer_bit_ptr[i >> 2] |= ENCODING_G << shift_kmer;
				break;
			case 'G':
				kmer_bit_ptr[i >> 2] |= ENCODING_C << shift_kmer;
				break;
			case 'g':
				kmer_bit_ptr[i >> 2] |= ENCODING_C << shift_kmer;
				break;
			case 'T':
				kmer_bit_ptr[i >> 2] |= ENCODING_A << shift_kmer;
				break;
			case 't':
				kmer_bit_ptr[i >> 2] |= ENCODING_A << shift_kmer;
				break;
			default:
				throw std::invalid_argument("invalid DNA base");
			}
		}
	}
	else{	// storing the DNA string as it is

		for (size_t i = 0; i < k; ++i){
			uint8_t shift_kmer = (i & ((1 << 2) - 1)) << 1;

			switch (kmer_str_ptr[i]){
			case 'A':
				kmer_bit_ptr[i >> 2] |= ENCODING_A << shift_kmer;
				break;
			case 'a':
				kmer_bit_ptr[i >> 2] |= ENCODING_A << shift_kmer;
				break;
			case 'C':
				kmer_bit_ptr[i >> 2] |= ENCODING_C << shift_kmer;
				break;
			case 'c':
				kmer_bit_ptr[i >> 2] |= ENCODING_C << shift_kmer;
				break;
			case 'G':
				kmer_bit_ptr[i >> 2] |= ENCODING_G << shift_kmer;
				break;
			case 'g':
				kmer_bit_ptr[i >> 2] |= ENCODING_G << shift_kmer;
				break;
			case 'T':
				kmer_bit_ptr[i >> 2] |= ENCODING_T << shift_kmer;
				break;
			case 't':
				kmer_bit_ptr[i >> 2] |= ENCODING_T << shift_kmer;
				break;
			default:
				throw std::invalid_argument("invalid DNA base");
			}
		}
	}

	return Nessie::check_kmer_bit(kmer_bit_ptr, k, start, end);
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	check_mirror_symmetry -- check the BIT_ARRAY encoding the kmer for the mirror symmetry
//
//	parameters:
//		monomer_bitarray_ptr - ptr to a *BIT_ARRAY[4]
//		k - length of the kmer
//		max_mm - max number of mismatch allowed
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
bool Nessie::check_mirror_symmetry(BIT_ARRAY **monomer_bitarray_ptr, size_t k, size_t max_mm){

	// Variables
	bool mirror = false;

	// Checking simmetry
	size_t mm_c = 0;	// mismatch counter
	for(size_t i = 0; i < 4; ++i){
		if (bit_array_get(monomer_bitarray_ptr[i], 0) && !bit_array_get(monomer_bitarray_ptr[i], k - 1)) {
			return mirror;
		}
		else{
			for(size_t j = 1; j < (k >> 1); ++j){
				if (bit_array_get(monomer_bitarray_ptr[i], j) && !bit_array_get(monomer_bitarray_ptr[i], k - 1 - j)) {
					++mm_c;
				}
			}
		}
	}

	if (mm_c <= max_mm){
		mirror = true;
	}

	return mirror;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	check_palindrome_symmetry -- check the BIT_ARRAY encoding the kmer for the palindrome symmetry
//
//	parameters:
//		monomer_bitarray_ptr - ptr to a *BIT_ARRAY[4]
//		k - length of the kmer
//		max_mm - max number of mismatch allowed
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
bool Nessie::check_palindrome_symmetry(BIT_ARRAY **monomer_bitarray_ptr, size_t k, size_t max_mm){

	// Variables
	bool palindrome = false;

	// Checking simmetry
	size_t mm_c = 0;	// mismatch counter
	for(size_t i = 0; i < 4; ++i){
		if (bit_array_get(monomer_bitarray_ptr[i], 0) && !bit_array_get(monomer_bitarray_ptr[3 - i], k - 1)) {
			return palindrome;
		}
		else{
			for(size_t j = 1; j < (k >> 1); ++j){
				if (bit_array_get(monomer_bitarray_ptr[i], j) && !bit_array_get(monomer_bitarray_ptr[3 - i], k - 1 - j)) {
					++mm_c;
				}
			}
		}
	}

	if (mm_c <= max_mm){
		palindrome = true;
	}

	return palindrome;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_check_mirror_symmetry_interval -- check a sub-interval of the BIT_ARRAY encoding the kmer for the mirror symmetry
//
//	parameters:
//		monomer_bitarray_ptr - ptr to a *BIT_ARRAY[4]
//		k - length of the kmer
//		max_mm - max number of mismatch allowed
//		start - starting index of the interval to check
//		end - ending index of the interval to check
//		mm_c - counter for mismatch found
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
bool Nessie::routine_check_mirror_symmetry_interval(BIT_ARRAY **monomer_bitarray_ptr, size_t k, size_t max_mm, size_t start, size_t end, size_t mm_c){

	if (start > end){ throw std::invalid_argument("Mirror symmetry interval: starting index is larger than ending index"); }
	if ((end + 1) > k){ throw std::invalid_argument("Mirror symmetry interval: ending index is larger thank sequence end"); }

	// Variables
	bool mirror = false;

	// Checking simmetry
	for(size_t i = 0; i < 4; ++i){
		if (bit_array_get(monomer_bitarray_ptr[i], start) && !bit_array_get(monomer_bitarray_ptr[i], end)) {
			return mirror;
		}
		else{
			for(size_t j = 1; j < ((end - start + 1) >> 1); ++j){
				if (bit_array_get(monomer_bitarray_ptr[i], start + j) && !bit_array_get(monomer_bitarray_ptr[i], end - j)) {
					++mm_c;
				}
			}
		}
	}

	if (mm_c <= max_mm){
		mirror = true;
	}

	return mirror;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_check_palindrome_symmetry_interval -- check a sub-interval of the BIT_ARRAY encoding the kmer for the palindrome symmetry
//
//	parameters:
//		monomer_bitarray_ptr - ptr to a *BIT_ARRAY[4]
//		k - length of the kmer
//		max_mm - max number of mismatch allowed
//		start - starting index of the interval to check
//		end - ending index of the interval to check
//		mm_c - counter for mismatch found
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
bool Nessie::routine_check_palindrome_symmetry_interval(BIT_ARRAY **monomer_bitarray_ptr, size_t k, size_t max_mm, size_t start, size_t end, size_t mm_c){

	if (start > end){ throw std::invalid_argument("Palindrome symmetry interval: starting index is larger than ending index"); }
	if ((end + 1) > k){ throw std::invalid_argument("Palindrome symmetry interval: ending index is larger thank sequence end"); }

	// Variables
	bool palindrome = false;

	// Checking simmetry
	for(size_t i = 0; i < 4; ++i){
		if (bit_array_get(monomer_bitarray_ptr[i], start) && !bit_array_get(monomer_bitarray_ptr[3 - i], end)) {
			return palindrome;
		}
		else{
			for(size_t j = 1; j < ((end - start + 1) >> 1); ++j){
				if (bit_array_get(monomer_bitarray_ptr[i], start + j) && !bit_array_get(monomer_bitarray_ptr[3 - i], end - j)) {
					++mm_c;
				}
			}
		}
	}

	if (mm_c <= max_mm){
		palindrome = true;
	}

	return palindrome;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_check_global_alignment
//
//	parameters:
//		mask_kmer_ptr - ptr to a mask (uint8_t array) encoding the kmer
//		k - length of the kmer
//		max_mm - max number of mismatch allowed
//		max_gap - max number of gaps allowed
//		type - symmetry checked: 0 = mirror, 1 = palindrome
//
//	note: implementation of a modified Needleman-Wunsch algorithm for global alignment
//
/////////////////////////////////////////////////////////////////////////////////////
bool Nessie::routine_check_global_alignment(uint8_t *mask_kmer_ptr, std::vector<bool> *align_vector_ptr, size_t k, size_t max_mm, size_t max_gap, size_t max_gapmm, int type){

	// Variables
	bool alignment = false;
	long scoring_matrix[k + 1][k + 1];	// scoring matrix
	int m = 1, mm = -1, indls = -1; // m = match, mm = mismatch, indls = gap
	size_t mm_c = 0, gap_c = 0;	// counters for mm and gap
	align_vector_ptr->clear();

	// Other variables
	size_t max_i, max_j;
	long score_diag, score_up, score_left, max;

	// Check first and last bases symmetry
	if ((0 == type) && (mm == Nessie::routine_compare_bases(mask_kmer_ptr, 0, k - 1, m, mm))){
			return alignment;
	}
	else if ((1 == type) && (mm == Nessie::routine_compare_bases_complement(mask_kmer_ptr, 0, k - 1, m, mm))){
			return alignment;
	}

	// Initializing scoring_matrix
	for(size_t i = 0; i <= k; ++i){
		scoring_matrix[i][0] = indls * i;
	}

	for(size_t j = 0; j <= k; ++j){
		scoring_matrix[0][j] = indls * j;
	}

	// Filling matrix
	bool max_defined = false;
	for(size_t i = 1; i <= k; ++i){
		for(size_t j = 1; j <= k; ++j){
			if ((i + j) > k){ break; }
			if (0 == type){ score_diag = scoring_matrix[i - 1][j - 1] + Nessie::routine_compare_bases(mask_kmer_ptr, i - 1, k - j, m, mm); }	// mirror
			else if (1 == type){ score_diag = scoring_matrix[i - 1][j - 1] + Nessie::routine_compare_bases_complement(mask_kmer_ptr, i - 1, k - j, m, mm); }	// palindrome
			else{ throw std::invalid_argument("Alignment: type is not valid"); }
			score_up = scoring_matrix[i - 1][j]	+ indls;
			score_left = scoring_matrix[i][j - 1] + indls;
			scoring_matrix[i][j] = std::max(std::max(score_diag, score_up), score_left);
			if (k == (i + j)){
				if (!max_defined){
					max_defined = true;
					max = scoring_matrix[i][j];
					max_i = i;
					max_j = j;
				}
				else if (max < scoring_matrix[i][j]){
					max = scoring_matrix[i][j];
					max_i = i;
					max_j = j;
				}
			}
		}
	}

	// Check if possible alignment have been found
	long min_score = (max_gapmm) ? (std::max((max_i), (max_j)) - (max_gapmm << 1)) : (std::max((max_i), (max_j)) - ((max_mm + max_gap) << 1));
	if (max < min_score){ return alignment; }

	// Retrieving best alignment
	size_t i = max_i, j = max_j;
	int comp_score;
	while ((i > 0) || (j > 0)){
		if ((i > 0) && (j > 0)){
			comp_score = (0 == type) ? Nessie::routine_compare_bases(mask_kmer_ptr, i - 1, k - j, m, mm)
								 	 : Nessie::routine_compare_bases_complement(mask_kmer_ptr, i - 1, k - j, m, mm);
		}

		if ((i > 0) && (j > 0) && (scoring_matrix[i][j] == (scoring_matrix[i - 1][j - 1] + comp_score))){
			if (m == comp_score) {
				--i;
				--j;
				align_vector_ptr->push_back(true); // M = 01
				align_vector_ptr->push_back(false);
			}
			else{
				++mm_c;
				--i;
				--j;
				align_vector_ptr->push_back(false); // m = 10
				align_vector_ptr->push_back(true);
			}
		}
		else if ((i > 0) && (scoring_matrix[i][j] == (scoring_matrix[i - 1][j] + indls))){
			++gap_c;
			--i;
			align_vector_ptr->push_back(false);	// u = 00
			align_vector_ptr->push_back(false);
		}
		else{
			++gap_c;
			--j;
			align_vector_ptr->push_back(true); // l = 11
			align_vector_ptr->push_back(true);
		}
	}

	// Check gaps and mismatches
	if (max_gapmm){
		if (!max_gap && !max_mm){
			if ((gap_c + mm_c) <= max_gapmm){ alignment = true; }
		}
		else if (!max_gap && max_mm){
			if (((gap_c + mm_c) <= max_gapmm) && (mm_c <= max_mm)){ alignment = true; }
		}
		else if (max_gap && !max_mm){
			if (((gap_c + mm_c) <= max_gapmm) && (gap_c <= max_gap)){ alignment = true; }
		}
		else {
			if (((gap_c + mm_c) <= max_gapmm) && (gap_c <= max_gap) && (mm_c <= max_mm)){ alignment = true; }
		}
	}
	else if ((gap_c <= max_gap) && (mm_c <= max_mm)){ alignment = true; }

	return alignment;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_check_global_alignment_interval
//
//	parameters:
//		mask_kmer_ptr - ptr to a mask (uint8_t array) encoding the kmer
//		k - length of the kmer
//		max_mm - max number of mismatch allowed
//		max_gap - max number of gaps allowed
//		type - symmetry checked: 0 = mirror, 1 = palindrome
//		start - starting index of the interval to check
//		end - ending index of the interval to check
//
//	note: implementation of a modified Needleman-Wunsch algorithm for an interval
//
/////////////////////////////////////////////////////////////////////////////////////
bool Nessie::routine_check_global_alignment_interval(uint8_t *mask_kmer_ptr, std::vector<bool> *align_vector_ptr, size_t k, size_t max_mm, size_t max_gap, size_t max_gapmm, int type, size_t start, size_t end){

	if (start > end){ throw std::invalid_argument("Alignment interval: starting index is larger than ending index"); }
	if ((end + 1) > k){ throw std::invalid_argument("Alignment interval: ending index is larger than sequence end"); }

	// Variables
	bool alignment = false;
	size_t l = end - start + 1;;
	long scoring_matrix[l + 1][l + 1];	// scoring matrix
	int m = 1, mm = -1, indls = -1; // m = match, mm = mismatch, indls = gap
	size_t mm_c = 0, gap_c = 0;	// counters for mm and gap
	align_vector_ptr->clear();

	// Other variables
	size_t max_i, max_j;
	long score_diag, score_up, score_left, max;

	// Check first and last bases symmetry
	if ((0 == type) && (mm == Nessie::routine_compare_bases(mask_kmer_ptr, start, end, m, mm))){
			return alignment;
	}
	else if ((1 == type) && (mm == Nessie::routine_compare_bases_complement(mask_kmer_ptr, start, end, m, mm))){
			return alignment;
	}

	// Initializing scoring_matrix
	for(size_t i = 0; i <= l; ++i){
		scoring_matrix[i][0] = indls * i;
	}

	for(size_t j = 0; j <= l; ++j){
		scoring_matrix[0][j] = indls * j;
	}

	// Filling matrix
	bool max_defined = false;
	for(size_t i = 1; i <= l; ++i){
		for(size_t j = 1; j <= l; ++j){
			if ((i + j) > l){ break; }
			if (0 == type){ score_diag = scoring_matrix[i - 1][j - 1] + Nessie::routine_compare_bases(mask_kmer_ptr, start + i - 1, end + 1 - j, m, mm); }	// mirror
			else if (1 == type){ score_diag = scoring_matrix[i - 1][j - 1] + Nessie::routine_compare_bases_complement(mask_kmer_ptr, start + i - 1, end + 1 - j, m, mm); }	// palindrome
			else{ throw std::invalid_argument("Alignment: type is not valid"); }
			score_up = scoring_matrix[i - 1][j]	+ indls;
			score_left = scoring_matrix[i][j - 1] + indls;
			scoring_matrix[i][j] = std::max(std::max(score_diag, score_up), score_left);
			if (l == (i + j)){
				if (!max_defined){
					max_defined = true;
					max = scoring_matrix[i][j];
					max_i = i;
					max_j = j;
				}
				else if (max < scoring_matrix[i][j]){
					max = scoring_matrix[i][j];
					max_i = i;
					max_j = j;
				}
			}
		}
	}

	// Check if possible alignment have been found
	long min_score = (max_gapmm) ? (std::max((max_i), (max_j)) - (max_gapmm << 1)) : (std::max((max_i), (max_j)) - ((max_mm + max_gap) << 1));
	if (max < min_score){ return alignment; }

	// Retrieving best alignment
	size_t i = max_i, j = max_j;
	int comp_score;
	while ((i > 0) || (j > 0)){
		if ((i > 0) && (j > 0)){
			comp_score = (0 == type) ? Nessie::routine_compare_bases(mask_kmer_ptr, start + i - 1, end + 1 - j, m, mm)
								 	 : Nessie::routine_compare_bases_complement(mask_kmer_ptr, start + i - 1, end + 1 - j, m, mm);
		}

		if ((i > 0) && (j > 0) && (scoring_matrix[i][j] == (scoring_matrix[i - 1][j - 1] + comp_score))){
			if (m == comp_score) {
				--i;
				--j;
				align_vector_ptr->push_back(true); // M = 01
				align_vector_ptr->push_back(false);
			}
			else{
				++mm_c;
				--i;
				--j;
				align_vector_ptr->push_back(false); // m = 10
				align_vector_ptr->push_back(true);
			}
		}
		else if ((i > 0) && (scoring_matrix[i][j] == (scoring_matrix[i - 1][j] + indls))){
			++gap_c;
			--i;
			align_vector_ptr->push_back(false);	// u = 00
			align_vector_ptr->push_back(false);
		}
		else{
			++gap_c;
			--j;
			align_vector_ptr->push_back(true);	// l = 11
			align_vector_ptr->push_back(true);
		}
	}

	// Check gaps and mismatches
	if (max_gapmm){
		if (!max_gap && !max_mm){
			if ((gap_c + mm_c) <= max_gapmm){ alignment = true; }
		}
		else if (!max_gap && max_mm){
			if (((gap_c + mm_c) <= max_gapmm) && (mm_c <= max_mm)){ alignment = true; }
		}
		else if (max_gap && !max_mm){
			if (((gap_c + mm_c) <= max_gapmm) && (gap_c <= max_gap)){ alignment = true; }
		}
		else {
			if (((gap_c + mm_c) <= max_gapmm) && (gap_c <= max_gap) && (mm_c <= max_mm)){ alignment = true; }
		}
	}
	else if ((gap_c <= max_gap) && (mm_c <= max_mm)){ alignment = true; }

	return alignment;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_compare_bases -- compare i-th and j-th bases of mask_kmer_ptr returning score associated to a match or a mismatch
//
//	parameters:
//		mask_kmer_ptr - ptr to a mask (uint8_t array) encoding the kmer
//		i - index i
//		j -	index j
//		m - match score
//		mm - mismatch score
//
/////////////////////////////////////////////////////////////////////////////////////
int Nessie::routine_compare_bases(uint8_t *mask_kmer_ptr, size_t i, size_t j, int m, int mm){

	uint8_t shift_i = (i & ((1 << 2) - 1)) << 1;
	uint8_t shift_j = (j & ((1 << 2) - 1)) << 1;
	uint8_t base_i = (mask_kmer_ptr[i >> 2] & (BASE_MASK << shift_i)) >> shift_i;
	uint8_t base_j = (mask_kmer_ptr[j >> 2] & (BASE_MASK << shift_j)) >> shift_j;
	if (base_i == base_j){
		return m;
	}
	return mm;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_compare_bases_complement -- compare i-th and j-th bases (in complement) of mask_kmer_ptr returning score associated to a match or a mismatch
//
//	parameters:
//		mask_kmer_ptr - ptr to a mask (uint8_t array) encoding the kmer
//		i - index i
//		j -	index j
//		m - match score
//		mm - mismatch score
//
/////////////////////////////////////////////////////////////////////////////////////
int Nessie::routine_compare_bases_complement(uint8_t *mask_kmer_ptr, size_t i, size_t j, int m, int mm){

	uint8_t shift_i = (i & ((1 << 2) - 1)) << 1;
	uint8_t shift_j = (j & ((1 << 2) - 1)) << 1;
	uint8_t base_i = (mask_kmer_ptr[i >> 2] & (BASE_MASK << shift_i)) >> shift_i;
	uint8_t base_j = (mask_kmer_ptr[j >> 2] & (BASE_MASK << shift_j)) >> shift_j;
	if (base_i == (3 - base_j)){
		return m;
	}
	return mm;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_get_kmers_k_mirror -- returns a ptr to an HashTable containing Kmer objects (kmers of length k with a mirror symmetry in the interval)
//
//	parameters:
//		k - length of the kmers searched
//		max_mm - maximum number of mismatch allowed
//		start - starting index of the interval to search
//		end - ending index of the interval to search
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
HashTable *Nessie::routine_get_kmers_k_mirror(size_t k, size_t max_mm, size_t start, size_t end){	//size_t check = 0;

	// Some variables
	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (k > (end - start + 1)){ throw std::invalid_argument("Get mirrors k: k is longer than the sequence interval"); }

	// Initializing HashTable
	HashTable *hash_table_ptr = new HashTable(false);

	// Creating BIT_ARRAY to store monomers indexes for the mask
	BIT_ARRAY *monomer_bitarray[4];
	BIT_ARRAY **monomer_bitarray_ptr = monomer_bitarray;
	for (size_t i = 0; i < 4; ++i){
		monomer_bitarray_ptr[i] = bit_array_create(k);
	}

	// Defining bytes necessary to store the kmer
	size_t dna_bytes = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

	// Defining the mask to encode the kmer
	uint8_t mask_kmer[dna_bytes];
	uint8_t *mask_kmer_ptr = mask_kmer;
	std::memset(mask_kmer, 0, dna_bytes);

	// Initializing bit_array and mask for the first kmer of length k in the interval
	Nessie::routine_init_bitarray_and_mask(monomer_bitarray_ptr, mask_kmer_ptr, k, start);

	// Check mirror simmetry for the first kmer
	if (Nessie::check_mirror_symmetry(monomer_bitarray_ptr, k, max_mm)){
		Kmer *kmer_ptr = new Kmer(k);	// defining ptr to new Kmer object
		copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_ptr, dna_bytes);
		kmer_ptr->indexes.push_back(start);
		kmer_ptr->counts += 1;
		hash_table_ptr->insert_kmer(kmer_ptr);	// adding Kmer to the HashTable
	}

	// Sliding by one base at each iteration to get successive kmers and checking their simmetry
	for (size_t i = (start + 1); i <= (end - k + 1); ++i){
		Nessie::routine_shift_bitarray_and_mask(monomer_bitarray_ptr, mask_kmer_ptr, dna_bytes, k, i);

		if (Nessie::check_mirror_symmetry(monomer_bitarray_ptr, k, max_mm)){		//++check;
			Kmer *kmer_ptr = new Kmer(k);
			copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_ptr, dna_bytes);
			kmer_ptr->indexes.push_back(i);
			kmer_ptr->counts += 1;
			hash_table_ptr->insert_kmer(kmer_ptr);
		}
	}
	//std::cout << "CHECK " << check << std::endl;
	return hash_table_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_get_max_kmer_mirror -- returns a ptr to an HashTable containing Kmer objects (longest kmers of maximum length max_k and minimum length min_k with a mirror symmetry in the interval)
//
//	parameters:
//		max_k - max length of the kmers searched
//		min_k - min length of the kmers searched
//		modulo - max mismatch allowed are calculated as integer division k / modulo
//		start - starting index of the interval to search
//		end - ending index of the interval to search
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
HashTable *Nessie::routine_get_max_kmer_mirror(size_t k_max, size_t k_min, size_t modulo, size_t start, size_t end){	//size_t check = 0;

	// Some variables
	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Get max mirror: max_k is longer than the sequence interval"); }

	// Variables
	bool added;
	size_t added_end = 0;

	// Initializing HashTable
	HashTable *hash_table_ptr = new HashTable(false);

	// Creating BIT_ARRAY to store monomers indexes for the mask
	BIT_ARRAY *monomer_bitarray[4];
	BIT_ARRAY **monomer_bitarray_ptr = monomer_bitarray;
	for (size_t i = 0; i < 4; ++i){
		monomer_bitarray_ptr[i] = bit_array_create(k_max);
	}

	// Defining bytes necessary to store the kmer
	size_t dna_bytes = (k_max >> 2) + (0 != (k_max & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

	// Defining the mask to encode the kmer
	uint8_t mask_kmer[dna_bytes];
	uint8_t *mask_kmer_ptr = mask_kmer;
	std::memset(mask_kmer, 0, dna_bytes);

	// Initializing mask for the first kmer of length max_k in the interval
	Nessie::routine_init_bitarray_and_mask(monomer_bitarray_ptr, mask_kmer_ptr, k_max, start);

	size_t k = k_max;
	size_t end_i = start + k - 1;
	while (k >= k_min){
		size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
		if (Nessie::routine_check_mirror_symmetry_interval(monomer_bitarray_ptr, k_max, max_mm, 0, k - 1)){	//++check;
			Kmer *kmer_ptr = new Kmer(k);
			added_end = end_i;

			// Defining bytes necessary to store the kmer k
			size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

			// Defining the mask to encode the kmer k
			uint8_t mask_kmer_k[dna_bytes_k];
			uint8_t *mask_kmer_k_ptr = mask_kmer_k;
			std::memset(mask_kmer_k, 0, dna_bytes_k);

			// Initializing mask for the first kmer of length k in the interval
			Nessie::routine_init_mask(mask_kmer_k_ptr, k, start);

			copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
			kmer_ptr->indexes.push_back(start);
			kmer_ptr->counts += 1;
			hash_table_ptr->insert_kmer_var_len(kmer_ptr);
			break;
		}
	--k;
	--end_i;
	}

	// Sliding by one base at each iteration to get successive kmers and checking their simmetry
	size_t last_i = start + 1;
	for (size_t i = (start + 1); i <= (end - k_max + 1); ++i){
		Nessie::routine_shift_bitarray_and_mask(monomer_bitarray_ptr, mask_kmer_ptr, dna_bytes, k_max, i);

		size_t k = k_max;
		size_t end_i = i + k - 1;
		while ((k >= k_min) && (end_i > added_end)){
			size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
			if (Nessie::routine_check_mirror_symmetry_interval(monomer_bitarray_ptr, k_max, max_mm, 0, k - 1)){	//++check;
				Kmer *kmer_ptr = new Kmer(k);
				added_end = end_i;

				// Defining bytes necessary to store the kmer k
				size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

				// Defining the mask to encode the kmer k
				uint8_t mask_kmer_k[dna_bytes_k];
				uint8_t *mask_kmer_k_ptr = mask_kmer_k;
				std::memset(mask_kmer_k, 0, dna_bytes_k);

				// Initializing mask for the first kmer of length k in the interval
				Nessie::routine_init_mask(mask_kmer_k_ptr, k, i);

				copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
				kmer_ptr->indexes.push_back(i);
				kmer_ptr->counts += 1;
				hash_table_ptr->insert_kmer_var_len(kmer_ptr);
				break;
			}
		--k;
		--end_i;
		}
		last_i = i;
	}

	// Check last interval that may be shorter than k_max and skipped in previous for loop
	size_t last_interval_length = end - last_i;
	if (last_interval_length >= k_min){

		// Defining bytes necessary to store the last interval kmer
		size_t dna_bytes_l = (last_interval_length >> 2) + (0 != (last_interval_length & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

		// Creating BIT_ARRAY to store monomers indexes for the mask
		BIT_ARRAY *monomer_bitarray_l[4];
		BIT_ARRAY **monomer_bitarray_l_ptr = monomer_bitarray_l;
		for (size_t i = 0; i < 4; ++i){
			monomer_bitarray_l_ptr[i] = bit_array_create(last_interval_length);
		}

		// Defining the mask to encode the last interval kmer
		uint8_t mask_kmer_l[dna_bytes_l];
		uint8_t *mask_kmer_l_ptr = mask_kmer_l;
		std::memset(mask_kmer_l, 0, dna_bytes_l);

		// Initializing mask for the last interval kmer
		Nessie::routine_init_bitarray_and_mask(monomer_bitarray_l_ptr, mask_kmer_l_ptr, last_interval_length, last_i + 1);

		for (size_t i = 0; i <= (last_interval_length - k_min); ++i){

			size_t k = last_interval_length - i;
			size_t end_i = last_interval_length - 1; //std::cout << i << " " << k << std::endl;
			while ((k >= k_min) && ((last_i + 1 + i + k - 1) > added_end)){ //std::cout << "- " << i << " " << end_i << " " << k << std::endl;
				size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
				if (Nessie::routine_check_mirror_symmetry_interval(monomer_bitarray_l_ptr, last_interval_length, max_mm, i, end_i)){	//++check;
					Kmer *kmer_ptr = new Kmer(k);
					added_end = last_i + 1 + i + k - 1;

					// Defining bytes necessary to store the kmer k
					size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

					// Defining the mask to encode the kmer k
					uint8_t mask_kmer_k[dna_bytes_k];
					uint8_t *mask_kmer_k_ptr = mask_kmer_k;
					std::memset(mask_kmer_k, 0, dna_bytes_k);

					// Initializing mask for the first kmer of length k in the interval
					Nessie::routine_init_mask(mask_kmer_k_ptr, k, last_i + 1 + i);

					copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
					kmer_ptr->indexes.push_back(last_i + 1 + i);
					kmer_ptr->counts += 1;
					hash_table_ptr->insert_kmer_var_len(kmer_ptr);
					break;
				}
			--k;
			--end_i;
			}
		}
	}

	//std::cout << "CHECK " << check << std::endl;
	return hash_table_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_get_kmers_k_gap -- returns a ptr to an HashTable containing Kmer objects (kmers of length k with a specific symmetry in the interval),
//							   allows gap presence
//
//	parameters:
//		k - length of the kmers searched
//		max_mm - maximum number of mismatch allowed
//		max_gap - max number of gaps allowed
//		start - starting index of the interval to search
//		end - ending index of the interval to search
//		type - symmetry checked: 0 = mirror, 1 = palindrome
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
HashTable *Nessie::routine_get_kmers_k_gap(size_t k, size_t max_mm, size_t max_gap, size_t max_gapmm, size_t start, size_t end, int type){	//size_t check = 0;

	// Some variables
	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (k > (end - start + 1)){ throw std::invalid_argument("Get kmers k gap: k is longer than the sequence interval"); }

	// Initializing HashTable
	HashTable *hash_table_ptr = new HashTable(false);

	// Defining bytes necessary to store the kmer
	size_t dna_bytes = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

	// Defining the mask to encode the kmer
	uint8_t mask_kmer[dna_bytes];
	uint8_t *mask_kmer_ptr = mask_kmer;
	std::memset(mask_kmer, 0, dna_bytes);

	// Initializing mask for the first kmer of length k in the interval
	Nessie::routine_init_mask(mask_kmer_ptr, k, start);
	std::vector<bool> *align_vector_ptr = new std::vector<bool>;

	// Check mirror simmetry for the first kmer
	if (Nessie::routine_check_global_alignment(mask_kmer_ptr, align_vector_ptr, k, max_mm, max_gap, max_gapmm, type)){
		Kmer *kmer_ptr = new Kmer(k);	// defining ptr to new Kmer object
		copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_ptr, dna_bytes);
		kmer_ptr->indexes.push_back(start);
		kmer_ptr->counts += 1;
		kmer_ptr->alignment_ptr = align_vector_ptr;
		hash_table_ptr->insert_kmer(kmer_ptr);	// adding Kmer to the HashTable
	}
	else{
		delete align_vector_ptr;
	}

	// Sliding by one base at each iteration to get successive kmers and checking their simmetry
	for (size_t i = (start + 1); i <= (end - k + 1); ++i){
		Nessie::routine_shift_mask(mask_kmer_ptr, dna_bytes, k, i);
		std::vector<bool> *align_vector_ptr_i = new std::vector<bool>;

		if (Nessie::routine_check_global_alignment(mask_kmer_ptr, align_vector_ptr_i, k, max_mm, max_gap, max_gapmm, type)){		//++check;
			Kmer *kmer_ptr = new Kmer(k);
			copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_ptr, dna_bytes);
			kmer_ptr->indexes.push_back(i);
			kmer_ptr->counts += 1;
			kmer_ptr->alignment_ptr = align_vector_ptr_i;
			hash_table_ptr->insert_kmer(kmer_ptr);
		}
		else{
			delete align_vector_ptr_i;
		}
	}
	//std::cout << "CHECK " << check << std::endl;
	return hash_table_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_get_max_kmer_gap -- returns a ptr to an HashTable containing Kmer objects (longest kmers of maximum length max_k and minimum length min_k with a specific symmetry in the interval),
//							    allows gap presence
//
//	parameters:
//		max_k - max length of the kmers searched
//		min_k - min length of the kmers searched
//		modulo - max mismatch allowed are calculated as integer division k / modulo
//		modulo_gap - max gaps allowed are calculated as integer division k / modulo_gap
//		start - starting index of the interval to search
//		end - ending index of the interval to search
//		type - symmetry checked: 0 = mirror, 1 = palindrome
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
HashTable *Nessie::routine_get_max_kmer_gap(size_t k_max, size_t k_min, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t start, size_t end, int type){	//size_t check = 0;

	// Some variables
	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Get max kmer gap: max_k is longer than the sequence interval"); }

	// Variables
	bool added;
	size_t added_end = 0;

	// Initializing HashTable
	HashTable *hash_table_ptr = new HashTable(false);

	// Defining bytes necessary to store the kmer
	size_t dna_bytes = (k_max >> 2) + (0 != (k_max & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

	// Defining the mask to encode the kmer
	uint8_t mask_kmer[dna_bytes];
	uint8_t *mask_kmer_ptr = mask_kmer;
	std::memset(mask_kmer, 0, dna_bytes);

	// Initializing mask for the first kmer of length max_k in the interval
	Nessie::routine_init_mask(mask_kmer_ptr, k_max, start);
	std::vector<bool> *align_vector_ptr = new std::vector<bool>;
	added = false;

	size_t k = k_max;
	size_t end_i = start + k - 1;
	while (k >= k_min){
		size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
		size_t max_gap = (modulo_gap) ? (k * modulo_gap) / 100 : 0;
		size_t max_gapmm = (modulo_gapmm) ? (k * modulo_gapmm) / 100 : 0;
		if (Nessie::routine_check_global_alignment_interval(mask_kmer_ptr, align_vector_ptr, k_max, max_mm, max_gap, max_gapmm, type, 0, k - 1)){	//++check;
			Kmer *kmer_ptr = new Kmer(k);
			added_end = end_i;

			// Defining bytes necessary to store the kmer k
			size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

			// Defining the mask to encode the kmer k
			uint8_t mask_kmer_k[dna_bytes_k];
			uint8_t *mask_kmer_k_ptr = mask_kmer_k;
			std::memset(mask_kmer_k, 0, dna_bytes_k);

			// Initializing mask for the first kmer of length k in the interval
			Nessie::routine_init_mask(mask_kmer_k_ptr, k, start);

			copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
			kmer_ptr->indexes.push_back(start);
			kmer_ptr->counts += 1;
			kmer_ptr->alignment_ptr = align_vector_ptr;
			added = true;
			hash_table_ptr->insert_kmer_var_len(kmer_ptr);
			break;
		}
	--k;
	--end_i;
	}

	if (!added){
		delete align_vector_ptr;
	}

	// Sliding by one base at each iteration to get successive kmers and checking their simmetry
	size_t last_i = start + 1;
	for (size_t i = (start + 1); i <= (end - k_max + 1); ++i){
		Nessie::routine_shift_mask(mask_kmer_ptr, dna_bytes, k_max, i);
		std::vector<bool> *align_vector_ptr_i = new std::vector<bool>;
		added = false;

		size_t k = k_max;
		size_t end_i = i + k - 1;
		while ((k >= k_min) && (end_i > added_end)){
			size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
			size_t max_gap = (modulo_gap) ? (k * modulo_gap) / 100 : 0;
			size_t max_gapmm = (modulo_gapmm) ? (k * modulo_gapmm) / 100 : 0;
			if (Nessie::routine_check_global_alignment_interval(mask_kmer_ptr, align_vector_ptr_i, k_max, max_mm, max_gap, max_gapmm, type, 0, k - 1)){	//++check;
				Kmer *kmer_ptr = new Kmer(k);
				added_end = end_i;

				// Defining bytes necessary to store the kmer k
				size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

				// Defining the mask to encode the kmer k
				uint8_t mask_kmer_k[dna_bytes_k];
				uint8_t *mask_kmer_k_ptr = mask_kmer_k;
				std::memset(mask_kmer_k, 0, dna_bytes_k);

				// Initializing mask for the first kmer of length k in the interval
				Nessie::routine_init_mask(mask_kmer_k_ptr, k, i);

				copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
				kmer_ptr->indexes.push_back(i);
				kmer_ptr->counts += 1;
				kmer_ptr->alignment_ptr = align_vector_ptr_i;
				added = true;
				hash_table_ptr->insert_kmer_var_len(kmer_ptr);
				break;
			}
		--k;
		--end_i;
		}
		last_i = i;

		if (!added){
			delete align_vector_ptr_i;
		}
	}

	// Check last interval that may be shorter than k_max and skipped in previous for loop
	size_t last_interval_length = end - last_i;
	if (last_interval_length >= k_min){

		// Defining bytes necessary to store the last interval kmer
		size_t dna_bytes_l = (last_interval_length >> 2) + (0 != (last_interval_length & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

		// Defining the mask to encode the last interval kmer
		uint8_t mask_kmer_l[dna_bytes_l];
		uint8_t *mask_kmer_l_ptr = mask_kmer_l;
		std::memset(mask_kmer_l, 0, dna_bytes_l);

		// Initializing mask for the last interval kmer
		Nessie::routine_init_mask(mask_kmer_l_ptr, last_interval_length, last_i + 1);

		for (size_t i = 0; i <= (last_interval_length - k_min); ++i){
			std::vector<bool> *align_vector_ptr_l = new std::vector<bool>;
			added = false;

			size_t k = last_interval_length - i;
			size_t end_i = last_interval_length - 1; //std::cout << i << " " << k << std::endl;
			while ((k >= k_min) && ((last_i + 1 + i + k - 1) > added_end)){ //std::cout << "- " << i << " " << end_i << " " << k << std::endl;
				size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
				size_t max_gap = (modulo_gap) ? (k * modulo_gap) / 100 : 0;
				size_t max_gapmm = (modulo_gapmm) ? (k * modulo_gapmm) / 100 : 0;
				if (Nessie::routine_check_global_alignment_interval(mask_kmer_l_ptr, align_vector_ptr_l, last_interval_length, max_mm, max_gap, max_gapmm, type, i, end_i)){	//++check;
					Kmer *kmer_ptr = new Kmer(k);
					added_end = last_i + 1 + i + k - 1;

					// Defining bytes necessary to store the kmer k
					size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

					// Defining the mask to encode the kmer k
					uint8_t mask_kmer_k[dna_bytes_k];
					uint8_t *mask_kmer_k_ptr = mask_kmer_k;
					std::memset(mask_kmer_k, 0, dna_bytes_k);

					// Initializing mask for the first kmer of length k in the interval
					Nessie::routine_init_mask(mask_kmer_k_ptr, k, last_i + 1 + i);

					copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
					kmer_ptr->indexes.push_back(last_i + 1 + i);
					kmer_ptr->counts += 1;
					kmer_ptr->alignment_ptr = align_vector_ptr_l;
					added = true;
					hash_table_ptr->insert_kmer_var_len(kmer_ptr);
					break;
				}
			--k;
			--end_i;
			}

			if (!added){
				delete align_vector_ptr_l;
			}
		}
	}

	//std::cout << "CHECK " << check << std::endl;
	return hash_table_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_get_kmers_k_palindrome -- returns a ptr to an HashTable containing Kmer objects (kmers of length k with a palindrome symmetry in the interval)
//
//	parameters:
//		k - length of the kmers searched
//		max_mm - maximum number of mismatch allowed
//		start - starting index of the interval to search
//		end - ending index of the interval to search
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
HashTable *Nessie::routine_get_kmers_k_palindrome(size_t k, size_t max_mm, size_t start, size_t end){	//size_t check = 0;

	// Some variables
	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (k > (end - start + 1)){ throw std::invalid_argument("Get palindromes k: k is longer than the sequence interval"); }

	// Initializing HashTable
	HashTable *hash_table_ptr = new HashTable(false);

	// Creating BIT_ARRAY to store monomers indexes for the mask
	BIT_ARRAY *monomer_bitarray[4];
	BIT_ARRAY **monomer_bitarray_ptr = monomer_bitarray;
	for (size_t i = 0; i < 4; ++i){
		monomer_bitarray_ptr[i] = bit_array_create(k);
	}

	// Defining bytes necessary to store the kmer
	size_t dna_bytes = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

	// Defining the mask to encode the kmer
	uint8_t mask_kmer[dna_bytes];
	uint8_t *mask_kmer_ptr = mask_kmer;
	std::memset(mask_kmer, 0, dna_bytes);

	// Initializing bit_array and mask for the first kmer of length k in the interval
	Nessie::routine_init_bitarray_and_mask(monomer_bitarray_ptr, mask_kmer_ptr, k, start);

	// Check palindrome simmetry for the first kmer
	if (Nessie::check_palindrome_symmetry(monomer_bitarray_ptr, k, max_mm)){
		Kmer *kmer_ptr = new Kmer(k);	// defining ptr to new Kmer object
		copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_ptr, dna_bytes);
		kmer_ptr->indexes.push_back(start);
		kmer_ptr->counts += 1;
		hash_table_ptr->insert_kmer(kmer_ptr);	// adding Kmer to the HashTable
	}

	// Sliding by one base at each iteration to get successive kmers and checking their simmetry
	for (size_t i = (start + 1); i <= (end - k + 1); ++i){
		Nessie::routine_shift_bitarray_and_mask(monomer_bitarray_ptr, mask_kmer_ptr, dna_bytes, k, i);

		if (Nessie::check_palindrome_symmetry(monomer_bitarray_ptr, k, max_mm)){	//++check;
			Kmer *kmer_ptr = new Kmer(k);
			copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_ptr, dna_bytes);
			kmer_ptr->indexes.push_back(i);
			kmer_ptr->counts += 1;
			hash_table_ptr->insert_kmer(kmer_ptr);
		}
	}
	//std::cout << "CHECK " << check << std::endl;
	return hash_table_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_get_max_kmer_palindrome -- returns a ptr to an HashTable containing Kmer objects (longest kmers of maximum length max_k and minimum length min_k with a palindrome symmetry in the interval)
//
//	parameters:
//		max_k - max length of the kmers searched
//		min_k - min length of the kmers searched
//		modulo - max mismatch allowed are calculated as integer division k / modulo
//		start - starting index of the interval to search
//		end - ending index of the interval to search
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
HashTable *Nessie::routine_get_max_kmer_palindrome(size_t k_max, size_t k_min, size_t modulo, size_t start, size_t end){	//size_t check = 0;

	// Some variables
	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Get max palindrome: max_k is longer than the sequence interval"); }

	// Variables
	bool added;
	size_t added_end = 0;

	// Initializing HashTable
	HashTable *hash_table_ptr = new HashTable(false);

	// Creating BIT_ARRAY to store monomers indexes for the mask
	BIT_ARRAY *monomer_bitarray[4];
	BIT_ARRAY **monomer_bitarray_ptr = monomer_bitarray;
	for (size_t i = 0; i < 4; ++i){
		monomer_bitarray_ptr[i] = bit_array_create(k_max);
	}

	// Defining bytes necessary to store the kmer
	size_t dna_bytes = (k_max >> 2) + (0 != (k_max & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

	// Defining the mask to encode the kmer
	uint8_t mask_kmer[dna_bytes];
	uint8_t *mask_kmer_ptr = mask_kmer;
	std::memset(mask_kmer, 0, dna_bytes);

	// Initializing mask for the first kmer of length max_k in the interval
	Nessie::routine_init_bitarray_and_mask(monomer_bitarray_ptr, mask_kmer_ptr, k_max, start);

	size_t k = k_max;
	size_t end_i = start + k - 1;
	while (k >= k_min){
		size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
		if (Nessie::routine_check_palindrome_symmetry_interval(monomer_bitarray_ptr, k_max, max_mm, 0, k - 1)){	//++check;
			Kmer *kmer_ptr = new Kmer(k);
			added_end = end_i;

			// Defining bytes necessary to store the kmer k
			size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

			// Defining the mask to encode the kmer k
			uint8_t mask_kmer_k[dna_bytes_k];
			uint8_t *mask_kmer_k_ptr = mask_kmer_k;
			std::memset(mask_kmer_k, 0, dna_bytes_k);

			// Initializing mask for the first kmer of length k in the interval
			Nessie::routine_init_mask(mask_kmer_k_ptr, k, start);

			copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
			kmer_ptr->indexes.push_back(start);
			kmer_ptr->counts += 1;
			hash_table_ptr->insert_kmer_var_len(kmer_ptr);
			break;
		}
	--k;
	--end_i;
	}

	// Sliding by one base at each iteration to get successive kmers and checking their simmetry
	size_t last_i = start + 1;
	for (size_t i = (start + 1); i <= (end - k_max + 1); ++i){
		Nessie::routine_shift_bitarray_and_mask(monomer_bitarray_ptr, mask_kmer_ptr, dna_bytes, k_max, i);

		size_t k = k_max;
		size_t end_i = i + k - 1;
		while ((k >= k_min) && (end_i > added_end)){
			size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
			if (Nessie::routine_check_palindrome_symmetry_interval(monomer_bitarray_ptr, k_max, max_mm, 0, k - 1)){	//++check;
				Kmer *kmer_ptr = new Kmer(k);
				added_end = end_i;

				// Defining bytes necessary to store the kmer k
				size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

				// Defining the mask to encode the kmer k
				uint8_t mask_kmer_k[dna_bytes_k];
				uint8_t *mask_kmer_k_ptr = mask_kmer_k;
				std::memset(mask_kmer_k, 0, dna_bytes_k);

				// Initializing mask for the first kmer of length k in the interval
				Nessie::routine_init_mask(mask_kmer_k_ptr, k, i);

				copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
				kmer_ptr->indexes.push_back(i);
				kmer_ptr->counts += 1;
				hash_table_ptr->insert_kmer_var_len(kmer_ptr);
				break;
			}
		--k;
		--end_i;
		}
		last_i = i;
	}

	// Check last interval that may be shorter than k_max and skipped in previous for loop
	size_t last_interval_length = end - last_i;
	if (last_interval_length >= k_min){

		// Defining bytes necessary to store the last interval kmer
		size_t dna_bytes_l = (last_interval_length >> 2) + (0 != (last_interval_length & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

		// Creating BIT_ARRAY to store monomers indexes for the mask
		BIT_ARRAY *monomer_bitarray_l[4];
		BIT_ARRAY **monomer_bitarray_l_ptr = monomer_bitarray_l;
		for (size_t i = 0; i < 4; ++i){
			monomer_bitarray_l_ptr[i] = bit_array_create(last_interval_length);
		}

		// Defining the mask to encode the last interval kmer
		uint8_t mask_kmer_l[dna_bytes_l];
		uint8_t *mask_kmer_l_ptr = mask_kmer_l;
		std::memset(mask_kmer_l, 0, dna_bytes_l);

		// Initializing mask for the last interval kmer
		Nessie::routine_init_bitarray_and_mask(monomer_bitarray_l_ptr, mask_kmer_l_ptr, last_interval_length, last_i + 1);

		for (size_t i = 0; i <= (last_interval_length - k_min); ++i){

			size_t k = last_interval_length - i;
			size_t end_i = last_interval_length - 1; //std::cout << i << " " << k << std::endl;
			while ((k >= k_min) && ((last_i + 1 + i + k - 1) > added_end)){ //std::cout << "- " << i << " " << end_i << " " << k << std::endl;
				size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
				if (Nessie::routine_check_palindrome_symmetry_interval(monomer_bitarray_l_ptr, last_interval_length, max_mm, i, end_i)){ //++check;
					Kmer *kmer_ptr = new Kmer(k);
					added_end = last_i + 1 + i + k - 1;

					// Defining bytes necessary to store the kmer k
					size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

					// Defining the mask to encode the kmer k
					uint8_t mask_kmer_k[dna_bytes_k];
					uint8_t *mask_kmer_k_ptr = mask_kmer_k;
					std::memset(mask_kmer_k, 0, dna_bytes_k);

					// Initializing mask for the first kmer of length k in the interval
					Nessie::routine_init_mask(mask_kmer_k_ptr, k, last_i + 1 + i);

					copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
					kmer_ptr->indexes.push_back(last_i + 1 + i);
					kmer_ptr->counts += 1;
					hash_table_ptr->insert_kmer_var_len(kmer_ptr);
					break;
				}
			--k;
			--end_i;
			}
		}
	}

	//std::cout << "CHECK " << check << std::endl;
	return hash_table_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_init_bitarray_and_mask -- initializes the BIT_ARRAY and the mask for the first kmer of length k in the interval
//
//	parameters:
//		monomer_bitarray_ptr - ptr to a *BIT_ARRAY[4]
//		mask_kmer_ptr - ptr to a uint8_t array that encodes a kmer
//		k - length of the kmer
//		start - starting index of the interval (starting index of the first kmer to build the mask for)
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::routine_init_bitarray_and_mask(BIT_ARRAY **monomer_bitarray_ptr, uint8_t *mask_kmer_ptr, size_t k, size_t start){

	size_t c = 0;	// counter for the mask_kmer creation
	for (size_t i = start; i < (start + k); ++i){
		uint8_t shift_DNA = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
		uint8_t shift_mask = (c & ((1 << 2) - 1)) << 1;	// shift to move from one base to the next one in the mask_kmer (uint8_t array)
		uint8_t base = (string_bit_ptr->data_ptr[i >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING
		mask_kmer_ptr[c >> 2] |= (base << shift_mask);	// insert base in the mask_kmer
		bit_array_set(monomer_bitarray_ptr[base], c);	// set bit corresponding to the base index in the BIT_ARRAY
		++c;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_init_mask -- initializes the mask for the first kmer of length k in the interval
//
//	parameters:
//		mask_kmer_ptr - ptr to a uint8_t array that encodes a kmer
//		k - length of the kmer
//		start - starting index of the interval (starting index of the first kmer to build the mask for)
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::routine_init_mask(uint8_t *mask_kmer_ptr, size_t k, size_t start){

	size_t c = 0;	// counter for the mask_kmer creation
	for (size_t i = start; i < (start + k); ++i){
		uint8_t shift_DNA = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
		uint8_t shift_mask = (c & ((1 << 2) - 1)) << 1;	// shift to move from one base to the next one in the kmer_mask (uint8_t array)
		uint8_t base = (string_bit_ptr->data_ptr[i >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING
		mask_kmer_ptr[c >> 2] |= (base << shift_mask);	// insert base in the mask_kmer
		++c;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_shift_bitarray_and_mask -- shifts the BIT_ARRAY and the mask by one base
//
//	parameters:
//		monomer_bitarray_ptr - ptr to a *BIT_ARRAY[4]
//		mask_kmer_ptr - ptr to a uint8_t array that encodes a kmer
//		mask_kmer_len - length of the uint8_t array that encodes the kmer
//		k - length of the kmer searched
//		i - index from the loop step
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::routine_shift_bitarray_and_mask(BIT_ARRAY **monomer_bitarray_ptr, uint8_t *mask_kmer_ptr, size_t mask_kmer_len, size_t k, size_t i){

	// Variables
	size_t last_index = k - 1;
	uint8_t shift_mask = (last_index & ((1 << 2) - 1)) << 1;

	// Shift
	uint8_t shift_DNA = ((i + last_index) & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
	uint8_t base = (string_bit_ptr->data_ptr[(i + last_index) >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING
	shift_2_right(mask_kmer_ptr, mask_kmer_len);	// shift the uint8_t array encoding the kmer (removes the first base)
	mask_kmer_ptr[last_index >> 2] |= base << shift_mask;	// adding the new base to mask_kmer
	for(size_t i = 0; i < 4; ++i){	//shift monomer_bitarray by one
		if (i == base) {
			bit_array_shift_right(monomer_bitarray_ptr[i], 1, 1);	// set bit in the monomer_bitarray corresponding to the new base index
		}
		else{
			bit_array_shift_right(monomer_bitarray_ptr[i], 1, 0);
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_shift_mask -- shifts the mask by one base
//
//	parameters:
//		mask_kmer_ptr - ptr to a uint8_t array that encodes a kmer
//		mask_kmer_len - length of the uint8_t array that encodes the kmer
//		k - length of the kmer searched
//		i - index from the loop step
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::routine_shift_mask(uint8_t *mask_kmer_ptr, size_t mask_kmer_len, size_t k, size_t i){

	// Variables
	size_t last_index = k - 1;
	uint8_t shift_mask = (last_index & ((1 << 2) - 1)) << 1;

	// Shift
	uint8_t shift_DNA = ((i + last_index) & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
	uint8_t base = (string_bit_ptr->data_ptr[(i + last_index) >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING
	shift_2_right(mask_kmer_ptr, mask_kmer_len);	// shift the uint8_t array encoding the kmer (removes the first base)
	mask_kmer_ptr[last_index >> 2] |= base << shift_mask;	// adding the new base to mask_kmer
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	get_kmers_mirror -- returns a ptr to a LinkedlistKmer that stores all the Kmers with mirror symmetry of length [k_min..k_max] in the interval
//
//	parameters:
//		k_min - minimum length of the kmers searched
//		k_max - maximum length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
LinkedlistKmer *Nessie::get_kmers_mirror(size_t k_min, size_t k_max, size_t modulo, size_t start, size_t end){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Get mirrors: starting index is larger than ending index"); }
	if (!k_max){ k_max = k_min; }	// if k_max is not defined it is set to default as k_min, only kmers of length k_min are searched
	if (k_min > k_max){ throw std::invalid_argument("Get mirrors: k_min is larger than k_max"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Get mirrors: k_max is longer than the sequence interval"); }

	// Initializing LinkedlistKmer
	LinkedlistKmer *ll_kmer_ptr = new LinkedlistKmer;

	// Check range k_min..k_max
	for (size_t i = k_min; i <= k_max; ++i){
		size_t max_mm;
		if (!modulo){
			max_mm = 0;
		}
		else{
			max_mm = (i * modulo) / 100;
		}
		HashTable *hash_table_ptr = Nessie::routine_get_kmers_k_mirror(i, max_mm, start, end);
		hash_table_ptr->append_to_LinkedlistKmer(ll_kmer_ptr);
	}

	return ll_kmer_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	get_kmers_mirror_gap -- returns a ptr to a LinkedlistKmer that stores all the Kmers with mirror symmetry of length [k_min..k_max] in the interval,
//							allows for gaps
//
//	parameters:
//		k_min - minimum length of the kmers searched
//		k_max - maximum length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		modulo_gap - max gaps allowed are calculated as integer division k / modulo_gap [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
LinkedlistKmer *Nessie::get_kmers_mirror_gap(size_t k_min, size_t k_max, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t start, size_t end){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Get mirrors: starting index is larger than ending index"); }
	if (!k_max){ k_max = k_min; }	// if k_max is not defined it is set to default as k_min, only kmers of length k_min are searched
	if (k_min > k_max){ throw std::invalid_argument("Get mirrors: k_min is larger than k_max"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Get mirrors: k_max is longer than the sequence interval"); }

	// Initializing LinkedlistKmer
	LinkedlistKmer *ll_kmer_ptr = new LinkedlistKmer;
	HashTable *hash_table_ptr;

	// Check range k_min..k_max
	for (size_t i = k_min; i <= k_max; ++i){
		size_t max_mm = (modulo) ? (i * modulo) / 100 : 0;
		size_t max_gap = (modulo_gap) ? (i * modulo_gap) / 100 : 0;
		size_t max_gapmm = (modulo_gapmm) ? (i * modulo_gapmm) / 100 : 0;

		// Check if gaps allowed or not
		if (max_gap || max_gapmm){
			hash_table_ptr = Nessie::routine_get_kmers_k_gap(i, max_mm, max_gap, max_gapmm, start, end, 0);
		}
		else{
			hash_table_ptr = Nessie::routine_get_kmers_k_mirror(i, max_mm, start, end);
		}

		// Add to LinkedlistKmer
		hash_table_ptr->append_to_LinkedlistKmer(ll_kmer_ptr);
	}

	return ll_kmer_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	get_max_kmers_mirror_gap -- returns a ptr to a LinkedlistKmer that stores all the maximum Kmers with mirror symmetry (maximum length max_k and minimum length min_k) in the interval,
//								allows for gaps
//
//	parameters:
//		max_k - max length of the kmers searched
//		min_k - min length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		modulo_gap - max gaps allowed are calculated as integer division k / modulo_gap [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
LinkedlistKmer *Nessie::get_max_kmers_mirror_gap(size_t k_max, size_t k_min, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t start, size_t end){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Get max mirrors: starting index is larger than ending index"); }
	if (!k_min){ k_min = k_max; }	// if min_k is not defined it is set to default as max_k, only kmers of length max_k are searched
	if (k_min > k_max){ throw std::invalid_argument("Get max mirrors: min_k is larger than max_k"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Get max mirrors: max_k is longer than the sequence interval"); }

	// Initializing LinkedlistKmer
	LinkedlistKmer *ll_kmer_ptr = new LinkedlistKmer;
	HashTable *hash_table_ptr;

	// Check if gaps allowed or not
	if (modulo_gap || modulo_gapmm){
		hash_table_ptr = Nessie::routine_get_max_kmer_gap(k_max, k_min, modulo, modulo_gap, modulo_gapmm, start, end, 0);
	}
	else{
		hash_table_ptr = Nessie::routine_get_max_kmer_mirror(k_max, k_min, modulo, start, end);
	}

	// Add to LinkedlistKmer
	hash_table_ptr->append_to_LinkedlistKmer(ll_kmer_ptr);

	return ll_kmer_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	get_kmers_palindrome -- returns a ptr to a LinkedlistKmer that stores all the Kmers with palindrome symmetry of length [k_min..k_max] in the interval
//
//	parameters:
//		k_min - minimum length of the kmers searched
//		k_max - maximum length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
LinkedlistKmer *Nessie::get_kmers_palindrome(size_t k_min, size_t k_max, size_t modulo, size_t start, size_t end){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Get palindromes: starting index is larger than ending index"); }
	if (!k_max){ k_max = k_min; }	// if k_max is not defined it is set to default as k_min, only kmers of length k_min are searched
	if (k_min > k_max){ throw std::invalid_argument("Get palindromes: k_min is larger than k_max"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Get palindromes: k_max is longer than the sequence interval"); }

	// Initializing LinkedlistKmer
	LinkedlistKmer *ll_kmer_ptr = new LinkedlistKmer;

	// Check range k_min..k_max
	for (size_t i = k_min; i <= k_max; ++i){
		size_t max_mm;
		if (!modulo){
			max_mm = 0;
		}
		else{
			max_mm = (i * modulo) / 100;
		}
		HashTable *hash_table_ptr = Nessie::routine_get_kmers_k_palindrome(i, max_mm, start, end);
		hash_table_ptr->append_to_LinkedlistKmer(ll_kmer_ptr);
	}

	return ll_kmer_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	get_kmers_palindrome_gap -- returns a ptr to a LinkedlistKmer that stores all the Kmers with palindrome symmetry of length [k_min..k_max] in the interval,
//								allows for gaps
//
//	parameters:
//		k_min - minimum length of the kmers searched
//		k_max - maximum length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		modulo_gap - max gaps allowed are calculated as integer division k / modulo_gap [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
LinkedlistKmer *Nessie::get_kmers_palindrome_gap(size_t k_min, size_t k_max, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t start, size_t end){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Get palindromes: starting index is larger than ending index"); }
	if (!k_max){ k_max = k_min; }	// if k_max is not defined it is set to default as k_min, only kmers of length k_min are searched
	if (k_min > k_max){ throw std::invalid_argument("Get palindromes: k_min is larger than k_max"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Get palindromes: k_max is longer than the sequence interval"); }

	// Initializing LinkedlistKmer
	LinkedlistKmer *ll_kmer_ptr = new LinkedlistKmer;
	HashTable *hash_table_ptr;

	// Check range k_min..k_max
	for (size_t i = k_min; i <= k_max; ++i){
		size_t max_mm = (modulo) ? (i * modulo) / 100 : 0;
		size_t max_gap = (modulo_gap) ? (i * modulo_gap) / 100 : 0;
		size_t max_gapmm = (modulo_gapmm) ? (i * modulo_gapmm) / 100 : 0;

		// Check if gaps allowed or not
		if (max_gap || max_gapmm){
			hash_table_ptr = Nessie::routine_get_kmers_k_gap(i, max_mm, max_gap, max_gapmm, start, end, 1);
		}
		else{
			hash_table_ptr = Nessie::routine_get_kmers_k_palindrome(i, max_mm, start, end);
		}

		// Add to LinkedlistKmer
		hash_table_ptr->append_to_LinkedlistKmer(ll_kmer_ptr);
	}

	return ll_kmer_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	get_max_kmers_palindrome_gap -- returns a ptr to a LinkedlistKmer that stores all the maximum Kmers with palindrome symmetry (maximum length max_k and minimum length min_k) in the interval,
//									allows for gaps
//
//	parameters:
//		max_k - max length of the kmers searched
//		min_k - min length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		modulo_gap - max gaps allowed are calculated as integer division k / modulo_gap [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
LinkedlistKmer *Nessie::get_max_kmers_palindrome_gap(size_t k_max, size_t k_min, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t start, size_t end){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Get max palindromes: starting index is larger than ending index"); }
	if (!k_min){ k_min = k_max; }	// if min_k is not defined it is set to default as max_k, only kmers of length max_k are searched
	if (k_min > k_max){ throw std::invalid_argument("Get max palindromes: min_k is larger than max_k"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Get max palindromes: max_k is longer than the sequence interval"); }

	// Initializing LinkedlistKmer
	LinkedlistKmer *ll_kmer_ptr = new LinkedlistKmer;
	HashTable *hash_table_ptr;

	// Check if gaps allowed or not
	if (modulo_gap || modulo_gapmm){
		hash_table_ptr = Nessie::routine_get_max_kmer_gap(k_max, k_min, modulo, modulo_gap, modulo_gapmm, start, end, 1);
	}
	else{
		hash_table_ptr = Nessie::routine_get_max_kmer_palindrome(k_max, k_min, modulo, start, end);
	}

	// Add to LinkedlistKmer
	hash_table_ptr->append_to_LinkedlistKmer(ll_kmer_ptr);

	return ll_kmer_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_kmers_mirror -- prints all the Kmers with mirror symmetry of length [k_min..k_max] in the interval
//
//	parameters:
//
//		k_min - minimum length of the kmers searched
//		k_max - maximum length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//		fout - ostream element to be used for printing [cout]
//		counts - bool value, if true print the counts information for the Kmer [true]
//		indexes - bool value, if true print the indexes list for the Kmer [true]
//		start_idx - starting index used to shift all other printed indexes [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::print_kmers_mirror(size_t k_min, size_t k_max, size_t modulo, size_t start, size_t end, std::ostream &fout, bool counts, bool indexes, size_t start_idx){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Print mirrors: starting index is larger than ending index"); }
	if (!k_max){ k_max = k_min; }	// if k_max is not defined it is set to default as k_min, only kmers of length k_min are searched
	if (k_min > k_max){ throw std::invalid_argument("Print mirrors: k_min is larger than k_max"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Print mirrors: k_max is longer than the sequence interval"); }

	// Check range k_min..k_max
	for (size_t i = k_min; i <= k_max; ++i){
		size_t max_mm;
		if (!modulo){
			max_mm = 0;
		}
		else{
			max_mm = (i * modulo) / 100;
		}
		HashTable *hash_table_ptr = Nessie::routine_get_kmers_k_mirror(i, max_mm, start, end);
		if (!start_idx){
			hash_table_ptr->print_table(fout, counts, indexes);
		}
		else{
			hash_table_ptr->print_table_shifted_indexes(start_idx, fout, counts, indexes);
		}
		delete hash_table_ptr;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_kmers_mirror_gap -- prints all the Kmers with mirror symmetry of length [k_min..k_max] in the interval,
//							  allows for gaps
//
//	parameters:
//		k_min - minimum length of the kmers searched
//		k_max - maximum length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		modulo_gap - max gaps allowed are calculated as integer division k / modulo_gap [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//		fout - ostream element to be used for printing [cout]
//		counts - bool value, if true print the counts information for the Kmer [true]
//		indexes - bool value, if true print the indexes list for the Kmer [true]
//		start_idx - starting index used to shift all other printed indexes [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::print_kmers_mirror_gap(size_t k_min, size_t k_max, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t start, size_t end, std::ostream &fout, bool counts, bool indexes, size_t start_idx){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Print mirrors gap: starting index is larger than ending index"); }
	if (!k_max){ k_max = k_min; }	// if k_max is not defined it is set to default as k_min, only kmers of length k_min are searched
	if (k_min > k_max){ throw std::invalid_argument("Print mirrors gap: k_min is larger than k_max"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Print mirrors gap: k_max is longer than the sequence interval"); }

	// Variables
	HashTable *hash_table_ptr;

	// Check range k_min..k_max
	for (size_t i = k_min; i <= k_max; ++i){
		size_t max_mm = (modulo) ? (i * modulo) / 100 : 0;
		size_t max_gap = (modulo_gap) ? (i * modulo_gap) / 100 : 0;
		size_t max_gapmm = (modulo_gapmm) ? (i * modulo_gapmm) / 100 : 0;

		// Check if gaps allowed or not
		if (max_gap || max_gapmm){
			hash_table_ptr = Nessie::routine_get_kmers_k_gap(i, max_mm, max_gap, max_gapmm, start, end, 0);
		}
		else{
			hash_table_ptr = Nessie::routine_get_kmers_k_mirror(i, max_mm, start, end);
		}

		// Printing
		if (!start_idx){
			hash_table_ptr->print_table(fout, counts, indexes);
		}
		else{
			hash_table_ptr->print_table_shifted_indexes(start_idx, fout, counts, indexes);
		}
		delete hash_table_ptr;	// need to delete the HashTable at each iteration otherwise too much memory is used
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_max_kmers_mirror_gap -- prints all the maximum Kmers with mirror symmetry (maximum length max_k and minimum length min_k) in the interval,
//								  allows for gaps
//
//	parameters:
//		max_k - max length of the kmers searched
//		min_k - min length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		modulo_gap - max gaps allowed are calculated as integer division k / modulo_gap [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//		fout - ostream element to be used for printing [cout]
//		counts - bool value, if true print the counts information for the Kmer [true]
//		indexes - bool value, if true print the indexes list for the Kmer [true]
//		start_idx - starting index used to shift all other printed indexes [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::print_max_kmers_mirror_gap(size_t k_max, size_t k_min, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t start, size_t end, std::ostream &fout, bool counts, bool indexes, size_t start_idx){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Print max mirrors gap: starting index is larger than ending index"); }
	if (!k_min){ k_min = k_max; }	// if min_k is not defined it is set to default as max_k, only kmers of length max_k are searched
	if (k_min > k_max){ throw std::invalid_argument("Print max mirrors gap: min_k is larger than max_k"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Print max mirrors gap: max_k is longer than the sequence interval"); }

	// Variables
	HashTable *hash_table_ptr;

	// Check if gaps allowed or not
	if (modulo_gap || modulo_gapmm){
		hash_table_ptr = Nessie::routine_get_max_kmer_gap(k_max, k_min, modulo, modulo_gap, modulo_gapmm, start, end, 0);
	}
	else{
		hash_table_ptr = Nessie::routine_get_max_kmer_mirror(k_max, k_min, modulo, start, end);
	}

	// Printing
	if (!start_idx){
		hash_table_ptr->print_table(fout, counts, indexes);
	}
	else{
		hash_table_ptr->print_table_shifted_indexes(start_idx, fout, counts, indexes);
	}
	delete hash_table_ptr;	// need to delete the HashTable at each iteration otherwise too much memory is used
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_kmers_palindrome -- prints all the Kmers with palindrome symmetry of length [k_min..k_max] in the interval
//
//	parameters:
//		k_min - minimum length of the kmers searched
//		k_max - maximum length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//		fout - ostream element to be used for printing [cout]
//		counts - bool value, if true print the counts information for the Kmer [true]
//		indexes - bool value, if true print the indexes list for the Kmer [true]
//		start_idx - starting index used to shift all other printed indexes [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::print_kmers_palindrome(size_t k_min, size_t k_max, size_t modulo, size_t start, size_t end, std::ostream &fout, bool counts, bool indexes, size_t start_idx){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Print palindromes: starting index is larger than ending index"); }
	if (!k_max){ k_max = k_min; }	// if k_max is not defined it is set to default as k_min, only kmers of length k_min are searched
	if (k_min > k_max){ throw std::invalid_argument("Print palindromes: k_min is larger than k_max"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Print palindromes: k_max is longer than the sequence interval"); }

	// Check range k_min..k_max
	for (size_t i = k_min; i <= k_max; ++i){
		size_t max_mm;
		if (!modulo){
			max_mm = 0;
		}
		else{
			max_mm = (i * modulo) / 100;
		}
		HashTable *hash_table_ptr = Nessie::routine_get_kmers_k_palindrome(i, max_mm, start, end);
		if (!start_idx){
			hash_table_ptr->print_table(fout, counts, indexes);
		}
		else{
			hash_table_ptr->print_table_shifted_indexes(start_idx, fout, counts, indexes);
		}
		delete hash_table_ptr;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_kmers_palindrome_gap -- prints all the Kmers with palindrome symmetry of length [k_min..k_max] in the interval,
//							  	  allows for gaps
//
//	parameters:
//		k_min - minimum length of the kmers searched
//		k_max - maximum length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		modulo_gap - max gaps allowed are calculated as integer division k / modulo_gap [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//		fout - ostream element to be used for printing [cout]
//		counts - bool value, if true print the counts information for the Kmer [true]
//		indexes - bool value, if true print the indexes list for the Kmer [true]
//		start_idx - starting index used to shift all other printed indexes [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::print_kmers_palindrome_gap(size_t k_min, size_t k_max, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t start, size_t end, std::ostream &fout, bool counts, bool indexes, size_t start_idx){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Print palindromes gap: starting index is larger than ending index"); }
	if (!k_max){ k_max = k_min; }	// if k_max is not defined it is set to default as k_min, only kmers of length k_min are searched
	if (k_min > k_max){ throw std::invalid_argument("Print palindromes gap: k_min is larger than k_max"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Print palindromes gap: k_max is longer than the sequence interval"); }

	// Variables
	HashTable *hash_table_ptr;

	// Check range k_min..k_max
	for (size_t i = k_min; i <= k_max; ++i){
		size_t max_mm = (modulo) ? (i * modulo) / 100 : 0;
		size_t max_gap = (modulo_gap) ? (i * modulo_gap) / 100 : 0;
		size_t max_gapmm = (modulo_gapmm) ? (i * modulo_gapmm) / 100 : 0;

		// Check if gaps allowed or not
		if (max_gap || max_gapmm){
			hash_table_ptr = Nessie::routine_get_kmers_k_gap(i, max_mm, max_gap, max_gapmm, start, end, 1);
		}
		else{
			hash_table_ptr = Nessie::routine_get_kmers_k_palindrome(i, max_mm, start, end);
		}

		// Printing
		if (!start_idx){
			hash_table_ptr->print_table(fout, counts, indexes);
		}
		else{
			hash_table_ptr->print_table_shifted_indexes(start_idx, fout, counts, indexes);
		}
		delete hash_table_ptr;	// need to delete the HashTable at each iteration otherwise too much memory is used
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_max_kmers_palindrome_gap -- prints all the maximum Kmers with palindrome symmetry (maximum length max_k and minimum length min_k) in the interval,
//								  allows for gaps
//
//	parameters:
//		max_k - max length of the kmers searched
//		min_k - min length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		modulo_gap - max gaps allowed are calculated as integer division k / modulo_gap [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//		fout - ostream element to be used for printing [cout]
//		counts - bool value, if true print the counts information for the Kmer [true]
//		indexes - bool value, if true print the indexes list for the Kmer [true]
//		start_idx - starting index used to shift all other printed indexes [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::print_max_kmers_palindrome_gap(size_t k_max, size_t k_min, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t start, size_t end, std::ostream &fout, bool counts, bool indexes, size_t start_idx){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Print max palindromes gap: starting index is larger than ending index"); }
	if (!k_min){ k_min = k_max; }	// if min_k is not defined it is set to default as max_k, only kmers of length max_k are searched
	if (k_min > k_max){ throw std::invalid_argument("Print max palindromes gap: min_k is larger than max_k"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Print max palindromes gap: max_k is longer than the sequence interval"); }

	// Variables
	HashTable *hash_table_ptr;

	// Check if gaps allowed or not
	if (modulo_gap || modulo_gapmm){
		hash_table_ptr = Nessie::routine_get_max_kmer_gap(k_max, k_min, modulo, modulo_gap, modulo_gapmm, start, end, 1);
	}
	else{
		hash_table_ptr = Nessie::routine_get_max_kmer_palindrome(k_max, k_min, modulo, start, end);
	}

	// Printing
	if (!start_idx){
		hash_table_ptr->print_table(fout, counts, indexes);
	}
	else{
		hash_table_ptr->print_table_shifted_indexes(start_idx, fout, counts, indexes);
	}
	delete hash_table_ptr;	// need to delete the HashTable at each iteration otherwise too much memory is used
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_get_kmers_k -- returns a ptr to an HashTable containing Kmer objects (kmers of length k in the interval)
//
//	parameters:
//		k - length of the kmers searched
//		start - starting index of the interval to search
//		end - ending index of the interval to search
//
//	note: this function works well with sequences of length up to ten thousands of bp,
//		  but it becomes inefficient (both in space and time) with strings of hundred thousands of bp
//
/////////////////////////////////////////////////////////////////////////////////////
HashTable *Nessie::routine_get_kmers_k(size_t k, size_t start, size_t end){

	// Some variables
	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (k > (end - start + 1)){ throw std::invalid_argument("Get kmers k: k is longer than the sequence interval"); }

	// Initializing HashTable
	HashTable *hash_table_ptr = new HashTable(false);	// bigger HashTable is used for faster look up

	// Defining bytes necessary to store the kmer
	size_t dna_bytes = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

	// Defining the mask to encode the kmer
	uint8_t mask_kmer[dna_bytes];
	uint8_t *mask_kmer_ptr = mask_kmer;
	std::memset(mask_kmer, 0, dna_bytes);

	// Initializing mask for the first kmer of length k in the interval
	size_t c = 0;	// counter for the mask_kmer creation
	for (size_t i = start; i < (start + k); ++i){
		uint8_t shift_DNA = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
		uint8_t shift_mask = (c & ((1 << 2) - 1)) << 1;	// shift to move from one base to the next one in the kmer_mask (uint8_t array)
		uint8_t base = (string_bit_ptr->data_ptr[i >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING
		mask_kmer_ptr[c >> 2] |= (base << shift_mask);	// insert base in the mask_kmer
		++c;
	}

	// Creating and adding first Kmer to the HashTable
	Kmer *kmer_ptr = new Kmer(k);	// defining ptr to new Kmer object
	copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_ptr, dna_bytes);
	kmer_ptr->indexes.push_back(start);
	kmer_ptr->counts += 1;
	hash_table_ptr->insert_kmer(kmer_ptr);	// adding Kmer to the HashTable

	// Sliding by one base at each iteration to get successive kmers
	size_t last_index = k - 1;
	uint8_t shift_mask = (last_index & ((1 << 2) - 1)) << 1;
	for (size_t i = (start + 1); i <= (end - k + 1); ++i){	//std::cout << i << std::endl;
		uint8_t shift_DNA = ((i + last_index) & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
		uint8_t base = (string_bit_ptr->data_ptr[(i + last_index) >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING
		shift_2_right(mask_kmer_ptr, dna_bytes);	// shift the uint8_t array encoding the kmer (removes the first base)
		mask_kmer_ptr[last_index >> 2] |= base << shift_mask;	// adding the new base to mask_kmer

		// Creating and adding i-th Kmer to the HashTable
		Kmer *kmer_ptr = new Kmer(k);
		copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_ptr, dna_bytes);
		kmer_ptr->indexes.push_back(i);
		kmer_ptr->counts += 1;
		hash_table_ptr->insert_kmer(kmer_ptr);
	}

	return hash_table_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_kmers -- prints all the Kmers of length [k_min..k_max] in the interval
//
//	parameters:
//		k_min - minimum length of the kmers searched
//		k_max - maximum length of the kmers searched [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//		fout - ostream element to be used for printing [cout]
//		counts - bool value, if true print the counts information for the Kmer [true]
//		indexes - bool value, if true print the indexes list for the Kmer [true]
//		start_idx - starting index used to shift all other printed indexes [0]
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::print_kmers(size_t k_min, size_t k_max, size_t start, size_t end, std::ostream &fout, bool counts, bool indexes, size_t start_idx){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Print kmers: starting index is larger than ending index"); }
	if (!k_max){ k_max = k_min; }	// if k_max is not defined it is set to default as k_min, only kmers of length k_min are searched
	if (k_min > k_max){ throw std::invalid_argument("Print kmers: k_min is larger than k_max"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Print kmers: k_max is longer than the sequence interval"); }

	// Check range k_min..k_max
	for (size_t i = k_min; i <= k_max; ++i){
		HashTable *hash_table_ptr = Nessie::routine_get_kmers_k(i, start, end);
		if (!start_idx){
			hash_table_ptr->print_table(fout, counts, indexes);
		}
		else{
			hash_table_ptr->print_table_shifted_indexes(start_idx, fout, counts, indexes);
		}
		delete hash_table_ptr;	// need to delete the HashTable at each iteration otherwise too much memory is used
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_init_counts
//
//	parameters:
// 		array_counts_ptr - ptr to a size_t array that stores the counts for the bases in the sequence
//		start - starting index of the interval
//		end - ending index of the interval
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::routine_init_counts(size_t *array_counts_ptr, size_t start, size_t end){

	for (size_t i = start; i <= end; ++i){
		uint8_t shift_DNA = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
														// 2 * (i % 4);
		uint8_t base = (string_bit_ptr->data_ptr[i >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING
		++array_counts_ptr[base];
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_shift_counts
//
//	parameters:
// 		array_counts_ptr - ptr to a size_t array that stores the counts for the bases in the sequence
//		interval len - length of the interval
//		shift - shift of the interval
//		start - starting index of the interval
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::routine_shift_counts(size_t *array_counts_ptr, size_t interval_len, size_t shift, size_t start){

	// Variables
	size_t end = start + interval_len;
	size_t end_shifted = end + shift;
	size_t start_shifted = start + shift;

	// Decreasing counts for bases lost after the shift
	for (size_t i = start; i < start_shifted; ++i){
		uint8_t shift_DNA = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
														// 2 * (i % 4);
		uint8_t base = (string_bit_ptr->data_ptr[i >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING
		--array_counts_ptr[base];
	}

	// Increasing counts for bases added after the shift
	for (size_t i = end; i < end_shifted; ++i){
		uint8_t shift_DNA = (i & ((1 << 2) - 1)) << 1;
		uint8_t base = (string_bit_ptr->data_ptr[i >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING
		++array_counts_ptr[base];
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_shannon_entropy -- returns the Shannon entropy score for a sequence of sequence_len
//
//	parameters:
//		array_counts_ptr - ptr to a size_t array that stores the counts for the bases in the sequence
//		sequence_len - length of the sequence
//
//	note: -SUM[A..T](p * log2(p))
//		  p is the frequency for the base
//
////////////////////////////////////////////////////////////////////////////////////
double Nessie::routine_shannon_entropy(size_t *array_counts_ptr, size_t sequence_len){

	// Variables
	double base_frequency;	// base frequency, p = base count / string length
	double base_probability;	// base probability, s = p * log2(p)
	double probabilities_sum = 0; //sum of the probabilities for each base
	double entropy; // Shannon entropy score

	for (size_t i = 0; i < 4; ++i) {
		if (array_counts_ptr[i]){
			base_frequency = (double) array_counts_ptr[i] / sequence_len;
			base_probability = base_frequency * log2(base_frequency);
		}
		else{
			base_probability = 0;
		}
		probabilities_sum += base_probability;
	}

	/* Calculating Shannon entropy and normalizing for 4 bases alphabet DNA */
	entropy = - probabilities_sum * 0.5;

	return entropy;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	shannon_entropy_sliding -- returns a ptr to a std::vector<double> storing the Shannon entropy scores for a sliding interval
//
//	parameters:
//		interval_len - length of the interval
//		shift - shift of the interval
//		start - starting index of the interval [0]
//		end - ending index of the interval [0]
//
//	note:
//
////////////////////////////////////////////////////////////////////////////////////
std::vector<double> *Nessie::shannon_entropy_sliding(size_t interval_len, size_t shift, size_t start, size_t end){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }
	if (start > end){ throw std::invalid_argument("Shannon sliding: starting index is larger than ending index"); }
	if (interval_len > (end - start + 1)){ throw std::invalid_argument("Shannon sliding: sliding interval is longer than the sequence"); }
	if (!shift){ throw std::invalid_argument("Shannon sliding: shift is not declared"); }

	// Variables
	size_t array_counts[4] = {0};
	size_t *array_counts_ptr = array_counts;
	double entropy;
	std::vector<double> *entropy_vector_ptr = new std::vector<double>;

	// Initializing array_counts for the first interval
	size_t end_0 = start + interval_len - 1;
	Nessie::routine_init_counts(array_counts_ptr, start, end_0);

	// Shifting interval and calculating entropy
	for (size_t i = start; i <= end - interval_len + 1; i += shift){	//std::cout << i << " - " << end - interval_len + 1 << std::endl;
		entropy = Nessie::routine_shannon_entropy(array_counts_ptr, interval_len);
		entropy_vector_ptr->push_back(entropy);
		Nessie::routine_shift_counts(array_counts_ptr, interval_len, shift, i);
	}

	return entropy_vector_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_shannon_entropy_sliding
//
//	parameters:
//		fout - ostream element to be used for printing [cout]
//		interval_len - length of the interval
//		shift - shift of the interval
//		start - starting index of the interval [0]
//		end - ending index of the interval [0]
//
//	note:
//
////////////////////////////////////////////////////////////////////////////////////
void Nessie::print_shannon_entropy_sliding(size_t interval_len, size_t shift, std::ostream &fout, size_t start, size_t end){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }
	if (start > end){ throw std::invalid_argument("Print Shannon sliding: starting index is larger than ending index"); }
	if (interval_len > (end - start + 1)){ throw std::invalid_argument("Print Shannon sliding: sliding interval is longer than the sequence"); }
	if (!shift){ throw std::invalid_argument("Print Shannon sliding: shift is not declared"); }

	// Variables
	size_t array_counts[4] = {0};
	size_t *array_counts_ptr = array_counts;
	double entropy;

	// Initializing array_counts for the first interval
	size_t end_0 = start + interval_len - 1;
	Nessie::routine_init_counts(array_counts_ptr, start, end_0);

	// Shifting interval and calculating entropy
	for (size_t i = start; i <= end - interval_len + 1; i += shift){	//std::cout << i << " - " << end - interval_len + 1 << std::endl;
		entropy = Nessie::routine_shannon_entropy(array_counts_ptr, interval_len);
		fout << i << "\t" << std::setprecision (16) << entropy;
		fout << "\tA:" << array_counts_ptr[0];
		fout << "\tC:" << array_counts_ptr[1];
		fout << "\tG:" << array_counts_ptr[2];
		fout << "\tT:" << array_counts_ptr[3];
		fout << std::endl;
		Nessie::routine_shift_counts(array_counts_ptr, interval_len, shift, i);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	shannon_entropy_interval -- returns the Shannon entropy score for an interval
//
//	parameters:
//		start - starting index of the interval [0]
//		end - ending index of the interval [0]
//
//	note:
//
////////////////////////////////////////////////////////////////////////////////////
double Nessie::shannon_entropy_interval(size_t start, size_t end){

	// Variables
	size_t array_counts[4] = {0};
	size_t *array_counts_ptr = array_counts;
	size_t sequence_len;

	// Calculating Shannon
	if (!end && !start){	// full string
		for (size_t i = 0; i < 4; ++i){
			array_counts_ptr[i] = string_bit_ptr->array_counts_UP_ptr[i] + string_bit_ptr->array_counts_LOW_ptr[i];
			sequence_len = string_bit_ptr->data_len;
		}
	}
	else{	// interval
		if (!end){ end = (string_bit_ptr->data_len) - 1; }
		if (start > end){ throw std::invalid_argument("Shannon interval: starting index is larger than ending index"); }
		Nessie::routine_init_counts(array_counts_ptr, start, end);
		sequence_len = end - start + 1;
	}

	return Nessie::routine_shannon_entropy(array_counts_ptr, sequence_len);
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_shannon_entropy_interval
//
//	parameters:
//		fout - ostream element to be used for printing [cout]
//		start - starting index of the interval [0]
//		end - ending index of the interval [0]
//
//	note:
//
////////////////////////////////////////////////////////////////////////////////////
void Nessie::print_shannon_entropy_interval(std::ostream &fout, size_t start, size_t end){

	// Variables
	size_t array_counts[4] = {0};
	size_t *array_counts_ptr = array_counts;
	size_t sequence_len;

	// Calculating Shannon
	if (!end && !start){	// full string
		for (size_t i = 0; i < 4; ++i){
			array_counts_ptr[i] = string_bit_ptr->array_counts_UP_ptr[i] + string_bit_ptr->array_counts_LOW_ptr[i];
			sequence_len = string_bit_ptr->data_len;
		}
	}
	else{	// interval
		if (!end){ end = (string_bit_ptr->data_len) - 1; }
		if (start > end){ throw std::invalid_argument("Print Shannon interval: starting index is larger than ending index"); }
		Nessie::routine_init_counts(array_counts_ptr, start, end);
		sequence_len = end - start + 1;
	}

	fout << std::setprecision (16) << Nessie::routine_shannon_entropy(array_counts_ptr, sequence_len);
//	fout << "\tA:" << array_counts_ptr[0];
//	fout << "\tC:" << array_counts_ptr[1];
//	fout << "\tG:" << array_counts_ptr[2];
//	fout << "\tT:" << array_counts_ptr[3];
	fout << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	linguistic_complexity_interval -- returns the Linguistic complexity for an interval
//
//	parameters:
//		start - starting index of the interval [0]
//		end - ending index of the interval [0]
//		k_min - minimum length of the kmers considered [0]
//		k_max - maximum length of the kmers considered [0]
//
//	note: (SUM[k_min..k_max] Vk) / (SUM[k_min..k_max] Vmaxk)
//		Vk is the number of different kmers of length k in the sequence
//		Vmaxk is the maximum possible number of kmers of the length k in the sequence
//
////////////////////////////////////////////////////////////////////////////////////
double Nessie::linguistic_complexity_interval(size_t start, size_t end, size_t k_min, size_t k_max){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }
	if (!k_min){ k_min = 1; }
	if (!k_max){ k_max = (string_bit_ptr->data_len < 20) ? string_bit_ptr->data_len : 20; }
	if (start > end){ throw std::invalid_argument("Linguistic interval: starting index is larger than ending index"); }
	if (k_min > k_max){ throw std::invalid_argument("Linguistic interval: k_min is larger than k_max"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Linguistic interval: k_max is longer than the sequence interval"); }

	// Variables
	double complexity;
	size_t Vmaxk_tot = 0;
	size_t Vk_tot = 0;

	// Calculating Vmaxk total
	for (size_t i = k_min; i <= k_max; ++i){ //std::cout << i << std::endl;

		// Possible kmers counts
		size_t possible_kmers = (uint64_t) 2 << ((i << 1) - 1);
		size_t possible_kmers_sequence = end - start - i + 2;
		Vmaxk_tot += (possible_kmers < possible_kmers_sequence) ? possible_kmers : possible_kmers_sequence;

		// Real kmers counts
		std::list<uint64_t> *list_ptr = routine_get_kmers_k_unique(i, start, end);
		Vk_tot += list_ptr->size();
		delete list_ptr;
	}

	// Calculating complexity
	complexity =  (double) Vk_tot / Vmaxk_tot;

	return complexity;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_linguistic_complexity_interval -- returns the Linguistic complexity for an interval
//
//	parameters:
//		fout - ostream element to be used for printing [cout]
//		start - starting index of the interval [0]
//		end - ending index of the interval [0]
//		k_min - minimum length of the kmers considered [0]
//		k_max - maximum length of the kmers considered [0]
//
//	note: (SUM[k_min..k_max] Vk) / (SUM[k_min..k_max] Vmaxk)
//		Vk is the number of different kmers of length k in the sequence
//		Vmaxk is the maximum possible number of kmers of the length k in the sequence
//
////////////////////////////////////////////////////////////////////////////////////
void Nessie::print_linguistic_complexity_interval(std::ostream &fout, size_t start, size_t end, size_t k_min, size_t k_max){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }
	if (!k_min){ k_min = 1; }
	if (!k_max){ k_max = (string_bit_ptr->data_len < 20) ? string_bit_ptr->data_len : 20; }
	if (start > end){ throw std::invalid_argument("Print Linguistic interval: starting index is larger than ending index"); }
	if (k_min > k_max){ throw std::invalid_argument("Print Linguistic interval: k_min is larger than k_max"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Print Linguistic interval: k_max is longer than the sequence interval"); }

	// Variables
	double complexity;
	size_t Vmaxk_tot = 0;
	size_t Vk_tot = 0;

	// Calculating Vmaxk total
	for (size_t i = k_min; i <= k_max; ++i){ //std::cout << i << std::endl;

		// Possible kmers counts
		size_t possible_kmers = (uint64_t) 2 << ((i << 1) - 1);
		size_t possible_kmers_sequence = end - start - i + 2;
		Vmaxk_tot += (possible_kmers < possible_kmers_sequence) ? possible_kmers : possible_kmers_sequence;

		// Real kmers counts
		std::list<uint64_t> *list_ptr = routine_get_kmers_k_unique(i, start, end);
		Vk_tot += list_ptr->size();
		delete list_ptr;
	}

	// Calculating complexity
	complexity =  (double) Vk_tot / Vmaxk_tot;

	fout << std::setprecision (16) << complexity << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	linguistic_complexity_sliding -- returns a ptr to a std::vector<double> storing the Linguistic complexity for a sliding interval
//
//	parameters:
//		interval_len - length of the interval [0]
//		shift - shift of the interval [0]
//		start - starting index of the interval [0]
//		end - ending index of the interval [0]
//		k_min - minimum length of the kmers considered [0]
//		k_max - maximum length of the kmers considered [0]
//
//	note:
//
////////////////////////////////////////////////////////////////////////////////////
std::vector<double> *Nessie::linguistic_complexity_sliding(size_t interval_len, size_t shift, size_t start, size_t end, size_t k_min, size_t k_max){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }
	if (!interval_len){ interval_len = 20; }
	if (!shift){ shift = 10; }
	if (!k_min){ k_min = 1; }
	if (!k_max){ k_max = (interval_len < 20) ? interval_len : 20; }
	if (start > end){ throw std::invalid_argument("Linguistic sliding: starting index is larger than ending index"); }
	if (interval_len > (end - start + 1)){ throw std::invalid_argument("Linguistic sliding: sliding interval is longer than the sequence"); }
	if (k_min > k_max){ throw std::invalid_argument("Linguistic sliding: k_min is larger than k_max"); }
	if (k_max > interval_len){ throw std::invalid_argument("Linguistic sliding: k_max is longer than the sequence interval"); }

	// Variables
	double complexity;
	std::vector<double> *linguistic_vector_ptr = new std::vector<double>;

	// Calculating complexity
	for (size_t i = start; i <= end - interval_len + 1; i += shift){	//std::cout << i << " - " << end - interval_len + 1 << std::endl;
		size_t end_i = i + interval_len - 1;
		complexity = Nessie::linguistic_complexity_interval(i, end_i, k_min, k_max);
		linguistic_vector_ptr->push_back(complexity);
	}

	return linguistic_vector_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_linguistic_complexity_sliding -- returns a ptr to a std::vector<double> storing the Linguistic complexity for a sliding interval
//
//	parameters:
//		fout - ostream element to be used for printing [cout]
//		interval_len - length of the interval [0]
//		shift - shift of the interval [0]
//		start - starting index of the interval [0]
//		end - ending index of the interval [0]
//		k_min - minimum length of the kmers considered [0]
//		k_max - maximum length of the kmers considered [0]
//
//	note:
//
////////////////////////////////////////////////////////////////////////////////////
void Nessie::print_linguistic_complexity_sliding(size_t interval_len, size_t shift, std::ostream &fout, size_t start, size_t end, size_t k_min, size_t k_max){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }
	if (!interval_len){ interval_len = 20; }
	if (!shift){ shift = 10; }
	if (!k_min){ k_min = 1; }
	if (!k_max){ k_max = (interval_len < 20) ? interval_len : 20; }
	if (start > end){ throw std::invalid_argument("Print Linguistic sliding: starting index is larger than ending index"); }
	if (interval_len > (end - start + 1)){ throw std::invalid_argument("Print Linguistic sliding: sliding interval is longer than the sequence"); }
	if (k_min > k_max){ throw std::invalid_argument("Print Linguistic sliding: k_min is larger than k_max"); }
	if (k_max > interval_len){ throw std::invalid_argument("Print Linguistic sliding: k_max is longer than the sequence interval"); }

	// Variables
	double complexity;

	// Calculating complexity
	for (size_t i = start; i <= end - interval_len + 1; i += shift){	//std::cout << i << " - " << end - interval_len + 1 << std::endl;
		size_t end_i = i + interval_len - 1;
		complexity = Nessie::linguistic_complexity_interval(i, end_i, k_min, k_max);
		fout << i << "\t" << std::setprecision (16) << complexity << std::endl;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_convert_to_uint64
//
//	parameters:
//		mask_kmer_ptr - ptr to a uint8_t array that encodes a kmer
//		mask_kmer_len - length of the uint8_t array that encodes the kmer
//
/////////////////////////////////////////////////////////////////////////////////////
uint64_t Nessie::routine_convert_to_uint64(uint8_t *mask_kmer_ptr, size_t mask_kmer_len){

	uint64_t idx = 0;

	for (size_t i = 0; i < mask_kmer_len; ++i){
		idx |= mask_kmer_ptr[i] << (i * 8);
	}

	return idx;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_get_kmers_k_unique
//
//	parameters:
//		k - length of the kmers searched
//		start - starting index of the interval to search
//		end - ending index of the interval to search
//
//	note: this function works best with short sequences
//
/////////////////////////////////////////////////////////////////////////////////////
std::list<uint64_t> *Nessie::routine_get_kmers_k_unique(size_t k, size_t start, size_t end){

	// Some variables
	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (k > (end - start + 1)){ throw std::invalid_argument("Get kmers k unique: k is longer than the sequence interval"); }

	// Initializing list
	uint64_t idx;
	std::list<uint64_t> *list_ptr = new std::list<uint64_t>;

	// Defining bytes necessary to store the kmer
	size_t dna_bytes = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

	// Defining the mask to encode the kmer
	uint8_t mask_kmer[dna_bytes];
	uint8_t *mask_kmer_ptr = mask_kmer;
	std::memset(mask_kmer, 0, dna_bytes);

	// Initializing mask for the first kmer of length k in the interval
	size_t c = 0;	// counter for the mask_kmer creation
	for (size_t i = start; i < (start + k); ++i){
		uint8_t shift_DNA = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
		uint8_t shift_mask = (c & ((1 << 2) - 1)) << 1;	// shift to move from one base to the next one in the mask_kmer (uint8_t array)
		uint8_t base = (string_bit_ptr->data_ptr[i >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING
		mask_kmer_ptr[c >> 2] |= (base << shift_mask);	// insert base in the mask_kmer
		++c;
	}

	// Adding first kmer
	idx = routine_convert_to_uint64(mask_kmer_ptr, dna_bytes);
	list_ptr->push_back(idx);

	// Sliding by one base at each iteration to get successive kmers
	size_t last_index = k - 1;
	uint8_t shift_mask = (last_index & ((1 << 2) - 1)) << 1;
	for (size_t i = (start + 1); i <= (end - k + 1); ++i){	//std::cout << i << std::endl;
		uint8_t shift_DNA = ((i + last_index) & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
		uint8_t base = (string_bit_ptr->data_ptr[(i + last_index) >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING
		shift_2_right(mask_kmer_ptr, dna_bytes);	// shift the uint8_t array encoding the kmer (removes the first base)
		mask_kmer_ptr[last_index >> 2] |= base << shift_mask;	// adding the new base to mask_kmer

		// Adding i-th kmer
		idx = routine_convert_to_uint64(mask_kmer_ptr, dna_bytes);
		list_ptr->push_back(idx);
	}

	// Removing duplicates
	list_ptr->sort();
	list_ptr->unique();

	return list_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_check_triplex_forming -- check the BIT_ARRAY encoding the kmer for the triplex forming potential
//
//	parameters:
//		monomer_bitarray_ptr - ptr to a *BIT_ARRAY[4]
//		k - length of the kmer
//		max_mm - max number of mismatch allowed
//		max_purine - max number of non-purines allowed
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
bool Nessie::routine_check_triplex_forming(BIT_ARRAY **monomer_bitarray_ptr, size_t k, size_t max_mm, size_t max_purine){

	// Variables
	size_t array_counts[4] = {0};
	bool triplex = false;

	// Checking composition
	for(size_t i = 0; i < 4; ++i){
			array_counts[i] += bit_array_num_bits_set(monomer_bitarray_ptr[i]);
	}

	if (!((array_counts[0] + array_counts[2]) <= max_purine || (array_counts[1] + array_counts[3]) <= max_purine)){
		return triplex;
	}

	return Nessie::check_mirror_symmetry(monomer_bitarray_ptr, k, max_mm);
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_check_triplex_forming_interval -- check a sub-interval of the BIT_ARRAY encoding the kmer for the triplex forming potential
//
//	parameters:
//		monomer_bitarray_ptr - ptr to a *BIT_ARRAY[4]
//		k - length of the kmer
//		max_mm - max number of mismatch allowed
//		max_purine - max number of non-purines allowed
//		start - starting index of the interval to check
//		end - ending index of the interval to check
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
bool Nessie::routine_check_triplex_forming_interval(BIT_ARRAY **monomer_bitarray_ptr, size_t k, size_t max_mm, size_t max_purine, size_t start, size_t end){

	if (start > end){ throw std::invalid_argument("Triplex forming interval: starting index is larger than ending index"); }
	if ((end + 1) > k){ throw std::invalid_argument("Triplex forming interval: ending index is larger thank sequence end"); }

	// Variables
	size_t array_counts[4] = {0};
	bool triplex = false;

	// Checking composition
	for(size_t i = 0; i < 4; ++i){
		for(size_t j = start; j <= end; ++j){
			array_counts[i] += bit_array_get_bit(monomer_bitarray_ptr[i], j);
		}
	}

	if (!((array_counts[0] + array_counts[2]) <= max_purine || (array_counts[1] + array_counts[3]) <= max_purine)){
		return triplex;
	}

	return Nessie::routine_check_mirror_symmetry_interval(monomer_bitarray_ptr, k, max_mm, start, end);
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_check_triplex_forming_gap
//
//	parameters:
//		mask_kmer_ptr - ptr to a mask (uint8_t array) encoding the kmer
//		k - length of the kmer
//		max_mm - max number of mismatch allowed
//		max_gap - max number of gaps allowed
//		max_purine - max number of non-purines allowed
//
/////////////////////////////////////////////////////////////////////////////////////
bool Nessie::routine_check_triplex_forming_gap(uint8_t *mask_kmer_ptr, std::vector<bool> *align_vector_ptr, size_t k, size_t max_mm, size_t max_gap, size_t max_gapmm, size_t max_purine){

	// Variables
	size_t array_counts[4] = {0};
	bool triplex = false;

	for (size_t i = 0; i < k; ++i){
		uint8_t shift_DNA = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
														// 2 * (i % 4);
		uint8_t base = (mask_kmer_ptr[i >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING
		array_counts[base] += 1;
	}

	if (!((array_counts[0] + array_counts[2]) <= max_purine || (array_counts[1] + array_counts[3]) <= max_purine)){
		return triplex;
	}

	return routine_check_global_alignment(mask_kmer_ptr, align_vector_ptr, k, max_mm, max_gap, max_gapmm, 0);
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_check_triplex_forming_gap_interval
//
//	parameters:
//		mask_kmer_ptr - ptr to a mask (uint8_t array) encoding the kmer
//		k - length of the kmer
//		max_mm - max number of mismatch allowed
//		max_gap - max number of gaps allowed
//		max_purine - max number of non-purines allowed
//		start - starting index of the interval to check
//		end - ending index of the interval to check
//
/////////////////////////////////////////////////////////////////////////////////////
bool Nessie::routine_check_triplex_forming_gap_interval(uint8_t *mask_kmer_ptr, std::vector<bool> *align_vector_ptr, size_t k, size_t max_mm, size_t max_gap, size_t max_gapmm, size_t max_purine, size_t start, size_t end){

	if (start > end){ throw std::invalid_argument("Triplex forming gap interval: starting index is larger than ending index"); }
	if ((end + 1) > k){ throw std::invalid_argument("Triplex forming gap interval: ending index is larger than sequence end"); }

	// Variables
	size_t array_counts[4] = {0};
	bool triplex = false;

	for (size_t i = start; i <= end; ++i){
		uint8_t shift_DNA = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
														// 2 * (i % 4);
		uint8_t base = (mask_kmer_ptr[i >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING
		array_counts[base] += 1;
	}

	if (!((array_counts[0] + array_counts[2]) <= max_purine || (array_counts[1] + array_counts[3]) <= max_purine)){
		return triplex;
	}

	return Nessie::routine_check_global_alignment_interval(mask_kmer_ptr, align_vector_ptr, k, max_mm, max_gap, max_gapmm, 0, start, end);
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_get_kmers_k_triplex -- returns a ptr to an HashTable containing Kmer objects (kmers of length k with triplex forming potential in the interval)
//
//	parameters:
//		k - length of the kmers searched
//		max_mm - maximum number of mismatch allowed
//		max_purine - max number of non-purines allowed
//		start - starting index of the interval to search
//		end - ending index of the interval to search
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
HashTable *Nessie::routine_get_kmers_k_triplex(size_t k, size_t max_mm, size_t max_purine, size_t start, size_t end){	//size_t check = 0;

	// Some variables
	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (k > (end - start + 1)){ throw std::invalid_argument("Get triplex k: k is longer than the sequence interval"); }

	// Initializing HashTable
	HashTable *hash_table_ptr = new HashTable(false);

	// Creating BIT_ARRAY to store monomers indexes for the mask
	BIT_ARRAY *monomer_bitarray[4];
	BIT_ARRAY **monomer_bitarray_ptr = monomer_bitarray;
	for (size_t i = 0; i < 4; ++i){
		monomer_bitarray_ptr[i] = bit_array_create(k);
	}

	// Defining bytes necessary to store the kmer
	size_t dna_bytes = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

	// Defining the mask to encode the kmer
	uint8_t mask_kmer[dna_bytes];
	uint8_t *mask_kmer_ptr = mask_kmer;
	std::memset(mask_kmer, 0, dna_bytes);

	// Initializing bit_array and mask for the first kmer of length k in the interval
	Nessie::routine_init_bitarray_and_mask(monomer_bitarray_ptr, mask_kmer_ptr, k, start);

	// Check mirror simmetry for the first kmer
	if (Nessie::routine_check_triplex_forming(monomer_bitarray_ptr, k, max_mm, max_purine)){
		Kmer *kmer_ptr = new Kmer(k);	// defining ptr to new Kmer object
		copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_ptr, dna_bytes);
		kmer_ptr->indexes.push_back(start);
		kmer_ptr->counts += 1;
		hash_table_ptr->insert_kmer(kmer_ptr);	// adding Kmer to the HashTable
	}

	// Sliding by one base at each iteration to get successive kmers and checking their simmetry
	for (size_t i = (start + 1); i <= (end - k + 1); ++i){
		Nessie::routine_shift_bitarray_and_mask(monomer_bitarray_ptr, mask_kmer_ptr, dna_bytes, k, i);

		if (Nessie::routine_check_triplex_forming(monomer_bitarray_ptr, k, max_mm, max_purine)){		//++check;
			Kmer *kmer_ptr = new Kmer(k);
			copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_ptr, dna_bytes);
			kmer_ptr->indexes.push_back(i);
			kmer_ptr->counts += 1;
			hash_table_ptr->insert_kmer(kmer_ptr);
		}
	}
	//std::cout << "CHECK " << check << std::endl;
	return hash_table_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_get_kmers_k_triplex_gap -- returns a ptr to an HashTable containing Kmer objects (kmers of length k with with triplex forming potential in the interval),
//							   		   allows gap presence
//
//	parameters:
//		k - length of the kmers searched
//		max_mm - maximum number of mismatch allowed
//		max_gap - max number of gaps allowed
//		max_purine - max number of non-purines allowed
//		start - starting index of the interval to search
//		end - ending index of the interval to search
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
HashTable *Nessie::routine_get_kmers_k_triplex_gap(size_t k, size_t max_mm, size_t max_gap, size_t max_gapmm, size_t max_purine, size_t start, size_t end){	//size_t check = 0;

	// Some variables
	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (k > (end - start + 1)){ throw std::invalid_argument("Get triplex k gap: k is longer than the sequence interval"); }

	// Initializing HashTable
	HashTable *hash_table_ptr = new HashTable(false);

	// Defining bytes necessary to store the kmer
	size_t dna_bytes = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

	// Defining the mask to encode the kmer
	uint8_t mask_kmer[dna_bytes];
	uint8_t *mask_kmer_ptr = mask_kmer;
	std::memset(mask_kmer, 0, dna_bytes);

	// Initializing mask for the first kmer of length k in the interval
	Nessie::routine_init_mask(mask_kmer_ptr, k, start);
	std::vector<bool> *align_vector_ptr = new std::vector<bool>;

	// Check mirror simmetry for the first kmer
	if (Nessie::routine_check_triplex_forming_gap(mask_kmer_ptr, align_vector_ptr, k, max_mm, max_gap, max_gapmm, max_purine)){
		Kmer *kmer_ptr = new Kmer(k);	// defining ptr to new Kmer object
		copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_ptr, dna_bytes);
		kmer_ptr->indexes.push_back(start);
		kmer_ptr->counts += 1;
		kmer_ptr->alignment_ptr = align_vector_ptr;
		hash_table_ptr->insert_kmer(kmer_ptr);	// adding Kmer to the HashTable
	}
	else{
		delete align_vector_ptr;
	}

	// Sliding by one base at each iteration to get successive kmers and checking their simmetry
	for (size_t i = (start + 1); i <= (end - k + 1); ++i){
		Nessie::routine_shift_mask(mask_kmer_ptr, dna_bytes, k, i);
		std::vector<bool> *align_vector_ptr_i = new std::vector<bool>;

		if (Nessie::routine_check_triplex_forming_gap(mask_kmer_ptr, align_vector_ptr_i, k, max_mm, max_gap, max_gapmm, max_purine)){		//++check;
			Kmer *kmer_ptr = new Kmer(k);
			copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_ptr, dna_bytes);
			kmer_ptr->indexes.push_back(i);
			kmer_ptr->counts += 1;
			kmer_ptr->alignment_ptr = align_vector_ptr_i;
			hash_table_ptr->insert_kmer(kmer_ptr);
		}
		else {
			delete align_vector_ptr_i;
		}
	}
	//std::cout << "CHECK " << check << std::endl;
	return hash_table_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_get_max_kmer_triplex -- returns a ptr to an HashTable containing Kmer objects (longest kmers of maximum length max_k and minimum length min_k with triplex
//									forming potential in the interval)
//
//	parameters:
//		max_k - max length of the kmers searched
//		min_k - min length of the kmers searched
//		modulo - max mismatch allowed are calculated as integer division k / modulo
//		modulo_purine - max number on non-purines allowed are calculated as integer division k / modulo_purine
//		start - starting index of the interval to search
//		end - ending index of the interval to search
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
HashTable *Nessie::routine_get_max_kmer_triplex(size_t k_max, size_t k_min, size_t modulo, size_t modulo_purine, size_t start, size_t end){	//size_t check = 0;

	// Some variables
	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Get max triplex: max_k is longer than the sequence interval"); }

	// Variables
	bool added;
	size_t added_end = 0;

	// Initializing HashTable
	HashTable *hash_table_ptr = new HashTable(false);

	// Creating BIT_ARRAY to store monomers indexes for the mask
	BIT_ARRAY *monomer_bitarray[4];
	BIT_ARRAY **monomer_bitarray_ptr = monomer_bitarray;
	for (size_t i = 0; i < 4; ++i){
		monomer_bitarray_ptr[i] = bit_array_create(k_max);
	}

	// Defining bytes necessary to store the kmer
	size_t dna_bytes = (k_max >> 2) + (0 != (k_max & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

	// Defining the mask to encode the kmer
	uint8_t mask_kmer[dna_bytes];
	uint8_t *mask_kmer_ptr = mask_kmer;
	std::memset(mask_kmer, 0, dna_bytes);

	// Initializing mask for the first kmer of length max_k in the interval
	Nessie::routine_init_bitarray_and_mask(monomer_bitarray_ptr, mask_kmer_ptr, k_max, start);

	size_t k = k_max;
	size_t end_i = start + k - 1;
	while (k >= k_min){
		size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
		size_t max_purine = (modulo_purine) ? (k * modulo_purine) / 100 : 0;
		if (Nessie::routine_check_triplex_forming_interval(monomer_bitarray_ptr, k_max, max_mm, max_purine, 0, k - 1)){	//++check;
			Kmer *kmer_ptr = new Kmer(k);
			added_end = end_i;

			// Defining bytes necessary to store the kmer k
			size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

			// Defining the mask to encode the kmer k
			uint8_t mask_kmer_k[dna_bytes_k];
			uint8_t *mask_kmer_k_ptr = mask_kmer_k;
			std::memset(mask_kmer_k, 0, dna_bytes_k);

			// Initializing mask for the first kmer of length k in the interval
			Nessie::routine_init_mask(mask_kmer_k_ptr, k, start);

			copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
			kmer_ptr->indexes.push_back(start);
			kmer_ptr->counts += 1;
			hash_table_ptr->insert_kmer_var_len(kmer_ptr);
			break;
		}
	--k;
	--end_i;
	}

	// Sliding by one base at each iteration to get successive kmers and checking their simmetry
	size_t last_i = start + 1;
	for (size_t i = (start + 1); i <= (end - k_max + 1); ++i){
		Nessie::routine_shift_bitarray_and_mask(monomer_bitarray_ptr, mask_kmer_ptr, dna_bytes, k_max, i);

		size_t k = k_max;
		size_t end_i = i + k - 1;
		while ((k >= k_min) && (end_i > added_end)){
			size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
			size_t max_purine = (modulo_purine) ? (k * modulo_purine) / 100 : 0;
			if (Nessie::routine_check_triplex_forming_interval(monomer_bitarray_ptr, k_max, max_mm, max_purine, 0, k - 1)){	//++check;
				Kmer *kmer_ptr = new Kmer(k);
				added_end = end_i;

				// Defining bytes necessary to store the kmer k
				size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

				// Defining the mask to encode the kmer k
				uint8_t mask_kmer_k[dna_bytes_k];
				uint8_t *mask_kmer_k_ptr = mask_kmer_k;
				std::memset(mask_kmer_k, 0, dna_bytes_k);

				// Initializing mask for the first kmer of length k in the interval
				Nessie::routine_init_mask(mask_kmer_k_ptr, k, i);

				copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
				kmer_ptr->indexes.push_back(i);
				kmer_ptr->counts += 1;
				hash_table_ptr->insert_kmer_var_len(kmer_ptr);
				break;
			}
		--k;
		--end_i;
		}
		last_i = i;
	}

	// Check last interval that may be shorter than k_max and skipped in previous for loop
	size_t last_interval_length = end - last_i;
	if (last_interval_length >= k_min){

		// Defining bytes necessary to store the last interval kmer
		size_t dna_bytes_l = (last_interval_length >> 2) + (0 != (last_interval_length & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

		// Creating BIT_ARRAY to store monomers indexes for the mask
		BIT_ARRAY *monomer_bitarray_l[4];
		BIT_ARRAY **monomer_bitarray_l_ptr = monomer_bitarray_l;
		for (size_t i = 0; i < 4; ++i){
			monomer_bitarray_l_ptr[i] = bit_array_create(last_interval_length);
		}

		// Defining the mask to encode the last interval kmer
		uint8_t mask_kmer_l[dna_bytes_l];
		uint8_t *mask_kmer_l_ptr = mask_kmer_l;
		std::memset(mask_kmer_l, 0, dna_bytes_l);

		// Initializing mask for the last interval kmer
		Nessie::routine_init_bitarray_and_mask(monomer_bitarray_l_ptr, mask_kmer_l_ptr, last_interval_length, last_i + 1);

		for (size_t i = 0; i <= (last_interval_length - k_min); ++i){

			size_t k = last_interval_length - i;
			size_t end_i = last_interval_length - 1; //std::cout << i << " " << k << std::endl;
			while ((k >= k_min) && ((last_i + 1 + i + k - 1) > added_end)){ //std::cout << "- " << i << " " << end_i << " " << k << std::endl;
				size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
				size_t max_purine = (modulo_purine) ? (k * modulo_purine) / 100 : 0;
				if (Nessie::routine_check_triplex_forming_interval(monomer_bitarray_l_ptr, last_interval_length, max_mm, max_purine, i, end_i)){	//++check;
					Kmer *kmer_ptr = new Kmer(k);
					added_end = last_i + 1 + i + k - 1;

					// Defining bytes necessary to store the kmer k
					size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

					// Defining the mask to encode the kmer k
					uint8_t mask_kmer_k[dna_bytes_k];
					uint8_t *mask_kmer_k_ptr = mask_kmer_k;
					std::memset(mask_kmer_k, 0, dna_bytes_k);

					// Initializing mask for the first kmer of length k in the interval
					Nessie::routine_init_mask(mask_kmer_k_ptr, k, last_i + 1 + i);

					copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
					kmer_ptr->indexes.push_back(last_i + 1 + i);
					kmer_ptr->counts += 1;
					hash_table_ptr->insert_kmer_var_len(kmer_ptr);
					break;
				}
			--k;
			--end_i;
			}
		}
	}

	//std::cout << "CHECK " << check << std::endl;
	return hash_table_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	routine_get_max_kmer_triplex_gap -- returns a ptr to an HashTable containing Kmer objects (longest kmers of maximum length max_k and minimum length min_k with triplex
//										forming potential in the interval), allows gap presence
//
//	parameters:
//		max_k - max length of the kmers searched
//		min_k - min length of the kmers searched
//		modulo - max mismatch allowed are calculated as integer division k / modulo
//		modulo_gap - max gaps allowed are calculated as integer division k / modulo_gap
//		modulo_purine - max number on non-purines allowed are calculated as integer division k / modulo_purine
//		start - starting index of the interval to search
//		end - ending index of the interval to search
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
HashTable *Nessie::routine_get_max_kmer_triplex_gap(size_t k_max, size_t k_min, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t modulo_purine, size_t start, size_t end){	//size_t check = 0;

	// Some variables
	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Get max triplex gap: max_k is longer than the sequence interval"); }

	// Variables
	bool added;
	size_t added_end = 0;

	// Initializing HashTable
	HashTable *hash_table_ptr = new HashTable(false);

	// Defining bytes necessary to store the kmer
	size_t dna_bytes = (k_max >> 2) + (0 != (k_max & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

	// Defining the mask to encode the kmer
	uint8_t mask_kmer[dna_bytes];
	uint8_t *mask_kmer_ptr = mask_kmer;
	std::memset(mask_kmer, 0, dna_bytes);

	// Initializing mask for the first kmer of length max_k in the interval
	Nessie::routine_init_mask(mask_kmer_ptr, k_max, start);
	std::vector<bool> *align_vector_ptr = new std::vector<bool>;
	added = false;

	size_t k = k_max;
	size_t end_i = start + k - 1;
	while (k >= k_min){
		size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
		size_t max_gap = (modulo_gap) ? (k * modulo_gap) / 100 : 0;
		size_t max_gapmm = (modulo_gapmm) ? (k * modulo_gapmm) / 100 : 0;
		size_t max_purine = (modulo_purine) ? (k * modulo_purine) / 100 : 0;
		if (Nessie::routine_check_triplex_forming_gap_interval(mask_kmer_ptr, align_vector_ptr, k_max, max_mm, max_gap, max_gapmm, max_purine, 0, k - 1)){	//++check;
			Kmer *kmer_ptr = new Kmer(k);
			added_end = end_i;

			// Defining bytes necessary to store the kmer k
			size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

			// Defining the mask to encode the kmer k
			uint8_t mask_kmer_k[dna_bytes_k];
			uint8_t *mask_kmer_k_ptr = mask_kmer_k;
			std::memset(mask_kmer_k, 0, dna_bytes_k);

			// Initializing mask for the first kmer of length k in the interval
			Nessie::routine_init_mask(mask_kmer_k_ptr, k, start);

			copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
			kmer_ptr->indexes.push_back(start);
			kmer_ptr->counts += 1;
			kmer_ptr->alignment_ptr = align_vector_ptr;
			added = true;
			hash_table_ptr->insert_kmer_var_len(kmer_ptr);
			break;
		}
	--k;
	--end_i;
	}

	if (!added){
		delete align_vector_ptr;
	}

	// Sliding by one base at each iteration to get successive kmers and checking their simmetry
	size_t last_i = start + 1;
	for (size_t i = (start + 1); i <= (end - k_max + 1); ++i){
		Nessie::routine_shift_mask(mask_kmer_ptr, dna_bytes, k_max, i);
		std::vector<bool> *align_vector_ptr_i = new std::vector<bool>;
		added = false;

		size_t k = k_max;
		size_t end_i = i + k - 1;
		while ((k >= k_min) && (end_i > added_end)){
			size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
			size_t max_gap = (modulo_gap) ? (k * modulo_gap) / 100 : 0;
			size_t max_gapmm = (modulo_gapmm) ? (k * modulo_gapmm) / 100 : 0;
			size_t max_purine = (modulo_purine) ? (k * modulo_purine) / 100 : 0;
			if (Nessie::routine_check_triplex_forming_gap_interval(mask_kmer_ptr, align_vector_ptr_i, k_max, max_mm, max_gap, max_gapmm, max_purine, 0, k - 1)){	//++check;
				Kmer *kmer_ptr = new Kmer(k);
				added_end = end_i;

				// Defining bytes necessary to store the kmer k
				size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

				// Defining the mask to encode the kmer k
				uint8_t mask_kmer_k[dna_bytes_k];
				uint8_t *mask_kmer_k_ptr = mask_kmer_k;
				std::memset(mask_kmer_k, 0, dna_bytes_k);

				// Initializing mask for the first kmer of length k in the interval
				Nessie::routine_init_mask(mask_kmer_k_ptr, k, i);

				copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
				kmer_ptr->indexes.push_back(i);
				kmer_ptr->counts += 1;
				kmer_ptr->alignment_ptr = align_vector_ptr_i;
				added = true;
				hash_table_ptr->insert_kmer_var_len(kmer_ptr);
				break;
			}
		--k;
		--end_i;
		}
		last_i = i;

		if (!added){
			delete align_vector_ptr_i;
		}
	}

	// Check last interval that may be shorter than k_max and skipped in previous for loop
	size_t last_interval_length = end - last_i;
	if (last_interval_length >= k_min){

		// Defining bytes necessary to store the last interval kmer
		size_t dna_bytes_l = (last_interval_length >> 2) + (0 != (last_interval_length & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

		// Defining the mask to encode the last interval kmer
		uint8_t mask_kmer_l[dna_bytes_l];
		uint8_t *mask_kmer_l_ptr = mask_kmer_l;
		std::memset(mask_kmer_l, 0, dna_bytes_l);

		// Initializing mask for the last interval kmer
		Nessie::routine_init_mask(mask_kmer_l_ptr, last_interval_length, last_i + 1);

		for (size_t i = 0; i <= (last_interval_length - k_min); ++i){
			std::vector<bool> *align_vector_ptr_l = new std::vector<bool>;
			added = false;

			size_t k = last_interval_length - i;
			size_t end_i = last_interval_length - 1; //std::cout << i << " " << k << std::endl;
			while ((k >= k_min) && ((last_i + 1 + i + k - 1) > added_end)){ //std::cout << "- " << i << " " << end_i << " " << k << std::endl;
					size_t max_mm = (modulo) ? (k * modulo) / 100 : 0;
					size_t max_gap = (modulo_gap) ? (k * modulo_gap) / 100 : 0;
					size_t max_gapmm = (modulo_gapmm) ? (k * modulo_gapmm) / 100 : 0;
					size_t max_purine = (modulo_purine) ? (k * modulo_purine) / 100 : 0;
				if (Nessie::routine_check_triplex_forming_gap_interval(mask_kmer_l_ptr, align_vector_ptr_l, last_interval_length, max_mm, max_gap, max_gapmm, max_purine, i, end_i)){	//++check;
					Kmer *kmer_ptr = new Kmer(k);
					added_end = last_i + 1 + i + k - 1;

					// Defining bytes necessary to store the kmer k
					size_t dna_bytes_k = (k >> 2) + (0 != (k & ((1 << 2) - 1)));	// (k / 4) + (0 != (k % 4))

					// Defining the mask to encode the kmer k
					uint8_t mask_kmer_k[dna_bytes_k];
					uint8_t *mask_kmer_k_ptr = mask_kmer_k;
					std::memset(mask_kmer_k, 0, dna_bytes_k);

					// Initializing mask for the first kmer of length k in the interval
					Nessie::routine_init_mask(mask_kmer_k_ptr, k, last_i + 1 + i);

					copy_uint8_t_arry(kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len, mask_kmer_k_ptr, dna_bytes_k);
					kmer_ptr->indexes.push_back(last_i + 1 + i);
					kmer_ptr->counts += 1;
					kmer_ptr->alignment_ptr = align_vector_ptr_l;
					added = true;
					hash_table_ptr->insert_kmer_var_len(kmer_ptr);
					break;
				}
			--k;
			--end_i;
			}

			if (!added){
				delete align_vector_ptr_l;
			}
		}
	}

	//std::cout << "CHECK " << check << std::endl;
	return hash_table_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	get_kmers_triplex_gap -- returns a ptr to a LinkedlistKmer that stores all the Kmers with triplex forming potential of length [k_min..k_max] in the interval,
//							 allows for gaps
//
//	parameters:
//		k_min - minimum length of the kmers searched
//		k_max - maximum length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		modulo_gap - max gaps allowed are calculated as integer division k / modulo_gap [0]
//		modulo_purine - max number on non-purines allowed are calculated as integer division k / modulo_purine [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
LinkedlistKmer *Nessie::get_kmers_triplex_gap(size_t k_min, size_t k_max, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t modulo_purine, size_t start, size_t end){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Get triplexes: starting index is larger than ending index"); }
	if (!k_max){ k_max = k_min; }	// if k_max is not defined it is set to default as k_min, only kmers of length k_min are searched
	if (k_min > k_max){ throw std::invalid_argument("Get triplexes: k_min is larger than k_max"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Get triplexes: k_max is longer than the sequence interval"); }

	// Initializing LinkedlistKmer
	LinkedlistKmer *ll_kmer_ptr = new LinkedlistKmer;
	HashTable *hash_table_ptr;

	// Check range k_min..k_max
	for (size_t i = k_min; i <= k_max; ++i){
		size_t max_mm = (modulo) ? (i * modulo) / 100 : 0;
		size_t max_gap = (modulo_gap) ? (i * modulo_gap) / 100 : 0;
		size_t max_gapmm = (modulo_gapmm) ? (i * modulo_gapmm) / 100 : 0;
		size_t max_purine = (modulo_purine) ? (i * modulo_purine) / 100 : 0;

		// Check if gaps allowed or not
		if (max_gap || max_gapmm){
			hash_table_ptr = Nessie::routine_get_kmers_k_triplex_gap(i, max_mm, max_gap, max_gapmm, max_purine, start, end);
		}
		else{
			hash_table_ptr = Nessie::routine_get_kmers_k_triplex(i, max_mm, max_purine, start, end);
		}

		// Add to LinkedlistKmer
		hash_table_ptr->append_to_LinkedlistKmer(ll_kmer_ptr);
	}

	return ll_kmer_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	get_max_kmers_triplex_gap -- returns a ptr to a LinkedlistKmer that stores all the maximum Kmers with with triplex forming potential
//								 (maximum length max_k and minimum length min_k) in the interval, allows for gaps
//
//	parameters:
//		max_k - max length of the kmers searched
//		min_k - min length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		modulo_gap - max gaps allowed are calculated as integer division k / modulo_gap [0]
//		modulo_purine - max number on non-purines allowed are calculated as integer division k / modulo_purine [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
LinkedlistKmer *Nessie::get_max_kmers_triplex_gap(size_t k_max, size_t k_min, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t modulo_purine, size_t start, size_t end){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Get max triplexes: starting index is larger than ending index"); }
	if (!k_min){ k_min = k_max; }	// if min_k is not defined it is set to default as max_k, only kmers of length max_k are searched
	if (k_min > k_max){ throw std::invalid_argument("Get max triplexes: min_k is larger than max_k"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Get max triplexes: max_k is longer than the sequence interval"); }

	// Initializing LinkedlistKmer
	LinkedlistKmer *ll_kmer_ptr = new LinkedlistKmer;
	HashTable *hash_table_ptr;

	// Check if gaps allowed or not
	if (modulo_gap || modulo_gapmm){
		hash_table_ptr = Nessie::routine_get_max_kmer_triplex_gap(k_max, k_min, modulo, modulo_gap, modulo_gapmm, modulo_purine, start, end);
	}
	else{
		hash_table_ptr = Nessie::routine_get_max_kmer_triplex(k_max, k_min, modulo, modulo_purine, start, end);
	}

	// Add to LinkedlistKmer
	hash_table_ptr->append_to_LinkedlistKmer(ll_kmer_ptr);

	return ll_kmer_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_kmers_triplex_gap -- prints all the Kmers with triplex forming potential of length [k_min..k_max] in the interval,
//							  allows for gaps
//
//	parameters:
//		k_min - minimum length of the kmers searched
//		k_max - maximum length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		modulo_gap - max gaps allowed are calculated as integer division k / modulo_gap [0]
//		modulo_purine - max number on non-purines allowed are calculated as integer division k / modulo_purine [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//		fout - ostream element to be used for printing [cout]
//		counts - bool value, if true print the counts information for the Kmer [true]
//		indexes - bool value, if true print the indexes list for the Kmer [true]
//		start_idx - starting index used to shift all other printed indexes [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::print_kmers_triplex_gap(size_t k_min, size_t k_max, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t modulo_purine, size_t start, size_t end, std::ostream &fout, bool counts, bool indexes, size_t start_idx){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Print triplexes: starting index is larger than ending index"); }
	if (!k_max){ k_max = k_min; }	// if k_max is not defined it is set to default as k_min, only kmers of length k_min are searched
	if (k_min > k_max){ throw std::invalid_argument("Print triplexes: k_min is larger than k_max"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Print triplexes: k_max is longer than the sequence interval"); }

	// Variables
	HashTable *hash_table_ptr;

	// Check range k_min..k_max
	for (size_t i = k_min; i <= k_max; ++i){
		size_t max_mm = (modulo) ? (i * modulo) / 100 : 0;
		size_t max_gap = (modulo_gap) ? (i * modulo_gap) / 100 : 0;
		size_t max_gapmm = (modulo_gapmm) ? (i * modulo_gapmm) / 100 : 0;
		size_t max_purine = (modulo_purine) ? (i * modulo_purine) / 100 : 0;

		// Check if gaps allowed or not
		if (max_gap || max_gapmm){
			hash_table_ptr = Nessie::routine_get_kmers_k_triplex_gap(i, max_mm, max_gap, max_gapmm, max_purine, start, end);
		}
		else{
			hash_table_ptr = Nessie::routine_get_kmers_k_triplex(i, max_mm, max_purine, start, end);
		}

		// Printing
		if (!start_idx){
			hash_table_ptr->print_table(fout, counts, indexes);
		}
		else{
			hash_table_ptr->print_table_shifted_indexes(start_idx, fout, counts, indexes);
		}
		delete hash_table_ptr;	// need to delete the HashTable at each iteration otherwise too much memory is used
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_max_kmers_triplex_gap -- prints all the maximum Kmers with triplex forming potential (maximum length max_k and minimum length min_k) in the interval,
//								   allows for gaps
//
//	parameters:
//		max_k - max length of the kmers searched
//		min_k - min length of the kmers searched [0]
//		modulo - max mismatch are calculated as integer division k / modulo, if 0 then max_mm is set to 0 [0]
//		modulo_gap - max gaps allowed are calculated as integer division k / modulo_gap [0]
//		modulo_purine - max number on non-purines allowed are calculated as integer division k / modulo_purine [0]
//		start - starting index of the interval to search [0]
//		end - ending index of the interval to search [0]
//		fout - ostream element to be used for printing [cout]
//		counts - bool value, if true print the counts information for the Kmer [true]
//		indexes - bool value, if true print the indexes list for the Kmer [true]
//		start_idx - starting index used to shift all other printed indexes [0]
//
//	note:
//
/////////////////////////////////////////////////////////////////////////////////////
void Nessie::print_max_kmers_triplex_gap(size_t k_max, size_t k_min, size_t modulo, size_t modulo_gap, size_t modulo_gapmm, size_t modulo_purine, size_t start, size_t end, std::ostream &fout, bool counts, bool indexes, size_t start_idx ){

	if (!end){ end = (string_bit_ptr->data_len) - 1; }	// if end is not defined it is set to default as the end of the string
	if (start > end){ throw std::invalid_argument("Print max mirrors gap: starting index is larger than ending index"); }
	if (!k_min){ k_min = k_max; }	// if min_k is not defined it is set to default as max_k, only kmers of length max_k are searched
	if (k_min > k_max){ throw std::invalid_argument("Print max mirrors gap: min_k is larger than max_k"); }
	if (k_max > (end - start + 1)){ throw std::invalid_argument("Print max mirrors gap: max_k is longer than the sequence interval"); }

	// Variables
	HashTable *hash_table_ptr;

	// Check if gaps allowed or not
	if (modulo_gap || modulo_gapmm){
		hash_table_ptr = Nessie::routine_get_max_kmer_triplex_gap(k_max, k_min, modulo, modulo_gap, modulo_gapmm, modulo_purine, start, end);
	}
	else{
		hash_table_ptr = Nessie::routine_get_max_kmer_triplex(k_max, k_min, modulo, modulo_purine, start, end);
	}

	// Printing
	if (!start_idx){
		hash_table_ptr->print_table(fout, counts, indexes);
	}
	else{
		hash_table_ptr->print_table_shifted_indexes(start_idx, fout, counts, indexes);
	}
	delete hash_table_ptr;	// need to delete the HashTable at each iteration otherwise too much memory is used
}
