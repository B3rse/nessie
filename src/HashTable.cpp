/**************************************************************************************
*
**	FUNCTIONS (HashTable.cpp)
*		Implements the functions of the HashTable class.
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


#include "HashTable.h"

/////////////////////////////////////////////////////////////////////////////////////
//
//	HashTable (constructor): initializes the HashTable
//
//	parameters:
//		short_array - bool, if true an array of length 256 (2^8) is used / if false a longer array of 65536 (2^16) is used [true]
//
/////////////////////////////////////////////////////////////////////////////////////
HashTable::HashTable(bool short_array){
	if (short_array){	//std::cout << "SHORT" << std::endl;
		array_short_ptr = new LinkedlistKmer[256];
		array_long_ptr = NULL;
	}
	else{	//std::cout << "LONG" << std::endl;
		array_short_ptr = NULL;
		array_long_ptr = new LinkedlistKmer[65536];
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	~HashTable (destructor): destructs the HashTable
//
/////////////////////////////////////////////////////////////////////////////////////
HashTable::~HashTable() {
	if (array_short_ptr){
		delete[] array_short_ptr;
		array_short_ptr = NULL;
	}
	else{
		delete[] array_long_ptr;
		array_long_ptr = NULL;
	}

}

/////////////////////////////////////////////////////////////////////////////////////
//
//	insert_kmer: inserts the Kmer into the HashTable if not already present,
//				 otherwise increases the counter and add the new indexes to the one that is already present
//
//	parameters:
//		kmer_ptr - a ptr to a Kmer object
//
//	note: works only with kmers of the same length
//
/////////////////////////////////////////////////////////////////////////////////////
void HashTable::insert_kmer(Kmer *kmer_ptr){

	// Variables
	LinkedlistKmer *array_ptr;
	uint16_t idx;

	if (array_short_ptr){	// this is the default
		array_ptr = array_short_ptr;
		idx = kmer_ptr->kmer_mask_ptr[0];
	}
	else{
		array_ptr = array_long_ptr;
		idx = (kmer_ptr->kmer_mask_ptr[1] << 8) | kmer_ptr->kmer_mask_ptr[0];
	}

	array_ptr[idx].update_insert_kmer_end(kmer_ptr);
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	insert_kmer_var_len: inserts the Kmer into the HashTable if not already present,
//				 		 otherwise increases the counter and add the new indexes to the one that is already present,
//						 handles kmers of variable length
//
//	parameters:
//		kmer_ptr - a ptr to a Kmer object
//
/////////////////////////////////////////////////////////////////////////////////////
void HashTable::insert_kmer_var_len(Kmer *kmer_ptr){

	// Variables
	LinkedlistKmer *array_ptr;
	uint16_t idx;

	if (array_short_ptr){	// this is the default
		array_ptr = array_short_ptr;
		idx = kmer_ptr->kmer_mask_ptr[0];
	}
	else{
		array_ptr = array_long_ptr;
		idx = (kmer_ptr->kmer_mask_ptr[1] << 8) | kmer_ptr->kmer_mask_ptr[0];
	}

	array_ptr[idx].update_insert_kmer_end_var_len(kmer_ptr);
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	count_kmers: returns the number of Kmers stored in the HashTable
//
/////////////////////////////////////////////////////////////////////////////////////
size_t HashTable::count_kmers(){

	// Variables
	size_t count_kmers = 0;
	LinkedlistKmer *array_ptr;
	size_t array_len;

	if (array_short_ptr){
		array_ptr = array_short_ptr;
		array_len = 256;
	}
	else{
		array_ptr = array_long_ptr;
		array_len = 65536;
	}

	for (size_t i = 0; i < array_len; ++i){
		count_kmers += array_ptr[i].get_len();
	}

	return count_kmers;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	append_to_LinkedlistKmer: appends Kmers in the HashTable to a LinkedlistKmer
//
//	parameters:
//		ll_kmer_ptr - ptr to a LinkedlistKmer to append Kmers to
//
/////////////////////////////////////////////////////////////////////////////////////
void HashTable::append_to_LinkedlistKmer(LinkedlistKmer *ll_kmer_ptr){

	LinkedlistKmer *array_ptr;
	size_t array_len;

	if (array_short_ptr){
		array_ptr = array_short_ptr;
		array_len = 256;
	}
	else{
		array_ptr = array_long_ptr;
		array_len = 65536;
	}

	for (size_t i = 0; i < array_len; ++i){
		ll_kmer_ptr->concatenate(&array_ptr[i]);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_table: prints every Kmer in the HashTable
//
//	parameters:
//		fout - ostream element to be used for printing
//		counts - bool value, if true print the counts information for the Kmer [true]
//		indexes - bool value, if true print the indexes list for the Kmer [true]
//
/////////////////////////////////////////////////////////////////////////////////////
void HashTable::print_table(std::ostream &fout, bool counts, bool indexes){

	LinkedlistKmer *array_ptr;
	size_t array_len;

	if (array_short_ptr){
		array_ptr = array_short_ptr;
		array_len = 256;
	}
	else{
		array_ptr = array_long_ptr;
		array_len = 65536;
	}

	for (size_t i = 0; i < array_len; ++i){
		array_ptr[i].print_list(fout, counts, indexes);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_table_shifted_indexes: prints every Kmer in the HashTable,
//								 indexes are printed shifted by start_idx
//
//	parameters:
//		start_idx - starting index used to shift all other printed indexes
//		fout - ostream element to be used for printing
//		counts - bool value, if true print the counts information for the Kmer [true]
//		indexes - bool value, if true print the indexes list for the Kmer [true]
//
/////////////////////////////////////////////////////////////////////////////////////
void HashTable::print_table_shifted_indexes(size_t start_idx, std::ostream &fout, bool counts, bool indexes){

	LinkedlistKmer *array_ptr;
	size_t array_len;

	if (array_short_ptr){
		array_ptr = array_short_ptr;
		array_len = 256;
	}
	else{
		array_ptr = array_long_ptr;
		array_len = 65536;
	}

	for (size_t i = 0; i < array_len; ++i){
		array_ptr[i].print_list_shifted_indexes(start_idx, fout, counts, indexes);
	}
}

