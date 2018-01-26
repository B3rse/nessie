/**************************************************************************************
*
**	FUNCTIONS (LinkedlistKmer.cpp)
*		Implements the functions of the LinkedlistKmer class.
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


#include "LinkedlistKmer.h"

/////////////////////////////////////////////////////////////////////////////////////
// 								CLASS Kmer								   		   //
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//
//	print
//
/////////////////////////////////////////////////////////////////////////////////////
void Kmer::print(std::ostream &fout, bool counts, bool indexes){

	// Variables
	std::vector<size_t>::iterator IT;

	// Print
	fout << "$|" << k << '|';
	print_to_string(kmer_mask_ptr, k, fout);	// print kmer as string
	if (counts){
		fout << "@counts: " << this->counts << std::endl;
	}
	if (indexes){
		fout << "@indexes: ";
		for (IT = this->indexes.begin(); IT != this->indexes.end();){
			fout << *IT << '|';
			++IT;
		}
		fout << std::endl;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_full
//
/////////////////////////////////////////////////////////////////////////////////////
void Kmer::print_full(std::ostream &fout, bool counts, bool indexes){

	// Variables
	std::vector<size_t>::iterator IT;
	std::vector<bool>::reverse_iterator it;

	// Print
	fout << "$|" << k << '|';
	print_to_string_no_nl(kmer_mask_ptr, k, fout);	// print kmer as string
	if (alignment_ptr){
		fout << '|';
		for (it = alignment_ptr->rbegin(); it != alignment_ptr->rend();){
			fout << *it;
			++it;
		}
	}
	fout << std::endl;
	if (counts){
		fout << "@counts: " << this->counts << std::endl;
	}
	if (indexes){
		fout << "@indexes: ";
		for (IT = this->indexes.begin(); IT != this->indexes.end();){
			fout << *IT << '|';
			++IT;
		}
		fout << std::endl;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_full_shifted_indexes
//
/////////////////////////////////////////////////////////////////////////////////////
void Kmer::print_full_shifted_indexes(size_t start_idx, std::ostream &fout, bool counts, bool indexes){

	// Variables
	std::vector<size_t>::iterator IT;
	std::vector<bool>::reverse_iterator it;

	// Print
	fout << "$|" << k << '|';
	print_to_string_no_nl(kmer_mask_ptr, k, fout);	// print kmer as string
	if (alignment_ptr){
		fout << '|';
		for (it = alignment_ptr->rbegin(); it != alignment_ptr->rend();){
			fout << *it;
			++it;
		}
	}
	fout << std::endl;
	if (counts){
		fout << "@counts: " << this->counts << std::endl;
	}
	if (indexes){
		fout << "@indexes: ";
		for (IT = this->indexes.begin(); IT != this->indexes.end();){
			fout << *IT + start_idx << '|';
			++IT;
		}
		fout << std::endl;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
// 								CLASS LinkedlistKmer							   //
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//
//	LinkedlistKmer (constructor): initializes the LinkedlistKmer object
//
/////////////////////////////////////////////////////////////////////////////////////
LinkedlistKmer::LinkedlistKmer(){

	// Defining class variables
	first_ptr = NULL;	// no elements yet to point to
	len = 0;	// list is empty
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	~LinkedlistKmer (destructor): destructs the LinkedlistKmer object
//
/////////////////////////////////////////////////////////////////////////////////////
LinkedlistKmer::~LinkedlistKmer(){

	Kmer *ptr = first_ptr;
	Kmer *next_ptr = first_ptr;
	while (next_ptr){
		ptr = next_ptr;
		next_ptr = ptr->next_kmer_ptr;
		delete ptr;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	get_len: returns the length of the LinkedlistKmer
//
/////////////////////////////////////////////////////////////////////////////////////
size_t LinkedlistKmer::get_len(){

	return len;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	insert_kmer_front: inserts a Kmer object at the beginning of the LinkedlistKmer
//
//	parameters:
//		kmer_ptr -  ptr to a Kmer object to insert at the beginning of the LinkedlistKmer
//
/////////////////////////////////////////////////////////////////////////////////////
void LinkedlistKmer::insert_kmer_front(Kmer *kmer_ptr){

	len += 1;	// increasing list length
	kmer_ptr->next_kmer_ptr = first_ptr;	// new Kmer object is inserted at the beginning of the list
	first_ptr = kmer_ptr;	// first ptr is reassigned to point to the new object inserted
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	update_insert_kmer_end: updates the Kmer information if present or insert the Kmer at the end of the LinkedlistKmer if not present
//
//	parameters:
//		kmer_ptr - ptr to a Kmer object
//
//	note: works for kmer of the same length only
//
/////////////////////////////////////////////////////////////////////////////////////
void LinkedlistKmer::update_insert_kmer_end(Kmer *kmer_ptr){

	bool updated = false;

	if (first_ptr){
		Kmer *ptr = first_ptr;
		Kmer *next_ptr = first_ptr;
		while (next_ptr){
			ptr = next_ptr;

			if (hamming_distance_0(ptr->kmer_mask_ptr, ptr->kmer_mask_len, kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len)){
				ptr->indexes.insert(ptr->indexes.end(), kmer_ptr->indexes.begin(), kmer_ptr->indexes.end());	// non ho capito perchè l'editor da errore, sembra però andare bene
				ptr->counts += kmer_ptr->counts;
				updated = true;
				delete kmer_ptr;
				break;
			}
			next_ptr = ptr->next_kmer_ptr;
		}

		if (!updated){	// adding the Kmer at the end since it is not already in the LinkedlistKmer
			len += 1;
			ptr->next_kmer_ptr = kmer_ptr;
		}
	}
	else{
		this->insert_kmer_front(kmer_ptr);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	update_insert_kmer_end_var_len: updates the Kmer information if present or insert the Kmer at the end of the LinkedlistKmer if not present,
//									handles kmers of variable length
//
//	parameters:
//		kmer_ptr - ptr to a Kmer object
//
/////////////////////////////////////////////////////////////////////////////////////
void LinkedlistKmer::update_insert_kmer_end_var_len(Kmer *kmer_ptr){

	bool updated = false;

	if (first_ptr){
		Kmer *ptr = first_ptr;
		Kmer *next_ptr = first_ptr;
		while (next_ptr){
			ptr = next_ptr;

			if (ptr->k == kmer_ptr->k){	// compared kmers are of the same length
				if (hamming_distance_0(ptr->kmer_mask_ptr, ptr->kmer_mask_len, kmer_ptr->kmer_mask_ptr, kmer_ptr->kmer_mask_len)){
					ptr->indexes.insert(ptr->indexes.end(), kmer_ptr->indexes.begin(), kmer_ptr->indexes.end());	// non ho capito perchè l'editor da errore, sembra però andare bene
					ptr->counts += kmer_ptr->counts;
					updated = true;
					delete kmer_ptr;
					break;
				}
			}
			next_ptr = ptr->next_kmer_ptr;
		}

		if (!updated){	// adding the Kmer at the end since it is not already in the LinkedlistKmer
			len += 1;
			ptr->next_kmer_ptr = kmer_ptr;
		}
	}
	else{
		this->insert_kmer_front(kmer_ptr);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	concatenate: concatenates a LinkedlistKmer at the end of the LinkedlistKmer
//
//	parameters:
//		ll_kmer_ptr - ptr to a LinkedlistKmer to concatenate
//
/////////////////////////////////////////////////////////////////////////////////////
void LinkedlistKmer::concatenate(LinkedlistKmer *ll_kmer_ptr){

	len += ll_kmer_ptr->get_len();

	// Retrieving last element next_kmer_ptr
	if (first_ptr){	// there is already a non empty LinkedlistKmer, need to find the last element ptr (NULL) and set it to the first element of the concatenated LinkedlistKmer
		Kmer *ptr = first_ptr;
		Kmer *next_ptr = first_ptr;
		while (next_ptr){
			ptr = next_ptr;
			next_ptr = ptr->next_kmer_ptr;
		}
		ptr->next_kmer_ptr = ll_kmer_ptr->first_ptr;
	}
	else{	// LinkedlistKmer is empty so the first ptr simply need to be set to point to the first element of the concatenated LinkedlistKmer
		first_ptr = ll_kmer_ptr->first_ptr;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_list: prints each node in the LinkedlistKmer in a consecutive order
//
//	parameters:
//		fout - ostream element to be used for printing
//		counts - bool value, if true print the counts information for the Kmer [true]
//		indexes - bool value, if true print the indexes list for the Kmer [true]
//
/////////////////////////////////////////////////////////////////////////////////////
void LinkedlistKmer::print_list(std::ostream &fout, bool counts, bool indexes){

	// Defining variables
	Kmer *ptr = first_ptr;
	Kmer *next_ptr = first_ptr;
	std::vector<size_t>::iterator IT;

	// Printing kmers information
	while (next_ptr){
		ptr = next_ptr;
		ptr->print_full(fout, counts, indexes);
		next_ptr = ptr->next_kmer_ptr;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_list_shifted_indexes: prints each node in the LinkedlistKmer in a consecutive order,
//								printed indexes are shifted by start_idx
//
//	parameters:
//		start_idx - starting index used to shift all other printed indexes
//		fout - ostream element to be used for printing
//		counts - bool value, if true print the counts information for the Kmer [true]
//		indexes - bool value, if true print the indexes list for the Kmer [true]
//
/////////////////////////////////////////////////////////////////////////////////////
void LinkedlistKmer::print_list_shifted_indexes(size_t start_idx, std::ostream &fout, bool counts, bool indexes){

	// Defining variables
	Kmer *ptr = first_ptr;
	Kmer *next_ptr = first_ptr;
	std::vector<size_t>::iterator IT;

	// Printing kmers information
	while (next_ptr){
		ptr = next_ptr;
		ptr->print_full_shifted_indexes(start_idx, fout, counts, indexes);
		next_ptr = ptr->next_kmer_ptr;
	}
}
