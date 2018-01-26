/**************************************************************************************
*
**	CLASS (LinkedlistKmer.h)
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

#ifndef __FUNCTIONS_H_INCLUDED
#define __FUNCTIONS_H_INCLUDED
#include "Functions.h"
#endif /* __FUNCTIONS_H_INCLUDED */


// CLASS
#ifndef LINKEDLISTKMER_H
#define LINKEDLISTKMER_H

/////////////////////////////////////////////////////////////////////////////////////
//
//	CLASS Kmer DEFINITION
//		Class to store the information on a retrieved kmer
//
/////////////////////////////////////////////////////////////////////////////////////
class Kmer{

public:
	// Variables
	uint8_t *kmer_mask_ptr;	// ptr to the uint8_t array encoding the kmer as bit set
	size_t kmer_mask_len;	// length of the kmer_mask_ptr array
	std::vector<bool> *alignment_ptr;
    size_t k;	// length of the kmer
    size_t counts;	// counts of all the copies of the kmer in the string
    std::vector<size_t> indexes;	// std::vector<size_t> containing the indexes for all the copies of the kmer in the string
    Kmer *next_kmer_ptr;	// ptr to the next Kmer in a linked list, allows Kmer to be used to build a linked list

    // Functions
    inline Kmer(size_t kmer_len){
    	// Bytes necessary to store the kmer as a bit set
    	size_t dna_bytes = (kmer_len >> 2) + (0 != (kmer_len & ((1 << 2) - 1)));
    	kmer_mask_len = dna_bytes;

    	// Defining variables
    	kmer_mask_ptr = new uint8_t[dna_bytes];	// defining ptr to kmer_mask_ptr
    	std::memset(kmer_mask_ptr, 0, dna_bytes);	// initializing every bit of the kmer_mask_ptr array to 0
    	alignment_ptr = NULL;
    	k = kmer_len;
    	counts = 0;
    	next_kmer_ptr = NULL;
    }

    inline ~Kmer(){
    	delete[] kmer_mask_ptr;
    	kmer_mask_ptr = NULL;
    	if (alignment_ptr){
    		delete alignment_ptr;
    		alignment_ptr = NULL;
    	}
    	next_kmer_ptr = NULL;
    }

    void print(std::ostream &fout, bool counts = true, bool indexes = true);
    void print_full(std::ostream &fout, bool counts = true, bool indexes = true);
    void print_full_shifted_indexes(size_t start_idx, std::ostream &fout, bool counts = true, bool indexes = true);
};

/////////////////////////////////////////////////////////////////////////////////////
//
//	CLASS LinkedlistKmer DEFINITION
//		LinkedlistKmer -- class constructor
//		~LinkedlistKmer -- class destructor
//
//		get_len -- returns the length of the LinkedlistKmer
//		insert_kmer_front -- inserts a Kmer object at the beginning of the LinkedlistKmer
//		concatenate -- concatenates a LinkedlistKmer at the end of the LinkedlistKmer
//		update_insert_kmer_end -- updates information for the Kmer if present or insert the Kmer at the end of the LinkedlistKmer if not present
//		update_insert_kmer_end_var_len -- updates information for the Kmer if present or insert the Kmer at the end of the LinkedlistKmer if not present, handles kmers of variable length
//		update_insert_from_list_shifted_index -- updates information for the Kmer if present or insert the Kmer at the end of the LinkedlistKmer if not present,
//												 indexes are shifted by start_idx
//		print_list -- prints each node in the LinkedlistKmer in a consecutive order
//		print_list_shifted_indexes -- prints each node in the LinkedlistKmer in a consecutive order, printed indexes are shifted by start_idx
//
/////////////////////////////////////////////////////////////////////////////////////
class LinkedlistKmer{

private:
	Kmer *first_ptr;	// ptr to the first element of the LinkedlistKmer, initialized to NULL
	size_t len;	// length of the LinkedlistKmer

public:
	// Basic functions
	LinkedlistKmer();
	~LinkedlistKmer();

	size_t get_len();
	void insert_kmer_front(Kmer *kmer_ptr);
	void update_insert_kmer_end(Kmer *kmer_ptr);
	void update_insert_kmer_end_var_len(Kmer *kmer_ptr);
	void concatenate(LinkedlistKmer *ll_kmer_ptr);
	void print_list(std::ostream &fout = std::cout, bool counts = true, bool indexes = true);
	void print_list_shifted_indexes(size_t start_idx, std::ostream &fout = std::cout, bool counts = true, bool indexes = true);
};

#endif /* LINKEDLISTKMER_H */
