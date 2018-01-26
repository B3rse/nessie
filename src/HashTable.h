/**************************************************************************************
*
**	CLASS (HashTable.h)
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

#ifndef __LINKEDLISTKMER_H_INCLUDED
#define __LINKEDLISTKMER_H_INCLUDED
#include "LinkedlistKmer.h"
#endif /* __LINKEDLISTKMER_H_INCLUDED */


// CLASS
#ifndef HASHTABLE_H
#define HASHTABLE_H

/////////////////////////////////////////////////////////////////////////////////////
//
//	CLASS HashTable DEFINITION
//		HashTable -- class constructor
//		~HashTable -- class destructor
//
//		insert_kmer -- inserts the Kmer into the HashTable if not already present, otherwise increases the counter and add the new indexes
//		insert_kmer_var_len -- inserts the Kmer into the HashTable if not already present, otherwise increases the counter and add the new indexes, handles kmers of variable length
//		count_kmers -- returns the number of Kmers stored in the HashTable
//		append_to_LinkedlistKmer -- appends Kmers in the HashTable to a LinkedlistKmer
//		print_table -- prints every Kmer in the HashTable
//		print_table_shifted_indexes -- prints every Kmer in the HashTable, indexes are printed shifted by start_idx
//
/////////////////////////////////////////////////////////////////////////////////////
class HashTable {

private:
	// Variables
	LinkedlistKmer *array_short_ptr;
	LinkedlistKmer *array_long_ptr;

public:
	//Functions
	HashTable(bool short_array = true);
	~HashTable();

	void insert_kmer(Kmer *kmer_ptr);
	void insert_kmer_var_len(Kmer *kmer_ptr);
	size_t count_kmers();
	void append_to_LinkedlistKmer(LinkedlistKmer *ll_kmer_ptr);
	void print_table(std::ostream &fout = std::cout, bool counts = true, bool indexes = true);
	void print_table_shifted_indexes(size_t start_idx, std::ostream &fout = std::cout, bool counts = true, bool indexes = true);
};

#endif /* HASHTABLE_H */
