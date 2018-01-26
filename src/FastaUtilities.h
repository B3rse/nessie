/**************************************************************************************
*
**	HEADER (FastaUtilities.h)
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

#ifndef __LIST_H_INCLUDED
#define __LIST_H_INCLUDED
#include <list>
#endif /*__LIST_H_INCLUDED */


//CLASS
#ifndef FASTAUTILITIES_H
#define FASTAUTILITIES_H

/////////////////////////////////////////////////////////////////////////////////////
//
//	CLASS Fasta DEFINITION
//		Class to store the information on a fasta file
//
/////////////////////////////////////////////////////////////////////////////////////
class Fasta{
private:
	std::string id;	//stores the id ">..." for the sequence
	std::string sequence;	//stores the input sequence from ifstream
	std::vector<size_t> idx_unknown;	//stores the indexes of bases that differ from A, C, G, T, a, c, g, t

public:
	Fasta(std::string &seq_id, std::string &seq_sequence, std::vector<size_t> &seq_idx_unknown);
	std::string &get_id();
	std::string &get_sequence();
	std::vector<size_t> &get_idx_unknown();
	friend std::ostream &operator<<(std::ostream &out, Fasta &fasta);
};

/////////////////////////////////////////////////////////////////////////////////////
//
//	CLASS MultiFasta DEFINITION
//		Class to store the information from a multi-fasta file
//
/////////////////////////////////////////////////////////////////////////////////////
class MultiFasta{
private:
	std::vector<Fasta> sequences_vector;

public:
//	MultiFasta(std::istream &inFile);
	void get_data(std::istream &inFile);
	std::vector<Fasta> &get_sequences_vector();
	friend std::ostream &operator<<(std::ostream &out, MultiFasta &multifasta);
};

#endif /* FASTAUTILITIES_H */
