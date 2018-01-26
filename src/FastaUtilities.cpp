/**************************************************************************************
*
**	FUNCTIONS (FastaUtilities.cpp)
*		Implements the functions of the FastaUtilities header.
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

#include "FastaUtilities.h"

/////////////////////////////////////////////////////////////////////////////////////
// 									CLASS Fasta									   //
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//
//	Fasta (constructor)
//
//	parameters:
//		seq_id - sequence id ">..."
//		seq_sequence - sequence
//		seq_idx - std::vector<size_t> containing indexes of the unknown bases
//
/////////////////////////////////////////////////////////////////////////////////////
Fasta::Fasta(std::string &seq_id, std::string &seq_sequence, std::vector<size_t> &seq_idx_unknown){

	id = seq_id;
	sequence = seq_sequence;
	idx_unknown = seq_idx_unknown;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	get_id
//
/////////////////////////////////////////////////////////////////////////////////////
std::string &Fasta::get_id(){

	return id;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	get_sequence
//
/////////////////////////////////////////////////////////////////////////////////////
std::string &Fasta::get_sequence(){

	return sequence;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	get_idx_unknown
//
/////////////////////////////////////////////////////////////////////////////////////
std::vector<size_t> &Fasta::get_idx_unknown(){

	return idx_unknown;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	overloaded <<
//
/////////////////////////////////////////////////////////////////////////////////////
std::ostream &operator<<(std::ostream &out, Fasta &fasta){

	//Variables
	std::vector<size_t>::iterator it;

	out << fasta.id << std::endl;
	out << fasta.sequence << std::endl;

	for (it = fasta.idx_unknown.begin(); it != fasta.idx_unknown.end();){
		out << *it << "|";
		++it;
	} out << std::endl;

	return out;
}


/////////////////////////////////////////////////////////////////////////////////////
// 								CLASS MultiFasta								   //
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//
//	get_data
//
//	parameters:
//		inFile - istream element to be used for reading
//
/////////////////////////////////////////////////////////////////////////////////////
void MultiFasta::get_data(std::istream &inFile){

	// Variables
	char c;	//current read char
	bool first = true;
	bool add_id = false;
	size_t curr_base_idx = 0;

	//Variable Fasta
	std::string tmp_id;
	std::string tmp_sequence;
	std::vector<size_t> tmp_idx_unknown;

	// Reding inFile by char
	while (inFile.get(c)){
		if (('>' == c) && first){
			add_id = true;
			first = false;
		}
		else if (('>' == c) && !first){
			Fasta fasta_seq(tmp_id, tmp_sequence, tmp_idx_unknown);
			sequences_vector.push_back(fasta_seq);
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
	Fasta fasta_seq(tmp_id, tmp_sequence, tmp_idx_unknown);
	sequences_vector.push_back(fasta_seq);
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	get_sequences_vector
//
/////////////////////////////////////////////////////////////////////////////////////
std::vector<Fasta> &MultiFasta::get_sequences_vector(){

	return sequences_vector;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	overloaded <<
//
/////////////////////////////////////////////////////////////////////////////////////
std::ostream &operator<<(std::ostream &out, MultiFasta &multifasta){

	//Variables
	std::vector<Fasta>::iterator it;

	for (it = multifasta.sequences_vector.begin(); it != multifasta.sequences_vector.end();){
		out << *it;
		++it;
	}

	return out;
}
