/**************************************************************************************
*
**	FUNCTIONS (Functions.h)
*		Implements few basic functions
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


//FUNCTIONS
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

/////////////////////////////////////////////////////////////////////////////////////
//
//	shift_2_right -- shift a uint8_t array by 2 bits to the right
//
//	parameters:
//		array_ptr - ptr to a uint8_t array
//		len_array - length of the array
//
/////////////////////////////////////////////////////////////////////////////////////
inline void shift_2_right(uint8_t *array_ptr, size_t len_array){

	array_ptr[0] >>= 2;

	for (size_t i = 1; i < len_array; ++i){
		uint8_t carry = 0x0;
		carry |= (BASE_MASK & array_ptr[i]);
		array_ptr[i - 1] |= carry << 6;
		array_ptr[i] >>= 2;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	hamming_distance_dna: calculate the Hamming distance between two uint8_t arrays encoding DNA strings
//						  the DNA strings has to be of equal length or an error is thrown
//
//	parameters:
//		str_1_ptr - ptr to a uint8_t array containing a DNA string encoded as bit
//		str_1_len - length of the encoded string
//		str_2_ptr - ptr to another uint8_t array containing a DNA string encoded as bit
//		str_2_len - length of the encoded string
//
/////////////////////////////////////////////////////////////////////////////////////
inline size_t hamming_distance_dna(uint8_t *str_1_ptr, size_t str_1_len, uint8_t *str_2_ptr, size_t str_2_len){

	// Defining variables
	size_t h_dist = 0;	// Hamming distance score

	// Check objects length
	if (str_1_len != str_2_len){
		throw std::invalid_argument("Hamming distance dna: the two strings are not of the same length");
	}
	else{
		// Number of bytes used to store the string in the uint8_t array
		size_t dna_bytes = (str_1_len >> 2) + (0 != (str_1_len & ((1 << 2) - 1)));

		// Iterating trough the uint8_t array bytes
		for (size_t i = 0; i < dna_bytes; ++i){
			h_dist += __builtin_popcount(str_1_ptr[i] ^ str_2_ptr[i]); // compiler function that works on int and counts the number of set bits
		}
	}

	return h_dist;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	hamming_distance: check if the Hamming distance is 0 between two uint8_t arrays,
//					  the arrays has to be of equal length or an error is thrown
//
//	parameters:
//		array_1_ptr - ptr to a uint8_t array
//		array_1_len - length of the array
//		array_2_ptr - ptr to another uint8_t array
//		array_2_len - length of the array
//
/////////////////////////////////////////////////////////////////////////////////////
inline bool hamming_distance_0(uint8_t *array_1_ptr, size_t array_1_len, uint8_t *array_2_ptr, size_t array_2_len){

	// Defining variables
	bool h_dist_0 = true;	// Hamming distance is 0

	// Check objects length
	if (array_1_len != array_2_len){
		throw std::invalid_argument("Hamming distance 0: the two arrays are not of the same length");
	}
	else{
		// Iterating trough the uint8_t array bytes
		for (size_t i = 0; i < array_1_len; ++i){
			if (__builtin_popcount(array_1_ptr[i] ^ array_2_ptr[i])){	// compiler function that works on int and counts the number of set bits
				h_dist_0 = false;
				break;
			}
		}
	}

	return h_dist_0;
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	print_to_string: print an uint8_t array encoding a DNA sequence back as a string
//
//	parameters:
//		sequence_bit_ptr - ptr to a uint8_t array containing a DNA string encoded as bit
//		sequence_len - length of the DNA string
//
/////////////////////////////////////////////////////////////////////////////////////
inline void print_to_string(uint8_t *sequence_bit_ptr, size_t sequence_len, std::ostream &fout = std::cout){

	for (size_t i = 0; i < sequence_len; ++i){
		uint8_t shift_DNA = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
														// 2 * (i % 4);
		uint8_t base = (sequence_bit_ptr[i >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING

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
//	print_to_string_no_nl: print an uint8_t array encoding a DNA sequence back as a string, doesn't add e new line at the end
//
//	parameters:
//		sequence_bit_ptr - ptr to a uint8_t array containing a DNA string encoded as bit
//		sequence_len - length of the DNA string
//
/////////////////////////////////////////////////////////////////////////////////////
inline void print_to_string_no_nl(uint8_t *sequence_bit_ptr, size_t sequence_len, std::ostream &fout = std::cout){

	for (size_t i = 0; i < sequence_len; ++i){
		uint8_t shift_DNA = (i & ((1 << 2) - 1)) << 1;	// shift from 0 to 6 with a step of two to move from one base to the next one in the uint8_t array
														// 2 * (i % 4);
		uint8_t base = (sequence_bit_ptr[i >> 2] & (BASE_MASK << shift_DNA)) >> shift_DNA;	// retrieving the ENCODING

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
}

/////////////////////////////////////////////////////////////////////////////////////
//
//	copy_uint8_t_arry: copy an uint8_t array into an empty one
//
//	parameters:
//		empty_array_ptr -- ptr to the empty array
// 		empty_array_len -- length of the empty array
//		array_to_copy_ptr -- ptr to the array to be copied into the empty one
//		array_to_copy_len -- length of the array to be copied into the empty one
//
/////////////////////////////////////////////////////////////////////////////////////
inline void copy_uint8_t_arry(uint8_t *empty_array_ptr, size_t empty_array_len, uint8_t *array_to_copy_ptr, size_t array_to_copy_len){

	// Check objects length
	if (empty_array_len != array_to_copy_len){
		throw std::invalid_argument("Copy array: the two arrays are not of the same length");
	}
	else{
		// Iterating trough the uint8_t array bytes
		for (size_t i = 0; i < empty_array_len; ++i){
			empty_array_ptr[i] |= array_to_copy_ptr[i];
		}
	}
}

#endif /* FUNCTIONS_H */
