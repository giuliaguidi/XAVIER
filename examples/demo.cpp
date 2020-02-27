//===========================================================================
// Title:  Xavier: High-Performance X-Drop Adaptive Banded Pairwise Alignment
// Author: G. Guidi, E. Younis
// Update:   24 February 2020
//===========================================================================

#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include <chrono>
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>
#include <iterator>
#include <x86intrin.h>
#include "xavier.h"
#include "../benchmark/utils.h"

//======================================================================
// DEMO
//======================================================================

int main(int argc, char const *argv[])
{
	srand(time(NULL));
	std::string seq1, seq2;

	/* Sequences length */
	// Just an example, XAVIER works best on long sequences
	int slen = 150;

	/* X-drop value */
	unsigned short x = 50;

	/* Seed/k-mer length */
	unsigned short k = 0;

	/* Bandwidth (the alignment path of the input sequence and the result does not go out of the band) */
	short bw = 64;

	/* Error rate composition */
	double pgap = 0.10; // indels probability
	double pmis = 0.05; // substitution probability

	/* Penalties (XAVIER supports only linear gap penalty) */
	short match    =  1;
	short mismatch = -1;
	short gap 	   = -1;

	/* Initialize scoring scheme */
	xavier::ScoringScheme scoringScheme(match, mismatch, gap);

	/* Generate pair of sequences */
	generate_random_sequence(seq1, slen);
	generate_random_sequence(seq2, slen);
	// seq1 = generate_mutated_sequence(seq1, slen, pmis, pgap, bw);
	// seq2 = generate_mutated_sequence(seq1, slen, pmis, pgap, bw);

	// seq1 = "GGGGGCGCAATTTTTCAGTTCCTGCCGGCAGTAGGGGACTCCGTTCTGATGAAGCTAACGTCCGTATCAGCAGCCCCCCAATGTTTGACACTTCTGCCAGGAGGCGGCGCTGGTTAAGTGCGCGTCATTCGATGCGTGAGAGGCAAGAAA";
	// seq2 = seq1;

	/* Seed starting position on seq1, seed starting position on seq2, k-mer length */
	xavier::Seed seed(0, 0, k);

	//===================================================================
	// XAVIER: High-Performance X-Drop Adaptive Banded Pairwise Alignment
	//===================================================================

	std::cout << "seq1	" 	<< seq1 << std::endl;
	std::cout << "seq2	" 	<< seq2 << std::endl;

	xavier::AlignmentResult result = xavier::seed_and_extend_right( seq1, seq2, scoringScheme, x, seed );

	std::cout << "result.bestScore	" 	<< result.bestScore << std::endl;
	std::cout << "result.exitScore	" 	<< result.exitScore << std::endl;

	std::cout << "result.begH	" 		<< result.begH << std::endl;
	std::cout << "result.endH	" 		<< result.endH << std::endl;
	std::cout << "result.begV	" 		<< result.begV << std::endl;
	std::cout << "result.endV	" 		<< result.endV << std::endl;
	
	std::cout << "result.matches " 		<< result.matched_pair.matches << std::endl;
	std::cout << "result.cigar "   		<< result.matched_pair.cigar   << std::endl;

	return 0;
}
