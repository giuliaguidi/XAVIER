//==============================================================================
// Title:  Xavier: High-Performance X-Drop Adaptive Banded Pairwise Alignment
// Author: G. Guidi, E. Younis
// Date:   29 April 2019
//==============================================================================

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <vector>
#include <iostream>
#include <string.h>
#include <omp.h>
#include <algorithm>
#include <chrono>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include <iterator>
#include <x86intrin.h>
#include <random>

#include "xavier.h"
#include "ksw2/ksw2.h"
#include "ksw2/ksw2_extz2_sse.c"

#include "parasail/parasail.h"
#include "parasail/parasail/io.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/stats.h"
#include "ssw.h"
#include "ssw_cpp.h"

// #include <edlib.h>

// #include <seqan/align.h>
// #include <seqan/seeds/seeds_extension.h>
// #include <seqan/sequence.h>
// #include <seqan/align.h>
// #include <seqan/seeds.h>
// #include <seqan/score.h>
// #include <seqan/modifier.h>
// #include <seqan/basic.h>
// #include <seqan/stream.h>

#ifdef __cplusplus
extern "C" {
#endif
#include "libgaba/gaba.h" 	 // sometimes the forefront vector will not reach the end
							 // of the sequences. It is more likely to occur when the input
							 // sequence lengths greatly differ
#ifdef __cplusplus
}
#endif

//======================================================================================
// GLOBAL VARIABLE DECLARATION
//======================================================================================

#define MAT		( 1)		// match score
#define MIS		(-1)		// mismatch score
#define GAP		(-1)		// gap score
#define PMIS 	(0.03)		// substitution probability
#define PGAP 	(0.13)		// insertion/deletion probability
#define BW 		(128)		// bandwidth (the alignment path of the input sequence and the result does not go out of the band)

#define KSW2
#define GABA
#define PARASAIL1
// #define EDLIB
#define PARASAIL2
#define SSW

//======================================================================================
// READ SIMULATOR
//======================================================================================

// Functions to generate sequences from https://github.com/ocxtal/libgaba
char random_base(void)
{
	char const table[4] = {'A', 'C', 'G', 'T'};
	return(table[rand() % 4]);
}

void generate_random_sequence(std::string& seq, int& len1)
{
	for(int i = 0; i < len1; i++)
		seq.append(1, random_base());
}

std::string generate_mutated_sequence(const std::string& seq, int& len2)
{
	int i, j, wave = 0;	// wave is q-coordinate of the alignment path
	std::string mutated;

	for(i = 0, j = 0; i < len2; i++)
	{
		if(((double)rand()/(double)RAND_MAX) < PMIS)
		{
			mutated.append(1, random_base());
			j++; // mismatch
		}
		else if(((double)rand()/(double)RAND_MAX) < PGAP)
		{
			if(rand() & 0x01 && wave > (-BW + 1))
			{
				char tmp = (j < len2) ? seq[j++] : random_base();
				mutated.append(1, tmp);
				j++; wave--; // deletion
			}
			else if(wave < (BW-2))
			{
				mutated.append(1, random_base());
				wave++; // insertion
			}
			else
			{
				char tmp = (j < len2) ? seq[j++] : random_base();
				mutated.append(1, tmp);
			}
		}
		else
		{
			char tmp = (j < len2) ? seq[j++] : random_base();
			mutated.append(1, tmp);
		}
	}

	return mutated;
}

//======================================================================================
// BENCHMARK CODE
//======================================================================================

int main(int argc, char const *argv[])
{
	/* Simulate read pair */
	std::string targetSeg;

	// std::default_random_engine generator;
	// std::normal_distribution<float> distribution(11000.0, 2000.0);

    // int len1 = (int)distribution(generator);
	// int len2 = (int)distribution(generator);
    int len1 = 10000;
	int len2 = 12000;

	std::cout << len1 << "	" << len2 << std::endl;

	generate_random_sequence(targetSeg, len1);
	std::string querySeg = generate_mutated_sequence(targetSeg, len2);

	int xdrop = std::stoi(argv[1]);
	int bw = std::stoi(argv[2]);

	//======================================================================================
	// XAVIER (vectorized SSE2 and AVX2, banded, (not yet) x-drop)
	//======================================================================================

	/* Initialize scoring scheme */
	xavier::ScoringScheme penalties ( MAT, MIS, GAP );

	/* Seed starting position on seq1, seed starting position on seq2, k-mer lenght */
	xavier::Seed seed(0, 0, 17);

	// 1st prototype without seed and x-drop termination
	std::chrono::duration<double> diff1;
	auto start1 = std::chrono::high_resolution_clock::now();

	xavier::AlignmentResult result = xavier::seed_and_extend_right( targetSeg, querySeg, penalties, xdrop, seed );

	auto end1 = std::chrono::high_resolution_clock::now();
	diff1 = end1-start1;

	std::cout << std::endl;
	std::cout << "X	" << result.bestScore << "	" << diff1.count() << "	" <<  (double)len1 / diff1.count() << std::endl;

	//======================================================================================
	// KSW2 GLOBAL AND EXTENSION (vectorized SSE4.1 and x-drop, not banded)
	//======================================================================================

#ifdef KSW2

	int8_t a = MAT, b = MIS < 0? MIS : -MIS; // a>0 and b<0
	int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
	int tl = strlen(targetSeg.c_str()), ql = strlen(querySeg.c_str());
	uint8_t *ts, *qs, c[256];
	ksw_extz_t ez;

	memset(&ez, 0, sizeof(ksw_extz_t));
	memset(c, 4, 256);

	// build the encoding table
	c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
	c['G'] = c['g'] = 2; c['T'] = c['t'] = 3;
	ts = (uint8_t*)malloc(tl);
	qs = (uint8_t*)malloc(ql);

	// encode to 0/1/2/3
	for (int i = 0; i < tl; ++i)
	{
		ts[i] = c[(uint8_t)targetSeg[i]];
	}
	for (int i = 0; i < ql; ++i)
	{
		qs[i] = c[(uint8_t)querySeg[i]];
	}

	std::chrono::duration<double> diff2;
	auto start2 = std::chrono::high_resolution_clock::now();

	ksw_extz2_sse(0, ql, qs, tl, ts, 5, mat, 0, -GAP, bw, xdrop, 0, KSW_EZ_SCORE_ONLY, &ez);

	auto end2 = std::chrono::high_resolution_clock::now();
	diff2 = end2-start2;

	free(ts); free(qs);

	std::cout << "K	" << ez.score << "	" << diff2.count() << "	" << (double)len1 / diff2.count() << std::endl;
#endif

	//======================================================================================
	// LIBGABA SEMI-GLOBAL ALIGNMENT (SegFault, vetorized SSE2 and AVX2, x-drop, banded)
	//======================================================================================

#ifdef GABA
	gaba_t *ctx = gaba_init(GABA_PARAMS(
		// match award, mismatch penalty, gap open penalty (G_i), and gap extension penalty (G_e)
		GABA_SCORE_SIMPLE(1, 1, 0, 1),
		gfa : 0,
		gfb : 0,
		xdrop : (int8_t)xdrop,
		filter_thresh : 0,
	));

	std::chrono::duration<double> diff3;
	auto start3 = std::chrono::high_resolution_clock::now();

	char const t[64] = {0};	// tail array
	gaba_section_t asec = gaba_build_section(0, (uint8_t const *)targetSeg.c_str(), (uint32_t)targetSeg.size());
	gaba_section_t bsec = gaba_build_section(2, (uint8_t const *)querySeg.c_str(),  (uint32_t)querySeg.size());
	gaba_section_t tail = gaba_build_section(4, t, 64);

	// create thread-local object
	gaba_dp_t *dp = gaba_dp_init(ctx);	// dp[0] holds a 64-cell-wide context

	// init section pointers
	gaba_section_t const *ap = &asec, *bp = &bsec;
	gaba_fill_t const *f = gaba_dp_fill_root(dp, // dp -> &dp[_dp_ctx_index(band_width)] makes the band width selectable
		ap, 0,					// a-side (reference side) sequence and start position
		bp, 0,					// b-side (query)
		UINT32_MAX				// max extension length
	);

	// until x-drop condition is detected
	gaba_fill_t const *m = f;

	while((f->status & GABA_TERM) == 0) {
		// substitute the pointer by the tail section's if it reached the end
		if(f->status & GABA_UPDATE_A) { ap = &tail; }
		if(f->status & GABA_UPDATE_B) { bp = &tail; }

		f = gaba_dp_fill(dp, f, ap, bp, UINT32_MAX);	// extend the banded matrix
		m = f->max > m->max ? f : m;					// swap if maximum score was updated
	}

	// alignment path
	gaba_alignment_t *r = gaba_dp_trace(dp,
		m,		// section with the max
		NULL	// custom allocator: see struct gaba_alloc_s in gaba.h
	);

	auto end3 = std::chrono::high_resolution_clock::now();
	diff3 = end3-start3;

	std::cout << "G	" << r->score << "	" << diff3.count() << "	" << (double)len1 / diff3.count() << std::endl;

	// clean up
	gaba_dp_res_free(dp, r); gaba_dp_clean(dp);
	gaba_clean(ctx);

#endif

	//======================================================================================
	// SSW LOCAL ALIGNMENT (SSE2 vectorized, not banded)
	//======================================================================================

#ifdef SSW
	int32_t maskLen = strlen(querySeg.c_str())/2;
	maskLen = maskLen < 15 ? 15 : maskLen;

	StripedSmithWaterman::Aligner aligner(MAT, MIS, 0, GAP);
	StripedSmithWaterman::Alignment alignment;

	StripedSmithWaterman::Filter filter;
	filter.report_begin_position = false;
	filter.report_cigar = false;

	std::chrono::duration<double> diff5;
	auto start5 = std::chrono::high_resolution_clock::now();

	aligner.Align(querySeg.c_str(), targetSeg.c_str(), targetSeg.size(), filter, &alignment, maskLen);

	auto end5 = std::chrono::high_resolution_clock::now();
	diff5 = end5-start5;

	std::cout << "S	" << alignment.sw_score << "	" << diff5.count() << "	" << (double)len1 / diff5.count() << std::endl;
#endif

	//======================================================================================
	// EDLIB GLOBAL ALIGNMENT
	//======================================================================================

#ifdef EDLIB
	std::chrono::duration<double> diff7;
	auto start7 = std::chrono::high_resolution_clock::now();

	EdlibAlignResult edresult = edlibAlign(targetSeg.c_str(), targetSeg.size(), querySeg.c_str(), querySeg.size(),
		edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
	auto end7 = std::chrono::high_resolution_clock::now();
	diff7 = end7-start7;

	std::cout << "E	" << edresult.editDistance << "	" << diff7.count() << "	" << (double)len1 / diff7.count() << std::endl;
	edlibFreeAlignResult(edresult);
#endif

	//======================================================================================
	// PARASAIL GLOBAL ALIGNMENT (banded, not vectorized, SegFault)
	//======================================================================================

#ifdef PARASAIL1
	int s1Len = (int)strlen(targetSeg.c_str());
	int s2Len = (int)strlen(querySeg.c_str());
	parasail_result_t *result1 = NULL;

	const char alphabet[4] = {'A','T','C', 'G'};
	const parasail_matrix_t *matrix = parasail_matrix_create(alphabet, MAT, MIS);

	std::chrono::duration<double> diff6;
	auto start6 = std::chrono::high_resolution_clock::now();

	result1 = parasail_sg(targetSeg.c_str(), s1Len, querySeg.c_str(), s2Len, 0, 1, matrix);

	auto end6 = std::chrono::high_resolution_clock::now();
	diff6 = end6-start6;

	std::cout << "P2	" << parasail_result_get_score(result1) << "	" << diff6.count() << "	" << (double)len1 / diff6.count() << std::endl;
	std::cout << std::endl;
	parasail_result_free(result1);
#endif

	//======================================================================================
	// PARASAIL GLOBAL ALIGNMENT (not banded, vectorized)
	//======================================================================================

#ifdef PARASAIL2
	parasail_result_t *result2 = NULL;

	std::chrono::duration<double> diff8;
	auto start8 = std::chrono::high_resolution_clock::now();

	result2 = parasail_nw_banded(targetSeg.c_str(), s1Len, querySeg.c_str(), s2Len, 0, 1, bw, matrix); // check if properly triggered intrinsics

	auto end8 = std::chrono::high_resolution_clock::now();
	diff8 = end8-start8;

	std::cout << "P2	" << parasail_result_get_score(result2) << "	" << diff8.count() << "	" << (double)len1 / diff8.count() << std::endl;
	std::cout << std::endl;
	parasail_result_free(result2);
#endif
	return 0;
}