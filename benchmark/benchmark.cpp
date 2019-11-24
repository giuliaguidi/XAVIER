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

#include "xavier.h"
#include "ksw2/ksw2.h"
#include "ksw2/ksw2_extz2_sse.c"

// #include "parasail/parasail.h"
// #include "parasail/parasail/io.h"
// #include "parasail/parasail/memory.h"
// #include "parasail/parasail/stats.h"
// #include "Complete-Striped-Smith-Waterman-Library/src/ssw_cpp.h"
// #include <edlib.h>

// #include <seqan/align.h>
// #include <seqan/sequence.h>
// #include <seqan/align.h>
// #include <seqan/seeds.h>
// #include <seqan/score.h>
// #include <seqan/modifier.h>
// #include <seqan/basic.h>
// #include <seqan/stream.h>

// #ifdef __cplusplus
// extern "C" {
// #endif
// #include "libgaba/gaba.h" 		 // sometimes the forefront vector will not reach the end
// 								 // of the sequences. It is more likely to occur when the input
// 								 // sequence lengths greatly differ
// #ifdef __cplusplus
// }
// #endif

//======================================================================================
// GLOBAL VARIABLE DECLARATION
//======================================================================================

#define LEN1 	(10000)		// read length (this is going to be a distribution of length in the adaptive version)
#define LEN2 	(10000)		// 2nd read length
#define MAT		( 1)		// match score
#define MIS		(-1)		// mismatch score
#define GAP		(-1)		// gap score
#define PMIS 	(0.03)		// substitution probability
#define PGAP 	(0.13)		// insertion/deletion probability
#define BW 		(128)		// bandwidth (the alignment path of the input sequence and the result does not go out of the band)

#define KSW2
// #define GABA
// #define NOSIMD
// #define SEQAN
// #define SSW
// #define PARASAIL1
// #define PARASAIL2
// #define EDLIB

//======================================================================================
// READ SIMULATOR
//======================================================================================

// Functions to generate sequences from https://github.com/ocxtal/libgaba
char random_base(void)
{
	char const table[4] = {'A', 'C', 'G', 'T'};
	return(table[rand() % 4]);
}

void generate_random_sequence(std::string& seq)
{
	for(int i = 0; i < LEN1; i++)
		seq.append(1, random_base());
}

std::string generate_mutated_sequence(const std::string& seq)
{
	int i, j, wave = 0;	// wave is q-coordinate of the alignment path
	std::string mutated;

	for(i = 0, j = 0; i < LEN2; i++)
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
				char tmp = (j < LEN2) ? seq[j++] : random_base();
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
				char tmp = (j < LEN2) ? seq[j++] : random_base();
				mutated.append(1, tmp);
			}
		}
		else
		{
			char tmp = (j < LEN2) ? seq[j++] : random_base();
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
	generate_random_sequence(targetSeg);
	std::string querySeg = generate_mutated_sequence(targetSeg);

	int xdrop = std::stoi(argv[1]);

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
	std::cout << "result.bestScore	" << result.bestScore << std::endl;
	std::cout << "result.exitScore	" << result.exitScore << std::endl;
	std::cout << "result.begH	" << result.begH << std::endl;
	std::cout << "result.endH	" << result.endH << std::endl;
	std::cout << "result.begV	" << result.begV << std::endl;
	std::cout << "result.endV	" << result.endV << std::endl;

	std::cout << "time  " << diff1.count() << "\t" << (double)LEN1 / diff1.count() << "\tbases aligned per second" << std::endl;

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

	ksw_extz2_sse(0, ql, qs, tl, ts, 5, mat, 0, -GAP, 64, xdrop, 0, KSW_EZ_SCORE_ONLY, &ez);

	auto end2 = std::chrono::high_resolution_clock::now();
	diff2 = end2-start2;

	free(ts); free(qs);

	std::cout << std::endl;
	std::cout << "result.bestScore	" << ez.score << std::endl;
	std::cout << "time  " << diff2.count() << "\t" << (double)LEN1 / diff2.count() << "\tbases aligned per second" << std::endl;
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
	// x-drop has to be within [-127, 127] in libagaba
	gaba_fill_t const *m = f;
	// track max

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

	std::cout << xdrop << "\t" << r->score << "\t" << diff3.count() << "\t" << (double)LEN1 / diff3.count() << "\tbases aligned per second" << std::endl;

	// clean up
	gaba_dp_res_free(dp, r); gaba_dp_clean(dp);
	gaba_clean(ctx);

#endif

	//======================================================================================
	// SEQAN SEED EXTENSION (not vectorized, not banded, x-drop)
	//======================================================================================

#ifdef NOSIMD
	// SeqAn
	std::cout << "SeqAn" << std::endl;
	seqan::Score<int, seqan::Simple> scoringSchemeSeqAn(MAT, MIS, GAP);
	seqan::Seed<seqan::Simple> seed1(0, 0, 0);
	std::chrono::duration<double> diff4;
	auto start4 = std::chrono::high_resolution_clock::now();
	int score = seqan::extendSeed(seed1, targetSeg, querySeg, seqan::EXTEND_RIGHT,
		scoringSchemeSeqAn, xdrop, seqan::GappedXDrop(), 0);
	auto end4 = std::chrono::high_resolution_clock::now();
	diff4 = end4-start4;

	std::cout << xdrop << "\t" << score << "\t" << diff4.count() << "\t" << (double)LEN1 / diff4.count() << "\tbases aligned per second" << std::endl;
#endif
	std::cout << std::endl;
	//======================================================================================
	// SEQAN BANDED GLOBAL (vectorized, banded)
	//======================================================================================

#ifdef SEQAN
	seqan::Score<int16_t, seqan::Simple> scoringSchemeSeqAn(MAT, MIS, GAP);
	using TSequence    = seqan::String<seqan::Dna>;
	using TThreadModel = seqan::WavefrontAlignment<seqan::BlockOffsetOptimization>;
	using TVectorSpec  = seqan::Vectorial;
	using TExecPolicy  = seqan::ExecutionPolicy<TThreadModel, TVectorSpec>;

	seqan::StringSet<TSequence> seqs1;
	seqan::StringSet<TSequence> seqs2;

	appendValue(seqs1, TSequence{targetSeg.c_str()});
	appendValue(seqs2, TSequence{querySeg.c_str()});

	TExecPolicy execPolicy;
	setNumThreads(execPolicy, 1);
	std::chrono::duration<double> diff9;
	auto start9 = std::chrono::high_resolution_clock::now();
	seqan::String<int16_t> scores = seqan::globalAlignmentScore(execPolicy, seqs1, seqs2, scoringSchemeSeqAn);
	auto end9 = std::chrono::high_resolution_clock::now();
	diff9 = end9-start9;

	std::cout << "SeqAn's best (not banded, vectorized) " << scores[0] << " in " << diff9.count() << " sec " << std::endl;
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

	// SeqAn is doing more computation
	std::cout << "SSW's best (not banded, vectorized) " << alignment.sw_score << " in " << diff5.count() << " sec " << std::endl;
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

	std::cout << "Edlib's edit distance (banded, not vectorized) " << edresult.editDistance << " in " << diff7.count() << " sec " << std::endl;
	edlibFreeAlignResult(edresult);
#endif

	//======================================================================================
	// PARASAIL GLOBAL ALIGNMENT (banded, not vectorized, SegFault)
	//======================================================================================

#ifdef PARASAIL1
	int s1Len = (int)strlen(targetSeg.c_str());
	int s2Len = (int)strlen(querySeg.c_str());
	parasail_result_t *result = NULL;

	const char alphabet[4] = {'A','T','C', 'G'};
	const parasail_matrix_t *matrix = parasail_matrix_create(alphabet, MAT, MIS);

	std::chrono::duration<double> diff6;
	auto start6 = std::chrono::high_resolution_clock::now();

	result = parasail_nw_banded(targetSeg.c_str(), s1Len, querySeg.c_str(), s2Len, 0, 1, 16, matrix);
	parasail_result_free(result);

	auto end6 = std::chrono::high_resolution_clock::now();
	diff6 = end6-start6;

	// SeqAn is doing more computation
	std::cout << "Parasail's best (banded, not vectorized) " << parasail_result_get_score(result) << " in " << diff6.count() << " sec " << std::endl;
#endif

	//======================================================================================
	// PARASAIL GLOBAL ALIGNMENT (not banded, vectorized)
	//======================================================================================

#ifdef PARASAIL2
	parasail_result_t *result2 = NULL;

	std::chrono::duration<double> diff8;
	auto start8 = std::chrono::high_resolution_clock::now();

	result2 = parasail_nw(targetSeg.c_str(), s1Len, querySeg.c_str(), s2Len, 0, 1, matrix); // check if properly triggered intrinsics
	parasail_result_free(result2);

	auto end8 = std::chrono::high_resolution_clock::now();
	diff8 = end8-start8;

	// SeqAn is doing more computation
	std::cout << "Parasail's best (not banded, vectorized) " << parasail_result_get_score(result2) << " in " << diff8.count() << " sec " << std::endl;
#endif
	return 0;
}