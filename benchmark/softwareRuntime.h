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

#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/seeds.h>
#include <seqan/score.h>
#include <seqan/modifier.h>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/seeds/seeds_extension.h>

#include "utils.h"

#include "xavier.h"
#include "ksw2/ksw2.h"
#include "ksw2/ksw2_extz2_sse.c"

#ifdef __cplusplus
extern "C" {
#endif
#include "libgaba/gaba.h" 	 // sometimes the forefront vector will not reach the end
							 // of the sequences. It is more likely to occur when the input
							 // sequence lengths greatly differ
#ifdef __cplusplus
}
#endif

double seqanAlign(int mat, int mis, int gap, int k, int xdrop, 
    std::string& targetSeg, std::string& querySeg)
{

	seqan::Score<int, seqan::Simple> scoringSchemeSeqAn(mat, mis, gap);
	seqan::Seed<seqan::Simple> seed(0, 0, k);

	std::chrono::duration<double> diff4;
	auto start4 = std::chrono::high_resolution_clock::now();
	int score = seqan::extendSeed(seed, targetSeg, querySeg, seqan::EXTEND_RIGHT, scoringSchemeSeqAn, xdrop, seqan::GappedXDrop(), k);
	auto end4 = std::chrono::high_resolution_clock::now();
	diff4 = end4-start4;

	return diff4.count();
}

double xavireAlign(int mat, int mis, int gap, int k, int xdrop, 
    std::string& targetSeg, std::string& querySeg)
{
	// init scoring scheme
	xavier::ScoringScheme penalties (mat, mis, gap);
	// seed starting position on seq1, seed starting position on seq2, k-mer length
	xavier::Seed seed(0, 0, k);

	std::chrono::duration<double> diff1;
	auto start1 = std::chrono::high_resolution_clock::now();

	xavier::AlignmentResult result = xavier::seed_and_extend_right(targetSeg, querySeg, penalties, xdrop, seed);

	auto end1 = std::chrono::high_resolution_clock::now();
	diff1 = end1-start1;

    return diff1.count();
}

double ksw2Align(int mat, int mis, int gap, int k, int xdrop, 
    std::string& targetSeg, std::string& querySeg, int bw)
{
    // init
	int8_t a = mat, b = mis < 0? mis : -mis; // a>0 and b<0
	int8_t matrix[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
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

	// CHECK: KSW_EZ_SCORE_ONLY
	std::chrono::duration<double> diff2;
	auto start2 = std::chrono::high_resolution_clock::now();

  ksw_extz2_sse(0, ql, qs, tl, ts, 5, matrix, 0, -gap, bw, xdrop, 0, 0, &ez);
//   ksw_extz2_sse(0, ql, qs, tl, ts, 5, matrix, 0, -gap, bw, xdrop, 0, KSW_EZ_SCORE_ONLY, &ez);

	auto end2 = std::chrono::high_resolution_clock::now();
	diff2 = end2-start2;

	free(ts); free(qs);

    return diff2.count();
}

double gabaAlign(int mat, int mis, int gap, int k, int xdrop, 
    std::string& targetSeg, std::string& querySeg)
{
gaba_t *ctx = gaba_init(GABA_PARAMS(
		// match award, mismatch penalty, gap open penalty (G_i), and gap extension penalty (G_e)
		GABA_SCORE_SIMPLE(mat, -mis, 0, -gap),
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

	// clean up
	gaba_dp_res_free(dp, r); gaba_dp_clean(dp);
	gaba_clean(ctx);    

    return diff3.count();
}