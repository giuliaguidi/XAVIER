//==============================================================================
// Title:   Xavier: High-Performance X-Drop Adaptive Banded Pairwise Alignment
// Author:  G. Guidi, E. Younis
// Date:    7 February 2020
// Test:    
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
#include <bits/stdc++.h> 
#include <omp.h>

#include "utils.h"
#include "software2.h"

int main(int argc, char const *argv[])
{
	// variables declaration
	int xdrop = std::stoi(argv[1]);
	int bw	  = std::stoi(argv[2]);

	double pmis = std::stoi(argv[3]);
	double pgap = std::stoi(argv[4]);

	int mat = std::stoi(argv[5]);
	int mis = std::stoi(argv[6]);
	int gap = std::stoi(argv[7]);

	int k   = std::stoi(argv[8]);
	int min = std::stoi(argv[9]);
	int max = std::stoi(argv[10]);

	int maxt = 1;
	#pragma omp parallel 
	{
		maxt = omp_get_num_threads();
	}

	std::default_random_engine generator;

	std::vector<std::tuple<int, double, double, double, double>> runtime(max-min);
	std::vector<std::vector<std::tuple<int, double, double, double, double>>> vruntime(maxt);

#pragma omp parallel for
	for(int i = min; i < max; i++)
	{
		int tid = omp_get_thread_num();

		std::string randomSeg;
		generate_random_sequence(randomSeg, i);

		std::string querySeg  = generate_mutated_sequence(randomSeg, i, pmis, pgap, bw);
		std::string targetSeg = generate_mutated_sequence(randomSeg, i, pmis, pgap, bw);

		// aligner
		double stime = seqanAlign (mat, mis, gap, k, xdrop, targetSeg, querySeg);
		double xtime = xavireAlign(mat, mis, gap, k, xdrop, targetSeg, querySeg);
		double gtime = gabaAlign  (mat, mis, gap, k, xdrop, targetSeg, querySeg);
		double ktime = ksw2Align  (mat, mis, gap, k, xdrop, targetSeg, querySeg, bw);

		vruntime[tid].push_back(std::make_tuple(i, stime, xtime, ktime, gtime));
	}

	unsigned int sofar = 0;

	for(int t = 0; t < maxt; ++t)
	{
		std::copy(vruntime[t].begin(), vruntime[t].end(), runtime.begin() + sofar);
		sofar += vruntime[t].size();
	}
	
	for(int i = 0; i < max-min; i++)
	{
		std::cout << std::get<0>(runtime[i]) << "\t" << std::get<1>(runtime[i]) << "\t" << std::get<2>(runtime[i]) 
			<< "\t" << std::get<3>(runtime[i]) << "\t" << std::get<4>(runtime[i]) << std::endl;
	}

	return 0;
}