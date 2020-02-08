//==============================================================================
// Title:   Xavier: High-Performance xy[i].first-Drop Adaptive Banded Pairwise Alignment
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

#include "utils.h"
#include "software.h"

// pearson correlation coefficient
float pc(std::vector<std::pair<int, int>> xy) 
{ 
	int n = xy.length();

    int sumx = 0, sumy = 0, sumxy = 0; 
    int sum2x = 0, sum2y = 0; 
  
    for (int i = 0; i < n; i++) 
    { 
        sumx 	= sumx 	+ xy[i].first; 
        sumy 	= sumy 	+ xy[i].second; 
        sumxy 	= sumxy + xy[i].first * xy[i].second; 
  
        sum2x 	= sum2x + xy[i].first 	* xy[i].first; 
        sum2y 	= sum2y + xy[i].second 	* xy[i].second; 
    } 
  
    float pcorr = (float)(n * sumxy - sumx * sumy) / std::sqrt((n * sum2x - sumx * sumx) * (n * sum2y - sumy * sumy)); 

    return pcorr; 
}

int main(int argc, char const *argv[])
{
	// variables declaration
	int xdrop = std::stoi(argv[1]);
	int bw	  = std::stoi(argv[2]);
	int k	  = std::stoi(argv[8]);
	int iter  = std::stoi(argv[9]);

	int mat = std::stoi(argv[5]);
	int mis = std::stoi(argv[6]);
	int gap = std::stoi(argv[7]);

	double pmis = std::stoi(argv[3]);
	double pgap = std::stoi(argv[4]);

	std::default_random_engine generator;
	std::vector<std::pair<int, int>> xpc, kpc, gpc;
	std::vector<std::vector<std::pair<int, int>>> vxpc, vkpc, vgpc;

#pragma omp parallel for
	for(int i = 0; i < iter; i++)
	{
		int tid = omp_get_num_threds();

		std::string randomSeg;
		std::normal_distribution<float> distribution(1200.0, 3000.0);

		int len1 = (int)distribution(generator);
		int len2 = (int)distribution(generator);
		int len3 = (int)distribution(generator);

		generate_random_sequence(randomSeg, len1);
		std::string querySeg  = generate_mutated_sequence(randomSeg, len2, pmis, pgap, bw);
		std::string targetSeg = generate_mutated_sequence(randomSeg, len3, pmis, pgap, bw);

		// reference
		// int rscore = truthAlign(mat, mis, gap, k, xdrop, targetSeg, querySeg);

		// aligner
		int xscore = xavireAlign(mat, mis, gap, k, xdrop, targetSeg, querySeg);
		int gscore = gabaAlign  (mat, mis, gap, k, xdrop, targetSeg, querySeg);
		int kscore = ksw2Align  (mat, mis, gap, k, xdrop, targetSeg, querySeg, bw);

		vxpc[tid].push_back(std::make_pair(rscore, xscore));
		vkpc[tid].push_back(std::make_pair(rscore, kscore));
		vgpc[tid].push_back(std::make_pair(rscore, gscore));
	}

	unsigned int xpcsofar = 0;
	unsigned int kpcsofar = 0;
	unsigned int gpcsofar = 0;

	for(int t = 0; t < MAXTHREADS; ++t)
	{
		copy(vxpc[t].begin(), vxpc[t].end(), xpc.begin() + xpcsofar);
		xpcsofar += vxpc[t].size();

		copy(vkpc[t].begin(), vkpc[t].end(), kpc.begin() + kpcsofar);
		kpcsofar += vkpc[t].size();

		copy(vgpc[t].begin(), vgpc[t].end(), gpc.begin() + gpcsofar);
		gpcsofar += vgpc[t].size();
	}

	for(int i = 0; i < iter; i++)
		std::cout << xpc[i].first << "\t" << xpc[i].second << std::endl;
	float pcx = pc(xpc);
	std::cout << pcx << "\n" << std::endl;

	for(int i = 0; i < iter; i++)
		std::cout << kpc[i].first << "\t" << kpc[i].second << std::endl;
	float pck = pc(kpc);
	std::cout << pck << "\n" << std::endl;

	for(int i = 0; i < iter; i++)
		std::cout << gpc[i].first << "\t" << gpc[i].second << std::endl;
	float pcg = pc(gpc);
	std::cout << pcg << "\n" << std::endl;

	return 0;
}