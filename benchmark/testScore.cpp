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
#include <random>
#include <omp.h>

#include "utils.h"
#include "softwareScore.h"

#define MAT ( 1)
#define MIS (-1)
#define GAP (-1)
#define BWD (64)
#define KML (15)

int main(int argc, const char *argv[])
{
	// Command line input
	int x = std::stoi(argv[1]);

	std::ifstream seqs1(argv[2]);
	std::ifstream seqs2(argv[3]);

	// Set up OMP thread
	int maxt = 1;
	#pragma omp parallel 
	{
		maxt = omp_get_num_threads();
	}

    int n = std::count(std::istreambuf_iterator<char>(seqs1), std::istreambuf_iterator<char>(), '\n');
    seqs1.seekg(0, std::ios_base::beg);

    std::vector<std::pair<std::string, std::string>> entries;
    std::vector<std::stringstream> local(maxt);  

	// Read seqs1 and seqs2
    if(seqs1 && seqs2)
        for (int i = 0; i < n; ++i)
        {
            std::string seq1, seq2;
            std::getline(seqs1, seq1);
			std::getline(seqs2, seq2);
            entries.push_back(std::make_pair(seq1, seq2));
        }
    seqs1.close(); 
	seqs2.close(); 

	std::vector<std::tuple<int, int, int, int>> scores(n);
	std::vector<std::vector<std::tuple<int, int, int, int>>> vscores(maxt);

#pragma omp parallel for
	for(int i = 0; i < n; i++)
	{
		int tid = omp_get_thread_num();

		std::string seq1 = entries[i].first;
		std::string seq2 = entries[i].second;

		int sscore = seqanAlign (MAT, MIS, GAP, KML, x, seq1, seq2);
		int xscore = xavireAlign(MAT, MIS, GAP, KML, x, seq1, seq2);
		int gscore = gabaAlign  (MAT, MIS, GAP, KML, x, seq1, seq2);
		int kscore = ksw2Align  (MAT, MIS, GAP, KML, x, seq1, seq2, BWD);

		vscores[tid].push_back(std::make_tuple(sscore, xscore, gscore, kscore));
	}

	unsigned int sofar = 0;

	for(int t = 0; t < maxt; ++t)
	{
		std::copy(vscores[t].begin(), vscores[t].end(), scores.begin() + sofar);
		sofar += vscores[t].size();
	}
	std::cout << "seqan\txavier\tlibgaba\tksw2" << std::endl;
	for(int i = 0; i < n; i++)
	{
		std::cout << std::get<0>(scores[i]) << "\t" << std::get<1>(scores[i]) << "\t" << std::get<2>(scores[i]) 
			<< "\t" << std::get<3>(scores[i]) << std::endl;
	}

	return 0;
}