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
#include "softwareRuntime.h"

#define MAT ( 1)
#define MIS (-1)
#define GAP (-1)
#define BWD (128)
#define KML (15)
#define PMIS (0.05)
#define PGAP (0.10)

int main(int argc, const char *argv[])
{
	// Command line input
	int x = std::stoi(argv[1]);
	std::vector<std::tuple<int, double, double, double, double>> runtime;

	std::cout << "length\tseqan\txavier\tlibgaba\tksw2" << std::endl;

	int init  = 2500;
	for(int i = init; i < (init*20+1); i+init)
	{
		int tid = omp_get_thread_num();

		std::string seq1;
		std::string seq2;

		generate_random_sequence(seq1, i);
		seq2 = generate_mutated_sequence(seq1, i, PMIS, PGAP, BWD);
		seq1 = seq2;

		double stime = seqanAlign (MAT, MIS, GAP, KML, x, seq1, seq2);
		double xtime = xavireAlign(MAT, MIS, GAP, KML, x, seq1, seq2);
		double gtime = gabaAlign  (MAT, MIS, GAP, KML, x, seq1, seq2);
		double ktime = ksw2Align  (MAT, MIS, GAP, KML, x, seq1, seq2, std::stoi(argv[2]));

		std::cout << seq1.length() << "\t" << stime << "\t" << xtime
			<< "\t" << gtime << "\t" << ktime << std::endl;

		i = i+init;
	}

	return 0;
}