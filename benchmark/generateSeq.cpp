//==============================================================================
// Title:   Xavier: High-Performance xy[i].first-Drop Adaptive Banded Pairwise Alignment
// Author:  G. Guidi, E. Younis
// Date:    10 February 2020
// Test:    Generate sequences for test purposes
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

#define BWD (256)

int main(int argc, char const *argv[])
{
    // TODO: move to optlist
    int n = std::stoi(argv[1]);

	double pmis = std::stod(argv[2]);
	double pgap = std::stod(argv[3]);

    bool random       = std::stoi(argv[4]); 
	bool constLenPair = std::stoi(argv[5]); 

    int mean = std::stoi(argv[6]);
    int sd   = std::stoi(argv[7]);

	std::default_random_engine generator;
	std::vector<std::pair<std::string, std::string>> seqs;

	for(int i = 0; i < n; i++)
	{
		std::normal_distribution<double> distribution(mean, sd);
	
		int len1 = (int)distribution(generator);
		int len2 = (int)distribution(generator);
        int len3 = (int)distribution(generator);

        std::string randomSeg;
		std::string targetSeg;
        std::string querySeg;
    
        if(random)
        {
			if(constLenPair)
			{
			    generate_random_sequence(querySeg , len1);
		    	generate_random_sequence(targetSeg, len1);			
			}
			else
			{	
		   		generate_random_sequence(querySeg , len1);
		    	generate_random_sequence(targetSeg, len2);
			}

        }
        else
        {
            generate_random_sequence(randomSeg, len1);

			if(constLenPair)
			{
		    	querySeg  = generate_mutated_sequence(randomSeg, len2, pmis, pgap, BWD);
		    	targetSeg = generate_mutated_sequence(randomSeg, len2, pmis, pgap, BWD);
			}
			else
			{
			    querySeg  = generate_mutated_sequence(randomSeg, len2, pmis, pgap, BWD);
		    	targetSeg = generate_mutated_sequence(randomSeg, len3, pmis, pgap, BWD);			
			}
        }
        
		seqs.push_back(std::make_pair(querySeg, targetSeg));
	}

    std::string name1 = "S1BWD" + std::to_string(BWD) + "N" + std::to_string(n) + "PMIS" 
							+ std::to_string(pmis) + "PGAP" + std::to_string(pgap) + "R" + std::to_string(random) 
								+ "C" + std::to_string(constLenPair) + "AV" + std::to_string(mean) + "SD" + std::to_string(sd) + ".seq";  

    std::string name2 = "S2BWD" + std::to_string(BWD) + "N" + std::to_string(n) + "PMIS" 
							+ std::to_string(pmis) + "PGAP" + std::to_string(pgap) + "R" + std::to_string(random) 
								+ "C" + std::to_string(constLenPair) + "AV" + std::to_string(mean) + "SD" + std::to_string(sd) + ".seq"; 

    std::ofstream ofs1(name1, std::ofstream::out);
    std::ofstream ofs2(name2, std::ofstream::out);

  	for(int i = 0; i < n; ++i)
	{
		ofs1 << seqs[i].first  << std::endl;
        ofs2 << seqs[i].second << std::endl;
	}  

    ofs1.close();
    ofs2.close();

	return 0;
}