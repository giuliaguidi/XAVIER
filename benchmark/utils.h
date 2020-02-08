#ifndef UTILS_H
#define UTILS_H

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

/* Functions to generate sequences from https://github.com/ocxtal/libgaba */
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

std::string generate_mutated_sequence(const std::string& seq, int& len2, double& pmis, double& pgap, int& bw)
{
	int i, j, wave = 0;	// wave is q-coordinate of the alignment path
	std::string mutated;

	for(i = 0, j = 0; i < len2; i++)
	{
		if(((double)rand()/(double)RAND_MAX) < pmis)
		{
			mutated.append(1, random_base());
			j++; // mismatch
		}
		else if(((double)rand()/(double)RAND_MAX) < pgap)
		{
			if(rand() & 0x01 && wave > (-bw + 1))
			{
				char tmp = (j < len2) ? seq[j++] : random_base();
				mutated.append(1, tmp);
				j++; wave--; // deletion
			}
			else if(wave < (bw-2))
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

#endif