/**
 * File: trace.h
 * Author: G. Guidi, E. Younis
 * Description: Xavier Trace Type Header.
 *
 * Xavier: High-Performance X-Drop Adaptive Banded Pairwise Alignment (Xavier)
 * Copyright (c) 2019, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * (1) Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * (2) Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * (3) Neither the name of the University of California, Lawrence Berkeley
 * National Laboratory, U.S. Dept. of Energy nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * You are under no obligation whatsoever to provide any bug fixes, patches,
 * or upgrades to the features, functionality or performance of the source
 * code ("Enhancements") to anyone; however, if you choose to make your
 * Enhancements available either publicly, or directly to Lawrence Berkeley
 * National Laboratory, without imposing a separate written license agreement
 * for such Enhancements, then you hereby grant the following license: a
 * non-exclusive, royalty-free perpetual license to install, use, modify,
 * prepare derivative works, incorporate into other computer software,
 * distribute, and sublicense such enhancements or derivative works thereof,
 * in binary and source code form.
 */

#ifndef XAVIER_TYPES_TRACE_H
#define XAVIER_TYPES_TRACE_H

#include <list>
#include <vector>
#include "vectors.h"
#include "score.h"
#include <algorithm>

namespace xavier
{

	class TraceEntry
	{
	public:
		TraceEntry():
		antiDiag1(0),
		antiDiag2(0),
		antiDiag3(0),
		vqueryh(0),
		vqueryv(0),
		scoreOffset(0),
		lastMove(0)
		{}

		/**
		 * Constructors
		 */
		TraceEntry ( const VectorRegister& _ad1, const VectorRegister& _ad2,
		             const VectorRegister& _ad3, const VectorRegister& _vqh,
		             const VectorRegister& _vqv, const int64_t offset,
		             const int _lastMove );

		TraceEntry ( const TraceEntry& copy );

		VectorRegister antiDiag1;
		VectorRegister antiDiag2;
		VectorRegister antiDiag3;

		VectorRegister vqueryh;
		VectorRegister vqueryv;
		int64_t 	   scoreOffset;
		int            lastMove;
	};

	class Trace
	{
	public:

		struct AlignmentPair
		{
			// std::string alignH;
			std::string cigar;
			size_t matches;
		};

		/**
		 * Default Constructor
		 */
		Trace( const ScoringScheme& score ):
		scoringScheme( score ) {}

		/**
		 * Store another trace entry corresponding to the current state.
		 */
		void pushbackState ( const VectorRegister& _ad1, const VectorRegister& _ad2,
		                     const VectorRegister& _ad3, const VectorRegister& _vqh,
		                     const VectorRegister& _vqv, const int64_t offset,
		                     const int _lastMove );

		/**
		 * Backtrace Algorithm
		 */
		AlignmentPair getAlignment();

		/**
		 * Record Keeping Functions
		 */
		void recordGlobalMaxPos();
		void saveOpeningPhaseDPMatrix ( std::vector< std::vector<int> > _DPMatrix, int8_t* _queryh, int8_t* _queryv );
		std::string compression(const std::string& str);

	private:
		std::vector<TraceEntry> trace;
		ScoringScheme scoringScheme;
		std::vector< std::vector<int> > DPMatrix;
		int8_t* queryh;
		int8_t* queryv;
		size_t  maxPos;
	};
}

#endif /* XAVIER_TYPES_TRACE_H */