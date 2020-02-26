/**
 * File: trace.cpp
 * Author: G. Guidi, E. Younis
 * Description: Xavier Trace Type Source.
 */

#include <iostream>
#include "trace.h"

namespace xavier
{

	TraceEntry::TraceEntry ( const VectorRegister& _ad1, const VectorRegister& _ad2,
	             			 const VectorRegister& _ad3, const VectorRegister& _vqh,
	             			 const VectorRegister& _vqv, const int64_t offset,
	             			 const int _lastMove ):
		antiDiag1( _ad1 ),
		antiDiag2( _ad2 ),
		antiDiag3( _ad3 ),
		vqueryh( _vqh ),
		vqueryv( _vqv ),
		scoreOffset( offset ),
		lastMove( _lastMove )
		{ }

	TraceEntry::TraceEntry ( const TraceEntry& copy ):
		antiDiag1( copy.antiDiag1 ),
		antiDiag2( copy.antiDiag2 ),
		antiDiag3( copy.antiDiag3 ),
		vqueryh( copy.vqueryh ),
		vqueryv( copy.vqueryv ),
		scoreOffset( copy.scoreOffset ),
		lastMove( copy.lastMove )
		{ }

	void Trace::pushbackState ( const VectorRegister& _ad1, const VectorRegister& _ad2,
	   	                        const VectorRegister& _ad3, const VectorRegister& _vqh,
		                        const VectorRegister& _vqv, const int64_t offset,
		                        const int _lastMove )
	{
		trace.emplace_back( _ad1, _ad2, _ad3, _vqh, _vqv, offset, _lastMove );
	}

	void Trace::recordGlobalMaxPos()
	{
		maxPos = trace.size() - 1;
	}

	// GG: back trace
	Trace::AlignmentPair Trace::getAlignment()
	{
		// Initialize return struct
		// GG: match horiz, match verti, matches
		AlignmentPair alignments = { "", "", 0 };

		// Find antiDiag3's max => This is exit score (need position)
		// GG: go to the last element of trace, find the max in the last trace
		// GG: but it might not be the global maximum
		auto itAtMax = trace.rbegin() + (trace.size() - 1 - maxPos);
		size_t dp_pos = itAtMax->antiDiag3.argmax();

		// Follow dp_pos back and track alignment
		size_t sq_left_pos  = 0;
		size_t sq_above_pos = 0;
		size_t sq_diag_pos  = 0;

		// GG: get antidiag index of antiDiag3 (32 elements-wide)
		// GG: keep track of the max
		// GG: find elements that created that max and so on and so for
		// GG: it's hard bc we have antidiags and not an actual DP matrix
		// GG: much work in converting indexes
		for ( auto it = itAtMax; std::next(it) != trace.rend(); ++it )
		{
			// Calculate necessary position in antiDiag1 and antiDiag2
			auto nit = std::next(it);

			sq_left_pos  = dp_pos;
			sq_above_pos = dp_pos + 1;
			sq_diag_pos  = dp_pos;

			// Calculate where the max value came from
			char queryHChar = it->vqueryh[dp_pos];
			char queryVChar = it->vqueryv[dp_pos];

			int st = it->antiDiag3[dp_pos] + it->scoreOffset;

			int offset = nit == trace.rend() ? 0 : nit->scoreOffset;
			int sa = sq_above_pos >= VectorRegister::VECTORWIDTH ? VectorRegister::NINF : it->antiDiag2[sq_above_pos] + scoringScheme.getGapScore() + offset;
			int sl = sq_left_pos  >= VectorRegister::VECTORWIDTH ? VectorRegister::NINF : it->antiDiag2[sq_left_pos]  + scoringScheme.getGapScore() + offset;
			int sd = sq_diag_pos  >= VectorRegister::VECTORWIDTH ? VectorRegister::NINF : it->antiDiag1[sq_diag_pos]  + scoringScheme.score( queryHChar, queryVChar ) + offset;

			if ( sd != VectorRegister::NINF && sd == st )
			{
				if ( queryHChar == queryVChar )
					alignments.matches++;

				alignments.alignH.push_back( queryHChar );
				alignments.alignV.push_back( queryVChar );

				dp_pos = sq_diag_pos + (it->lastMove == nit->lastMove ? (-2 * it->lastMove) + 1 : 0);
				++it;
			}
			else if ( sl != VectorRegister::NINF && sl == st )
			{
				alignments.alignH.push_back( queryHChar );
				alignments.alignV.push_back( '-' );

				dp_pos = sq_left_pos - it->lastMove;
			}
			else if ( sa != VectorRegister::NINF && sa == st )
			{
				alignments.alignH.push_back( '-' );
				alignments.alignV.push_back( queryVChar );

				dp_pos = sq_above_pos - it->lastMove;
			}
			else
			{
				std::cout << "ERROR: Failure to find backpath in part 1 of traceback" << std::endl;
			}
		}

		// Traceback in Opening phase

		// Find position in matrix
		// dp_pos is pos in antiDiag2 now, need to convert to DPMatrix coord
		int i = VectorRegister::LOGICALWIDTH + 1 - dp_pos;
		int j = dp_pos + 2;

		while (i > 0 && j > 0)
		{
			char queryHChar = queryh[j-1];
			char queryVChar = queryv[i-1];

			int st = DPMatrix[i][j];
			int sa = DPMatrix[i-1][j]   + scoringScheme.getGapScore();
			int sl = DPMatrix[i][j-1]   + scoringScheme.getGapScore();
			int sd = DPMatrix[i-1][j-1] + scoringScheme.score( queryHChar, queryVChar );

			// GG: none of these is satified at the first iteration than it seems just fine
			if (sd == st)
			{
				if (queryHChar == queryVChar)
					alignments.matches++;

				alignments.alignH.push_back( queryHChar );
				alignments.alignV.push_back( queryVChar );
				--i;
				--j;
			}
			else if ( sl == st )
			{
				alignments.alignH.push_back( queryHChar );
				alignments.alignV.push_back( '-' );
				--j;
			}
			else if ( sa != VectorRegister::NINF && sa == st )
			{
				alignments.alignH.push_back( '-' );
				alignments.alignV.push_back( queryVChar );
				--i;
			}
			else
			{
				std::cout << "ERROR: Failure to find backpath in part 2 of traceback" << std::endl;
			}
		}

		std::reverse( alignments.alignH.begin(), alignments.alignH.end() );
		std::reverse( alignments.alignV.begin(), alignments.alignV.end() );
		return alignments;
	}

	void Trace::saveOpeningPhaseDPMatrix ( std::vector< std::vector<int> > _DPMatrix, int8_t* _queryh, int8_t* _queryv )
	{
		DPMatrix = _DPMatrix;
		queryh = _queryh;
		queryv = _queryv;
	}

}