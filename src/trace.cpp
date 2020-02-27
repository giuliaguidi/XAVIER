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

	/**
	 * CIGAR Standard Format
	 * D : Deletion (gap in the target sequence).
	 * I : Insertion (gap in the query sequence).
	 * = : Alignment column containing two identical letters. USEARCH can read CIGAR strings using this operation, but does not generate them.
	 * X : Alignment column containing a mismatch, i.e. two different letters. USEARCH can read CIGAR strings using this operation, but does not generate them.
	 */ 

	std::string Trace::compression(const std::string& str)
	{
	    int i = str.size();
	    std::string cigar;
	
	    for (int j = 0; j < i; ++j)
		{
	        int count = 1;
	        while (str[j] == str[j+1])
			{
	            count++;
	            j++;
	        }
	        cigar += std::to_string(count);
	        cigar.push_back(str[j]);
	    }
	    return cigar;
	}

	Trace::AlignmentPair Trace::getAlignment()
	{
		// GG: cigar, matches
		AlignmentPair traceback = {"", 0};

		// Find antiDiag3's max => This is exit score (need position)
		auto itAtMax = trace.rbegin() + (trace.size() - 1 - maxPos);
		size_t dp_pos = itAtMax->antiDiag3.argmax();

		// Follow dp_pos back and track alignment
		size_t sq_left_pos  = 0;
		size_t sq_above_pos = 0;
		size_t sq_diag_pos  = 0;

		for (auto it = itAtMax; std::next(it) != trace.rend(); ++it)
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

			if (sd != VectorRegister::NINF && sd == st)
			{
				if (queryHChar == queryVChar)
				{
					// GG: match
					traceback.matches++;
					traceback.cigar.push_back('=');
				}
				else
				{
					// GG: mismatch
					traceback.cigar.push_back('X');	
				}
				
				dp_pos = sq_diag_pos + (it->lastMove == nit->lastMove ? (-2 * it->lastMove) + 1 : 0);
				++it;
			}
			else if ( sl != VectorRegister::NINF && sl == st )
			{
				// GG: gap in the query 
				traceback.cigar.push_back('I');
				dp_pos = sq_left_pos - it->lastMove;		
			}
			else if ( sa != VectorRegister::NINF && sa == st )
			{
				// GG: gap in the target
				traceback.cigar.push_back('D');
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

			if (sd == st)
			{
				if (queryHChar == queryVChar)
				{
					// GG: match
					traceback.matches++;
					traceback.cigar.push_back('=');
				}
				else
				{
					// GG: mismatch
					traceback.cigar.push_back('X');	
				}
				--i;
				--j;
			}
			else if ( sl == st )
			{
				traceback.cigar.push_back('I');
				--j;
			}
			else if ( sa != VectorRegister::NINF && sa == st )
			{
				traceback.cigar.push_back('D');
				--i;
			}
			else
			{
				std::cout << "ERROR: Failure to find backpath in part 2 of traceback" << std::endl;
			}
		}
		
		// TODO: double check how we handle this when extending both left and right
		// (no need to reverse if we are extending left)
		std::reverse(traceback.cigar.begin(), traceback.cigar.end());

		// GG: cigar compression
		traceback.cigar = compression(traceback.cigar);

		return traceback;
	}

	void Trace::saveOpeningPhaseDPMatrix ( std::vector< std::vector<int> > _DPMatrix, int8_t* _queryh, int8_t* _queryv )
	{
		DPMatrix = _DPMatrix;
		queryh = _queryh;
		queryv = _queryv;
	}

}
