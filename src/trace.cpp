/**
 * File: trace.cpp
 * Author: G. Guidi, E. Younis
 * Description: Xavier Trace Type Source.
 */

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

	void Trace::pushbackState ( const VectorRegister& _ad1, const VectorRegister& _ad2,
	   	                        const VectorRegister& _ad3, const VectorRegister& _vqh,
		                        const VectorRegister& _vqv, const int64_t offset,
		                        const int _lastMove )
	{
		trace.push_back( TraceEntry( _ad1, _ad2, _ad3, _vqh, _vqv, offset, _lastMove ) );
	}

	Trace::AlignmentPair Trace::getAlignment()
	{

		// Start on last trace entry
		auto it = trace.rbegin();
		TraceEntry& e = *it;

		// Find antiDiag3's max => This is exit score (need position)
		size_t dp_pos = e.antiDiag3.argmax();

		// Initialize return struct
		AlignmentPair alignments = { "", "", 0 };

		// Follow dp_pos back and track alignment
		for ( auto it = trace.rbegin(); it != trace.rend(); ++it )
		{
			e = *it;

			size_t sq_left_pos  = dp_pos;
			size_t sq_above_pos = dp_pos + 1;
			size_t sq_diag_pos  = dp_pos;

			char queryHChar = e.vqueryh[dp_pos];
			char queryVChar = e.vqueryv[dp_pos];

			if ( queryHChar == queryVChar )
				alignments.matches++;

			int sa = e.antiDiag2[sq_above_pos] + scoringScheme.getGapScore();
			int sl = e.antiDiag2[sq_left_pos]  + scoringScheme.getGapScore();
			int sd = e.antiDiag1[sq_diag_pos]  + scoringScheme.score( queryHChar, queryVChar );

			// EY: What should the default choice be??
			if ( sa > sl && sa > sd )
			{
				alignments.alignH.append( 1, '-' );
				dp_pos = sq_above_pos;
			}
			else if ( sl > sd )
			{
				alignments.alignV.append( 1, '-' );
				dp_pos = sq_left_pos;
			}
			else
			{
				alignments.alignH.append( 1, queryHChar );
				alignments.alignV.append( 1, queryVChar );
				dp_pos = sq_diag_pos;
				++it;
			}

		}

		// Handle Opening Phase Specially

		// Find position in matrix
		// dp_pos is pos in the 3rd anti diag, this is not in dp_matrix
		int i = VectorRegister::LOGICALWIDTH - dp_pos + 2;
		int j = dp_pos + 2;

		while ( i > 0 && j > 0 )
		{
			char queryHChar = e.vqueryh[i-2]; // Don't have access to full queries (this might be buggy)
			char queryVChar = e.vqueryv[j-2]; // Don't have access to full queries (this might be buggy)

			int sa = DPMatrix[i-1][j]   + scoringScheme.getGapScore();
			int sl = DPMatrix[i][j-1]   + scoringScheme.getGapScore();
			int sd = DPMatrix[i-1][j-1] + scoringScheme.score( queryHChar, queryVChar );

			if ( sa > sl && sa > sd )
			{
				alignments.alignH.append( 1, '-' );
				--i;
			}
			else if ( sl > sd )
			{
				alignments.alignV.append( 1, '-' );
				--j;
			}
			else
			{
				alignments.alignH.append( 1, queryHChar );
				alignments.alignV.append( 1, queryVChar );
				--i;
				--j;
			}

		}

		std::reverse( alignments.alignH.begin(), alignments.alignH.end() );
		std::reverse( alignments.alignV.begin(), alignments.alignV.end() );
		return alignments;
	}

	void Trace::saveOpeningPhaseDPMatrix ( std::vector< std::vector<int> > _DPMatrix )
	{
		DPMatrix = _DPMatrix;
	}

}