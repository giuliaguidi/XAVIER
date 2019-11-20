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
	             			 const VectorRegister& _vqv, const int64_t offset ):
		antiDiag1( _ad1 ),
		antiDiag2( _ad2 ),
		antiDiag3( _ad3 ),
		vqueryh( _vqh ),
		vqueryv( _vqv ),
		scoreOffset( offset )
		{ }

	Trace::Trace()
	{
	}

	void Trace::pushbackState ( const VectorRegister& _ad1, const VectorRegister& _ad2,
	   	                        const VectorRegister& _ad3, const VectorRegister& _vqh,
		                        const VectorRegister& _vqv, const int64_t offset)
	{
		static int i = 0;
		std::cout << "hello " << i++ << std::endl;
		trace.push_back( TraceEntry( _ad1, _ad2, _ad3, _vqh, _vqv, offset ) );
	}

}