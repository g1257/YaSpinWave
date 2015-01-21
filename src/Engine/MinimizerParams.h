#ifndef MINIMIZERPARAMS_H
#define MINIMIZERPARAMS_H
#include "Vector.h"

namespace yasw {

template<typename RealType>
struct MinimizerParams {

	enum EnumAlgo {SIMPLEX, CONJUGATE_GRADIENT};

	MinimizerParams(EnumAlgo algo_,
	                SizeType maxIter_,
	                RealType delta_,
	                RealType tol_,
	                bool verbose_)
	    : algo(algo_),
	      maxIter(maxIter_),
	      delta(delta_),
	      tol(tol_),
	      verbose(verbose_)
	{}

	EnumAlgo algo;
	SizeType maxIter;
	RealType delta;
	RealType tol;
	bool verbose;
};


template<typename RealType>
std::ostream& operator<<(std::ostream& os, const MinimizerParams<RealType>& m)
{
	os<<"maxIter= "<<m.maxIter<<"\n";
	os<<"delta= "<<m.delta<<"\n";
	os<<"tolerance= "<<m.tol<<"\n";
	os<<"verbose= "<<m.verbose<<"\n";
	return os;
}
} // namespace yasw

#endif // MINIMIZERPARAMS_H

