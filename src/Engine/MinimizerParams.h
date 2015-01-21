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
	                RealType delta2_,
	                RealType tol_,
	                bool verbose_)
	    : algo(algo_),
	      maxIter(maxIter_),
	      delta(delta_),
	      delta2(delta2_),
	      tol(tol_),
	      verbose(verbose_)
	{}

	EnumAlgo algo;
	SizeType maxIter;
	RealType delta;
	RealType delta2;
	RealType tol;
	bool verbose;
};


template<typename RealType>
std::ostream& operator<<(std::ostream& os, const MinimizerParams<RealType>& m)
{
	os<<"algo= "<<m.algo<<"\n";
	os<<"maxIter= "<<m.maxIter<<"\n";
	os<<"delta= "<<m.delta<<"\n";
	os<<"delta2= "<<m.delta2<<"\n";
	os<<"tolerance= "<<m.tol<<"\n";
	os<<"verbose= "<<m.verbose<<"\n";
	return os;
}
} // namespace yasw

#endif // MINIMIZERPARAMS_H

