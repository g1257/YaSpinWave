#ifndef MINIMIZERPARAMS_H
#define MINIMIZERPARAMS_H
#include "Vector.h"

namespace yasw {

template<typename RealType>
struct MinimizerParams {

	MinimizerParams(SizeType maxIter_, RealType delta_, RealType tol_, bool verbose_)
	    : maxIter(maxIter_),
	      delta(delta_),
	      tol(tol_),
	      verbose(verbose_)
	{}

	SizeType maxIter;
	RealType delta;
	RealType tol;
	bool verbose;
};

} // namespace yasw

#endif // MINIMIZERPARAMS_H

