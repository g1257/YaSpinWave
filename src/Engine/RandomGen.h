#ifndef RANDOMGEN_H
#define RANDOMGEN_H
#include "Vector.h"
#include "MersenneTwister.h"

namespace yasw {

template<typename ComplexOrRealType>
class RandomGen {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	RandomGen(int seed)
	: seed_(seed), rng_(0)
	{
		if (seed_ > 0)
			rng_ = new PsimagLite::MersenneTwister(seed);
	}

	~RandomGen()
	{
		delete rng_;
		rng_ = 0;
	}

	bool needsRandom() const { return (seed_ > 0); }

	void randomize(VectorRealType& data) const
	{
		if (!needsRandom())
			err("RandomGen: INTERNAL ERROR\n");

		SizeType lda = static_cast<SizeType>(data.size()*0.5);
		for (SizeType i = 0; i < lda; ++i) {
			data[2*i] = rng_->operator()() * M_PI;
			data[2*i+1] = rng_->operator()() * 2.0*M_PI;
		}
	}

	RealType random(RealType delta) const
	{
		return 2.0*rng_->operator()() * delta - delta;
	}

private:

	RandomGen(const RandomGen&) = delete;
	RandomGen& operator=(const RandomGen&) = delete;

	int seed_;
	PsimagLite::MersenneTwister* rng_;
};
}
#endif // RANDOMGEN_H
