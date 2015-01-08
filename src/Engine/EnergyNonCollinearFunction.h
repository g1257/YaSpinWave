#ifndef ENERGY_NC_FUNCTION_H
#define ENERGY_NC_FUNCTION_H
#include "SpaceConnectors.h"
#include "Vector.h"
#include "Minimizer.h"
#include "MersenneTwister.h"

namespace yasw {

template<typename RealType, typename ComplexOrRealType>
class EnergyNonCollinearFunction {

	typedef yasw::SpaceConnectors<RealType,ComplexOrRealType> SpaceConnectorsType;
	typedef EnergyNonCollinearFunction<RealType,ComplexOrRealType> ThisType;

public:

	typedef RealType FieldType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef VectorRealType ConfigurationType;

	EnergyNonCollinearFunction(PsimagLite::String jfile)
	    : sc_(jfile),data_(2*sc_.rows())
	{}

	void minimize(ConfigurationType& config)
	{
		PsimagLite::MersenneTwister rng(1234);
		SizeType lda = size();
		assert(data_.size() == 2*lda);

		for (SizeType i = 0; i < lda; ++i) {
			data_[2*i] = rng() * M_PI;
			data_[2*i+1] = rng() * 2.0*M_PI;
		}

		int maxIter = 100;
		PsimagLite::Minimizer<RealType,ThisType> min(*this, maxIter);
		int used = min.simplex(config);
		std::cerr<<"EnergyNonCollinearFunction::minimize(): done after ";
		std::cerr<<used<<" iterations.\n";

		std::cout<<"Angles\n";
		std::cout<<lda<<" 2\n";
		for (SizeType i = 0; i < lda; ++i) {
			std::cout<<data_[2*i]<<" "<<data_[2*i+1]<<"\n";
		}

		config = data_;
	}

	SizeType size() const { return sc_.rows(); }

	FieldType operator()(FieldType* data, SizeType n)
	{
		assert(n == size());
		FieldType sum = 0;
		for (SizeType i = 0; i < n; ++i)
			data_[i] = data[i];

		for (SizeType i = 0; i < sc_.size(); ++i) {
			sum += energyThisCell(i);
		}

		return sum;
	}

private:

	FieldType energyThisCell(SizeType ind) const
	{
		RealType sum = 0;
		SizeType lda = size();
		VectorRealType vi(3,0);
		VectorRealType vj(3,0);
		for (SizeType i = 0; i < lda; ++i) {
			buildVector(vi,data_,i);
			for (SizeType j = 0; j < lda; ++j) {
				buildVector(vj,data_,j);
				sum += std::real(sc_(i,j,ind))*scalarProduct(vi,vj);
			}
		}

		return sum;
	}

	void buildVector(VectorRealType& dst,
	                 const VectorRealType& src,
	                 int i) const
	{
		SizeType x = i*2;
		assert(dst.size() == 3);
		assert(src.size() > x+1);
		dst[0] = sin(src[x])*cos(src[1+x]);
		dst[1] = sin(src[x])*sin(src[1+x]);
		dst[2] = cos(src[x]);
	}

	SpaceConnectorsType sc_;
	VectorRealType data_;
};

}

#endif // ENERGY_NC_FUNCTION_H

