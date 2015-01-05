#ifndef ENERGY_C_FUNCTION_H
#define ENERGY_C_FUNCTION_H
#include "SpaceConnectors.h"
#include "Vector.h"

namespace yasw {

template<typename RealType, typename ComplexOrRealType>
class EnergyCollinearFunction {

	typedef yasw::SpaceConnectors<RealType,ComplexOrRealType> SpaceConnectorsType;
	typedef unsigned long long int  LongSizeType;

public:

	typedef LongSizeType ConfigurationType;

	EnergyCollinearFunction(PsimagLite::String jfile)
	    : sc_(jfile)
	{}

	void minimize(ConfigurationType& config) const
	{
		RealType minEnergy = 1e10;
		SizeType total = 1;
		SizeType lda = size();
		total <<= lda;
		for (SizeType ket = 0; ket < total; ++ket) {
			RealType energy = this->operator()(ket);
			if (energy < minEnergy) {
				config = ket;
				minEnergy = energy;
			}
		}

		std::cout<<"Angles\n";
		std::cout<<lda<<" 2\n";
		for (SizeType i = 0; i < lda; ++i) {
			int si = valueAt(config,i);
			if (si == 1)
				std::cout<<"0 0\n";
			else
				std::cout<<M_PI<<" 0\n";
		}
	}

private:

	SizeType size() const { return sc_.rows(); }

	RealType operator()(LongSizeType ket) const
	{
		RealType sum = 0.0;
		for (SizeType i = 0; i < sc_.size(); ++i) {
			sum += energyThisCell(i,ket);
		}

		return sum;
	}

	RealType energyThisCell(SizeType ind, LongSizeType ket) const
	{
		RealType sum = 0;
		SizeType lda = size();
		for (SizeType i = 0; i < lda; ++i) {
			int si = valueAt(ket,i);
			for (SizeType j = 0; j < lda; ++j) {
				sum += std::real(sc_(i,j,ind))*si*valueAt(ket,j);
			}
		}

		return sum;
	}

	int valueAt(LongSizeType ket, SizeType i) const
	{
		LongSizeType mask = 1;
		mask <<= i;
		return (mask & ket) ? -1 : 1;
	}

	SpaceConnectorsType sc_;
};

}

#endif // ENERGY_C_FUNCTION_H
