#ifndef ENERGY_C_FUNCTION_H
#define ENERGY_C_FUNCTION_H
#include "SpaceConnectors.h"
#include "Vector.h"

namespace yasw {

template<typename ComplexOrRealType>
class EnergyCollinearFunction {

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef yasw::SpaceConnectors<ComplexOrRealType> SpaceConnectorsType;
	typedef unsigned long long int  LongSizeType;

	class Configuration {

	public:

		Configuration(SizeType totalSpins,
		              SizeType fixedSpins,
		              bool,
		              PsimagLite::String afile = "",
		              SizeType seed = 0)
		    : data_(0),fixedSpins_(fixedSpins),mask_(0)
		{
			if (fixedSpins >= totalSpins) {
				PsimagLite::String str("Configuration: ");
				str += "Too many fixed spins\n";
				throw PsimagLite::RuntimeError(str);
			}

			LongSizeType tmp = 1;
			for (SizeType i = 0; i < fixedSpins_; ++i) {
				mask_ |= tmp;
				tmp <<= 1;
			}

			if (seed > 0) {
				PsimagLite::String str("Configuration: No seed support yet");
				str +=" for collinear configuration\n";
				std::cerr<<"WARNING "<<str;
			}

			if (totalSpins >= 256*sizeof(LongSizeType)) {
				PsimagLite::String str("Configuration: Configuration is");
				str += " too big\n";
				throw PsimagLite::RuntimeError(str);
			}

			if (afile != "") {
				PsimagLite::String str("Configuration: No angles support yet");
				str +=" for collinear configuration\n";
				std::cerr<<"WARNING "<<str;
			}
		}

		void fromRaw(LongSizeType x)
		{
			data_ = x;
		}

		bool isValid(LongSizeType x) const
		{
			return ((mask_&x) == mask_);
		}

		LongSizeType& operator()() { return data_; }

		const LongSizeType& operator()() const { return data_; }

	private:

		LongSizeType data_;
		SizeType fixedSpins_;
		LongSizeType mask_;
	};

public:

	typedef Configuration ConfigurationType;

	EnergyCollinearFunction(PsimagLite::String jfile, bool verbose)
	    : sc_(jfile, verbose)
	{}

	template<typename DummyType>
	void minimize(ConfigurationType& config, const DummyType&) const
	{
		RealType minEnergy = 1e10;
		SizeType total = 1;
		SizeType lda = size();
		total <<= lda;
		for (SizeType ket = 0; ket < total; ++ket) {
			if (!config.isValid(ket)) continue;
			RealType energy = this->operator()(ket);
			if (energy < minEnergy) {
				config.fromRaw(ket);
				minEnergy = energy;
			}
		}

		std::cout<<"Angles\n";
		std::cout<<lda<<" 2\n";
		for (SizeType i = 0; i < lda; ++i) {
			int si = valueAt(config(),i);
			if (si == 1)
				std::cout<<"0 0\n";
			else
				std::cout<<M_PI<<" 0\n";
		}
	}

	SizeType totalSpins() const { return sc_.rows(); }

	SizeType size() const { return sc_.rows(); }

private:

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
		return (mask & ket) ? 1 : -1;
	}

	SpaceConnectorsType sc_;
};

}

#endif // ENERGY_C_FUNCTION_H

