#ifndef ENERGY_C_FUNCTION_H
#define ENERGY_C_FUNCTION_H
#include "SpaceConnectors.h"
#include "Vector.h"
#include "InitConfig.h"

namespace yasw {

template<typename ComplexOrRealType_>
class EnergyCollinearFunction {

public:

	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef yasw::SpaceConnectors<ComplexOrRealType> SpaceConnectorsType;
	typedef unsigned long long int  LongSizeType;
	typedef InitConfig<ComplexOrRealType_, true> InitConfigType;

	class Configuration {

	public:

		Configuration(SizeType totalSpins,
		              SizeType fixedSpins,
		              bool,
		              const InitConfigType&)
		    : totalSpins_(totalSpins), fixedSpins_(fixedSpins), data_(0), mask_(0)
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

			if (totalSpins >= 256*sizeof(LongSizeType)) {
				PsimagLite::String str("Configuration: Configuration is");
				str += " too big\n";
				throw PsimagLite::RuntimeError(str);
			}
		}

		bool isValid(LongSizeType x) const
		{
			return ((mask_&x) == mask_);
		}

		void print(std::ostream& os) const
		{
			os<<totalSpins_<<" 2\n";
			for (SizeType i = 0; i < totalSpins_; ++i) {
				int si = valueAt(data_, i);
				if (si == 1)
					os<<"0 0\n";
				else
					os<<M_PI<<" 0\n";
			}
		}

		LongSizeType& operator()() { return data_; }

		const LongSizeType& operator()() const { return data_; }

		static int valueAt(LongSizeType ket, SizeType i)
		{
			LongSizeType mask = 1;
			mask <<= i;
			return (mask & ket) ? 1 : -1;
		}

	private:

		SizeType totalSpins_;
		SizeType fixedSpins_;
		LongSizeType data_;
		LongSizeType mask_;
	};

public:

	typedef  PsimagLite::Vector<double>::Type VectorRealType;
	typedef Configuration ConfigurationType;

	EnergyCollinearFunction(const SpaceConnectorsType& spaceConnectors,
	                        const VectorRealType&,
	                        const VectorRealType& moduli)
	    : moduli_(moduli), sc_(spaceConnectors)
	{}

	template<typename DummyType>
	RealType minimize(ConfigurationType& config, const DummyType&, SizeType) const
	{
		RealType minEnergy = 1e10;
		SizeType total = 1;
		SizeType lda = size();
		total <<= lda;
		for (LongSizeType ket = 0; ket < total; ++ket) {
			if (!config.isValid(ket)) continue;
			RealType energy = this->operator()(ket);
			if (energy < minEnergy) {
				config() = ket; // save this config
				minEnergy = energy;
			}
		}

		std::cout<<"Angles\n";
		std::cout<<lda<<" 2\n";
		for (SizeType i = 0; i < lda; ++i) {
			int si = Configuration::valueAt(config(), i);
			if (si == 1)
				std::cout<<"0 0\n";
			else
				std::cout<<M_PI<<" 0\n";
		}

		return minEnergy;
	}

	SizeType totalSpins() const { return sc_.rows(); }

	SizeType size() const { return sc_.rows(); }

private:

	RealType operator()(LongSizeType ket) const
	{
		RealType sum = 0.0;
		for (SizeType i = 0; i < sc_.size(); ++i) {
			sum += energyThisCell(i, ket);
		}

		return sum;
	}

	RealType energyThisCell(SizeType ind, LongSizeType ket) const
	{
		RealType sum = 0;
		SizeType lda = size();
		for (SizeType i = 0; i < lda; ++i) {
			const RealType si = Configuration::valueAt(ket, i)*moduli_[i];
			for (SizeType j = 0; j < lda; ++j) {
				sum += std::real(sc_(i,j,ind))*si*Configuration::valueAt(ket, j)*moduli_[j];
			}
		}

		return sum;
	}

	const VectorRealType& moduli_;
	const SpaceConnectorsType& sc_;
};
}

#endif // ENERGY_C_FUNCTION_H

