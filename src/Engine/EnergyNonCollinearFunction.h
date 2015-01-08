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

	class Configuration {

	public:

		typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

		Configuration(SizeType twiceTheSites,int seed = 0)
		    : data_(twiceTheSites,0.0)
		{
			if (seed > 0)
				randomize(seed);
		}

		void print(std::ostream& os) const
		{
			SizeType lda = static_cast<SizeType>(data_.size()*0.5);
			os<<lda<<" 2\n";
			os<<"0 0\n";
			for (SizeType i = 0; i < lda; ++i) {
				os<<data_[2*i]<<" "<<data_[2*i+1]<<"\n";
			}
		}

		void fromRaw(RealType* data, SizeType n)
		{
			assert(data_.size() == n);

			for (SizeType i = 0; i < n; ++i)
				data_[i] = data[i];
		}

		VectorRealType& operator()() { return data_; }

		const VectorRealType& operator()() const { return data_; }

		SizeType size() const { return data_.size(); }

	private:

		void randomize(int seed)
		{
			PsimagLite::MersenneTwister rng(seed);
			SizeType lda = static_cast<SizeType>(data_.size()*0.5);
			for (SizeType i = 0; i < lda; ++i) {
				data_[2*i] = rng() * M_PI;
				data_[2*i+1] = rng() * 2.0*M_PI;
			}
		}

		VectorRealType data_;
	};

public:

	typedef RealType FieldType;
	typedef typename Configuration::VectorRealType VectorRealType;
	typedef Configuration ConfigurationType;

	EnergyNonCollinearFunction(PsimagLite::String jfile)
	    : sc_(jfile),data_(0)
	{}

	void minimize(ConfigurationType& config)
	{
		data_ = config;
		assert(data_.size()+2 == 2*sc_.rows());

		int maxIter = 100;
		PsimagLite::Minimizer<RealType,ThisType> min(*this, maxIter);
		int used = min.simplex(config());
		std::cerr<<"EnergyNonCollinearFunction::minimize(): done after ";
		std::cerr<<used<<" iterations.\n";

		std::cout<<"Angles\n";
		data_.print(std::cout);

		config = data_;
	}

	SizeType size() const
	{
		assert(sc_.rows() > 1);
		return 2*(sc_.rows()-1);
	}

	FieldType operator()(FieldType* data, SizeType n)
	{
		data_.fromRaw(data, n);

		FieldType sum = 0;
		for (SizeType i = 0; i < sc_.size(); ++i) {
			sum += energyThisCell(i);
		}

		return sum;
	}

private:

	FieldType energyThisCell(SizeType ind) const
	{
		RealType sum = 0;
		SizeType lda = sc_.rows();
		VectorRealType vi(3,0);
		VectorRealType vj(3,0);
		for (SizeType i = 0; i < lda; ++i) {
			buildVector(vi,data_(),i);
			for (SizeType j = 0; j < lda; ++j) {
				buildVector(vj,data_(),j);
				sum += std::real(sc_(i,j,ind))*scalarProduct(vi,vj);
			}
		}

		return sum;
	}

	void buildVector(VectorRealType& dst,
	                 const VectorRealType& src,
	                 int i) const
	{
		if (i == 0) {
			dst[0] = dst[1] = 0;
			dst[2] = 1;
			return;
		}

		SizeType ii = i - 1;
		SizeType x = ii*2;
		assert(dst.size() == 3);
		assert(src.size() > x+1);
		dst[0] = sin(src[x])*cos(src[1+x]);
		dst[1] = sin(src[x])*sin(src[1+x]);
		dst[2] = cos(src[x]);
	}

	SpaceConnectorsType sc_;
	ConfigurationType data_;
};

}

#endif // ENERGY_NC_FUNCTION_H

