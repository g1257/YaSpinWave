#ifndef ENERGY_NC_FUNCTION_H
#define ENERGY_NC_FUNCTION_H
#include "SpaceConnectors.h"
#include "Vector.h"
#include "Minimizer.h"
#include "MersenneTwister.h"
#include "MinimizerParams.h"
#include "Angles.h"

namespace yasw {

template<typename RealType, typename ComplexOrRealType>
class EnergyNonCollinearFunction {

	typedef yasw::SpaceConnectors<RealType,ComplexOrRealType> SpaceConnectorsType;
	typedef EnergyNonCollinearFunction<RealType,ComplexOrRealType> ThisType;

	class Configuration {

	public:

		typedef Angles<RealType> AnglesType;
		typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

		Configuration(SizeType totalSpins,
		              SizeType fixedSpins,
		              PsimagLite::String afile = "" ,
		              int seed = 0)
		    : fixedSpins_(fixedSpins)
		{

			if (totalSpins <= fixedSpins) {
				PsimagLite::String msg("Configuration: ");
				msg += "Too many fixed spins\n";
				throw PsimagLite::RuntimeError(msg);
			}

			data_.resize(2*(totalSpins-fixedSpins),0.0);

			if (seed > 0 && afile != "") {
				PsimagLite::String msg("Configuration: Providing both ");
				msg += " seed and angles file is an error\n";
				throw PsimagLite::RuntimeError(msg);
			}

			if (afile != "") {
				AnglesType angles(afile);
				if (angles.size() != totalSpins) {
					PsimagLite::String msg("Configuration: angles file");
					msg += " has different size than j file\n";
					throw PsimagLite::RuntimeError(msg);
				}

				if (!isValidAngles(angles)) {
					PsimagLite::String msg("Configuration: FATAL: ");
					msg += "First angle of " + afile + " is not with theta=0\n";
					throw PsimagLite::RuntimeError(msg);
				}

				for (SizeType i = fixedSpins_; i < angles.size(); ++i) {
					SizeType twiceIndex = 2*(i - fixedSpins_);
					if (twiceIndex + 1 >= data_.size()) break;
					data_[twiceIndex] = angles.theta(i);
					data_[twiceIndex+1] = angles.phi(i);
				}

				return;
			}

			if (seed == 0) seed = 1234;
			randomize(seed);
		}

		void print(std::ostream& os) const
		{
			SizeType lda = static_cast<SizeType>(data_.size()*0.5);
			os<<(lda+fixedSpins_)<<" 2\n";
			for (SizeType i = 0; i < fixedSpins_; ++i)
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

		SizeType fixedSpins() const { return fixedSpins_; }

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

		bool isValidAngles(const AnglesType& angles) const
		{
			for (SizeType i = 0; i < fixedSpins_; ++i)
				if (fabs(angles.theta(i))>1e-6) return false;

			return true;
		}

		VectorRealType data_;
		SizeType fixedSpins_;
	};

public:

	typedef RealType FieldType;
	typedef typename Configuration::VectorRealType VectorRealType;
	typedef Configuration ConfigurationType;
	typedef PsimagLite::Minimizer<RealType,ThisType> MinimizerType;

	EnergyNonCollinearFunction(PsimagLite::String jfile)
	    : sc_(jfile),data_(1,0)
	{}

	void minimize(ConfigurationType& config,
	              const MinimizerParams<RealType>& minParams)
	{
		if (sc_.rows() <= config.fixedSpins()) {
			PsimagLite::String msg("Configuration: ");
			msg += "Too many fixed spins\n";
			throw PsimagLite::RuntimeError(msg);
		}

		data_ = config;
		assert(data_.size()+2*config.fixedSpins() == 2*sc_.rows());

		MinimizerType min(*this, minParams.maxIter,minParams.verbose);
		std::cerr<<"Initial config\n";
		config.print(std::cerr);
		int used = min.simplex(config(), minParams.delta, minParams.tol);
		data_ = config;
		std::cerr<<"Minimizer params\n";
		std::cerr<<minParams;
		std::cerr<<"EnergyNonCollinearFunction::minimize(): ";
		if (min.status() == MinimizerType::GSL_SUCCESS) {
			std::cerr<<" converged after ";
		} else {
			std::cerr<<"NOT CONVERGED after ";
		}

		std::cerr<<used<<" iterations.\n";
		std::cerr<<"Energy at Minimum= "<<operator()(&(data_()[0]),data_.size());
		std::cerr<<"\n";

		std::cout<<"Angles\n";
		data_.print(std::cout);

		config = data_;
	}

	SizeType totalSpins() const { return sc_.rows(); }

	SizeType size() const
	{
		assert(sc_.rows() > data_.fixedSpins());
		return 2*(sc_.rows() - data_.fixedSpins());
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
	                 SizeType i) const
	{
		SizeType fixedSpins = data_.fixedSpins();
		if (i < fixedSpins) {
			dst[0] = dst[1] = 0;
			dst[2] = 1;
			return;
		}

		SizeType ii = i - fixedSpins;
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

