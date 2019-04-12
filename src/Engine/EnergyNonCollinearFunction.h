#ifndef ENERGY_NC_FUNCTION_H
#define ENERGY_NC_FUNCTION_H
#include "SpaceConnectors.h"
#include "Vector.h"
#include "Minimizer.h"
#include "MersenneTwister.h"
#include "MinimizerParams.h"
#include "Angles.h"

namespace yasw {

template<typename ComplexOrRealType>
class EnergyNonCollinearFunction {

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef yasw::SpaceConnectors<ComplexOrRealType> SpaceConnectorsType;
	typedef EnergyNonCollinearFunction<ComplexOrRealType> ThisType;
	typedef std::complex<RealType> ComplexType;

	class Configuration {

	public:

		typedef Angles<RealType> AnglesType;
		typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

		Configuration(SizeType totalSpins,
		              SizeType fixedSpins,
		              bool verbose,
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
				AnglesType angles(afile, verbose);
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

		void fromRaw(const VectorRealType& data)
		{
			data_ = data;
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

	EnergyNonCollinearFunction(PsimagLite::String jfile,
	                           const VectorRealType& qvector,
	                           bool verbose)
	    : qvector_(qvector), sc_(jfile, verbose), data_(1,0,verbose)
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
		int used = 0;
		if (minParams.algo == MinimizerParams<RealType>::SIMPLEX) {
			used = min.simplex(config(),
			                   minParams.delta,
			                   minParams.tol);
		} else {
			used = min.conjugateGradient(config(),
			                             minParams.delta,
			                             minParams.delta2,
			                             minParams.tol,
			                             minParams.saveEvery);
		}

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
		std::cerr<<"Energy at Minimum= "<<operator()(data_());
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

	FieldType operator()(const VectorRealType& data)
	{
		data_.fromRaw(data);

		FieldType sum = 0;
		for (SizeType i = 0; i < sc_.size(); ++i) {
			sum += energyThisCell(i);
		}

		return sum;
	}

	void df(VectorRealType& df, const VectorRealType& data)
	{
		data_.fromRaw(data);

		const SizeType n = data.size();
		for (SizeType i = 0; i < n; ++i) {
			df[i] = 0.0;
			for (SizeType ind = 0; ind < sc_.size(); ++ind) {
				df[i] += derivativeEnergyThisCell(ind,i);
			}
		}

	}

private:

	FieldType energyThisCell(SizeType ind) const
	{
		ComplexType sum = 0;
		SizeType lda = sc_.rows();
		VectorRealType vi(3,0);
		VectorRealType vj(3,0);

		ComplexType phase = getUnitPhase(ind);
		for (SizeType i = 0; i < lda; ++i) {
			buildVector(vi,data_(),i);
			for (SizeType j = 0; j < lda; ++j) {
				buildVector(vj,data_(),j);
				sum += phase*std::real(sc_(i,j,ind))*scalarProduct(vi,vj);
			}
		}

		if (fabs(std::imag(sum)) > 1e-10) {
			throw PsimagLite::RuntimeError("Energy is complex\n");
		}

		return std::real(sum);
	}

	FieldType derivativeEnergyThisCell(SizeType ind, SizeType index) const
	{
		SizeType lda = sc_.rows();
		VectorRealType vi(3,0);
		VectorRealType vj(3,0);
		SizeType newIndex = index + 2*data_.fixedSpins();
		bool isPhi = (newIndex & 1);
		if (isPhi) newIndex--;
		newIndex = static_cast<SizeType>(0.5*newIndex);

		ComplexType phase = getUnitPhase(ind);

		buildDeltaVector(vi,data_(),newIndex,isPhi);
		ComplexType sum = 0;
		for (SizeType j = 0; j < lda; ++j) {
			buildVector(vj,data_(),j);
			sum += phase*std::real(sc_(newIndex,j,ind))*scalarProduct(vi,vj);
		}

		for (SizeType i = 0; i < lda; ++i) {
			buildVector(vi,data_(),i);
			buildDeltaVector(vj,data_(),newIndex,isPhi);
			sum += phase*std::real(sc_(i,newIndex,ind))*scalarProduct(vi,vj);
		}

		if (fabs(std::imag(sum)) > 1e-10) {
			throw PsimagLite::RuntimeError("Energy is complex\n");
		}

		return std::real(sum);
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

	void buildDeltaVector(VectorRealType& dst,
	                      const VectorRealType& src,
	                      SizeType i,
	                      SizeType isPhi) const
	{
		SizeType fixedSpins = data_.fixedSpins();
		if (i < fixedSpins) {
			dst[0] = dst[1] = dst[2] = 0;
			return;
		}

		SizeType ii = i - fixedSpins;
		SizeType x = ii*2;
		assert(dst.size() == 3);
		assert(src.size() > x+1);
		if (!isPhi) {
			dst[0] = cos(src[x])*cos(src[1+x]);
			dst[1] = cos(src[x])*sin(src[1+x]);
			dst[2] = -sin(src[x]);
			return;
		}

		dst[0] = -sin(src[x])*sin(src[1+x]);
		dst[1] = sin(src[x])*cos(src[1+x]);
		dst[2] = 0.0;
	}

	void canonicalAngles(FieldType* data, SizeType n) const
	{
		for (SizeType i = 0; i < n; ++i)
			data[i] = canonicalAngle(data[i],(i&1) > 0);
	}

	FieldType canonicalAngle(const FieldType& value,bool isPhi) const
	{
		FieldType x = value;
		if (isPhi) {
			while (x > 2*M_PI) x -= 2*M_PI;
		}

		while (x < 0) x += 2*M_PI;

		return x;
	}

	ComplexType getUnitPhase(SizeType ind) const
	{
		VectorRealType nvector(3,0.0);
		for (SizeType i = 0; i < 3; ++i) {
			nvector[i] = sc_.nmatrix(ind,i);
		}

		RealType tmp = scalarProduct(nvector,qvector_);
		return ComplexType(cos(tmp),sin(tmp));
	}

	const VectorRealType& qvector_;
	SpaceConnectorsType sc_;
	ConfigurationType data_;
};

}

#endif // ENERGY_NC_FUNCTION_H

