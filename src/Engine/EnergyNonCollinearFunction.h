#ifndef ENERGY_NC_FUNCTION_H
#define ENERGY_NC_FUNCTION_H
#include "SpaceConnectors.h"
#include "Vector.h"
#include "Minimizer.h"
#include "MinimizerParams.h"
#include "InitConfig.h"

namespace yasw {

template<typename ComplexOrRealType_>
class EnergyNonCollinearFunction {

public:

	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef yasw::SpaceConnectors<ComplexOrRealType> SpaceConnectorsType;
	typedef EnergyNonCollinearFunction<ComplexOrRealType> ThisType;
	typedef std::complex<RealType> ComplexType;
	typedef InitConfig<ComplexOrRealType_, false> InitConfigType;

	class Configuration {

	public:

		typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

		Configuration(SizeType totalSpins,
		              SizeType fixedSpins,
		              bool verbose,
		              const InitConfigType& initConfig)
		    : fixedSpins_(fixedSpins)
		{

			if (totalSpins <= fixedSpins) {
				PsimagLite::String msg("Configuration: ");
				msg += "Too many fixed spins\n";
				throw PsimagLite::RuntimeError(msg);
			}

			data_.resize(2*(totalSpins-fixedSpins),0.0);

			data_ <= initConfig;
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

		SizeType fixedSpins() const { return fixedSpins_; }

		VectorRealType& operator()() { return data_; }

		const VectorRealType& operator()() const { return data_; }

		SizeType size() const { return data_.size(); }

	private:

		VectorRealType data_;
		SizeType fixedSpins_;
	};

public:

	typedef RealType FieldType;
	typedef typename Configuration::VectorRealType VectorRealType;
	typedef Configuration ConfigurationType;
	typedef PsimagLite::Minimizer<RealType,ThisType> MinimizerType;

	EnergyNonCollinearFunction(const SpaceConnectorsType& spaceConnectors,
	                           const VectorRealType& qvector,
	                           const VectorRealType& moduli)
	    : qvector_(qvector), moduli_(moduli), sc_(spaceConnectors), fixedSpins_(0)
	{}

	int minimize(ConfigurationType& config,
	             RealType& energy,
	             const MinimizerParams<RealType>& minParams,
	             SizeType tryIndex)
	{
		fixedSpins_ = config.fixedSpins();

		if (sc_.rows() <= fixedSpins_) {
			PsimagLite::String msg("Configuration: ");
			msg += "Too many fixed spins\n";
			err(msg);
		}

		assert(config.size()+2*fixedSpins_ == 2*sc_.rows());

		MinimizerType min(*this, minParams.maxIter,minParams.verbose);

		const bool printInitConfig = (tryIndex == 0 || minParams.verbose);

		if (printInitConfig) {
			std::cerr<<"Initial config for try number "<<tryIndex<<"\n";
			config.print(std::cerr);
		}

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


		energy = operator ()(config());
		int returnStatus = (used > 0) ? 0 : 1;
		const bool printFooter = (tryIndex == 0 || minParams.verbose);

		if (!printFooter) return returnStatus; // <--- EARLY EXIT HERE

		std::cerr<<"Minimizer params\n";
		std::cerr<<minParams;
		std::cerr<<"EnergyNonCollinearFunction::minimize(): ";
		if (min.status() == MinimizerType::GSL_SUCCESS) {
			std::cerr<<" converged after ";
		} else {
			std::cerr<<"NOT CONVERGED after ";
		}

		++used;
		std::cerr<<used<<" iterations.\n";

		return returnStatus;
	}

	SizeType totalSpins() const { return sc_.rows(); }

	SizeType size() const
	{
		assert(sc_.rows() > fixedSpins_);
		return 2*(sc_.rows() - fixedSpins_);
	}

	FieldType operator()(const VectorRealType& data) const
	{
		FieldType sum = 0;
		for (SizeType i = 0; i < sc_.size(); ++i) {
			sum += energyThisCell(i, data);
		}

		return sum;
	}

	void df(VectorRealType& df, const VectorRealType& data) const
	{
		const SizeType n = data.size();
		for (SizeType i = 0; i < n; ++i) {
			df[i] = 0.0;
			for (SizeType ind = 0; ind < sc_.size(); ++ind) {
				df[i] += derivativeEnergyThisCell(ind, i, data);
			}
		}
	}

private:

	FieldType energyThisCell(SizeType ind, const VectorRealType& vdata) const
	{
		ComplexType sum = 0;
		const SizeType lda = sc_.rows();
		VectorRealType vi(3,0);
		VectorRealType vj(3,0);
		const SizeType pixelSize = sc_.pixelSize();

		ComplexType phase = getUnitPhase(ind);
		for (SizeType i = 0; i < lda; ++i) {

			buildVector(vi, vdata, i);

			for (SizeType j = 0; j < lda; ++j) {

				const RealType modulus = moduli_[i]*moduli_[j];
				buildVector(vj, vdata, j);

				if (pixelSize == 1)
					sum += phase*modulus*std::real(sc_(i, j, ind))*scalarProduct(vi, vj);
				else
					sum += phase*modulus*energyThisPixel(pixelSize, i, j, ind, vi, vj);
			}
		}

		if (fabs(std::imag(sum)) > 1e-10) {
			throw PsimagLite::RuntimeError("Energy is complex\n");
		}

		return std::real(sum);
	}

	FieldType energyThisPixel(SizeType pixelSize,
	                          SizeType row,
	                          SizeType col,
	                          SizeType ind,
	                          const VectorRealType& si,
	                          const VectorRealType& sj) const
	{
		FieldType sum = 0;
		for (SizeType x1 = 0; x1 < pixelSize; ++x1) {
			for (SizeType x2 = 0; x2 < pixelSize; ++x2) {
				sum += std::real(sc_(pixelSize*row + x1,
				                     pixelSize*col + x2,
				                     ind))
				        *si[x1]*sj[x2];
			}
		}

		return sum;
	}

	FieldType derivativeEnergyThisCell(SizeType ind,
	                                   SizeType index,
	                                   const VectorRealType& vdata) const
	{
		SizeType lda = sc_.rows();
		VectorRealType vi(3,0);
		VectorRealType vj(3,0);
		SizeType newIndex = index + 2*fixedSpins_;
		bool isPhi = (newIndex & 1);
		if (isPhi) newIndex--;
		newIndex = static_cast<SizeType>(0.5*newIndex);

		ComplexType phase = getUnitPhase(ind);

		buildDeltaVector(vi, vdata, newIndex, isPhi);
		const SizeType pixelSize = sc_.pixelSize();
		ComplexType sum = 0;
		for (SizeType j = 0; j < lda; ++j) {
			const RealType modulus = moduli_[newIndex]*moduli_[j];
			buildVector(vj, vdata, j);
			if (pixelSize == 1)
				sum += phase*modulus*std::real(sc_(newIndex,j,ind))*scalarProduct(vi,vj);
			else
				sum += phase*modulus*energyThisPixel(pixelSize, newIndex, j, ind, vi, vj);
		}

		for (SizeType i = 0; i < lda; ++i) {
			buildVector(vi, vdata, i);
			buildDeltaVector(vj, vdata, newIndex, isPhi);
			const RealType modulus = moduli_[newIndex]*moduli_[i];
			if (pixelSize == 1)
				sum += phase*modulus*std::real(sc_(i,newIndex,ind))*scalarProduct(vi,vj);
			else
				sum += phase*modulus*energyThisPixel(pixelSize, i, newIndex, ind, vi, vj);
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
		SizeType fixedSpins = fixedSpins_;
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
		const SizeType fixedSpins = fixedSpins_;
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
		RealType tmp = scalarProduct(sc_.nvector(ind), qvector_);
		return ComplexType(cos(tmp),sin(tmp));
	}

	const VectorRealType& qvector_;
	const VectorRealType& moduli_;
	const SpaceConnectorsType& sc_;
	SizeType fixedSpins_;
};

}

#endif // ENERGY_NC_FUNCTION_H

