#ifndef INIT_CONFIG_H
#define INIT_CONFIG_H
#include "Vector.h"
#include "PsimagLite.h"
#include "Angles.h"
#include "RandomGen.h"

namespace yasw {

template<typename ComplexOrRealType, bool>
class InitConfig {};

template<typename ComplexOrRealType>
class InitConfig<ComplexOrRealType, true> {

public:

	typedef RandomGen<ComplexOrRealType> RandomGenType;

	InitConfig(PsimagLite::String afile,
	           RandomGenType& randomGen,
	           SizeType totalSpins,
	           SizeType fixedSpins,
	           bool verbose)
	    : afile_(afile),
	      totalSpins_(totalSpins),
	      fixedSpins_(fixedSpins),
	      verbose_(verbose)
	{
		if (randomGen.needsRandom()) {
			PsimagLite::String str("InitConfig: No random support yet");
			str +=" for collinear configuration\n";
			std::cerr<<"WARNING "<<str;
		}


		if (afile != "") {
			PsimagLite::String str("InitConfig: No angles support yet");
			str +=" for collinear configuration\n";
			std::cerr<<"WARNING "<<str;
		}
	}

private:

	InitConfig(const InitConfig&) = delete;

	InitConfig& operator=(const InitConfig&) = delete;

	PsimagLite::String afile_;
	SizeType totalSpins_;
	SizeType fixedSpins_;
    bool verbose_;
};

template<typename ComplexOrRealType>
class InitConfig<ComplexOrRealType, false> {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef Angles<RealType> AnglesType;
	typedef RandomGen<ComplexOrRealType> RandomGenType;

	InitConfig(PsimagLite::String afile,
	           RandomGenType& randomGen,
	           SizeType totalSpins,
	           SizeType fixedSpins,
	           bool verbose)
	    : afile_(afile),
	      randomGen_(randomGen),
	      totalSpins_(totalSpins),
	      fixedSpins_(fixedSpins),
	      verbose_(verbose)
	{}

	friend void operator<=(VectorRealType& data,
	                       const InitConfig& initConfig)
	{
		const SizeType totalSpins = initConfig.totalSpins_;
        const SizeType fixedSpins = initConfig.fixedSpins_;
        const bool verbose = initConfig.verbose_;

		if (initConfig.randomGen_.needsRandom() && initConfig.afile_ != "") {
			PsimagLite::String msg("Configuration: Providing both ");
			msg += " seed and angles file is an error\n";
			err(msg);
		}

		if (initConfig.afile_ != "") {
			AnglesType angles(initConfig.afile_, verbose);
			if (angles.size() != totalSpins) {
				PsimagLite::String msg("Configuration: angles file");
				msg += " has different size than j file\n";
				err(msg);
			}

			if (!initConfig.isValidAngles(angles, fixedSpins)) {
				PsimagLite::String msg("Configuration: FATAL: ");
				msg += "First angle of " + initConfig.afile_ + " is not with theta=0\n";
				err(msg);
			}

			for (SizeType i = fixedSpins; i < angles.size(); ++i) {
				SizeType twiceIndex = 2*(i - fixedSpins);
				if (twiceIndex + 1 >= data.size()) break;
				data[twiceIndex] = angles.theta(i);
				data[twiceIndex+1] = angles.phi(i);
			}

			return;
		}

		initConfig.randomGen_.randomize(data);
	}

	InitConfig(const InitConfig&) = delete;

	InitConfig& operator=(const InitConfig&) = delete;

private:

	bool isValidAngles(const AnglesType& angles, SizeType fixedSpins) const
	{
		for (SizeType i = 0; i < fixedSpins; ++i)
			if (fabs(angles.theta(i))>1e-6) return false;

		return true;
	}

	PsimagLite::String afile_;
	const RandomGenType& randomGen_;
	SizeType totalSpins_;
	SizeType fixedSpins_;
    bool verbose_;
};
}
#endif // INIT_CONFIG_H
