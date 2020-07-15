#include <iostream>
#include <unistd.h>
#include "Minimizer.h"
#include "EnergyCollinearFunction.h"
#include "EnergyNonCollinearFunction.h"
#include "MinimizerParams.h"
#include "PsimagLite.h"
#include "InitConfig.h"
#include "RandomGen.h"
#include "SpinModulus.h"

void usage(const char *progName, const yasw::MinimizerParams<double>* minParams)
{
	std::cerr<<"Usage: "<<progName<<" [options] -j file \n";
	std::cerr<<"\t-v verbose\n";
	std::cerr<<"\t-p precision (default is 8)\n";
	std::cerr<<"\t-F spins (number of spins to fix, default is 1)\n";
	std::cerr<<"\t-c Use collinear\n";
	std::cerr<<"Below options only for non collinear\n";
	std::cerr<<"\t-a anglesFile (initial angles for minimizer)\n";
	std::cerr<<"\t-s seed (cannot be used with anglesFiles)\n";
	std::cerr<<"\t-N tries (number of random tries; cannot be used with anglesFiles)\n";
	std::cerr<<"\t-C Use conjugate gradient\n";
	std::cerr<<"\t-m maxIter (max. iterations for minimizer)\n";
	std::cerr<<"\t-d delta (x advancement for minimizer)\n";
	std::cerr<<"\t-D delta (advancement for gradient, ignored unless using -C)\n";
	std::cerr<<"\t-t tolerance (y tolerance for minimizer)\n";
	std::cerr<<"\t-S n (save work configuration each n steps)\n";
	std::cerr<<"\t-q q0,q1,q2 (q values to use in energy function)\n";
	if (!minParams) return;
	std::cerr<<"Defaults are\n";
	std::cerr<<(*minParams);
}

template<typename EnergyFunctionType>
void main2(const typename EnergyFunctionType::SpaceConnectorsType& spaceConnectors,
           SizeType fixedSpins,
           const typename EnergyFunctionType::VectorRealType& qvector,
           const yasw::SpinModulus<typename EnergyFunctionType::VectorRealType>& spinModulus,
           PsimagLite::String afile,
           int seed,
           SizeType randomTries,
           const yasw::MinimizerParams<double>& minParams)
{
	EnergyFunctionType energy(spaceConnectors, qvector, spinModulus());
	SizeType totalSpins = energy.totalSpins();
	typedef typename EnergyFunctionType::ComplexOrRealType ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef yasw::EnergyCollinearFunction<ComplexOrRealType> EnergyCollinearFunctionType;
	typedef yasw::RandomGen<ComplexOrRealType> RandomGenType;
	enum {isCollinear = PsimagLite::TypesEqual<EnergyCollinearFunctionType,
		  EnergyFunctionType>::True};
	typedef yasw::InitConfig<ComplexOrRealType, isCollinear> InitConfigType;
	typedef typename EnergyFunctionType::ConfigurationType ConfigurationType;

	if (seed == 0 && afile == "")
		seed = 1234;

	RandomGenType randomGen(seed);

	if (randomTries > 1 && isCollinear)
		err("findAngles: random tries > 1 cannot be used for collinear\n");

	RealType minEnergy = 0;
	ConfigurationType* savedConfig = 0;
	for (SizeType i = 0; i < randomTries; ++i) {
		InitConfigType initConfig(afile, randomGen, totalSpins, fixedSpins, minParams.verbose);

		ConfigurationType minConfig(totalSpins,
		                            fixedSpins,
		                            minParams.verbose,
		                            initConfig);

		const RealType energyValue = energy.minimize(minConfig, minParams, i);

		if (i == 0 || energyValue < minEnergy) {
			minEnergy = energyValue;
			if (!savedConfig)
				savedConfig = new ConfigurationType(minConfig);
			else
				*savedConfig = minConfig;
		}

		if (!isCollinear && minParams.verbose) {
			std::cerr<<"Energy at Minimum= "<<energyValue<<"\n";
			std::cerr<<"Angles\n";
			minConfig.print(std::cerr);
		}
	}

	std::cerr<<"FINAL Energy at Minimum= "<<minEnergy<<"\n";
	std::cout<<"Angles\n";
	savedConfig->print(std::cout);
	delete savedConfig;
	savedConfig = 0;
}

int main(int argc, char** argv)
{
	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef yasw::SpaceConnectors<ComplexType> SpaceConnectorsType;
	typedef yasw::EnergyCollinearFunction<ComplexType> EnergyCollinearFunctionType;
	typedef yasw::EnergyNonCollinearFunction<ComplexType> EnergyNonCollinearFunctionType;
	typedef EnergyNonCollinearFunctionType::VectorRealType VectorRealType;
	typedef yasw::MinimizerParams<RealType> MinimizerParamsType;
	typedef MinimizerParamsType::EnumAlgo EnumAlgo;

	int opt;
	PsimagLite::String jfile;
	bool collinear = false;
	int seed = 0;
	SizeType maxIter = 1000;
	RealType delta = 1e-1;
	RealType delta2 = 1e-1;
	RealType tol = 1e-4;
	bool verbose = false;
	RealType prec = 8;
	PsimagLite::String afile;
	PsimagLite::String spinModulusFile;
	SizeType fixedSpins = 1;
	EnumAlgo algo = MinimizerParamsType::SIMPLEX;
	SizeType saveEvery = 0;
	PsimagLite::Vector<PsimagLite::String>::Type tokens;
	PsimagLite::String delimiter = ",";
	SizeType randomTries = 1;
	SizeType pixel = 1;

	while ((opt = getopt(argc, argv,"j:s:m:d:D:t:p:P:a:F:S:q:N:M:cvC")) != -1) {
		switch (opt) {
		case 'j':
			jfile = optarg;
			break;
		case 's':
			seed = atoi(optarg);
			break;
		case 'c':
			collinear = true;
			break;
		case 'm':
			maxIter = atoi(optarg);
			break;
		case 'd':
			delta = atof(optarg);
			break;
		case 'D':
			delta2 = atof(optarg);
			break;
		case 't':
			tol = atof(optarg);
			break;
		case 'v':
			verbose = true;
			break;
		case 'p':
			prec = atoi(optarg);
			break;
		case 'P':
			pixel = atoi(optarg);
			break;
		case 'M':
			spinModulusFile = optarg;
			break;
		case 'a':
			afile = optarg;
			break;
		case 'F':
			fixedSpins = atoi(optarg);
			break;
		case 'N':
			randomTries = atoi(optarg);
			break;
		case 'C':
			algo = MinimizerParamsType::CONJUGATE_GRADIENT;
			break;
		case 'S':
			saveEvery = atoi(optarg);
			break;
		case 'q':
			PsimagLite::split(tokens, optarg, delimiter);
			break;
		default:
			usage(argv[0],0);
			return 1;
		}
	}

	std::cerr.precision(prec);
	std::cout.precision(prec);
	MinimizerParamsType minParams(algo,maxIter,delta,delta2,tol,saveEvery,verbose);

	if (jfile == "") {
		usage(argv[0],&minParams);
		return 1;
	}

	PsimagLite::Vector<RealType>::Type q(tokens.size(),0);
	for (SizeType i = 0; i < q.size(); ++i) {
		q[i] = atof(tokens[i].c_str());
		if (verbose) std::cerr<<"q["<<i<<"]= "<<q[i]<<"\n";
	}


	SpaceConnectorsType spaceConnectors(jfile, pixel, verbose);
	yasw::SpinModulus<VectorRealType> spinModulus(spinModulusFile,
	                                              spaceConnectors.rows());
	if (collinear) {
		main2<EnergyCollinearFunctionType>(spaceConnectors,
		                                   fixedSpins,
		                                   q,
		                                   spinModulus,
		                                   afile,
		                                   seed,
		                                   randomTries,
		                                   minParams);
	} else {
		main2<EnergyNonCollinearFunctionType>(spaceConnectors,
		                                      fixedSpins,
		                                      q,
		                                      spinModulus,
		                                      afile,
		                                      seed,
		                                      randomTries,
		                                      minParams);
	}
}

