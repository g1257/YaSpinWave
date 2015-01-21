#include <iostream>
#include <unistd.h>
#include "String.h"
#include "Minimizer.h"
#include "EnergyCollinearFunction.h"
#include "EnergyNonCollinearFunction.h"
#include "MinimizerParams.h"

void usage(const char *progName, const yasw::MinimizerParams<double>* minParams)
{
	std::cerr<<"Usage: "<<progName<<" [options] -j file \n";
	std::cerr<<"\t-v verbose\n";
	std::cerr<<"\t-p precision (default is 8)\n";
	std::cerr<<"\t-F spins (number of spins to fix, default is 1)\n";
	std::cerr<<"\t-c Use collinear\n";
	std::cerr<<"Below options only for non collinear\n";
	std::cerr<<"\t-a anglesFile (initial angles for minimizer)\n";
	std::cerr<<"\t-s seed\n";
	std::cerr<<"\t-C Use conjugate gradient\n";
	std::cerr<<"\t-m maxIter (max. iterations for minimizer)\n";
	std::cerr<<"\t-d delta (x advancement for minimizer)\n";
	std::cerr<<"\t-t tolerance (y tolerance for minimizer)\n";
	if (!minParams) return;
	std::cerr<<"Defaults are\n";
	std::cerr<<(*minParams);
}

template<typename EnergyFunctionType>
void main2(PsimagLite::String jfile,
           SizeType fixedSpins,
           PsimagLite::String afile,
           int seed,
           const yasw::MinimizerParams<double>& minParams)
{
	EnergyFunctionType energy(jfile, minParams.verbose);
	SizeType totalSpins = energy.totalSpins();
	typename EnergyFunctionType::ConfigurationType minConfig(totalSpins,
	                                                         fixedSpins,
	                                                         minParams.verbose,
	                                                         afile,
	                                                         seed);

	energy.minimize(minConfig,minParams);
}

int main(int argc, char** argv)
{
	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef yasw::EnergyCollinearFunction<ComplexType> EnergyCollinearFunctionType;
	typedef yasw::EnergyNonCollinearFunction<ComplexType> EnergyNonCollinearFunctionType;
	typedef yasw::MinimizerParams<RealType> MinimizerParamsType;
	typedef MinimizerParamsType::EnumAlgo EnumAlgo;

	int opt;
	PsimagLite::String jfile;
	bool collinear = false;
	int seed = 0;
	SizeType maxIter = 1000;
	RealType delta = 1e-1;
	RealType tol = 1e-4;
	bool verbose = false;
	RealType prec = 8;
	PsimagLite::String afile("");
	SizeType fixedSpins = 1;
	EnumAlgo algo = MinimizerParamsType::SIMPLEX;

	while ((opt = getopt(argc, argv,"j:s:m:d:t:p:a:F:cvC")) != -1) {
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
		case 't':
			tol = atof(optarg);
			break;
		case 'v':
			verbose = true;
			break;
		case 'p':
			prec = atoi(optarg);
			break;
		case 'a':
			afile = optarg;
			break;
		case 'F':
			fixedSpins = atoi(optarg);
			break;
		case 'C':
			algo = MinimizerParamsType::CONJUGATE_GRADIENT;
			break;
		default:
			usage(argv[0],0);
			return 1;
		}
	}

	std::cerr.precision(prec);
	std::cout.precision(prec);
	MinimizerParamsType minParams(algo,maxIter,delta,tol,verbose);

	if (jfile == "") {
		usage(argv[0],&minParams);
		return 1;
	}

	if (collinear) {
		main2<EnergyCollinearFunctionType>(jfile,fixedSpins,afile,seed,minParams);
	} else {
		main2<EnergyNonCollinearFunctionType>(jfile,fixedSpins,afile,seed,minParams);
	}
}

