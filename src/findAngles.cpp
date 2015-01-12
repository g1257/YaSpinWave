#include <iostream>
#include <unistd.h>
#include "String.h"
#include "Minimizer.h"
#include "EnergyCollinearFunction.h"
#include "EnergyNonCollinearFunction.h"
#include "MinimizerParams.h"

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" [options] -j file \n";
	std::cerr<<"\t-v verbose\n";
	std::cerr<<"\t-c Use collinear\n";
	std::cerr<<"\tBelow options only for non collinear\n";
	std::cerr<<"\t-s seed\n";
	std::cerr<<"\t-m maxIter (max. iterations for minimizer)\n";
	std::cerr<<"\t-d delta (x advancement for minimizer)\n";
	std::cerr<<"\t-t tolerance (y tolerance for minimizer)\n";
}

template<typename EnergyFunctionType>
void main2(PsimagLite::String jfile,
           int seed,
           const yasw::MinimizerParams<double>& minParams)
{
	EnergyFunctionType energy(jfile);
	typename EnergyFunctionType::ConfigurationType minConfig(energy.size(),seed);

	energy.minimize(minConfig,minParams);
}

int main(int argc, char** argv)
{
	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef yasw::EnergyCollinearFunction<RealType,ComplexType>
	        EnergyCollinearFunctionType;
	typedef yasw::EnergyNonCollinearFunction<RealType,ComplexType>
	        EnergyNonCollinearFunctionType;

	int opt;
	PsimagLite::String jfile;
	bool collinear = false;
	int seed = 0;
	SizeType maxIter = 1000;
	RealType delta = 1e-1;
	RealType tol = 1e-4;
	bool verbose = false;
	while ((opt = getopt(argc, argv,"j:s:m:d:t:cv")) != -1) {
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
		default:
			usage(argv[0]);
			return 1;
		}
	}

	if (jfile == "") {
		usage(argv[0]);
		return 1;
	}

	yasw::MinimizerParams<RealType> minParams(maxIter,delta,tol,verbose);

	if (collinear) {
		main2<EnergyCollinearFunctionType>(jfile,seed,minParams);
	} else {
		main2<EnergyNonCollinearFunctionType>(jfile,seed,minParams);
	}
}
