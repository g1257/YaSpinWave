#include <iostream>
#include <unistd.h>
#include "String.h"
#include "Minimizer.h"
#include "EnergyFunction.h"

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" -j file\n";
}

int main(int argc, char** argv)
{
	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef yasw::EnergyFunction<RealType,ComplexType> EnergyFunctionType;
	typedef EnergyFunctionType::LongSizeType LongSizeType;

	int opt;
	PsimagLite::String jfile;

	while ((opt = getopt(argc, argv,"j:")) != -1) {
		switch (opt) {
		case 'j':
			jfile = optarg;
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

	EnergyFunctionType energy(jfile);
	LongSizeType minConfig = 0;

	energy.minimize(minConfig);
	std::cerr<<"#Found config for minimization: "<<minConfig<<"\n";
}
