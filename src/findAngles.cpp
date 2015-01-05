#include <iostream>
#include <unistd.h>
#include "String.h"
#include "Minimizer.h"
#include "EnergyCollinearFunction.h"
#include "EnergyNonCollinearFunction.h"

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" -j file\n";
}

template<typename EnergyFunctionType>
void main2(PsimagLite::String jfile)
{
	EnergyFunctionType energy(jfile);
	typename EnergyFunctionType::ConfigurationType minConfig(0);

	energy.minimize(minConfig);
	std::cerr<<"#Found config for minimization: "<<minConfig<<"\n";
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

	while ((opt = getopt(argc, argv,"j:c")) != -1) {
		switch (opt) {
		case 'j':
			jfile = optarg;
			break;
		case 'c':
			collinear = true;
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

	if (collinear) {
		main2<EnergyCollinearFunctionType>(jfile);
	} else {
		main2<EnergyNonCollinearFunctionType>(jfile);
	}
}
