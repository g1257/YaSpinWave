#include <iostream>
#include <unistd.h>
#include "String.h"
#include "Minimizer.h"
#include "EnergyCollinearFunction.h"
#include "EnergyNonCollinearFunction.h"

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" -j file [-s seed | -c]\n";
}

template<typename EnergyFunctionType>
void main2(PsimagLite::String jfile, int seed)
{
	EnergyFunctionType energy(jfile);
	typename EnergyFunctionType::ConfigurationType minConfig(energy.size(),seed);

	energy.minimize(minConfig);
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
	int seed = 1234;
	while ((opt = getopt(argc, argv,"j:s:c")) != -1) {
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
		main2<EnergyCollinearFunctionType>(jfile,seed);
	} else {
		main2<EnergyNonCollinearFunctionType>(jfile,seed);
	}
}
