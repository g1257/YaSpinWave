#include <iostream>
#include <unistd.h>
#include "String.h"
#include "SpinWave.h"

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" -j file  -a file\n";
}

int main(int argc, char** argv)
{
	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef yasw::SpinWave<RealType,ComplexType> SpinWaveType;

	int opt;
	PsimagLite::String jfile;
	PsimagLite::String afile;

	while ((opt = getopt(argc, argv,"j:a:")) != -1) {
		switch (opt) {
		case 'j':
			jfile = optarg;
			break;
		case 'a':
			afile = optarg;
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}

	if (jfile == "" || afile == "") {
		usage(argv[0]);
		return 1;
	}

	SpinWaveType sw(jfile,afile);

	sw.printSpaceMatrices(std::cout);
}
