#include <iostream>
#include <unistd.h>
#include "String.h"
#include "MatrixSpaceCollinear.h"
#include "MatrixSpaceNonCollinear.h"
#include "SpinWave.h"

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" [options] -j file  -a file \n";
	std::cerr<<"Options\n";
	std::cerr<<"\t-v Verbose\n";
	std::cerr<<"\t-c Use algorithm for collinear\n";
	std::cerr<<"\t-H Compute Hamiltonian only\n";
	std::cerr<<"\t-A use alternative rotation\n";
}

template<typename MatrixSpaceType>
void main2(PsimagLite::String jfile,
           PsimagLite::String afile,
           bool verbose,
           bool hOnly)
{
	yasw::SpinWave<MatrixSpaceType> sw(jfile,afile,verbose);

	if (hOnly) {
		sw.printSpaceMatrices(std::cout);
		return;
	}

	if (verbose) sw.printSpaceMatrices(std::cerr);
	sw.printDynamicMatrix(std::cout);
}

int main(int argc, char** argv)
{
	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef yasw::MatrixSpaceCollinear<ComplexType> MatrixSpaceCollinearType;
	typedef yasw::MatrixSpaceNonCollinear<ComplexType> MatrixSpaceNonCollinearType;

	int opt;
	PsimagLite::String jfile;
	PsimagLite::String afile;
	bool collinear = false;
	bool hOnly = false;
	bool verbose = false;

	while ((opt = getopt(argc, argv,"j:a:cHv")) != -1) {
		switch (opt) {
		case 'j':
			jfile = optarg;
			break;
		case 'a':
			afile = optarg;
			break;
		case 'c':
			collinear = true;
			break;
		case 'H':
			hOnly = true;
			break;
		case 'v':
			verbose = true;
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

	if (collinear) {
		main2<MatrixSpaceCollinearType>(jfile,afile,verbose,hOnly);
	} else {
		main2<MatrixSpaceNonCollinearType>(jfile,afile,verbose,hOnly);
	}
}

