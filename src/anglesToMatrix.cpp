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
	std::cerr<<"\t-c Use algorithm for collinear\n";
	std::cerr<<"\t-H Compute Hamiltonian only\n";
	std::cerr<<"\t-A use alternative rotation\n";
}

template<typename MatrixSpaceType>
void main2(PsimagLite::String jfile,
           PsimagLite::String afile,
           bool hOnly,
           bool altRotation)
{
	yasw::SpinWave<MatrixSpaceType> sw(jfile,afile,altRotation);

	if (hOnly) {
		sw.printSpaceMatrices(std::cout);
		return;
	}

	sw.printSpaceMatrices(std::cerr);
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
	bool altRotation = false;

	while ((opt = getopt(argc, argv,"j:a:cHA")) != -1) {
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
		case 'A':
			altRotation = true;
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
		main2<MatrixSpaceCollinearType>(jfile,afile,hOnly,altRotation);
	} else {
		main2<MatrixSpaceNonCollinearType>(jfile,afile,hOnly,altRotation);
	}
}

