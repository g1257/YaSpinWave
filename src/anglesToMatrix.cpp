#include <iostream>
#include <unistd.h>
#include "MatrixSpaceCollinear.h"
#include "MatrixSpaceNonCollinear.h"
#include "SpinWave.h"

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" [options] -j file  -a file [-M modulusFile] ";
	std::cerr<<" [-P pixelSize=1] \n";
	std::cerr<<"Options\n";
	std::cerr<<"\t-v Verbose\n";
	std::cerr<<"\t-c Use algorithm for collinear\n";
	std::cerr<<"\t-H Compute Hamiltonian only\n";
}

template<typename MatrixSpaceType>
void main2(PsimagLite::String jfile,
           PsimagLite::String afile,
           PsimagLite::String spinModulusFile,
           SizeType pixelSize,
           bool verbose,
           bool hOnly)
{
	yasw::SpinWave<MatrixSpaceType> sw(jfile, afile, spinModulusFile, pixelSize, verbose);

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
	SizeType pixel = 1;
	PsimagLite::String spinModulusFile;

	while ((opt = getopt(argc, argv,"j:a:M:P:cHv")) != -1) {
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
		case 'M':
			spinModulusFile = optarg;
			break;
		case 'P':
			pixel = atoi(optarg);
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
		main2<MatrixSpaceCollinearType>(jfile, afile, spinModulusFile, pixel, verbose, hOnly);
	} else {
		main2<MatrixSpaceNonCollinearType>(jfile, afile, spinModulusFile, pixel, verbose, hOnly);
	}
}

