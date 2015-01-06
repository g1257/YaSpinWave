#include <iostream>
#include <unistd.h>
#include "String.h"
#include "MatrixSpaceCollinear.h"
#include "MatrixSpaceNonCollinear.h"
#include "SpinWave.h"

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" -j file  -a file [-c]\n";
}

template<typename MatrixSpaceType>
void main2(PsimagLite::String jfile, PsimagLite::String afile)
{
	yasw::SpinWave<MatrixSpaceType> sw(jfile,afile);

	sw.printSpaceMatrices(std::cerr);
	sw.printDynamicMatrix(std::cout);
}

int main(int argc, char** argv)
{
	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef yasw::MatrixSpaceCollinear<RealType,ComplexType> MatrixSpaceCollinearType;
	typedef yasw::MatrixSpaceNonCollinear<RealType,ComplexType>
	        MatrixSpaceNonCollinearType;

	int opt;
	PsimagLite::String jfile;
	PsimagLite::String afile;
	bool collinear = false;

	while ((opt = getopt(argc, argv,"j:a:c")) != -1) {
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
		main2<MatrixSpaceCollinearType>(jfile,afile);
	} else {
		main2<MatrixSpaceNonCollinearType>(jfile,afile);
	}
}
