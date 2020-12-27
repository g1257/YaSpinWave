#include "MatrixReciprocalSpace.h"
#include <iostream>
#include <unistd.h>
#include <complex>
#include "PsimagLite.h"

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" arguments\n";
	std::cerr<<"\tMandatory arguments\n";
	std::cerr<<"\t-m H.txt (Hamiltonian file, which is output of anglesToMatrix)\n";
	std::cerr<<"\t-c Case file without extension (?)\n";
	std::cerr<<"\t-A mapfile (?)\n";
	std::cerr<<"\t-s hsfile (?)\n";
	std::cerr<<"\t-k nk\n";
	std::cerr<<"\t-w width (?)\n";
	std::cerr<<"\t-O outputfile\n";
	std::cerr<<"\t-a anglesfile (which is the output of findAngles)\n";
	std::cerr<<"\t-f cutoff\n";
	std::cerr<<"\t-M modulus file\n";
	std::cerr<<"\n";
	std::cerr<<"\tOptions are [-P pixelSize=1] [-v]\n";

}

int main(int argc, char** argv)
{
	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef yasw::MatrixReciprocalSpace<ComplexType> MatrixReciprocalSpaceType;
//	typedef MatrixReciprocalSpaceType::VectorRealType VectorRealType;
	typedef MatrixReciprocalSpaceType::ReciprocalArgs ReciprocalArgsType;
	int opt;

	ReciprocalArgsType reciprocalArgs;
	bool verbose = false;
	SizeType pixel = 1;

	while ((opt = getopt(argc, argv,"m:c:A:s:k:w:O:a:f:M:P:v")) != -1) {
		switch (opt) {
		case 'm':
			reciprocalArgs.mfile = optarg;
			break;
		case 'c':
			reciprocalArgs.casefile = optarg;
			break;
		case 'A':
			reciprocalArgs.mapfile = optarg;
			break;
		case 's':
			reciprocalArgs.hsfile = optarg;
			break;
		case 'k':
			reciprocalArgs.nk = PsimagLite::atoi(optarg);
			break;
		case 'w':
			reciprocalArgs.width = PsimagLite::atoi(optarg);
			break;
		case 'O':
			reciprocalArgs.outputfile = optarg;
			break;
		case 'a':
			reciprocalArgs.anglesfile = optarg;
			break;
		case 'f':
			reciprocalArgs.cutoff = PsimagLite::atof(optarg);
			break;
		case 'M':
			reciprocalArgs.modulusfile = optarg;
			break;
		case 'P':
			pixel = PsimagLite::atoi(optarg);
			break;
		case 'v':
			verbose = true;
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}

	if (!reciprocalArgs.check()) {
		usage(argv[0]);
		return 1;
	}

	MatrixReciprocalSpaceType matrixReciprocal(reciprocalArgs, pixel, verbose);

	PsimagLite::String str;
	matrixReciprocal.mainLoop(str);
	return 0;
}

