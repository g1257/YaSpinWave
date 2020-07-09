#include "MatrixReciprocalSpace.h"
#include <iostream>
#include <unistd.h>
#include <complex>
#include "PsimagLite.h"

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" -m file -q q0,q1[,q2] [-P pixelSize=1] [-v]\n";
}

int main(int argc, char** argv)
{
	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef yasw::MatrixReciprocalSpace<ComplexType> MatrixReciprocalSpaceType;
	typedef MatrixReciprocalSpaceType::VectorRealType VectorRealType;

	int opt;
	PsimagLite::String mfile;
	bool verbose = false;
	PsimagLite::Vector<PsimagLite::String>::Type tokens;
	PsimagLite::String delimiter = ",";
	SizeType pixel = 1;

	while ((opt = getopt(argc, argv,"m:q:P:v")) != -1) {
		switch (opt) {
		case 'm':
			mfile = optarg;
			break;
		case 'v':
			verbose = true;
			break;
		case 'q':
			PsimagLite::split(tokens, optarg, delimiter);
			break;
		case 'P':
			pixel = atoi(optarg);
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}

	if (mfile == "" || tokens.size() < 2) {
		usage(argv[0]);
		return 1;
	}

	MatrixReciprocalSpaceType matrixReciprocal(mfile, pixel, verbose);

	MatrixReciprocalSpaceType::VectorRealType q(tokens.size());
	for (SizeType i = 0; i < q.size(); ++i) {
		q[i] = atof(tokens[i].c_str());
		if (verbose) std::cerr<<"q["<<i<<"]= "<<q[i]<<"\n";
	}

	const VectorRealType& dispersion = matrixReciprocal.dispersion(q);
	for (SizeType i = 0; i < dispersion.size(); ++i)
		std::cout<<dispersion[i]<<" ";
	std::cout<<"\n";

	return 0;
}

