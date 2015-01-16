#include "MatrixReciprocalSpace.h"
#include <iostream>
#include <unistd.h>
#include <complex>
#include "String.h"
#include "Tokenizer.h"

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" -m file -q q0,q1[,q2] [-v]\n";
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
	PsimagLite::String str;

	while ((opt = getopt(argc, argv,"m:q:v")) != -1) {
		switch (opt) {
		case 'm':
			mfile = optarg;
			break;
		case 'v':
			verbose = true;
			break;
		case 'q':
			str = optarg;
			PsimagLite::tokenizer(str,tokens,delimiter);
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

	MatrixReciprocalSpaceType matrixReciprocal(mfile,verbose);

	MatrixReciprocalSpaceType::VectorRealType q(tokens.size());
	for (SizeType i = 0; i < q.size(); ++i) {
		q[i] = atof(tokens[i].c_str());
		if (verbose) std::cerr<<"q["<<i<<"]= "<<q[i]<<"\n";
	}

	const VectorRealType& dispersion = matrixReciprocal.dispersion(q);
	if (dispersion.size() != 2) {
		PsimagLite::String msg(argv[0]);
		msg += ": Dispersion size is "+ ttos(dispersion.size()) + "for this q.\n";
		throw PsimagLite::RuntimeError(msg);
	}

	std::cout<<dispersion[0]<<" "<<dispersion[1]<<"\n";

	return 0;
}
