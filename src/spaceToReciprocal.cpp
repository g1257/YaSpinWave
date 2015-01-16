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
	typedef typename MatrixReciprocalSpaceType::MatrixType MatrixType;

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

	MatrixReciprocalSpaceType m(mfile,verbose);

	MatrixReciprocalSpaceType::VectorRealType q(tokens.size());
	for (SizeType i = 0; i < q.size(); ++i) {
		q[i] = atof(tokens[i].c_str());
		std::cerr<<"q["<<i<<"]= "<<q[i]<<"\n";
	}

	MatrixType a = m(q);
	std::cerr<<a;
	std::cerr<<"-------------\n";

	MatrixType vl(10,10), vr(10,10);
	typename PsimagLite::Vector<ComplexType>::Type eigenvalues(a.n_row());
	PsimagLite::geev('N','N',a,eigenvalues,vl,vr);

	for (SizeType i = 0; i < eigenvalues.size(); ++i)
		std::cout<<eigenvalues[i]<<"\n";

	return 0;
}

