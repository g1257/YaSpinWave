#ifndef QVECTORS_H
#define QVECTORS_H
#include "Vector.h"
#include <fstream>

namespace yasw {

template<typename MatrixRealType>
class Qvectors {

public:

	typedef typename MatrixRealType::value_type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	template<typename InputNgType>
	Qvectors(InputNgType& io)
	{
		bool hasQvector = false;
		bool hasQvector0 = false;

		VectorRealType v;
		SizeType rows = 0;

		try {
			io.read(v, "Qvector");
			hasQvector = true;
		} catch (std::exception&) {}

		try {
			io.readline(rows, "Qvectors");
			hasQvector0 = true;
		} catch (std::exception&) {}

		//bool valid = (hasQvector ^ hasQvector0);
		//if (!valid)
		//	err("Qvectors: If Qvector appears then Qvectors must not\n");

		if (hasQvector) {
			const SizeType cols = v.size();
			qmatrix_.resize(1, cols);
			for (SizeType i = 0; i < cols; ++i)
				qmatrix_(0, i) = v[i];
			return;
		}

		if (!hasQvector0 || rows == 0) return;

		SizeType cols = 0;
		for (SizeType row = 0; row < rows; ++row) {

			io.read(v, "qvector" + ttos(row));

			if (row == 0) {
				cols = v.size();
				qmatrix_.resize(rows, cols);
			} else {
				err("All spins must have the same length for qvector\n");
			}

			for (SizeType i = 0; i < cols; ++i)
				qmatrix_(0, i) = v[i];
		}

	}

	Qvectors(PsimagLite::Vector<PsimagLite::String>::Type& tokens,
	         PsimagLite::String qfile,
	         bool verbose)
	{
		const SizeType tsize = tokens.size();
		if (tsize > 0 && qfile != "")
			err("Qvectors::ctor(): FATAL: -q and -Q cannot be BOTH specified\n");

		if (tsize == 0 && qfile == "")
			return;

		if (tsize > 0) {
			qmatrix_.resize(1, tsize);
			if (verbose) std::cerr<<"q vectors independent of spin\n";
			for (SizeType i = 0; i < tsize; ++i) {
				qmatrix_(0, i) = atof(tokens[i].c_str());
				if (verbose) std::cerr<<"q["<<i<<"]= "<<qmatrix_(0, i)<<"\n";
			}
		} else {
			std::ifstream fin(qfile.c_str());
			if (!fin || fin.bad() || !fin.good())
				err("Cannot open file " + qfile + " for reading\n");

			SizeType nrows = 0;
			SizeType ncols = 0;

			PsimagLite::String temp;
			fin>>temp; // discard first word that must contain version number
			fin>>nrows;
			fin>>ncols;

			if (verbose) {
				std::cerr<<"Qvectors depend on spin: "<<nrows<<" "<<ncols<<"\n";
			}

			for (SizeType row = 0; row < nrows; ++row)
				for (SizeType col = 0; col < ncols; ++col)
					fin>>qmatrix_(row, col);
		}
	}

	const MatrixRealType& operator()() const
	{
		return qmatrix_;
	}

private:

	MatrixRealType qmatrix_;
};
}
#endif // QVECTORS_H
