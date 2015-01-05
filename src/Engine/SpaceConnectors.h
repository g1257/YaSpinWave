#ifndef SPACECONNECTORS_H
#define SPACECONNECTORS_H
#include "String.h"
#include "Matrix.h"
#include <fstream>

namespace yasw {

template<typename RealType, typename ComplexOrRealType>
class SpaceConnectors {

public:

	typedef PsimagLite::Matrix<RealType> MatrixRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixComplexOrRealType;
	typedef typename PsimagLite::Vector<MatrixComplexOrRealType>::Type
	VectorComplexOrRealMatrixType;

	SpaceConnectors(PsimagLite::String file, bool verbose = true)
	    : rows_(0),centralCellIndex_(0)
	{
		std::ifstream fin(file.c_str());
		if (!fin || fin.bad() || !fin.good()) {
			PsimagLite::String msg("SpaceConnectors: Problem reading file ");
			msg += file + "\n";
			throw PsimagLite::RuntimeError(msg);
		}

		PsimagLite::String str;
		getline(fin,str);
		std::cerr<<"SpaceConnectors::ctor(): Ignoring first line of ";
		std::cerr<<file<<"\n";
		SizeType n = 0;
		fin>>n;
		fin>>rows_;
		if (n>=100) {
			throw PsimagLite::RuntimeError("SpaceConnectors:: too many\n");
		}

		nmatrix_.resize(n,3);
		data_.resize(n);
		bool centralSeen = false;

		for (SizeType i = 0; i < n; ++i) {
			bool flag = true;
			for (SizeType j = 0; j < 3; ++j) {
				fin>>nmatrix_(i,j);
				if (nmatrix_(i,j) != 0) flag = false;
			}

			if (flag) {
				if (centralSeen) {
					PsimagLite::String msg("SpaceConnectors: Multiple cells ");
					msg += "are central and this is an error.\n";
					throw PsimagLite::RuntimeError(msg);
				}

				centralSeen = true;
				centralCellIndex_ = i;
				if (verbose)
					std::cerr<<"Found central cell index "<<centralCellIndex_<<"\n";
			}

			data_[i].resize(rows_,rows_);
			for (SizeType j = 0; j < rows_; ++j) {
				for (SizeType k = 0; k < rows_; ++k) {
					fin>>data_[i](j,k);
				}
			}
		}

		fin.close();

		if (verbose) printConnectors(std::cerr);
	}

	ComplexOrRealType operator()(SizeType i, SizeType j, SizeType ind) const
	{
		assert(ind < size());
		return data_[ind](i,j);
	}

	PsimagLite::String nvector(SizeType ind) const
	{
		assert(ind < nmatrix_.n_row());
		PsimagLite::String str("");
		for (SizeType i = 0; i < nmatrix_.n_col(); ++i)
			str += ttos(nmatrix_(ind,i)) + " ";
		return str;
	}

	bool isCentralCell(SizeType ind) const
	{
		return (ind == centralCellIndex_);
	}

	SizeType size() const { return data_.size(); }

	SizeType rows() const {return rows_; }

private:

	void printConnectors(std::ostream& os) const
	{
		os<<"#There are "<<data_.size()<<" connectors.\n";
		for (SizeType i = 0; i < data_.size(); ++i) {
			os<<data_[i];
		}
	}

	MatrixRealType nmatrix_;
	VectorComplexOrRealMatrixType data_;
	SizeType rows_;
	SizeType centralCellIndex_;
}; // SpaceConnectors

}
#endif // SPACECONNECTORS_H
