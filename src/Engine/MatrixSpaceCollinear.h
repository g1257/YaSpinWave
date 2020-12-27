#ifndef MATRIXSPACECOLLINEAR_H
#define MATRIXSPACECOLLINEAR_H
#include "MatrixSpaceCommon.h"

namespace yasw {

template<typename ComplexOrRealType_>
class MatrixSpaceCollinear {

public:

	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef MatrixSpaceCommon<ComplexOrRealType> MatrixSpaceCommonType;
	typedef typename MatrixSpaceCommonType::SpaceConnectorsType SpaceConnectorsType;
	typedef typename MatrixSpaceCommonType::AnglesType AnglesType;
	typedef typename MatrixSpaceCommonType::MatrixComplexOrRealType
	MatrixComplexOrRealType;
	typedef typename SpaceConnectorsType::VectorComplexOrRealMatrixType
	VectorComplexOrRealMatrixType;

	MatrixSpaceCollinear(PsimagLite::String jfile,
	                     PsimagLite::String afile,
	                     PsimagLite::String spinModulusFile,
	                     SizeType pixelSize,
	                     bool verbose)
	    : verbose_(verbose),
	      common_(jfile, afile, spinModulusFile, pixelSize, verbose),
	      data_(common_.size())
	{
		SizeType lda = common_.rows();
		for (SizeType i = 0; i < common_.size(); ++i) {
			MatrixComplexOrRealType m(2*lda,2*lda);
			fillThisMatrix(m,i);
			data_[i] = m;
		}
	}

	SizeType size() const { return common_.size(); }

	const MatrixComplexOrRealType& operator()(SizeType i) const
	{
		assert(i < data_.size());
		return data_[i];
	}

	PsimagLite::String nvectorToString(SizeType ind) const
	{
		return common_.nvectorToString(ind);
	}

private:

	void fillThisMatrix(MatrixComplexOrRealType& m, SizeType ind) const
	{
		SizeType lda = common_.rows();
		for (SizeType i = 0; i < lda; ++i) {
			for (SizeType j = 0; j < lda; ++j) {

				if (common_.a(i,j) < 0)
					m(i,j+lda) = m(i+lda,j) = common_.J(i,j,ind);

				if (common_.a(i,j) > 0)
					m(i,j) = m(i+lda,j+lda) = common_.J(i,j,ind);

				if (common_.isCentralCell(ind) && i == j)
					m(i,i) = m(i+lda,i+lda) = diagonal(i);

			}
		}
	}

	ComplexOrRealType diagonal(SizeType i) const
	{
		ComplexOrRealType sum = 0;
		SizeType lda = common_.rows();
		for (SizeType ind = 0; ind < common_.size(); ++ind) {
			for (SizeType j = 0; j < lda; ++j) {
				sum += common_.J(i,j,ind)*common_.a(i,j);
			}
		}

		return -sum;
	}

	bool verbose_;
	MatrixSpaceCommonType common_;
	VectorComplexOrRealMatrixType data_;
}; // MatrixSpaceCollinear

} // namespace yasw
#endif // MATRIXSPACECOLLINEAR_H

