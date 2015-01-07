#ifndef MATRIXSPACE_NON_COLLINEAR_H
#define MATRIXSPACE_NON_COLLINEAR_H
#include "MatrixSpaceCommon.h"

namespace yasw {

template<typename RealType_, typename ComplexOrRealType_>
class MatrixSpaceNonCollinear {

public:

	typedef RealType_ RealType;
	typedef std::complex<RealType> ComplexType;
	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef MatrixSpaceCommon<RealType,ComplexOrRealType> MatrixSpaceCommonType;
	typedef typename MatrixSpaceCommonType::SpaceConnectorsType SpaceConnectorsType;
	typedef typename MatrixSpaceCommonType::AnglesType AnglesType;
	typedef typename AnglesType::MatrixType MatrixRealType;
	typedef typename PsimagLite::Vector<MatrixRealType>::Type VectorMatrixRealType;
	typedef typename MatrixSpaceCommonType::MatrixComplexOrRealType
	MatrixComplexOrRealType;
	typedef typename SpaceConnectorsType::VectorComplexOrRealMatrixType
	VectorComplexOrRealMatrixType;

	MatrixSpaceNonCollinear(PsimagLite::String jfile,PsimagLite::String afile)
	    : common_(jfile,afile),data_(common_.size()),u_(common_.rows())
	{
		SizeType lda = common_.rows();
		for (SizeType i = 0; i < lda; ++i) {
			fillUmatrix(u_[i]);
		}

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

	PsimagLite::String nvector(SizeType ind) const
	{
		return common_.nvector(ind);
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
		RealType taui = u_[i](2,2);
		ComplexType xiistar = std::conj(xi(i));
		for (SizeType ind = 0; ind < common_.size(); ++ind) {
			for (SizeType j = 0; j < lda; ++j) {
				RealType factor = std::real(xiistar * xi(j)) + taui * u_[j](2,2);
				sum += common_.J(i,j,ind)*rotatedA(i,j)*factor;
			}
		}

		return -sum;
	}

	ComplexType xi(SizeType i) const
	{
		ComplexType x(u_[i](2,0),-u_[i](2,1));
		return 0.5*x;
	}

	RealType rotatedA(SizeType i, SizeType j) const
	{
		RealType x = fabs(common_.theta(i)-common_.theta(j));
		return (x < M_PI*0.5) ? 1.0 : -1.0;
	}

	void fillUmatrix(MatrixRealType& u) const
	{

	}

	MatrixSpaceCommonType common_;
	VectorComplexOrRealMatrixType data_;
	VectorMatrixRealType u_;
}; // MatrixSpaceNonCollinear

} // namespace yasw
#endif // MATRIXSPACE_NON_COLLINEAR_H

