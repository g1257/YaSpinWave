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
			u_[i].resize(3,3);
			fillUmatrix(i);
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
				RealType type = rotatedA(i,j);
				int aRotI = (common_.theta(i) < M_PI*0.5) ? 1 : -1;
				if (type < 0) {
					m(i,j+lda) = m(i+lda,j) = bminus(aRotI,i,j,ind);
					m(i,j) = m(i+lda,j+lda) = aminus(aRotI,i,j,ind);
				} else {
					m(i,j) = m(i+lda,j+lda) = aplus(aRotI,i,j,ind);
					m(i,j+lda) = m(i+lda,j) = bplus(aRotI,i,j,ind);
				}

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
				sum += common_.J(i,j,ind)*factor;
			}
		}

		return -sum;
	}

	ComplexType aplus(int aRotI, SizeType i, SizeType j, SizeType ind) const
	{
		if (aRotI < 0) return std::conj(aplus(1,i,j,ind));
		ComplexType tmp = std::conj(alpha(i)) * alpha(j);
		tmp += std::conj(beta(j))*beta(i);
		tmp += 2.0*std::conj(chi(i))*chi(j);
		return tmp * common_.J(i,j,ind);
	}

	ComplexType bplus(int aRotI, SizeType i, SizeType j, SizeType ind) const
	{
		if (aRotI < 0) return std::conj(bplus(1,i,j,ind));
		ComplexType tmp = std::conj(alpha(i)) * beta(j);
		tmp += std::conj(beta(i))*alpha(j);
		tmp += 2.0*chi(i)*chi(j);
		return 0.5*tmp * common_.J(i,j,ind);
	}

	ComplexType aminus(int aRotI, SizeType i, SizeType j, SizeType ind) const
	{
		if (aRotI < 0) return std::conj(aminus(1,i,j,ind));
		ComplexType tmp = bplus(1,i,j,ind);
		return tmp + std::conj(tmp);
	}

	ComplexType bminus(int aRotI, SizeType i, SizeType j, SizeType ind) const
	{
		return aplus(aRotI,i,j,ind);
	}

	ComplexType xi(SizeType i) const
	{
		return ComplexType(u_[i](0,2),u_[i](1,2));
	}

	RealType rotatedA(SizeType i, SizeType j) const
	{
		RealType x = fabs(common_.theta(i)-common_.theta(j));
		return (x < M_PI*0.5) ? 1.0 : -1.0;
	}

	void fillUmatrix(SizeType i)
	{
		assert(i < u_.size());
		MatrixRealType& u = u_[i];
		assert(u.n_row() == u.n_col());
		assert(u.n_row() == 3);
		RealType theta = common_.theta(i);
		RealType phi = common_.phi(i);
		u(0,0) = cos(theta)*cos(phi);
		u(0,1) = cos(theta)*sin(phi);
		u(0,2) = -sin(theta);
		u(1,0) = -sin(phi);
		u(1,1) = cos(phi);
		u(1,2) = 0;
		u(2,0) = sin(theta) * cos(phi);
		u(2,1) = sin(theta) * sin(phi);
		u(2,2) = cos(theta);
	}

	ComplexType alpha(SizeType i) const
	{
		return alphaOrBeta(i,1,-1);
	}

	ComplexType beta(SizeType i) const
	{
		return alphaOrBeta(i,-1,1);
	}

	ComplexType alphaOrBeta(SizeType i, int s1, int s2) const
	{
		RealType re = u_[i](0,0) + u_[i](1,1)*s1;
		RealType im = u_[i](1,0) + u_[i](0,1)*s2;
		return 0.5*ComplexType(re,im);
	}

	ComplexType chi(SizeType i) const
	{
		return 0.5*ComplexType(u_[i](2,0),-u_[i](2,1));
	}

	MatrixSpaceCommonType common_;
	VectorComplexOrRealMatrixType data_;
	VectorMatrixRealType u_;
}; // MatrixSpaceNonCollinear

} // namespace yasw
#endif // MATRIXSPACE_NON_COLLINEAR_H
