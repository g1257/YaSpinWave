#ifndef MATRIXSPACE_NON_COLLINEAR_H
#define MATRIXSPACE_NON_COLLINEAR_H
#include "MatrixSpaceCommon.h"

namespace yasw {

template<typename ComplexOrRealType_>
class MatrixSpaceNonCollinear {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType_>::Type RealType;
	typedef std::complex<RealType> ComplexType;
	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef MatrixSpaceCommon<ComplexOrRealType> MatrixSpaceCommonType;
	typedef typename MatrixSpaceCommonType::SpaceConnectorsType SpaceConnectorsType;
	typedef typename MatrixSpaceCommonType::AnglesType AnglesType;
	typedef typename AnglesType::MatrixType MatrixRealType;
	typedef typename PsimagLite::Vector<MatrixRealType>::Type VectorMatrixRealType;
	typedef typename MatrixSpaceCommonType::MatrixComplexOrRealType
	MatrixComplexOrRealType;
	typedef typename SpaceConnectorsType::VectorComplexOrRealMatrixType
	VectorComplexOrRealMatrixType;

	MatrixSpaceNonCollinear(PsimagLite::String jfile,
	                        PsimagLite::String afile,
	                        bool verbose,
	                        bool altRotation)
	    : common_(jfile,afile,verbose),
	      data_(common_.size()),
	      u_(common_.rows()),
	      alt_(altRotation)
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

				ComplexOrRealType tmp = (alt_) ? -gcoeff(i,j,ind,1) :  aplus(i,j,ind);
				if (std::norm(tmp)<1e-10) tmp = 0;
				m(i,j) = m(i+lda,j+lda)  = tmp;

				tmp = (alt_) ? -gcoeff(i,j,ind,-1) :  bplus(i,j,ind);
				ComplexOrRealType tmp2 = (alt_) ? -std::conj(gcoeff(j,i,ind,-1)) :  
				std::conj(bplus(j,i,ind));
				if (std::norm(tmp)<1e-10) tmp = 0;
				if (std::norm(tmp2)<1e-10) tmp2 = 0;
				m(i,j+lda) = tmp;
				m(i+lda,j) = tmp2;

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
				if (alt_) factor = diagonalFactor(i,j);
				sum += common_.J(i,j,ind)*factor;
			}
		}

		return -sum;
	}

	ComplexType aplus(SizeType i, SizeType j, SizeType ind) const
	{
		ComplexType tmp = std::conj(alpha(i)) * alpha(j);
		tmp += std::conj(beta(j))*beta(i);
		tmp += 2.0*std::conj(chi(i))*chi(j);
		return tmp * common_.J(i,j,ind);
	}

	ComplexType bplus(SizeType i, SizeType j, SizeType ind) const
	{
		ComplexType tmp = alpha(i) * std::conj(beta(j));
		tmp += std::conj(beta(i))*alpha(j);
		tmp += 2.0*std::conj(chi(i))*std::conj(chi(j));
		return tmp * common_.J(i,j,ind);
	}

	ComplexType xi(SizeType i) const
	{
		return ComplexType(u_[i](0,2),u_[i](1,2));
	}

	void fillUmatrix(SizeType i)
	{
		assert(i < u_.size());
		MatrixRealType& u = u_[i];
		assert(u.n_row() == u.n_col());
		assert(u.n_row() == 3);
		RealType theta = common_.theta(i);
		RealType phi = common_.phi(i);
		if (!alt_) {
			u(0,0) = cos(theta)*cos(phi);
			u(0,1) = -sin(phi);
			u(0,2) = -sin(theta)*cos(phi);
			u(1,0) = cos(theta)*sin(phi);
			u(1,1) = cos(phi);
			u(1,2) = -sin(theta)*sin(phi);
			u(2,0) = sin(theta);
			u(2,1) = 0;
			u(2,2) = cos(theta);
		} else {
			u(0,0) = cos(theta)*cos(phi);
			u(0,1) = cos(theta)*sin(phi);
			u(0,2) = -sin(theta);
			u(1,0) = -sin(phi);
			u(1,1) = cos(phi);
			u(1,2) = 0;
			u(2,0) = sin(theta)*cos(phi);
			u(2,1) = sin(theta)*sin(phi);
			u(2,2) = cos(theta);
		}
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

	ComplexType gcoeff(SizeType i, SizeType j, SizeType ind, int sign) const
	{
		MatrixRealType uinverse = u_[j];
		inverse(uinverse);
		MatrixRealType fmatrix = u_[i] * uinverse;
		RealType re = fmatrix(0,0) + fmatrix(1,1)*sign;
		RealType im = fmatrix(0,1) - fmatrix(1,0)*sign;
		ComplexOrRealType factor = -0.5*common_.J(i,j,ind);
		return factor*ComplexType(re,-im);
	}

	RealType diagonalFactor(SizeType i, SizeType j) const
	{
		assert(alt_);
		MatrixRealType uinverse = u_[j];
		inverse(uinverse);
		MatrixRealType fmatrix = u_[i] * uinverse;
		return fmatrix(2,2);
	}

	MatrixSpaceCommonType common_;
	VectorComplexOrRealMatrixType data_;
	VectorMatrixRealType u_;
	bool alt_;
}; // MatrixSpaceNonCollinear

} // namespace yasw
#endif // MATRIXSPACE_NON_COLLINEAR_H

