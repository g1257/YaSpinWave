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
	                        bool verbose)
	    : common_(jfile,afile,verbose),
	      data_(common_.size()),
	      u_(common_.rows())
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

				ComplexOrRealType tmp =  aplus(i,j,ind);
				if (std::norm(tmp)<1e-10) tmp = 0;
				m(i,j) = m(j+lda,i+lda)  = tmp;

				tmp =  bplus(i,j,ind);
				if (std::norm(tmp)<1e-10) tmp = 0;
				m(i,j+lda) = tmp;
				m(i+lda,j) = std::conj(tmp);

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

	ComplexType aplus(SizeType i, SizeType j, SizeType ind) const
	{
		RealType t1 = common_.theta(i);
		RealType p1 = common_.phi(i);
		RealType t2 = common_.theta(j);
		RealType p2 = common_.phi(j);
		RealType deltaPhi = p1-p2;
		RealType re = cos(deltaPhi)*(cos(t1)*cos(t2) + 1) + sin(t1)*sin(t2);
		RealType im = sin(deltaPhi)*(cos(t2) + cos(t1));
		ComplexType tmp(re,im);
		return 0.5*tmp*common_.J(i,j,ind);
	}

	ComplexType bplus(SizeType i, SizeType j, SizeType ind) const
	{
		RealType t1 = common_.theta(i);
		RealType p1 = common_.phi(i);
		RealType t2 = common_.theta(j);
		RealType p2 = common_.phi(j);
		RealType deltaPhi = p1-p2;
		RealType re = cos(deltaPhi)*(cos(t1)*cos(t2) - 1) + sin(t1)*sin(t2);
		RealType im = sin(deltaPhi)*(cos(t2) - cos(t1));
		ComplexType tmp(re,im);
		return 0.5*tmp*common_.J(i,j,ind);
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

		u(0,0) = cos(theta)*cos(phi);
		u(0,1) = -sin(phi);
		u(0,2) = -sin(theta)*cos(phi);
		u(1,0) = cos(theta)*sin(phi);
		u(1,1) = cos(phi);
		u(1,2) = -sin(theta)*sin(phi);
		u(2,0) = sin(theta);
		u(2,1) = 0;
		u(2,2) = cos(theta);
	}

	MatrixSpaceCommonType common_;
	VectorComplexOrRealMatrixType data_;
	VectorMatrixRealType u_;
}; // MatrixSpaceNonCollinear

} // namespace yasw
#endif // MATRIXSPACE_NON_COLLINEAR_H

