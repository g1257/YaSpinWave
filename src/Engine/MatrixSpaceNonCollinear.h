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
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

	MatrixSpaceNonCollinear(PsimagLite::String jfile,
	                        PsimagLite::String afile,
	                        PsimagLite::String spinModulusFile,
	                        SizeType pixelSize,
	                        bool verbose)
	    : common_(jfile, afile, spinModulusFile, pixelSize, verbose),
	      data_(common_.size()),
	      u_(common_.rows())
	{
		SizeType lda = common_.rows();
		for (SizeType i = 0; i < lda; ++i) {
			u_[i].resize(3,3);
			fillUmatrix(i);
		}

		VectorType centralCellAcc(2*lda);
		int centralCellIndex = -1;
		for (SizeType i = 0; i < common_.size(); ++i) {
			MatrixComplexOrRealType m(2*lda,2*lda);
			fillThisMatrix(m, centralCellAcc, i);
			data_[i] = m;
			if (common_.isCentralCell(i))
			        centralCellIndex = i;
		}

		if (centralCellIndex < 0)
			err("Could not find central cell index\n");

		assert(centralCellIndex >= 0);
		assert(static_cast<SizeType>(centralCellIndex) < data_.size());
		for (SizeType i = 0; i < lda; ++i) {
			data_[centralCellIndex](i, i) += centralCellAcc[i];
			data_[centralCellIndex](i + lda, i + lda) += centralCellAcc[i + lda];
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

	void fillThisMatrix(MatrixComplexOrRealType& m,
	                    VectorType& centralCellAcc,
	                    SizeType ind) const
	{
		switch (common_.pixelSize()) {
		case 1:
			fillThisMatrixPixel1(m, centralCellAcc, ind);
			break;
		case 3:
			fillThisMatrixPixel3(m, centralCellAcc, ind);
			break;
		default:
			err("fillThisMatrix: pixelSize must be one or three\n");
		}
	}

	void fillThisMatrixPixel3(MatrixComplexOrRealType& m,
	                          VectorType& centralCellAcc,
	                          SizeType ind) const
	{
		static const ComplexOrRealType sqrtMinusOne = ComplexOrRealType(0, 1);

		SizeType lda = common_.rows();
		for (SizeType i = 0; i < lda; ++i) {

			const RealType theta2 = common_.theta(i);
			const RealType phi2 = common_.phi(i);
			const RealType s2 = common_.modulus(i);

			for (SizeType j = 0; j < lda; ++j) {

				const RealType theta1 = common_.theta(j);
				const RealType phi1 = common_.phi(j);
				const RealType s1 = common_.modulus(j);

				const ComplexOrRealType Axx = common_.J(3*i + 0, 3*j + 0, ind);
				const ComplexOrRealType Axy = common_.J(3*i + 0, 3*j + 1, ind);
				const ComplexOrRealType Axz = common_.J(3*i + 0, 3*j + 2, ind);
				const ComplexOrRealType Ayx = common_.J(3*i + 1, 3*j + 0, ind);
				const ComplexOrRealType Ayy = common_.J(3*i + 1, 3*j + 1, ind);
				const ComplexOrRealType Ayz = common_.J(3*i + 1, 3*j + 2, ind);
				const ComplexOrRealType Azx = common_.J(3*i + 2, 3*j + 0, ind);
				const ComplexOrRealType Azy = common_.J(3*i + 2, 3*j + 1, ind);
				const ComplexOrRealType Azz = common_.J(3*i + 2, 3*j + 2, ind);

				const ComplexOrRealType a2daga1 = (sqrt(s2)*cos(phi2)*cos(theta2)/sqrt(2) -
				                                   sqrtMinusOne*sqrt(s2)*sin(phi2)/sqrt(2))*
				        (-sqrtMinusOne*Ayx*sqrt(s1)*cos(phi1)/sqrt(2) +
				         Axx*sqrt(s1)*cos(phi1)*cos(theta1)/sqrt(2) +
				         sqrtMinusOne*Axx*sqrt(s1)*sin(phi1)/sqrt(2) +
				         Ayx*sqrt(s1)*cos(theta1)*sin(phi1)/sqrt(2) -
				         Azx*sqrt(s1)*sin(theta1)/sqrt(2)) +
				        (sqrtMinusOne*sqrt(s2)*cos(phi2)/sqrt(2) +
				         sqrt(s2)*cos(theta2)*sin(phi2)/sqrt(2))*
				        (-sqrtMinusOne*Ayy*sqrt(s1)*cos(phi1)/sqrt(2) +
				         Axy*sqrt(s1)*cos(phi1)*cos(theta1)/sqrt(2) +
				         sqrtMinusOne*Axy*sqrt(s1)*sin(phi1)/sqrt(2) +
				         Ayy*sqrt(s1)*cos(theta1)*sin(phi1)/sqrt(2) -
				         Azy*sqrt(s1)*sin(theta1)/sqrt(2)) +
				        0.5*sqrtMinusOne*Ayz*sqrt(s1)*sqrt(s2)*cos(phi1)*sin(theta2) -
				        Axz*sqrt(s1)*sqrt(s2)*cos(phi1)*cos(theta1)*sin(theta2)/2. -
				        0.5*sqrtMinusOne*Axz*sqrt(s1)*sqrt(s2)*sin(phi1)*sin(theta2) -
				        Ayz*sqrt(s1)*sqrt(s2)*cos(theta1)*sin(phi1)*sin(theta2)/2. +
				        Azz*sqrt(s1)*sqrt(s2)*sin(theta1)*sin(theta2)/2.;

				const ComplexOrRealType a2a1= (sqrt(s2)*cos(phi2)*cos(theta2)/sqrt(2) +
				                               sqrtMinusOne*sqrt(s2)*sin(phi2)/sqrt(2))*
				        (-sqrtMinusOne*Ayx*sqrt(s1)*cos(phi1)/sqrt(2) +
				         Axx*sqrt(s1)*cos(phi1)*cos(theta1)/sqrt(2) +
				         sqrtMinusOne*Axx*sqrt(s1)*sin(phi1)/sqrt(2) +
				         Ayx*sqrt(s1)*cos(theta1)*sin(phi1)/sqrt(2) -
				         Azx*sqrt(s1)*sin(theta1)/sqrt(2)) +
				        (-sqrtMinusOne*sqrt(s2)*cos(phi2)/sqrt(2) +
				         sqrt(s2)*cos(theta2)*sin(phi2)/sqrt(2))*
				        (-sqrtMinusOne*Ayy*sqrt(s1)*cos(phi1)/sqrt(2) +
				         Axy*sqrt(s1)*cos(phi1)*cos(theta1)/sqrt(2) +
				         sqrtMinusOne*Axy*sqrt(s1)*sin(phi1)/sqrt(2) +
				         Ayy*sqrt(s1)*cos(theta1)*sin(phi1)/sqrt(2) -
				         Azy*sqrt(s1)*sin(theta1)/sqrt(2)) +
				        0.5*sqrtMinusOne*Ayz*sqrt(s1)*sqrt(s2)*cos(phi1)*sin(theta2) -
				        Axz*sqrt(s1)*sqrt(s2)*cos(phi1)*cos(theta1)*sin(theta2)/2. -
				        0.5*sqrtMinusOne*Axz*sqrt(s1)*sqrt(s2)*sin(phi1)*sin(theta2) -
				        Ayz*sqrt(s1)*sqrt(s2)*cos(theta1)*sin(phi1)*sin(theta2)/2. +
				        Azz*sqrt(s1)*sqrt(s2)*sin(theta1)*sin(theta2)/2.;

				const ComplexOrRealType  a1daga1 = s2*cos(theta2)*(-Azz*cos(theta1) -
				                                                   Axz*cos(phi1)*sin(theta1) -
				                                                   Ayz*sin(phi1)*sin(theta1)) +
				        s2*cos(phi2)*(-Azx*cos(theta1) -
				                      Axx*cos(phi1)*sin(theta1) -
				                      Ayx*sin(phi1)*sin(theta1))*sin(theta2) +
				        s2*sin(phi2)*(-Azy*cos(theta1) -
				                      Axy*cos(phi1)*sin(theta1) -
				                      Ayy*sin(phi1)*sin(theta1))*sin(theta2);

				const ComplexOrRealType adaga = Ayy*s1*pow(cos(phi1),2) -
				        2.0*Azz*s1*pow(cos(theta1),2) +
				        Axx*s1*pow(cos(phi1),2)*pow(cos(theta1),2) -
				        Axy*s1*cos(phi1)*sin(phi1) -
				        Ayx*s1*cos(phi1)*sin(phi1) +
				        Axy*s1*cos(phi1)*pow(cos(theta1),2)*sin(phi1) +
				        Ayx*s1*cos(phi1)*pow(cos(theta1),2)*sin(phi1) +
				        Axx*s1*pow(sin(phi1),2) +
				        Ayy*s1*pow(cos(theta1),2)*pow(sin(phi1),2) -
				        3.0*Axz*s1*cos(phi1)*cos(theta1)*sin(theta1) -
				        3.0*Azx*s1*cos(phi1)*cos(theta1)*sin(theta1) -
				        3.0*Ayz*s1*cos(theta1)*sin(phi1)*sin(theta1) -
				        3.0*Azy*s1*cos(theta1)*sin(phi1)*sin(theta1) +
				        Azz*s1*pow(sin(theta1),2) -
				        2.0*Axx*s1*pow(cos(phi1),2)*pow(sin(theta1),2) -
				        2.0*Axy*s1*cos(phi1)*sin(phi1)*pow(sin(theta1),2) -
				        2.0*Ayx*s1*cos(phi1)*sin(phi1)*pow(sin(theta1),2) -
				        2.0*Ayy*s1*pow(sin(phi1),2)*pow(sin(theta1),2);

				const ComplexOrRealType a2daga1dag = PsimagLite::conj(a2a1);

				const ComplexOrRealType a2a1dag = PsimagLite::conj(a2daga1);

				if (common_.isCentralCell(ind) && i == j) {
					m(i, i) += 0.5*adaga;
					m(i + lda, i) += a2a1;
					m(i, i + lda) += a2daga1dag;
					m(i + lda, i + lda) += 0.5*adaga;
				} else {
					m(i, j) += 0.5*a2daga1;
					m(i + lda, j) += 0.5*a2a1;
					m(i, j + lda) += 0.5*a2daga1dag;
					m(i + lda, j + lda)  += 0.5*a2a1dag;

					centralCellAcc[j] += 0.5*a1daga1;
					centralCellAcc[j + lda] +=0.5*a1daga1;
				}
			}
		}
	}

	void fillThisMatrixPixel1(MatrixComplexOrRealType& m,
	                          VectorType& centralCellAcc,
	                          SizeType ind) const
	{
		static const ComplexOrRealType sqrtMinusOne = ComplexOrRealType(0, 1);
		static const SizeType norbital = 1;

		SizeType lda = common_.rows();
		for (SizeType i = 0; i < lda; ++i) {

			const RealType theta2 = common_.theta(i);
			const RealType phi2 = common_.phi(i);

			for (SizeType j = 0; j < lda; ++j) {

				const RealType theta1 = common_.theta(j);
				const RealType phi1 = common_.phi(j);

				ComplexOrRealType a2daga1=0.5*cos(phi1)*cos(phi2)+
				        0.5*cos(theta1)*cos(theta2)*cos(phi1)*cos(phi2)+
				        0.5*sin(theta1)*sin(theta2)+
				        0.5*sqrtMinusOne*cos(theta1)*cos(phi2)*sin(phi1)+
				        0.5*sqrtMinusOne*cos(theta2)*cos(phi2)*sin(phi1)-
				        0.5*sqrtMinusOne*cos(theta1)*cos(phi1)*sin(phi2)
				        -0.5*sqrtMinusOne*cos(theta2)*cos(phi1)*sin(phi2)+
				        0.5*sin(phi1)*sin(phi2)+
				        0.5*cos(theta1)*cos(theta2)*sin(phi1)*sin(phi2);

				ComplexOrRealType a2a1 = -0.5*cos(phi1)*cos(phi2)+
				        0.5*cos(theta1)*cos(theta2)*cos(phi1)*cos(phi2)+
				        0.5*sin(theta1)*sin(theta2)-
				        0.5*sqrtMinusOne*cos(theta1)*cos(phi2)*sin(phi1)+
				        0.5*sqrtMinusOne*cos(theta2)*cos(phi2)*sin(phi1)+
				        0.5*sqrtMinusOne*cos(theta1)*cos(phi1)*sin(phi2)-
				        0.5*sqrtMinusOne*cos(theta2)*cos(phi1)*sin(phi2)-
				        0.5*sin(phi1)*sin(phi2)+
				        0.5*cos(theta1)*cos(theta2)*sin(phi1)*sin(phi2);

				ComplexOrRealType a1daga1 = -cos(theta1)*cos(theta2)-
				        cos(phi1)*cos(phi2)*sin(theta1)*sin(theta2)-
				        sin(theta1)*sin(theta2)*sin(phi1)*sin(phi2);

				ComplexOrRealType a2daga1dag = PsimagLite::conj(a2a1);

				ComplexOrRealType a2a1dag = PsimagLite::conj(a2daga1);

				m(i, j) += 0.5*common_.J(i, j, ind)*a2daga1;
				m(i + norbital, j) += 0.5*common_.J(i, j, ind)*a2a1;
				m(i, j + norbital) += 0.5*common_.J(i, j, ind)*a2daga1dag;
				m(i + norbital, j + norbital) += 0.5*common_.J(i, j, ind)*a2a1dag;

				centralCellAcc[j] += 0.5*common_.J(i, j, ind)*a1daga1;
				centralCellAcc[j + norbital] += 0.5*common_.J(i, j, ind)*a1daga1;
			}
		}
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

