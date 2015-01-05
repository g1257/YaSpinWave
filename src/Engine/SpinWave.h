#ifndef SPINWAVE_H
#define SPINWAVE_H
#include "String.h"
#include <iostream>
#include "SpaceConnectors.h"
#include "Angles.h"

namespace yasw {

template<typename RealType, typename ComplexOrRealType>
class SpinWave {

	typedef yasw::SpaceConnectors<RealType,ComplexOrRealType> SpaceConnectorsType;
	typedef yasw::Angles<RealType> AnglesType;

	typedef typename SpaceConnectorsType::MatrixComplexOrRealType
	MatrixComplexOrRealType;
	typedef typename SpaceConnectorsType::VectorComplexOrRealMatrixType
	VectorComplexOrRealMatrixType;

public:

	SpinWave(PsimagLite::String jfile,PsimagLite::String afile)
	    : sc_(jfile),a_(afile),data_(sc_.size())
	{
		SizeType lda = sc_.rows();
		for (SizeType i = 0; i < sc_.size(); ++i) {
			MatrixComplexOrRealType m(2*lda,2*lda);
			fillThisMatrix(m,i);
			data_[i] = m;
		}
	}

	void printSpaceMatrices(std::ostream& os)
	{
		std::cerr<<"#Here are the "<<data_.size()<< " Hamiltonian-factor matrices.\n";
		std::cerr<<"-------------------------------------------------------------\n";
		os<<"\n";
		assert(data_.size() > 0);
		os<<data_.size()<<" "<<data_[0].n_row()<<"\n";
		for (SizeType i = 0; i < data_.size(); ++i) {
			os<<sc_.nvector(i)<<"\n";
			printMatrix(os,data_[i],1);
		}
	}

	void printDynamicMatrix(std::ostream& os)
	{
		std::cerr<<"#Here are the "<<data_.size()<< "Dynamic-factor matrices.\n";
		std::cerr<<"-------------------------------------------------------------\n";
		os<<"\n";
		assert(data_.size() > 0);
		os<<data_.size()<<" "<<data_[0].n_row()<<"\n";
		for (SizeType i = 0; i < data_.size(); ++i) {
			os<<sc_.nvector(i)<<"\n";
			printMatrix(os,data_[i],-1);
		}
	}

private:

	void fillThisMatrix(MatrixComplexOrRealType& m, SizeType ind) const
	{
		SizeType lda = sc_.rows();
		for (SizeType i = 0; i < lda; ++i) {
			for (SizeType j = 0; j < lda; ++j) {

				if (a_(i,j) < 0)
					m(i,j+lda) = m(i+lda,j) = sc_(i,j,ind);

				if (a_(i,j) > 0)
					m(i,j) = m(i+lda,j+lda) = sc_(i,j,ind);

				if (sc_.isCentralCell(ind) && i == j)
					m(i,i) = m(i+lda,i+lda) = diagonal(i);

			}
		}
	}


	ComplexOrRealType diagonal(SizeType i) const
	{
		ComplexOrRealType sum = 0;
		SizeType lda = sc_.rows();
		for (SizeType ind = 0; ind < sc_.size(); ++ind) {
			for (SizeType j = 0; j < lda; ++j) {
				sum += sc_(i,j,ind)*a_(i,j);
			}
		}

		return sum;
	}

	void printMatrix(std::ostream& os,
	                 const MatrixComplexOrRealType& m,
	                 RealType sign) const
	{
		assert(!(m.n_row() & 1));
		SizeType nrowOver2 = static_cast<SizeType>(m.n_row()*0.5);
		for (SizeType i = 0; i < m.n_row(); ++i) {
			for (SizeType j=0; j < m.n_col(); ++j) {
				ComplexOrRealType tmp = m(i,j);
				if (i >= nrowOver2) tmp *= sign;
				os<<tmp<<" ";
			}
			os<<"\n";
		}
	}

	SpaceConnectorsType sc_;
	AnglesType a_;
	VectorComplexOrRealMatrixType data_;
}; // SpinWave

}
#endif // SPINWAVE_H
