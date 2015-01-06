#ifndef SPINWAVE_H
#define SPINWAVE_H
#include "String.h"
#include <iostream>
#include "MatrixSpaceCollinear.h"

namespace yasw {

template<typename RealType, typename ComplexOrRealType>
class SpinWave {

	typedef MatrixSpaceCollinear<RealType,ComplexOrRealType> MatrixSpaceType;
	typedef typename MatrixSpaceType::SpaceConnectorsType SpaceConnectorsType;
	typedef typename MatrixSpaceType::AnglesType AnglesType;
	typedef typename MatrixSpaceType::MatrixComplexOrRealType MatrixComplexOrRealType;

public:

	SpinWave(PsimagLite::String jfile,PsimagLite::String afile)
	    : matrixSpace_(jfile,afile)
	{}

	void printSpaceMatrices(std::ostream& os)
	{
		std::cerr<<"#Here are the "<<matrixSpace_.size();
		std::cerr<< " Hamiltonian-factor matrices.\n";
		std::cerr<<"-------------------------------------------------------------\n";
		os<<"\n";
		assert(matrixSpace_.size() > 0);
		os<<matrixSpace_.size()<<" "<<matrixSpace_(0).n_row()<<"\n";
		for (SizeType i = 0; i < matrixSpace_.size(); ++i) {
			os<<matrixSpace_.nvector(i)<<"\n";
			printMatrix(os,matrixSpace_(i),1);
		}
	}

	void printDynamicMatrix(std::ostream& os)
	{
		std::cerr<<"#Here are the "<<matrixSpace_.size();
		std::cerr<< "Dynamic-factor matrices.\n";
		std::cerr<<"-------------------------------------------------------------\n";
		os<<"\n";
		assert(matrixSpace_.size() > 0);
		os<<matrixSpace_.size()<<" "<<matrixSpace_(0).n_row()<<"\n";
		for (SizeType i = 0; i < matrixSpace_.size(); ++i) {
			os<<matrixSpace_.nvector(i)<<"\n";
			printMatrix(os,matrixSpace_(i),-1);
		}
	}

private:

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

	MatrixSpaceType matrixSpace_;
}; // SpinWave

}
#endif // SPINWAVE_H
