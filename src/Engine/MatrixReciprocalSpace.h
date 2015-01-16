#ifndef MATRIXRECIPROCALSPACE_H
#define MATRIXRECIPROCALSPACE_H
#include "String.h"
#include "Matrix.h"
#include "SpaceConnectors.h"

namespace yasw {

template<typename ComplexOrRealType>
class MatrixReciprocalSpace {

	typedef yasw::SpaceConnectors<ComplexOrRealType> SpaceConnectorsType;

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename SpaceConnectorsType::MatrixComplexOrRealType MatrixType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	MatrixReciprocalSpace(PsimagLite::String mfile, bool verbose)
	    : sc_(mfile,verbose)
	{}

	MatrixType operator()(const VectorRealType& q)
	{
		SizeType cells = sc_.size();
		SizeType rows = sc_.rows();
		MatrixType m(rows,rows);
		for (SizeType n = 0; n < cells; ++n) {
			procThisCell(m,n);
		}

		return m;
	}

private:

	void procThisCell(MatrixType& m,SizeType n)
	{
		ComplexOrRealType wqn = 1.0;
		m += wqn * sc_.getMatrix(n);
	}

	SpaceConnectorsType sc_;
	MatrixType data_;
};

} // namespace yasw

#endif // MATRIXRECIPROCALSPACE_H

