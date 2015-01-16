#ifndef MATRIXSPACE_COMMON_H
#define MATRIXSPACE_COMMON_H
#include "SpaceConnectors.h"
#include "Angles.h"

namespace yasw {

template<typename ComplexOrRealType>
class MatrixSpaceCommon {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef SpaceConnectors<ComplexOrRealType> SpaceConnectorsType;
	typedef Angles<RealType> AnglesType;
	typedef typename SpaceConnectorsType::MatrixComplexOrRealType
	MatrixComplexOrRealType;
	typedef typename SpaceConnectorsType::VectorComplexOrRealMatrixType
	VectorComplexOrRealMatrixType;

	MatrixSpaceCommon(PsimagLite::String jfile,
	                  PsimagLite::String afile,
	                  bool verbose)
	    : sc_(jfile,verbose),a_(afile,verbose)
	{}

	SizeType size() const { return sc_.size(); }

	PsimagLite::String nvector(SizeType ind) const
	{
		return sc_.nvector(ind);
	}

	SizeType rows() const { return sc_.rows(); }

	RealType a(SizeType i, SizeType j) const { return a_(i,j); }

	const RealType& theta(SizeType i) const {return a_.theta(i); }

	const RealType& phi(SizeType i) const {return a_.phi(i); }

	ComplexOrRealType J(SizeType i, SizeType j,SizeType ind) const
	{
		return sc_(i,j,ind);
	}

	bool isCentralCell(SizeType ind) const { return sc_.isCentralCell(ind); }

private:

	SpaceConnectorsType sc_;
	AnglesType a_;
}; // MatrixSpaceCommon

} // namespace yasw
#endif // MATRIXSPACE_COMMON_H

