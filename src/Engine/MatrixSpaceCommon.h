#ifndef MATRIXSPACE_COMMON_H
#define MATRIXSPACE_COMMON_H
#include "SpaceConnectors.h"
#include "Angles.h"

namespace yasw {

template<typename RealType, typename ComplexOrRealType>
class MatrixSpaceCommon {

public:

	typedef SpaceConnectors<RealType,ComplexOrRealType> SpaceConnectorsType;
	typedef Angles<RealType> AnglesType;
	typedef typename SpaceConnectorsType::MatrixComplexOrRealType
	MatrixComplexOrRealType;
	typedef typename SpaceConnectorsType::VectorComplexOrRealMatrixType
	VectorComplexOrRealMatrixType;

	MatrixSpaceCommon(PsimagLite::String jfile,PsimagLite::String afile)
	    : sc_(jfile),a_(afile)
	{}

	SizeType size() const { return sc_.size(); }

	PsimagLite::String nvector(SizeType ind) const
	{
		return sc_.nvector(ind);
	}

	SizeType rows() const { return sc_.rows(); }

	RealType a(SizeType i, SizeType j) const { return a_(i,j); }

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

