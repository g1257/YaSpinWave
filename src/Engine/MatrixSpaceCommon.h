#ifndef MATRIXSPACE_COMMON_H
#define MATRIXSPACE_COMMON_H
#include "SpaceConnectors.h"
#include "Angles.h"
#include "SpinModulus.h"

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
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	MatrixSpaceCommon(PsimagLite::String jfile,
	                  PsimagLite::String afile,
	                  PsimagLite::String spinModulusFile,
	                  SizeType pixelSize,
	                  bool verbose)
	    : sc_(jfile, pixelSize, verbose),
	      a_(afile, verbose),
	      spinModulus_(spinModulusFile, sc_.rows())
	{}

	SizeType size() const { return sc_.size(); }

	PsimagLite::String nvectorToString(SizeType ind) const
	{
		const VectorRealType& v = sc_.nvector(ind);
		const SizeType nsize = v.size();
		PsimagLite::String str("");
		for (SizeType i = 0; i < nsize; ++i)
			str += ttos(v[i]) + " ";
		return str;
	}

	SizeType rows() const { return sc_.rows(); }

	RealType a(SizeType i, SizeType j) const { return a_(i,j); }

	const RealType& theta(SizeType i) const {return a_.theta(i); }

	const RealType& phi(SizeType i) const {return a_.phi(i); }

	ComplexOrRealType J(SizeType i, SizeType j,SizeType ind) const
	{
		return sc_(i,j,ind);
	}

	SizeType pixelSize() const { return sc_.pixelSize(); }

	bool isCentralCell(SizeType ind) const { return sc_.isCentralCell(ind); }

	RealType modulus(SizeType x) const
	{
		assert( x < spinModulus_().size());
		return spinModulus_()[x];
	}

private:

	SpaceConnectorsType sc_;
	AnglesType a_;
	SpinModulus<VectorRealType> spinModulus_;
}; // MatrixSpaceCommon

} // namespace yasw
#endif // MATRIXSPACE_COMMON_H

