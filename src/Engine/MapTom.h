#ifndef MAPTOM_H
#define MAPTOM_H
#include "Matrix.h"

namespace yasw {

template<typename ComplexOrRealType>
class MapTom {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Vector<int>::Type VectorIntType;

	class RnIndex {

	};

	typedef RnIndex RnIndexType;
	typedef typename PsimagLite::Vector<RnIndex>::Type VectorRnIndexType;

	MapTom(PsimagLite::String file)
	{
		err("MapTom: Needs to read file\n");
	}

	RnIndexType operator()(const RnIndexType& RN) const
	{
		if(RN.second() >= vecRnIndex_.size())
			err("RnIndexType::operator(): N=" + ttos(RN.second()) + " is not part of the map.\n");

		VectorRealType vecTmp = vecRnIndex_[RN.second()].first()+alphaNalphaS_*RN.first();
		SizeType index = vecRnIndex_[RN.second()].second();
		return RnIndexType(vecTmp, index);
	}

private:

	VectorRnIndexType vecRnIndex_; // this is now contained instead of inherited
	MatrixRealType alphaNalphaS_;
	VectorRealType tauS_;

};
}
#endif // MAPTOM_H
