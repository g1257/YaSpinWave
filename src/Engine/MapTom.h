#ifndef MAPTOM_H
#define MAPTOM_H

#include "Matrix.h"

namespace yasw {

template<typename ComplexOrRealType>
class MapTom {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;

	MapTom(PsimagLite::String file)
	{
		err("MapTom: Needs to read file\n");
	}

};
}
#endif // MAPTOM_H
