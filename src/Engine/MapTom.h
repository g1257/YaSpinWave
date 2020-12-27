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

	public:

		RnIndex(const VectorIntType& v, SizeType second)
		    : v_(v), second_(second)
		{}

		const VectorIntType& first() const { return v_; }

		SizeType second() const { return second_; }

	private:

		VectorIntType v_;
		SizeType second_;
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

		const VectorIntType& v1 = vecRnIndex_[RN.second()].first();
		const VectorIntType& v2 = RN.first();
		const SizeType rows = alphaNalphaS_.rows();
		const SizeType cols = alphaNalphaS_.cols();
		assert(v2.size() == v1.size());
		assert(cols == v1.size());
		VectorIntType vecTmp(rows);
		for (SizeType i = 0; i < rows; ++i) {
			RealType sum = 0;
			for (SizeType j = 0; j < cols; ++j) {
				sum += v1[j] + alphaNalphaS_(i, j)*v2[j];
			}

			vecTmp[i] = static_cast<int>(sum);
		}

		SizeType index = vecRnIndex_[RN.second()].second();
		return RnIndexType(vecTmp, index);
	}

	const MatrixRealType& alphaNalphaS() const { return alphaNalphaS_; }

	int unitPerSuper()
	{
		throw PsimagLite::RuntimeError("unitPerSuper unimplemented\n");
		//return abs(det(alphaN_alphaS));
	}

private:

	VectorRnIndexType vecRnIndex_; // this is now contained instead of inherited
	MatrixRealType alphaNalphaS_;
	VectorRealType tauS_;

};
}
#endif // MAPTOM_H
