#ifndef ANGLES_H
#define ANGLES_H
#include "String.h"
#include "Matrix.h"
#include "InputCheck.h"
#include "InputNg.h"

namespace yasw {

template<typename RealType>
class Angles {

	typedef yasw::InputCheck InputCheckType;
	typedef PsimagLite::InputNg<InputCheckType> InputType;

public:

	typedef PsimagLite::Matrix<RealType> MatrixType;

	Angles(PsimagLite::String filename,bool verbose = true)
	{
		InputType::Writeable ioWriteable(filename,inputCheck_);
		InputType::Readable io(ioWriteable);
		io.readMatrix(data_,"Angles");
		if (verbose) {
			std::cerr<<"#Angles read from "<<filename<<"\n";
			std::cerr<<data_;
		}
	}

	RealType operator()(SizeType i, SizeType j) const
	{
		return cos(data_(i,0))*cos(data_(j,0));
	}

	const RealType& theta(SizeType i) const
	{
		return data_(i,0);
	}

	const RealType& phi(SizeType i) const
	{
		return data_(i,1);
	}

private:

	InputCheckType inputCheck_;
	MatrixType data_;
}; // Angles

}
#endif // ANGLES_H
