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
	typedef PsimagLite::Matrix<RealType> MatrixType;

public:

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
		return cos(data_(i,0));
	}

private:

	InputCheckType inputCheck_;
	MatrixType data_;
}; // Angles

}
#endif // ANGLES_H
