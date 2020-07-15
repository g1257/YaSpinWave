#ifndef SPIN_MODULUS_H
#define SPIN_MODULUS_H
#include "Vector.h"
#include <fstream>
#include "InputNg.h"
#include "InputCheck.h"

namespace yasw {

template<typename VectorRealType>
class SpinModulus {

public:

	SpinModulus(PsimagLite::String file, SizeType n)
	{
		if (file != "") {
			InputCheck inputCheck;
			PsimagLite::InputNg<InputCheck>::Writeable ioW(file, inputCheck);
			PsimagLite::InputNg<InputCheck>::Readable io(ioW);

			try {
				io.read(moduli_, "SpinModulus");
			} catch (std::exception&) {}
		}

		if (moduli_.size() == 0) moduli_.resize(n, 1.0);
		if (moduli_.size() != n)
			err("SpinModulus::ctor(): FATAL: Expecting spin moduli file " + file +
			    " to have + " + ttos(n) + " spins " +
			    "but found " + ttos(moduli_.size()) + " instead.\n");
	}

	const VectorRealType& operator()() const
	{
		return moduli_;
	}

private:

	VectorRealType moduli_;
};
}
#endif // SPIN_MODULUS_H
