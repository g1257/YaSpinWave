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
			std::ifstream fin(file.c_str());
			if (!fin || fin.bad() || !fin.good())
				err("SpinModulus::ctor(): Could not open file " + file + "\n");

			SizeType tmp = 0;
			fin>>tmp;
			if (tmp != n)
				err("SpinModulus::ctor(): FATAL: Expecting spin moduli file " + file +
				    " to have + " + ttos(n) + " spins " +
				    "but found " + ttos(tmp) + " instead.\n");

			moduli_.resize(n);
			for (SizeType i = 0; i < tmp; ++i)
				fin>>moduli_[i];

			fin.close();
		}

		if (moduli_.size() == 0) moduli_.resize(n, 1.0);
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
