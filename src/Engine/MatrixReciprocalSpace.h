#ifndef MATRIXRECIPROCALSPACE_H
#define MATRIXRECIPROCALSPACE_H
#include "Matrix.h"
#include "SpaceConnectors.h"

namespace yasw {

template<typename ComplexOrRealType>
class MatrixReciprocalSpace {

public:

	typedef yasw::SpaceConnectors<ComplexOrRealType> SpaceConnectorsType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

	struct ReciprocalArgs {
		PsimagLite::String mfile;
		PsimagLite::String casefile;
		PsimagLite::String mapfile;
		PsimagLite::String jfile;
		PsimagLite::String outputfile;
		PsimagLite::String anglesfile;
		PsimagLite::String modulusfile;
		SizeType nk = 400;
		SizeType width = 1;
		RealType cutoff = 0.0001;

		bool check() const
		{
			PsimagLite::String msg;
			if (mfile == "") msg = "Missing Hamiltonian file\n";
			if (casefile == "") msg = "Missing case file\n";
			if (mapfile == "") msg = "Missing map file\n";
			if (jfile == "") msg = "Missing jfile (couplings) file\n";
			if (outputfile == "") msg = "Missing outputfile file\n";
			if (anglesfile == "") msg = "Missing anglesfile file\n";
			if (modulusfile == "") msg = "Missing modulusfile file\n";
			if (msg != "") std::cerr<<msg;
			return (msg == "");
		}
	};

	typedef std::complex<RealType> ComplexType;
	typedef typename SpaceConnectorsType::MatrixComplexOrRealType MatrixType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	MatrixReciprocalSpace(const ReciprocalArgs& reciprocalArgs, SizeType pixelSize, bool verbose)
	    : sc_(reciprocalArgs.mfile, pixelSize, verbose), verbose_(verbose)
	{}

	VectorRealType dispersion(const VectorRealType& q) const
	{
		MatrixType a = getMatrix(q);
		if (verbose_) {
			std::cerr<<"Hamiltonian matrix in k-space:\n";
			std::cerr<<a;
			std::cerr<<"-------------\n";
		}

		if (!isHermitian(a,true)) {
			PsimagLite::String msg("MatrixSpaceNonCollinear: ");
			msg += "Hamiltonian matrix is not Hermitian\n";
//			throw PsimagLite::RuntimeError(msg);
		}

		multiplyByG(a);

		if (verbose_) {
			std::cerr<<"Dynamic matrix in k-space:\n";
			std::cerr<<a;
			std::cerr<<"-------------\n";
		}

		MatrixType vl(10,10), vr(10,10);
		typename PsimagLite::Vector<ComplexType>::Type eigenvalues(a.n_row());
		PsimagLite::geev('N','N',a,eigenvalues,vl,vr);
		VectorRealType ev;
		for (SizeType i = 0; i < eigenvalues.size(); ++i) {
			if (verbose_) std::cerr<<eigenvalues[i]<<"\n";
			RealType e = getRealPart(eigenvalues[i]);
			if (e < 0) continue;
			if (isInVector(ev,e)) continue;
			ev.push_back(e);
		}

		std::sort(ev.begin(),ev.end());
		return ev;
	}

private:

	void procThisCell(MatrixType& m,SizeType n, const VectorRealType& q) const
	{
		RealType wqn = 0.0;
		for (SizeType i = 0; i < q.size(); ++i) {
			wqn += sc_.nmatrix(n,i) * q[i];
		}

		wqn *= (2.0 * M_PI);

		m += ComplexType(cos(wqn),sin(wqn)) * sc_.getMatrix(n);
	}

	RealType getRealPart(ComplexType c) const
	{
		if (fabs(std::imag(c)) > 1e-6) {
			PsimagLite::String msg("MatrixReciprocalSpace::getRealPart(): ");
			msg += "The number " + ttos(c) + " is not real\n";
			throw PsimagLite::RuntimeError(msg);
		}

		return std::real(c);
	}

	bool isInVector(const VectorRealType& ev, RealType e) const
	{
		for (SizeType j = 0; j < ev.size(); ++j) {
			if (fabs(ev[j]-e) < 1e-4) return true;
		}

		return false;
	}

	void multiplyByG(MatrixType& a) const
	{
		SizeType n = a.n_row();
		SizeType nOver2 = static_cast<SizeType>(n*0.5);

		for (SizeType i = nOver2; i < n; ++i) {
			for (SizeType j = 0; j < n; ++j) {
				a(i,j) *= (-1.0);
			}
		}
	}

	MatrixType getMatrix(const VectorRealType& q) const
	{
		SizeType cells = sc_.size();
		SizeType rows = sc_.rows();
		MatrixType m(rows,rows);
		for (SizeType n = 0; n < cells; ++n) {
			procThisCell(m,n,q);
		}

		return m;
	}

	SpaceConnectorsType sc_;
	bool verbose_;
};

} // namespace yasw

#endif // MATRIXRECIPROCALSPACE_H

