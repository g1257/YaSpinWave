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
	typedef std::complex<RealType> ComplexType;
	typedef typename SpaceConnectorsType::MatrixComplexOrRealType MatrixType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<VectorRealType>::Type VectorVectorRealType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;

	struct ReciprocalArgs {
		PsimagLite::String mfile;
		PsimagLite::String casefile;
		PsimagLite::String mapfile;
		PsimagLite::String hsfile;
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
			if (hsfile == "") msg = "Missing hsfile (couplings) file\n";
			if (outputfile == "") msg = "Missing outputfile file\n";
			if (anglesfile == "") msg = "Missing anglesfile file\n";
			if (modulusfile == "") msg = "Missing modulusfile file\n";
			if (msg != "") std::cerr<<msg;
			return (msg == "");
		}
	};

	class CaseAux {

	public:

		CaseAux(PsimagLite::String file)
		{
			const SizeType rows = 3;
			const SizeType cols = 3;
			MatrixRealType xb(rows, cols);
			xbc_.resize(rows, cols);

			if (file.find(".aux") == PsimagLite::String::npos)
				file += ".aux";

			{
				std::ifstream fin(file.c_str());

				for (SizeType i = 0; i < rows; ++i)
					for (SizeType j = 0; j < cols; ++j)
						fin>>xbc_(i, j);

				for (SizeType i = 0; i < rows; ++i)
					for (SizeType j = 0; j < cols; ++j)
						fin>>xb(i, j);

				fin.close();
			}

			inverse(xb);
			bbc_ = xb*xbc_;
		}

	private:

		MatrixRealType xbc_;
		MatrixRealType bbc_;
	};

	MatrixReciprocalSpace(const ReciprocalArgs& reciprocalArgs, SizeType pixelSize, bool verbose)
	    : sc_(reciprocalArgs.hsfile, pixelSize, verbose),
	      verbose_(verbose),
	      reciprocalArgs_(reciprocalArgs),
	      caseAux_(reciprocalArgs.casefile)
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

	RealType computeKlength(const CaseAux& caseAux,
	                        const VectorRealType& kbegin,
	                        const VectorRealType& kend) const
	{
		const SizeType npanel = sc_.size();

		//klengthtot (from Tom B.)
		RealType klengthtot(0);
		for (SizeType i = 0; i < npanel; ++i)
			klengthtot += PsimagLite::norm(matMulVec(caseAux.xbc(), kend[i] - kbegin[i]));
		return klengthtot;
	}

	void createMesh(const CaseAux& caseAux,
	                const VectorRealType& kbegin,
	                const VectorRealType& kend)
	{
		const SizeType npanel = sc_.size();
		const RealType klengthtot = computeKlength(caseAux, kbegin, kend);

		//generate kmesh,klenght (taken from Tom B.)
		VectorVectorRealType kmesh;
		VectorRealType klength;
		RealType klengthtmp = 0;
		for (SizeType i = 0; i < npanel; ++i) {
			RealType panellenght = PsimagLite::norm(matMulVec(caseAux.xbc(), kend[i] - kbegin[i]));
			int panelnk = ceil(reciprocalArgs_.nk*panellenght/klengthtot);
			for (int j = 0; j < panelnk; ++j) {
				VectorRealType k = matMulVec(caseAux.bbc(),
				                             kbegin[i]+(1.*j/panelnk)*(kend[i] - kbegin[i]));
				kmesh.push_back(k);
				klength.push_back(klengthtmp+(1.*j/panelnk)*panellenght);
			}

			klengthtmp += panellenght;
		}

		kmesh.push_back(kend[npanel-1]);
		klength.push_back(klengthtot);
	}

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
	const ReciprocalArgs& reciprocalArgs_;
	CaseAux caseAux_;
};

} // namespace yasw

#endif // MATRIXRECIPROCALSPACE_H

