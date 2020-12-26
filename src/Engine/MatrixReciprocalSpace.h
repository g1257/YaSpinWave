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

	static VectorRealType matMulVec(const MatrixRealType& m, const VectorRealType& v)
	{
		const SizeType rows = m.rows();
		const SizeType cols = m.cols();
		assert(rows == v.size());
		VectorRealType w(rows);
		for (SizeType i = 0; i < rows; ++i) {
			for (SizeType j = 0; j < cols; ++j) {
				w[i] += m(i, j)*v[j];
			}
		}

		return w;
	}
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

		const MatrixRealType& xbc() const { return xbc_; }

		const MatrixRealType& bbc() const { return bbc_; }

	private:

		MatrixRealType xbc_;
		MatrixRealType bbc_;
	};

	class Hs {

	public:

		Hs(PsimagLite::String file, const CaseAux& caseAux, SizeType nk)
		{
			const SizeType cols = 3;
			int npanel = 0;

			{
				std::ifstream fin(file.c_str());
				fin>>npanel;
				assert(npanel > 0);
				kbegin_.resize(npanel);
				kend_.resize(npanel);
				kdiff_.resize(npanel);
				for (int i = 0; i < npanel; ++i) {
					kbegin_[i].resize(cols);
					for (SizeType j = 0; j < cols; ++j)
						fin>>kbegin_[i][j];

					kend_[i].resize(cols);
					kdiff_[i].resize(cols);
					for (SizeType j = 0; j < cols; ++j) {
						fin>>kend_[i][j];
						kdiff_[i][j] = kend_[i][j] - kbegin_[i][j];
					}
				}
			}

			createMesh(caseAux, nk);
		}

	private:

		RealType computeKlength(const CaseAux& caseAux) const
		{
			const SizeType npanel = kbegin_.size();

			//klengthtot (from Tom B.)
			RealType klengthtot(0);
			for (SizeType i = 0; i < npanel; ++i)
				klengthtot += PsimagLite::norm(matMulVec(caseAux.xbc(), kdiff(i)));
			return klengthtot;
		}

		const VectorRealType& kdiff(SizeType ind) const
		{
			assert(ind < kdiff_.size());
			return kdiff_[ind];
		}

		void createMesh(const CaseAux& caseAux, SizeType nk)
		{
			const SizeType npanel = kbegin_.size();
			const RealType klengthtot = computeKlength(caseAux);

			//generate kmesh,klenght (taken from Tom B.)
			RealType klengthtmp = 0;
			const SizeType cols = kbegin_.size();
			assert(cols == 3);
			VectorRealType tmpVec(cols);
			for (SizeType i = 0; i < npanel; ++i) {
				RealType panellenght = PsimagLite::norm(matMulVec(caseAux.xbc(), kdiff(i)));
				int panelnk = ceil(nk*panellenght/klengthtot);
				for (int j = 0; j < panelnk; ++j) {
					const RealType factor = (1.*j/panelnk);
					fillTmpVec(tmpVec, i, factor);
					VectorRealType k = matMulVec(caseAux.bbc(),
					                             tmpVec);
					kmesh_.push_back(k);
					klength_.push_back(klengthtmp+(1.*j/panelnk)*panellenght);
				}

				klengthtmp += panellenght;
			}

			kmesh_.push_back(kend_[npanel - 1]);
			klength_.push_back(klengthtot);
		}

		void fillTmpVec(VectorRealType& tmpVec, SizeType i, RealType factor) const
		{
			assert(i < kbegin_.size());
			const SizeType cols = kbegin_[i].size();
			assert(cols == 3);
			if (tmpVec.size() != cols) tmpVec.resize(cols);
			for (SizeType jj = 0; jj < cols; ++jj)
				tmpVec[jj] = kbegin_[i][jj] + factor*kdiff(i)[jj];
		}

		VectorVectorRealType kbegin_;
		VectorVectorRealType kend_;
		VectorVectorRealType kdiff_;
		VectorVectorRealType kmesh_;
		VectorRealType klength_;
	};

	MatrixReciprocalSpace(const ReciprocalArgs& reciprocalArgs, SizeType pixelSize, bool verbose)
	    : reciprocalArgs_(reciprocalArgs),
	      verbose_(verbose),
	      caseAux_(reciprocalArgs.casefile),
	      hs_(reciprocalArgs.hsfile, caseAux_, reciprocalArgs_.nk)
	{}

private:

	const ReciprocalArgs& reciprocalArgs_;
	bool verbose_;
	CaseAux caseAux_;
	Hs hs_;
};

} // namespace yasw

#endif // MATRIXRECIPROCALSPACE_H

