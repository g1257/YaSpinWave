#ifndef MATRIXRECIPROCALSPACE_H
#define MATRIXRECIPROCALSPACE_H
#include "Matrix.h"
#include "SpaceConnectors.h"
#include "MapTom.h"
#include "Angles.h"
#include "SpinModulus.h"

namespace yasw {

template<typename ComplexOrRealType>
class MatrixReciprocalSpace {

public:

	typedef yasw::SpaceConnectors<ComplexOrRealType> SpaceConnectorsType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef std::complex<RealType> ComplexType;
	typedef typename SpaceConnectorsType::MatrixComplexOrRealType MatrixType;
	typedef typename SpaceConnectorsType::VectorRealType VectorRealType;
	typedef typename SpaceConnectorsType::VectorVectorRealType VectorVectorRealType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename PsimagLite::Vector<MatrixType*>::Type VectorMatrixType;
	typedef MapTom<ComplexOrRealType> MapTomType;
	typedef typename MapTomType::VectorIntType VectorIntType;
	typedef typename MapTomType::RnIndexType RnIndexType;
	typedef Angles<RealType> AnglesType;

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
			xb_.resize(rows, cols);
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
						fin>>xb_(i, j);

				fin.close();
			}

			MatrixRealType bx = xb_;
			inverse(bx);
			bbc_ = bx*xbc_;
		}

		const MatrixRealType& xbc() const { return xbc_; }

		const MatrixRealType& bbc() const { return bbc_; }

		const MatrixRealType& xb() const { return xb_; }

	private:

		MatrixRealType xbc_;
		MatrixRealType bbc_;
		MatrixRealType xb_;
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

		SizeType kMeshSize() const { return kmesh_.size(); }

		const VectorRealType& kMesh(SizeType ik) const
		{
			assert(ik < kmesh_.size());
			return kmesh_[ik];
		}

		const RealType& klength(SizeType ind) const
		{
			assert(ind < klength_.size());
			return klength_[ind];
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
					const RealType factor = static_cast<RealType>(j)/panelnk;
					fillTmpVec(tmpVec, i, factor);
					VectorRealType k = matMulVec(caseAux.bbc(),
					                             tmpVec);
					kmesh_.push_back(k);
					klength_.push_back(klengthtmp + factor*panellenght);
				}

				klengthtmp += panellenght;
			}

			assert(npanel < kend_.size() + 1);
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
	      sc_(reciprocalArgs_.mfile, pixelSize, verbose),
	      caseAux_(reciprocalArgs.casefile),
	      mapTom_(reciprocalArgs_.mapfile),
	      hs_(reciprocalArgs.hsfile, caseAux_, reciprocalArgs_.nk),
	      a_(reciprocalArgs_.anglesfile, verbose),
	      spinModulus_(reciprocalArgs.modulusfile, sc_.rows())
	{}

	void mainLoop(PsimagLite::String& line)
	{
		const SizeType norbital = sc_.rows(); // may need to multiply by pixelSize FIXME TODO

		const SizeType nsite = a_.size();
		MatrixType Vm(nsite, 3);
		MatrixType Vp(nsite, 3);
		computeVminVmax(Vm, Vp);

		// original author: Tom B.
		RealType numImtot(0);

		const SizeType nkmesh = hs_.kMeshSize();
		VectorType E;
		MatrixType Ydag;
		MatrixType mvl;
		for (SizeType ik = 0; ik < nkmesh; ++ik) {
			PsimagLite::String kmeshStr = vectorToString(hs_.kMesh(ik));
			line += kmeshStr + " " + ttos(hs_.klength(ik));
			//std::cerr<<"q="<<kmesh[ik]<<'\n';

			//construct <kn1|H|kn2>
			MatrixType HK(norbital, norbital);
			fillHk(HK, ik);

			//diagonalize <kn1|H|kn2>
			geev('N', 'V', HK, E, mvl, Ydag); // mvl will be ignored apparently

			VectorBoolType deletedIndices(E.size(), false);
			const SizeType erased = cutoffImE(deletedIndices, E);
			numImtot += erased;

			printRe(line, deletedIndices, E);

			MatrixType Xdag(norbital, norbital - erased);
			orthogonalize(Xdag, deletedIndices, Ydag);
			normalize(Xdag, deletedIndices);

			VectorRealType xk = matMulVec(caseAux_.xb(), hs_.kMesh(ik));
			computeSq(Xdag, E, deletedIndices, xk, Vm, Vp);
		}
	}

private:

	void fillHk(MatrixType& HK, SizeType ik) const
	{
		const SizeType nsite = sc_.size();
		const SizeType norbital = HK.rows();
		const SizeType norbitalOver2 = norbital/2;
		assert(norbital == HK.cols());
		assert(norbital == sc_.rows()); // may need to multiply by pixelSize FIXME TODO

		// original author: Tom B.
		VectorIntType vzero = {0, 0, 0};
		for (SizeType ir = 0; ir < nsite; ++ir) {
			VectorRealType aR = matMulVec(mapTom_.alphaNalphaS(), sc_.nvector(ir)); // aR=aA*AR
			for (SizeType io1 = 0; io1 < norbital; ++io1) {
				for (SizeType io0 = 0; io0 < norbital; ++io0) {
					VectorIntType r0 = mapTom_(RnIndexType(vzero, io0)).first();
					VectorIntType r1 = mapTom_(RnIndexType(vzero, io1)).first();

					RealType arg = 0;
					const SizeType nrsize = r0.size();
					assert(nrsize == r1.size());
					assert(nrsize == aR.size());
					for (SizeType ii = 0; ii < nrsize; ++ii)
						arg += 2*M_PI*hs_.kMesh(ik)[ii]*(r0[ii] - r1[ii] - aR[ii]);

					ComplexOrRealType phase(cos(arg), sin(arg));
					const RealType sign = (io0 < norbitalOver2) ? 1 : -1;
					HK(io1, io0) += sign*phase*sc_(io1, io0, ir);
				}
			}
		}
	}

	SizeType cutoffImE(VectorBoolType& deletedIndices, const VectorType& E) const
	{
		// taken from Tom B.
		//rm |ImE|>cutImE
		// don't actually remove, just mark them
		const SizeType norbital = E.size();
		const RealType cutImE = reciprocalArgs_.cutoff;
		SizeType erased = 0;
		for (SizeType ie = 0; ie < norbital; ++ie) {
			if (fabs(PsimagLite::imag(E[ie])) <= cutImE) continue;
			++erased;
			deletedIndices[ie] = true;
		}

		//log numImE
		std::cerr<<"number of ImE="<<erased<<'\n';
		return erased;
	}

	void printRe(PsimagLite::String& line,
	             const VectorBoolType& deletedIndices,
	             const VectorType& E) const
	{
		// from Tom B.
		//print reE
		const SizeType norbital = E.size();
		for (SizeType ie = 0; ie < norbital; ++ie) {
			if (deletedIndices[ie]) continue;
			line += " " + ttos(PsimagLite::real(E[ie]));
		}

		for (SizeType ie = 0; ie < norbital; ++ie) {
			if (!deletedIndices[ie]) continue;
			line += " 123.456789";
		}
	}

	void orthogonalize(MatrixType& Xdag,
	                   const VectorBoolType& deletedIndices,
	                   const MatrixType& Ydag) const
	{
		// taken from Tom B.
		const SizeType norbital = Xdag.rows();
		const SizeType norbitalOver2 = norbital/2;
		const SizeType cols = Xdag.cols();
		VectorRealType X2X(cols);
		SizeType ie1Translated = 0;
		for (SizeType ie1 = 0; ie1 < norbital; ++ie1) {
			if (deletedIndices[ie1]) continue;
			for(SizeType io = 0; io < norbital; ++io)
				Xdag(io, ie1Translated) = Ydag(io, ie1);

			for(SizeType ie2Translated = 0; ie2Translated < ie1; ++ie2Translated) {
				ComplexOrRealType X2Y = 0;
				for (SizeType io = 0; io < norbital; ++io)  {
					const RealType sign = (io < norbitalOver2) ? 1 : -1;
					X2Y += sign*PsimagLite::conj(Xdag(io, ie2Translated))*Ydag(io, ie1);
				}

				const ComplexOrRealType factor = X2Y/X2X[ie2Translated];
				for (SizeType io = 0; io < norbital; ++io)
					Xdag(io, ie1Translated) -= Xdag(io, ie2Translated)*factor;
			}

			for (SizeType io = 0; io < norbital; ++io)  {
				const RealType sign = (io < norbitalOver2) ? 1 : -1;
				X2X[ie1Translated] += sign*PsimagLite::norm(Xdag(io, ie1Translated));
			}

			++ie1Translated;
		}

	}

	void normalize(MatrixType& Xdag, const VectorBoolType& deletedIndices) const
	{
		const SizeType norbital = Xdag.rows();
		const SizeType norbitalOver2 = norbital/2;

		for (SizeType ie = 0; ie < norbital; ++ie)  {
			if (deletedIndices[ie]) continue;
			ComplexOrRealType XNX = 0;
			for (SizeType io = 0; io < norbital; ++io)  {
				const RealType sign = (io < norbitalOver2) ? 1 : -1;
				XNX += sign * PsimagLite::norm(Xdag(io, ie));
			}

			ComplexOrRealType factor = 1.0/sqrt(XNX);
			for (SizeType io = 0; io < norbital; ++io)
				Xdag(io, ie) = Xdag(io, ie)*factor;;
		}
	}

	void computeSq(const MatrixType& Xdag,
	               const VectorType& E,
	               const VectorBoolType& deletedIndices,
	               const VectorRealType& xk,
	               const MatrixType& Vm,
	               const MatrixType& Vp)
	{
		static const ComplexOrRealType oneComplex = 1;
		const SizeType norbital = Xdag.rows();
		const SizeType norbitalOver2 = norbital/2;
		const SizeType cols = Xdag.cols();
		VectorMatrixType Sq(cols);
		VectorType Sqtot(cols);
		const RealType factor = 1.0/(2.0*mapTom_.unitPerSuper());
		const RealType normXk = PsimagLite::norm(xk);
		const RealType factor2 = 1.0/(normXk*normXk);

		SizeType ieTranslated = 0;
		for (SizeType ie = 0; ie < norbital; ++ie)  {
			if (deletedIndices[ie]) continue;
			RealType sign = (PsimagLite::real(E[ie]) < 0) ? -1 : 1;
			VectorType W(3);
			for (SizeType io = 0; io < norbitalOver2; ++io) {
				RealType spin = spinModulus_()[io];
				for (SizeType ii = 0; ii < 3; ++ii)
					W[ii] += sign*(Xdag(io, ieTranslated)*Vm(io, ii) -
					               Xdag(io + norbitalOver2, ieTranslated)*
					               Vp(io, ii))*sqrt(spin);
			}

			Sq[ieTranslated] = new MatrixType(3, 3);
			MatrixType& m = *(Sq[ieTranslated]);
			for (SizeType a = 0; a < 3; ++a) {
				for (SizeType b = 0; b < 3; ++b) {
					m(a, b) = W[a]*PsimagLite::conj(W[b])*factor;
					if (normXk == 0) continue;

					ComplexOrRealType tmp1 = xk[a]*xk[b]*factor2;
					ComplexOrRealType tmp2 = (a == b) ? (oneComplex - tmp1) : tmp1;

					Sqtot[ieTranslated] += tmp2*W[a]*PsimagLite::conj(W[b])*factor;
				}
			}

			++ieTranslated;
		}
	}

	void computeVminVmax(MatrixType& Vm, MatrixType& Vp)
	{
		const SizeType nsite = a_.size();

		for (SizeType io = 0; io < nsite; ++io) {
			MatrixRealType Uinv(3, 3);

			RealType theta = a_.theta(io);
			RealType phi = a_.phi(io);

			Uinv(0, 0) = cos(theta)*cos(phi);
			Uinv(1, 0) = cos(theta)*sin(phi);
			Uinv(2, 0) = -sin(theta);
			Uinv(0, 1) = -sin(phi);
			Uinv(1, 1) = cos(phi);
			Uinv(2,1) = 0;



			Vp(io, 0) = ComplexOrRealType(Uinv(0, 0), Uinv(0, 1));
			Vp(io, 1) = ComplexOrRealType(Uinv(1, 0), Uinv(1, 1));
			Vp(io, 2) = ComplexOrRealType(Uinv(2, 0), Uinv(2, 1));
			Vm(io, 0) = ComplexOrRealType(Uinv(0, 0), -Uinv(0, 1));
			Vm(io, 1) = ComplexOrRealType(Uinv(1, 0), -Uinv(1, 1));
			Vm(io, 2) = ComplexOrRealType(Uinv(2, 0), -Uinv(2, 1));
		}
	}

	template<typename T>
	static PsimagLite::String vectorToString(const std::vector<T>& v)
	{
		const SizeType n = v.size();
		PsimagLite::String str;
		for (SizeType i = 0; i < n; ++i)
			str += ttos(v[i]) + " ";
		return str;
	}

	const ReciprocalArgs& reciprocalArgs_;
	bool verbose_;
	SpaceConnectorsType sc_;
	CaseAux caseAux_;
	MapTomType mapTom_;
	Hs hs_;
	AnglesType a_;
	SpinModulus<VectorRealType> spinModulus_;
};

} // namespace yasw

#endif // MATRIXRECIPROCALSPACE_H

