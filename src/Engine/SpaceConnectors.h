#ifndef SPACECONNECTORS_H
#define SPACECONNECTORS_H
#include "Matrix.h"
#include <fstream>
#include "PsimagLite.h"

namespace yasw {

template<typename ComplexOrRealType>
class SpaceConnectors {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixComplexOrRealType;
	typedef typename PsimagLite::Vector<MatrixComplexOrRealType>::Type
	VectorComplexOrRealMatrixType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<VectorRealType>::Type VectorVectorRealType;

	SpaceConnectors(PsimagLite::String file, SizeType pixel, bool verbose)
	    : rows_(0), pixelSize_(pixel), centralCellIndex_(0)
	{
		std::ifstream fin(file.c_str());
		if (!fin || fin.bad() || !fin.good()) {
			PsimagLite::String msg("SpaceConnectors: Problem reading file ");
			msg += file + "\n";
			throw PsimagLite::RuntimeError(msg);
		}

		PsimagLite::String str;
		getline(fin,str);
		std::cerr<<"SpaceConnectors::ctor(): Ignoring first line of ";
		std::cerr<<file<<"\n";
		SizeType n = 0;
		fin>>n;

		PsimagLite::String tmp;
		fin>>tmp;
		SizeType fatRowSize = atoi(tmp.c_str());
		rows_ = fatRowSize/pixelSize_;

		// legacy reading
		size_t indForTmp = tmp.find("*");
		if (indForTmp != PsimagLite::String::npos) {

			std::cerr<<"WARNING: Legacy rows*pixelSize reading for line "<<tmp<<"\n";

			if (pixelSize_ != 1)
				err("Cannot combine legacy rows*pixelSize reading with -P option\n");

			PsimagLite::String sRows = tmp.substr(0, indForTmp);
			SizeType len = tmp.length() - indForTmp - 1;
			PsimagLite::String sPixelSize = tmp.substr(indForTmp + 1, len);
			pixelSize_ = atoi(sPixelSize.c_str());
			rows_ = atoi(sRows.c_str());
			fatRowSize = rows_*pixelSize_;
		}

		std::cerr<<__FILE__<<" : FAT Pixels, rows="<<rows_;
		std::cerr<<" pixelSize="<<pixelSize_<<"\n";

		if (n >= 100)
			err("SpaceConnectors:: too many\n");

		nmatrix_.resize(n);
		data_.resize(n);
		bool centralSeen = false;

		for (SizeType i = 0; i < n; ++i) {
			bool flag = true;
			nmatrix_[i].resize(3);
			for (SizeType j = 0; j < 3; ++j) {
				fin>>nmatrix_[i][j];
				if (nmatrix_[i][j] != 0) flag = false;
			}

			if (flag) {
				if (centralSeen) {
					PsimagLite::String msg("SpaceConnectors: Multiple cells ");
					msg += "are central and this is an error.\n";
					throw PsimagLite::RuntimeError(msg);
				}

				centralSeen = true;
				centralCellIndex_ = i;
				if (verbose)
					std::cerr<<"Found central cell index "<<centralCellIndex_<<"\n";
			}

			data_[i].resize(fatRowSize, fatRowSize);
			for (SizeType j = 0; j < fatRowSize; ++j) {
				for (SizeType k = 0; k < fatRowSize; ++k) {
					fin>>data_[i](j, k);
				}
			}
		}

		fin.close();

		if (!centralSeen) {
			data_.resize(n + 1);
			data_[n].resize(fatRowSize, fatRowSize);
			for (SizeType j = 0; j < fatRowSize; ++j) {
				for (SizeType k = 0; k < fatRowSize; ++k) {
					data_[n](j,k) = 0.0;
				}
			}

			VectorRealType cvector(3);

			nmatrix_.push_back(cvector);
			centralCellIndex_ = n;
		}

		if (verbose) printConnectors(std::cerr);
	}

	ComplexOrRealType operator()(SizeType i, SizeType j, SizeType ind) const
	{
		assert(ind < size());
		return data_[ind](i, j);
	}

	bool isCentralCell(SizeType ind) const
	{
		return (ind == centralCellIndex_);
	}

	SizeType size() const { return data_.size(); }

	SizeType rows() const {return rows_; }

	SizeType pixelSize() const { return pixelSize_; }

	const MatrixComplexOrRealType& getMatrix(SizeType ind) const
	{
		assert(ind < data_.size());
		return data_[ind];
	}

	const VectorRealType& nvector(SizeType ind) const
	{
		assert(ind < nmatrix_.size());
		return nmatrix_[ind];
	}

private:

	void printConnectors(std::ostream& os) const
	{
		os<<"#There are "<<data_.size()<<" connectors.\n";
		for (SizeType i = 0; i < data_.size(); ++i) {
			os<<data_[i];
		}
	}

	VectorVectorRealType nmatrix_;
	VectorComplexOrRealMatrixType data_;
	SizeType rows_;
	SizeType pixelSize_;
	SizeType centralCellIndex_;
}; // SpaceConnectors

}
#endif // SPACECONNECTORS_H
