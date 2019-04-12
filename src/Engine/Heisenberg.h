/*! \file Heisenberg
 *
 *  Heisenberg model
 *
 */
#ifndef YASW_HEISENBERG_H
#define YASW_HEISENBERG_H
#include "HeisenbergFields.h"
#include "../../../spf/v7/Engine/ModelBase.h"
#include "../../../spf/v7/ClassicalFields/Spin.h"
#include "CrsMatrix.h"
#include "EnergyNonCollinearFunction.h"

namespace yasw {
template<typename EngineParamsType, typename GeometryType>
class Heisenberg  : public Spf::ModelBase<Spf::Spin<typename EngineParamsType::RealType>,
        EngineParamsType,
        GeometryType> {

	typedef typename EngineParamsType::IoInType IoInType;
	typedef typename EngineParamsType::RealType RealType;
	typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;
	typedef std::complex<RealType> ComplexType;
	typedef EnergyNonCollinearFunction<ComplexType> EnergyFunctionType;
	typedef typename EnergyFunctionType::VectorRealType VectorRealType;

public:

	enum {OLDFIELDS,NEWFIELDS};

	typedef HeisenbergFields<RealType,GeometryType> DynVarsType;
	typedef typename DynVarsType::SpinOperationsType SpinOperationsType;
	typedef typename DynVarsType::SpinType SpinType;
	typedef PsimagLite::Matrix<RealType> MatrixType;

	Heisenberg(const EngineParamsType& engineParams,
	           const GeometryType& geometry,
	           IoInType& io)
	    : dynVars_(geometry.unitCellSize,engineParams),
	      spinOperations_(geometry,engineParams),
	      energy_(geometry.jfile, geometry.qvector, geometry.verbose)
	{}

	DynVarsType& dynVars() { return dynVars_; }

	SizeType hilbertSize() const { return 0; }

	void setOperation(SpinOperationsType** op,SizeType i)
	{
		assert(i == 0);
		*op = &spinOperations_;
	}

	void createHamiltonian(MatrixType& matrix, SizeType oldOrNewDynVars)
	{}

	void createHsparse(SparseMatrixType& sparseMatrix,SizeType oldOrNewDynVars)
	{}

	RealType deltaDirect(SizeType i,const SpinOperationsType& ops,int n)
	{
		SpinType* spin = 0;
		dynVars_.getField(&spin,0);
		return  calcEnergy(ops.dynVars2()) - calcEnergy(*spin);
	}

	template<typename GreenFunctionType,typename SomePackerType>
	void doMeasurements(GreenFunctionType& greenFunction,
	                    SizeType iter,
	                    SomePackerType& packer)
	{
		SpinType* dynVarsPtr = 0;
		dynVars_.getField(&dynVarsPtr,0);
		const SpinType& dynVars = *dynVarsPtr;

		packer.pack("iter=",iter);

		RealType temp2 = calcEnergy(dynVars);
		packer.pack("TotalEnergy=",temp2);
	}

	void setTpemThings(RealType& a,
	                   RealType& b,
	                   PsimagLite::Vector<SizeType>::Type& support) const
	{}

	template<typename SomeOutputType>
	void finalize(SomeOutputType& fout) {}

	RealType calcEnergy(const SpinType& spin)
	{
		if (data_.size() == 0) data_.resize(2*spin.size);
		for (SizeType i = 0; i < spin.size; ++i) {
			assert(2*i+1 < data_.size());
			data_[2*i] = spin.theta[i];
			data_[2*i+1] = spin.phi[i];
		}

		return energy_(data_);
	}

private:

	DynVarsType dynVars_;
	SpinOperationsType spinOperations_;
	EnergyFunctionType energy_;
	mutable VectorRealType data_;
}; // Heisenberg

template<typename EngineParamsType, typename GeometryType>
std::ostream& operator<<(std::ostream& os,
                         const Heisenberg<EngineParamsType,GeometryType>& model)
{
	os<<"ModelParameters\n";
	return os;
}
} // namespace yasw

/*@}*/
#endif // YASW_HEISENBERG_H

