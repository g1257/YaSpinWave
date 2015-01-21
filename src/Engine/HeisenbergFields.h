/*! \file HeisenbergFields.h
 *
 *  Wrapper around classical spin, to say that only spin is a MC variable
 *
 */
#ifndef YASW_HEISENBERG_FIELDS_H
#define YASW_HEISENBERG_FIELDS_H
#include "../../spf/v7/ClassicalFields/SpinOperations.h"
#include "loki/Typelist.h"

namespace yasw {
template<typename FieldType, typename GeometryType>
class HeisenbergFields {

public:

	typedef Spf::ClassicalSpinOperations<GeometryType,FieldType> SpinOperationsType;
	typedef typename SpinOperationsType::SpinType SpinType;

	typedef LOKI_TYPELIST_1(SpinOperationsType) OperationsList;

	template<typename SomeParamsType>
	HeisenbergFields(SizeType vol,
	                 const SomeParamsType& params)
	    : spin_(vol,params)
	{
		//spin_.modulus = modulus;
	}

	const PsimagLite::String& name(SizeType) const { return name_; }

	void getField(SpinType** field,SizeType i)
	{
		assert(i == 0);
		*field = &spin_;
	}

	void getField(SpinType const** field,SizeType i) const
	{
		assert(i == 0);
		*field = &spin_;
	}

	template<typename FieldType2,typename GeometryType2>
	friend std::ostream& operator<<(std::ostream& os,
	                                const HeisenbergFields<FieldType2,GeometryType2>& f);

private:
	static const PsimagLite::String name_;
	SpinType spin_;

}; // HeisenbergFields

template<typename FieldType, typename GeometryType>
std::ostream& operator<<(std::ostream& os,
                         const HeisenbergFields<FieldType,GeometryType>& f)
{
	os<<f.spin_;
	return os;
}

template<typename FieldType, typename GeometryType>
const PsimagLite::String HeisenbergFields<FieldType,GeometryType>::name_="spin";

} // namespace yasw

/*@}*/
#endif //YASW_HEISENBERG_FIELDS_H

