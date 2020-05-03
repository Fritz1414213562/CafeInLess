#ifndef CAFEINLESS_FITTING_HPP
#define CAFEINLESS_FITTING_HPP
#include<CafeInLess/analysis.dir/structure_fitting/Fitting_Base.hpp>
#include<coffee-makers/Containers.hpp>
#include<type_traits>
#include<vector>
#include<memory>



namespace CafeInLess::analysis {


template<
	template<typename> typename Shape, typename realT,
		bool isExtended = std::is_base_of<StructureFittingBase<Shape, realT>, Shape<realT>>::value>
class StructureFitting {
	static_assert(isExtended, "Structure Fitting Base is not extended.");
};


template<template<typename> typename Shape, typename realT>
class StructureFitting<Shape, realT, true> {

public:

	template<typename... Args>
	StructureFitting(Args... arguments) {_object = std::make_unique<Shape<realT>>(arguments...);}
	~StructureFitting() = default;

	void fit(const std::vector<makers::FixedVector<realT, 3>>& points) {_object->fit(points);}

	void set_DebugMode() {_object->set_DebugMode();}

private:
	std::unique_ptr<Shape<realT>> _object;

};

}


#endif /* CAFEINLESS_FITTING_HPP */
