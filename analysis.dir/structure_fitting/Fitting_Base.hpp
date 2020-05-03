#ifndef CAFEINLESS_FITTING_BASE_HPP
#define CAFEINLESS_FITTING_BASE_HPP
#include<coffee-makers/Containers.hpp>
#include<vector>

namespace CafeInLess::analysis {


template<template<typename> typename Shape, typename realT>
class StructureFittingBase {

public:

	void fit(const std::vector<makers::FixedVector<realT, 3>& points) {
		static_cast<Shape<realT>&>(this)->fit(points);
	}

	void set_DebugMode() {static_cast<Shape<realT>&>(this)->set_DebugMode();}

};

}



#endif /* CAFEINLESS_FITTING_BASE_HPP */
