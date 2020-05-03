#ifndef CAFEINLESS_VALIDATION_BASE_HPP
#define CAFEINLESS_VALIDATION_BASE_HPP
#include<coffee-makers/Containers/Containers.hpp>
#include<vector>

namespace CafeInLess::analysis {

template<class T>
class ClusteringValidationBase {

	template<typename scalarT>
	using MatX = makers::Matrix<scalarT, makers::Variable, makers::Variable>;

public:
	template<typename scalarT>
	scalarT run(const MatX<scalarT>& all_data, const std::vector<std::vector<int>>& cluster_indices_set) const {
		return static_cast<T&>(this)->run(all_data, cluster_indices_set);
	}

	void set_DebugMode() {static_cast<T&>(this)->set_DebugMode();}
};

}


#endif /* CAFEINLESS_VALIDATION_BASE_HPP */
