#ifndef CAFEINLESS_VALIDATION_INTERFACE_HPP
#define CAFEINLESS_VALIDATION_INTERFACE_HPP
#include<vector>

namespace CafeInLess::analysis {

template<class T>
class ClusteringValidationInterface {

public:
	const float run(const std::vector<std::vector<float>>& all_data, const std::vector<std::vector<std::size_t>>& cluster_indices_set) {
		return static_cast<T&>(this)->run(all_data, cluster_indices_set);
	}

	void set_DebugMode() {static_cast<T&>(this)->set_DebugMode();}
};

}


#endif /* CAFEINLESS_VALIDATION_INTERFACE_HPP */
