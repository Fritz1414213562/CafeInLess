#ifndef CAFEINLESS_VALIDATION_HPP
#define CAFEINLESS_VALIDATION_HPP
#include<CafeInLess/analysis.dir/clustering_validation/Validation_interface.hpp>
#include<type_traits>
#include<typeinfo>
#include<vector>
#include<memory>

namespace CafeInLess::analysis {

template<class T, bool is_extended = std::is_base_of<ClusteringValidationInterface<T>, T>::value>
class ClusteringValidation {
	static_assert(is_extended, "ClusteringValidation_Interface is not extended");
};


template<class T>
class ClusteringValidation<T, true> {
public:

	ClusteringValidation() {
		_object = std::make_unique<T>();
	}
	~ClusteringValidation() = default;

	const float run(const std::vector<std::vector<float>>& all_data, const std::vector<std::vector<std::size_t>>& cluster_indices_set) {
		return _object->run(all_data, cluster_indices_set);
	}

	void set_DebugMode() {_object->set_DebugMode();}

private:
	std::unique_ptr<T> _object;

};


}


#endif /* CAFEINLESS_VALIDATION_HPP */
