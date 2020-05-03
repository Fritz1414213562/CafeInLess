#ifndef CAFEINLESS_VALIDATION_HPP
#define CAFEINLESS_VALIDATION_HPP
#include<coffee-makers/Containers/Containers.hpp>
#include<CafeInLess/analysis.dir/clustering_validation/Validation_Base.hpp>
#include<type_traits>
#include<typeinfo>
#include<vector>
#include<memory>

namespace CafeInLess::analysis {

template<class T, bool is_extended = std::is_base_of<ClusteringValidationBase<T>, T>::value>
class ClusteringValidation {
	static_assert(is_extended, "ClusteringValidation_Base is not extended");
};


template<class T>
class ClusteringValidation<T, true> {

	template<typename scalarT>
	using MatX = makers::Matrix<scalarT, makers::Variable, makers::Variable>;

public:

	ClusteringValidation() {
		_object = std::make_unique<T>();
	}
	~ClusteringValidation() = default;

	template<typename scalarT>
	scalarT run(const MatX<scalarT>& all_data, const std::vector<std::vector<int>>& cluster_indices_set) const {
		return _object->run(all_data, cluster_indices_set);
	}

	void set_DebugMode() {_object->set_DebugMode();}

private:
	std::unique_ptr<T> _object;

};


}


#endif /* CAFEINLESS_VALIDATION_HPP */
