#ifndef CAFEINLESS_CLUSTERING_DEVICE_BASE_HPP
#define CAFEINLESS_CLUSTERING_DEVICE_BASE_HPP
#include<coffee-makers/Containers/Containers.hpp>
#include<vector>

namespace CafeInLess::analysis {

template<class T>
class ClusteringDeviceBase {

template<typename scalarT>
using MatX = makers::Matrix<scalarT, makers::Variable, makers::Variable>;

//using ScalarType = typename T::scalar_type;

public:

	template<typename scalarT>
	std::vector<std::vector<int>> run(const MatX<scalarT>& all_data) {
		return static_cast<T&>(this)->run(all_data);
	}

	void set_DebugMode() {static_cast<T&>(this)->set_DebugMode();}

};

}

#endif /* CAFEINLESS_CLUSTERING_DEVICE_BASE_HPP */
