#ifndef CAFEINLESS_CLUSTERING_DEVICE_INTERFACE_HPP
#define CAFEINLESS_CLUSTERING_DEVICE_INTERFACE_HPP
#include<vector>

namespace CafeInLess::analysis {

template<class T>
class ClusteringDeviceInterface {

public:
	const std::vector<std::vector<std::size_t>> run(const std::vector<std::vector<float>>& all_data) {
		return static_cast<T&>(this)->run(all_data);
	}

	void set_DebugMode() {static_cast<T&>(this)->set_DebugMode();}

};

}

#endif /* CAFEINLESS_CLUSTERING_DEVICE_INTERFACE_HPP */
