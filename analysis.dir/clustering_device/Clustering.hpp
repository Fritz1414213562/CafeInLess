#ifndef CAFEINLESS_CLUSTERING_HPP
#define CAFEINLESS_CLUSTERING_HPP
#include<CafeInLess/analysis.dir/clustering_device/ClusteringDevice_interface.hpp>
#include<type_traits>
#include<typeinfo>
#include<vector>
#include<memory>

namespace CafeInLess::analysis {

template<class T, bool is_extended = std::is_base_of<ClusteringDeviceInterface<T>, T>::value>
class ClusteringDevice {
	static_assert(is_extended, "ClusteringDevice_Interface  is not extended.");
};

template<class T>
class ClusteringDevice<T, true> {

public:
	
	template<typename... Args>
	ClusteringDevice(Args... arguments) {
		_object = std::make_unique<T>(arguments...);
	}
	~ClusteringDevice() = default;

	const std::vector<std::vector<std::size_t>> run(const std::vector<std::vector<float>>& all_data) {
		return _object->run(all_data);
	}

	void set_DebugMode() {_object->set_DebugMode();}

private:
	std::unique_ptr<T> _object;


};


}

#endif /* CAFEINLESS_CLUSTERING_HPP */
