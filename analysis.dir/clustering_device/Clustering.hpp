#ifndef CAFEINLESS_CLUSTERING_HPP
#define CAFEINLESS_CLUSTERING_HPP
#include<CafeInLess/analysis.dir/clustering_device/ClusteringDevice_Base.hpp>
#include<coffee-makers/Containers/Containers.hpp>
#include<type_traits>
#include<typeinfo>
#include<vector>
#include<memory>

namespace CafeInLess::analysis {

template<class T, bool is_extended = std::is_base_of<ClusteringDeviceBase<T>, T>::value>
class ClusteringDevice {
	static_assert(is_extended, "ClusteringDevice_Base is not extended.");
};

template<class T>
class ClusteringDevice<T, true> {

	template<typename scalarT>
	using MatX = makers::Matrix<scalarT, makers::Variable, makers::Variable>;

//	using ScalarType = typename T::scalar_type;

public:
	
	template<typename... Args>
	ClusteringDevice(Args... arguments) {
		_object = std::make_unique<T>(arguments...);
	}
	~ClusteringDevice() = default;

	template<typename scalarT>
	std::vector<std::vector<int>> run(const MatX<scalarT>& all_data) {
		return _object->run(all_data);
	}

	void set_DebugMode() {_object->set_DebugMode();}

private:
	std::unique_ptr<T> _object;


};


}

#endif /* CAFEINLESS_CLUSTERING_HPP */
