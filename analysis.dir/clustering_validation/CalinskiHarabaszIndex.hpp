#ifndef CAFEINLESS_CALINSKI_HARABASZ_INDEX_HPP
#define CAFEINLESS_CALINSKI_HARABASZ_INDEX_HPP
#include<CafeInLess/analysis.dir/clustering_validation/Validation_Base.hpp>
#include<CafeInLess/analysis.dir/clustering_validation/CalinskiHarabaszIndex_Flag.hpp>
#include<coffee-makers/Containers/Containers.hpp>
#include<CafeInLess/util.dir/arithmetic>
#include<vector>
#include<string>


namespace CafeInLess::analysis {

template<CH_Index_Flag Dist_Flag>
class CalinskiHarabaszIndex : ClusteringValidationBase<CalinskiHarabaszIndex<Dist_Flag>> {

	static_assert(
	(Dist_Flag == CH_DIST_L1) || (Dist_Flag == CH_DIST_L2),
	"Unknown distance calculation flag");

	template<typename scalarT>
	using MatX = makers::Matrix<scalarT, makers::Variable, makers::Variable>;
	template<typename scalarT>
	using VecX = makers::Vector<scalarT, makers::Variable>;

public:

	CalinskiHarabaszIndex() = default;
	~CalinskiHarabaszIndex() = default;

	template<typename scalarT>
	float run(const MatX<scalarT>& all_data, const std::vector<std::vector<int>>& cluster_indices_set) const {

		static_assert(std::is_floating_point<scalarT>::value,
		"The scalar type of input is not floating point");

		if ((Dist_Flag == CafeInLess::analysis::CH_DIST_L2) && isDebugMode)
			std::cout << "Distance Calculation Method: Euclidean distance" << std::endl;
		else if ((Dist_Flag == CafeInLess::analysis::CH_DIST_L1) && isDebugMode)
			std::cout << "Distance Calculation Method: Manhattan distance" << std::endl;

		if (cluster_indices_set.size() <= 1)
			throw std::invalid_argument(
			"The cluster number is smaller than 2.");

//		sout("Calculating BGSS");
		if (isDebugMode) std::cout << "Calculating BGSS" << std::endl;
		const scalarT& bgss = calc_BetweenGroupSquaredSum(all_data, cluster_indices_set);
		if (isDebugMode) std::cout << "BGSS = " << bgss << std::endl;
		if (isDebugMode) std::cout << ".... Done" << std::endl;

		if (isDebugMode) std::cout << "Calculating WGSS" << std::endl;
		const scalarT& wgss = calc_WithinGroupSquaredSum(all_data, cluster_indices_set);
		if (isDebugMode) std::cout << "WGSS = " << wgss << std::endl;
		if (isDebugMode) std::cout << ".... Done" << std::endl;

		if (isDebugMode) std::cout << "Calculating Index" << std::endl;
		const scalarT& cluster_number = static_cast<scalarT>(cluster_indices_set.size());
		const scalarT& data_size = static_cast<scalarT>(all_data.rows());

		const scalarT& result = ((data_size - cluster_number) / (cluster_number - 1)) *
								bgss / wgss;
		if (isDebugMode) std::cout << "Calinski-Harabasz index = " << result << std::endl;
		if (isDebugMode) std::cout << ".... Done" << std::endl;
		if (isDebugMode) std::cout << "Calculation ends." << std::endl;

		return result;

	}

	void set_DebugMode() {isDebugMode = true;}


private:


	bool isDebugMode = false;


	template<typename scalarT>
	scalarT calc_BetweenGroupSquaredSum(const MatX<scalarT>& all_data, const std::vector<std::vector<int>>& cluster_indices_set) const {
		const VecX<scalarT>& centroid_of_all = calc_CentroidOfAllData(all_data);

		scalarT result = 0.0;

		for (const std::vector<int>& cluster_indices : cluster_indices_set) {
			const VecX<scalarT>& centroid = calc_ClusterCentroid(all_data, cluster_indices);
			const scalarT& cluster_size = static_cast<scalarT>(cluster_indices.size());
			const scalarT& dist_centroids = calc_SquareDistance(centroid, centroid_of_all);
			result += cluster_size * dist_centroids;
		//	result = cluster_size * dist_centroids;
		}
		return result;
	}




	template<typename scalarT>
	scalarT calc_WithinGroupSquaredSum(const MatX<scalarT>& all_data, const std::vector<std::vector<int>>& cluster_indices_set) const {

		scalarT result = 0.0;

		for (const std::vector<int>& cluster_indices : cluster_indices_set) {
			const VecX<scalarT>& centroid = calc_ClusterCentroid(all_data, cluster_indices);
			for (const int& cluster_index : cluster_indices) {
				const scalarT& intracluster_dist = calc_SquareDistance(centroid, all_data.row(cluster_index));
				result += intracluster_dist;
			}
		}

		return result;
	}




//
//	template<typename scalarT>
//	std::vector<std::vector<float>> calc_Centroids(const std::vector<std::vector<float>>& all_data, const std::vector<std::vector<std::size_t>>& cluster_indices_set) {
//
//		std::vector<std::vector<float>> result(cluster_indices_set.size());
//
//		for (std::size_t i_cluster = 0; i_cluster < cluster_indices_set.size(); ++i_cluster) {
//			result[i_cluster] = CafeInLess::arithmetic::calc_PartialGeometricCentroid(all_data, cluster_indices_set[i_cluster]);
//		}
//
//		return result;
//	}
//


	template<typename scalarT>
	VecX<scalarT> calc_CentroidOfAllData(const MatX<scalarT>& all_data) const {

		VecX<scalarT> result(all_data.cols(), 1);

		for (int index = 0; index < all_data.rows(); ++index)
			result += all_data.row(index);

		result /= static_cast<scalarT>(all_data.rows());

		return result;
	}


	template<typename scalarT>
	VecX<scalarT> calc_ClusterCentroid(const MatX<scalarT>& all_data, const std::vector<int>& clustered_indices) const {

		VecX<scalarT> result(all_data.cols(), 1);

		for (const int& clustered_index : clustered_indices)
			result += all_data.row(clustered_index);

		result /= static_cast<scalarT>(clustered_indices.size());

		return result;
	}


	template<typename scalarT>
	scalarT calc_SquareDistance(const VecX<scalarT>& lhs, const VecX<scalarT>& rhs) const;

};



template<>
template<typename scalarT>
inline scalarT CalinskiHarabaszIndex<CH_DIST_L1>::calc_SquareDistance(const VecX<scalarT>& lhs, const VecX<scalarT>& rhs) const {
	const VecX<scalarT>& relative_diff = lhs - rhs;

	scalarT retval = 0.;
	for (int idx = 0; idx < relative_diff.size(); ++idx) {
		const scalarT& sign = (relative_diff[idx] > 0) - (relative_diff[idx] < 0);
		retval += sign * relative_diff[idx];
	}
	return retval;
}

template<>
template<typename scalarT>
inline scalarT CalinskiHarabaszIndex<CH_DIST_L2>::calc_SquareDistance(const VecX<scalarT>& lhs, const VecX<scalarT>& rhs) const {
	const VecX<scalarT>& relative_diff = lhs - rhs;
	const scalarT& retval = makers::dot(relative_diff, relative_diff);
	return retval;
}




}


#endif /* CAFEINLESS_CALINSKI_HARABASZ_INDEX_HPP */
