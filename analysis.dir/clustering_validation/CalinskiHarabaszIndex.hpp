#ifndef CAFEINLESS_CALINSKI_HARABASZ_INDEX_HPP
#define CAFEINLESS_CALINSKI_HARABASZ_INDEX_HPP
#include<CafeInLess/analysis.dir/clustering_validation/Validation_interface.hpp>
#include<CafeInLess/analysis.dir/clustering_validation/CalinskiHarabaszIndex_Flag.hpp>
#include<CafeInLess/IO.dir/StandardOutput.hpp>
#include<CafeInLess/IO.dir/ErrorMessage.hpp>
#include<CafeInLess/util.dir/arithmetic>
#include<vector>
#include<string>


namespace CafeInLess::analysis {

template<CH_Index_Flag Dist_Flag>
class CalinskiHarabaszIndex : ClusteringValidationInterface<CalinskiHarabaszIndex<Dist_Flag>> {

public:

	CalinskiHarabaszIndex() = default;
	~CalinskiHarabaszIndex() = default;


	const float run(const std::vector<std::vector<float>>& all_data, const std::vector<std::vector<std::size_t>>& cluster_indices_set) {

		sout[BLOCK_SIZE];
		sout("Calculation starts");

		if (Dist_Flag == CafeInLess::analysis::CH_DIST_L2) {
			sout("Distance Calculation Method: Euclidean distance");
		}
		else if (Dist_Flag == CafeInLess::analysis::CH_DIST_L1) {
			sout("Distance Calculation Method: Manhattan distance");
		}
		else {
			eout("Undefined Flag");
		}

		if (cluster_indices_set.size() <= 1) {
			eout("The cluster_number is smaller than 2.");
		}

		sout("Calculating BGSS");
		const float& bgss = calc_BetweenGroupSquaredSum(all_data, cluster_indices_set);
		sout("BGSS =", bgss);
		sout(".... Done");

		sout("Calculating WGSS");
		const float& wgss = calc_WithinGroupSquaredSum(all_data, cluster_indices_set);
		sout("WGSS =", wgss);
		sout(".... Done");

		sout("Calculating Index");
		const float& cluster_number = static_cast<float>(cluster_indices_set.size());
		const float& data_size = static_cast<float>(all_data.size());

		const float& result = ((data_size - cluster_number) / (cluster_number - 1)) *
								bgss / wgss;
		sout("Calinski-Harabasz index =", result);
		sout(".... Done");
		sout("Calculation ends.");

		return result;

	}



private:


	CafeInLess::IO::Standard_Output sout = CafeInLess::IO::Standard_Output();
	CafeInLess::IO::Error_Output eout = CafeInLess::IO::Error_Output();
	const std::size_t BLOCK_SIZE = 90;


	

	const float calc_BetweenGroupSquaredSum(const std::vector<std::vector<float>>& all_data, const std::vector<std::vector<std::size_t>>& cluster_indices_set) {
		const std::vector<float>& centroid_of_all = calc_CentroidOfAllData(all_data);

		float result = 0.0;

		for (const std::vector<std::size_t>& cluster_indices : cluster_indices_set) {
			const std::vector<float>& centroid = CafeInLess::arithmetic::calc_PartialGeometricCentroid(all_data, cluster_indices);
			const float& cluster_size = static_cast<float>(cluster_indices.size());
			const float& dist_centroids = calc_SquareDistance(centroid, centroid_of_all);
			result = cluster_size * dist_centroids;
		}
		return result;
	}




	const float calc_WithinGroupSquaredSum(const std::vector<std::vector<float>>& all_data, const std::vector<std::vector<std::size_t>>& cluster_indices_set) {

		float result = 0.0;

		for (const std::vector<std::size_t>& cluster_indices : cluster_indices_set) {
			const std::vector<float>& centroid = CafeInLess::arithmetic::calc_PartialGeometricCentroid(all_data, cluster_indices);
			for (const std::size_t& cluster_index : cluster_indices) {
				const float& intracluster_dist = calc_SquareDistance(centroid, all_data[cluster_index]);
				result += intracluster_dist;
			}
		}

		return result;
	}





	const std::vector<std::vector<float>> calc_Centroids(const std::vector<std::vector<float>>& all_data, const std::vector<std::vector<std::size_t>>& cluster_indices_set) {

		std::vector<std::vector<float>> result(cluster_indices_set.size());

		for (std::size_t i_cluster = 0; i_cluster < cluster_indices_set.size(); ++i_cluster) {
			result[i_cluster] = CafeInLess::arithmetic::calc_PartialGeometricCentroid(all_data, cluster_indices_set[i_cluster]);
		}

		return result;
	}



	const std::vector<float> calc_CentroidOfAllData(const std::vector<std::vector<float>>& all_data) {
		return CafeInLess::arithmetic::calc_GeometricCentroid(all_data);
	}



	const float calc_Distance(const std::vector<float>& lhs, const std::vector<float>& rhs);
	const float calc_SquareDistance(const std::vector<float>& lhs, const std::vector<float>& rhs);

};


template<>
inline const float CalinskiHarabaszIndex<CafeInLess::analysis::CH_DIST_L1>::calc_Distance(const std::vector<float>& lhs, const std::vector<float>& rhs) {
	return CafeInLess::arithmetic::calc_ManhattanDistance(lhs, rhs);
}


template<>
inline const float CalinskiHarabaszIndex<CafeInLess::analysis::CH_DIST_L2>::calc_Distance(const std::vector<float>& lhs, const std::vector<float>& rhs) {
	return CafeInLess::arithmetic::calc_EuclideanDistance(lhs, rhs);
}


template<>
inline const float CalinskiHarabaszIndex<CafeInLess::analysis::CH_DIST_L1>::calc_SquareDistance(const std::vector<float>& lhs, const std::vector<float>& rhs) {
	const float& manhattan_dist = CafeInLess::arithmetic::calc_ManhattanDistance(lhs, rhs);
	return manhattan_dist * manhattan_dist;
}

template<>
inline const float CalinskiHarabaszIndex<CafeInLess::analysis::CH_DIST_L2>::calc_SquareDistance(const std::vector<float>& lhs, const std::vector<float>& rhs) {
	return CafeInLess::arithmetic::calc_SquareEuclideanDistance(lhs, rhs);
}




}


#endif /* CAFEINLESS_CALINSKI_HARABASZ_INDEX_HPP */
