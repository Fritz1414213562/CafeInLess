#ifndef CAFEINLESS_S_DENS_BW_HPP
#define CAFEINLESS_S_DENS_BW_HPP
#include<CafeInLess/analysis.dir/clustering_validation/S_Dens_bw_Flag.hpp>
#include<CafeInLess/analysis.dir/clustering_validation/Validation_interface.hpp>
#include<CafeInLess/IO.dir/ErrorMessage.hpp>
#include<CafeInLess/IO.dir/StandardOutput.hpp>
#include<CafeInLess/util.dir/arithmetic>
#include<vector>
#include<string>
#include<cmath>


namespace CafeInLess::analysis {


template<S_Dens_bw_Flag Dist_Flag>
class S_Dens_bw : ClusteringValidationInterface<S_Dens_bw<Dist_Flag>> {

public:
	
	S_Dens_bw() = default;
	~S_Dens_bw() = default;





	const float run(const std::vector<std::vector<float>>& all_data, const std::vector<std::vector<std::size_t>>& cluster_indices_set) {

		sout[BLOCK_SIZE];
		sout("Calculation starts.");

		if (Dist_Flag == CafeInLess::analysis::S_DBW_DIST_L2) {
			sout("Distance Calculation Method: Euclidean distance");
		}

		if (cluster_indices_set.size() <= 1) {
			eout("The cluster number is smaller than 2.");
		}


		sout("Calculating Centroids.");

		const std::vector<std::vector<float>>& centroids = calc_Centroids(all_data, cluster_indices_set);
		sout("Calculating Variances");

		const std::vector<float>& variance_norms = calc_VarianceNorms(all_data, cluster_indices_set, centroids);

		const float& sum_of_squared_error_norm_of_all = calc_SumSquaredErrorOfAllData(all_data);

		sout("Calculating Standard Deviation");

		const float& average_standard_deviation = calc_AverageStandardDeviation(variance_norms);
		sout("SD =", average_standard_deviation);

		sout("Calculating Scatter Extent, Scat(k)");
		const float& scatter_extent = calc_ScatterExtent(variance_norms, sum_of_squared_error_norm_of_all, cluster_indices_set);
		sout("Scat(k) =", scatter_extent);

		sout("Calculating Density Extent, Dens(k)");
		const float& density_extent = calc_DensityExtent(all_data, cluster_indices_set, centroids, average_standard_deviation);

		sout("Dens(k) =", density_extent);

		const float& result = scatter_extent + density_extent;

		sout("Scat(k) + Dens(k) =", result);

		sout("Calculation ends.");

		return result;

	}























private:


	const std::size_t BLOCK_SIZE = 90;

	CafeInLess::IO::Standard_Output sout = CafeInLess::IO::Standard_Output();
	CafeInLess::IO::Error_Output eout = CafeInLess::IO::Error_Output();



// ------------------------------------------------------------------------------------





	const float calc_ScatterExtent(const std::vector<float>& variance_norms, const float& sum_of_squared_error_norm_of_all, const std::vector<std::vector<std::size_t>>& cluster_indices_set) {

		float result = 0.0;
		const std::size_t& cluster_number = variance_norms.size();

		for (std::size_t i_cluster = 0; i_cluster < cluster_number; ++i_cluster) {
			const float& cluster_size = cluster_indices_set[i_cluster].size();

			const float& variance_norm = variance_norms[i_cluster];
			result += variance_norm * cluster_size / sum_of_squared_error_norm_of_all;
		}

		result /= static_cast<float>(cluster_number);

		return result;
	}





	const float calc_DensityExtent(const std::vector<std::vector<float>>& all_data, const std::vector<std::vector<std::size_t>>& cluster_indices_set, const std::vector<std::vector<float>>& centroids, const float& average_standard_deviation) {

		const std::size_t& cluster_number = cluster_indices_set.size();

		float result = 0.0;

		for (std::size_t i_cluster = 0; i_cluster < cluster_number - 1; ++i_cluster) {
			const std::vector<std::size_t>& cluster_indices_of_i = cluster_indices_set[i_cluster];
			const std::vector<float>& centroid_of_i = centroids[i_cluster];

			const std::size_t& num_of_data_close2centroid_i = count_ClosestDatum2Centroid(all_data, cluster_indices_of_i, centroid_of_i, average_standard_deviation);

			for (std::size_t j_cluster = i_cluster + 1; j_cluster < cluster_number; ++j_cluster) {
				const std::vector<std::size_t>& cluster_indices_of_j = cluster_indices_set[j_cluster];
				const std::vector<float>& centroid_of_j = centroids[j_cluster];

				const std::size_t& num_of_data_close2centroid_j = count_ClosestDatum2Centroid(all_data, cluster_indices_of_j, centroid_of_j, average_standard_deviation);


				const std::vector<std::size_t>& cluster_indices_of_i_j = CafeInLess::arithmetic::combine_TwoVectors(cluster_indices_of_i, cluster_indices_of_j);
				const std::vector<float>& mid_of_ij_centroids = CafeInLess::arithmetic::calc_MidPoint(centroid_of_i, centroid_of_j);

				const float& num_of_data_close2mid_of_ij = static_cast<float>(count_ClosestDatum2Centroid(all_data, cluster_indices_of_i_j, mid_of_ij_centroids, average_standard_deviation));


				const float& max_num_between_ij = static_cast<float>(std::max(num_of_data_close2centroid_i, num_of_data_close2centroid_j));

				result += num_of_data_close2mid_of_ij / max_num_between_ij;
			}
		}

		const float& f_cluster_number = static_cast<float>(cluster_number);
		result *= (2 / (f_cluster_number * (f_cluster_number - 1)));

		return result;

	}







// -----------------------------------------------------------------------------------------
// calculating centroids


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





// ---------------------------------------------------------------------------------
// calculating variances, sse, and standard deviation





	const float calc_SumSquaredErrorOfAllData(const std::vector<std::vector<float>>& all_data) {

		const std::vector<float>& centroid_of_all = calc_CentroidOfAllData(all_data);

		std::vector<float> sum_of_squared_errors(centroid_of_all.size(), 0);

		for (const std::vector<float>& data_in_frame : all_data) {

			for (std::size_t i_datum = 0; i_datum < data_in_frame.size(); ++i_datum) {
				const float& datum_difference_from_center = data_in_frame[i_datum] - centroid_of_all[i_datum];
				sum_of_squared_errors[i_datum] += datum_difference_from_center * datum_difference_from_center;
			}
		}

		float result = 0.0;

		for (std::size_t i_datum = 0; i_datum < sum_of_squared_errors.size(); ++i_datum) {
			result += sum_of_squared_errors[i_datum] * sum_of_squared_errors[i_datum];
		}

		result = std::sqrt(result);

		return result;
	}




	const std::vector<float> calc_VarianceNorms(const std::vector<std::vector<float>>& all_data, const std::vector<std::vector<std::size_t>>& cluster_indices_set, const std::vector<std::vector<float>>& centroids) {

		const std::size_t& cluster_number = cluster_indices_set.size();

		std::vector<float> result(cluster_number, 0.0);

		for (std::size_t i_cluster = 0; i_cluster < cluster_number; ++i_cluster) {
			result[i_cluster] = calc_ClusterVarianceNorm(all_data, cluster_indices_set[i_cluster], centroids[i_cluster]);
		}

		return result;
	}




	const float calc_ClusterVarianceNorm(const std::vector<std::vector<float>>& all_data, const std::vector<std::size_t>& cluster_indices, const std::vector<float>& centroid) {

		const float& cluster_size = static_cast<float>(cluster_indices.size());

		std::vector<float> variances(centroid.size(), 0);

		for (const std::size_t& cluster_index : cluster_indices) {
			const std::vector<float>& data_in_frame = all_data[cluster_index];
			for (std::size_t i_datum = 0; i_datum < data_in_frame.size(); ++i_datum) {
				const float& datum_difference_from_center = data_in_frame[i_datum] - centroid[i_datum];
				variances[i_datum] += datum_difference_from_center * datum_difference_from_center;
			}
		}

		float result = 0.0;
		for (std::size_t i_datum = 0; i_datum < centroid.size(); ++i_datum) {
			variances[i_datum] /= cluster_size;
			result += variances[i_datum] * variances[i_datum];
		}

		result = std::sqrt(result);


		return result;
	}






	const float calc_AverageStandardDeviation(const std::vector<float>& variances_of_clusters) {

		float result = 0.0;

		const float& cluster_number = static_cast<float>(variances_of_clusters.size());

		for (const float& variances_of_cluster : variances_of_clusters) {
			result += variances_of_cluster;
		}

		result = std::sqrt(result) / cluster_number;

		return result;
	}



// ---------------------------------------------------------------------------------------
// count data number close to its centroid




	const std::size_t count_ClosestDatum2Centroid(const std::vector<std::vector<float>>& all_data, const std::vector<std::size_t>& cluster_indices, const std::vector<float>& centroid, const float& average_standard_deviation) {

		std::size_t result = 0;

		for (const std::size_t& cluster_index : cluster_indices) {
			const std::vector<float>& data_in_frame = all_data[cluster_index];
			if (is_near_Centroid(data_in_frame, centroid, average_standard_deviation)) {
				++result;
			}
		}

		return result;
	}




	const bool is_near_Centroid(const std::vector<float>& data_in_frame, const std::vector<float>& centroid, const float& average_standard_deviation) {
		const float& distance_from_centroid = calc_Distance(data_in_frame, centroid);
		const bool& result = (distance_from_centroid <= average_standard_deviation);

		return result;
	}




// -----------------------------------------------------------------------------------------


	const float calc_Distance(const std::vector<float>& lhs, const std::vector<float>& rhs);

};


template<>
inline const float S_Dens_bw<CafeInLess::analysis::S_DBW_DIST_L2>::calc_Distance(const std::vector<float>& lhs, const std::vector<float>& rhs) {
	return CafeInLess::arithmetic::calc_EuclideanDistance(lhs, rhs);
}

}


#endif /* CAFEINLESS_S_DENS_BW_HPP */
