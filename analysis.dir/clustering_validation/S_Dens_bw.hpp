#ifndef CAFEINLESS_S_DENS_BW_HPP
#define CAFEINLESS_S_DENS_BW_HPP
#include<CafeInLess/analysis.dir/clustering_validation/S_Dens_bw_Flag.hpp>
#include<CafeInLess/analysis.dir/clustering_validation/Validation_Base.hpp>
#include<CafeInLess/util.dir/arithmetic>
#include<coffee-makers/Containers/Containers.hpp>
#include<coffee-makers/utility/utility.hpp>
#include<vector>
#include<string>
#include<cmath>


namespace CafeInLess::analysis {


template<S_Dens_bw_Flag Dist_Flag>
class S_Dens_bw : ClusteringValidationBase<S_Dens_bw<Dist_Flag>> {

	static_assert(Dist_Flag == S_DBW_DIST_L2,
	"Unknown distance calculation method");

	template<typename scalarT>
	using MatX = makers::Matrix<scalarT, makers::Variable, makers::Variable>;
	template<typename scalarT>
	using VecX = makers::Vector<scalarT, makers::Variable>;


public:
	
	S_Dens_bw() = default;
	~S_Dens_bw() = default;


	template<typename scalarT>
	scalarT run(const MatX<scalarT>& all_data, const std::vector<std::vector<int>>& cluster_indices_set) const {

		if ((Dist_Flag == S_DBW_DIST_L2) || (isDebugMode))
			std::cout << "Distance Calculation Method: Euclidean distance" << std::endl;

		if (cluster_indices_set.size() <= 1)
			throw std::invalid_argument("The cluster number is smaller than 2.");


		if (isDebugMode) std::cout << "Calculating Centroids." << std::endl;

		const MatX<scalarT>& centroids = calc_Centroids(all_data, cluster_indices_set);
		if (isDebugMode) std::cout << "Calculating Variances" << std::endl;

		const VecX<scalarT>& variance_norms = calc_VarianceNorms(all_data, cluster_indices_set, centroids);

		const scalarT& sum_of_squared_error_norm_of_all = calc_SumSquaredErrorOfAllData(all_data);

		if (isDebugMode) std::cout << "Calculating Standard Deviation" << std::endl;

		const scalarT& average_standard_deviation = calc_AverageStandardDeviation(variance_norms);
		if (isDebugMode) std::cout << "SD = " << average_standard_deviation << std::endl;

		if (isDebugMode) std::cout << "Calculating Scatter Extent, Scat(k)" << std::endl;
		const scalarT& scatter_extent = calc_ScatterExtent(variance_norms, sum_of_squared_error_norm_of_all, cluster_indices_set);
		if (isDebugMode) std::cout << "Scat(k) = " << scatter_extent << std::endl;

		if (isDebugMode) std::cout << "Calculating Density Extent, Dens(k)" << std::endl;
		const scalarT& density_extent = calc_DensityExtent(all_data, cluster_indices_set, centroids, average_standard_deviation);

		if (isDebugMode) std::cout << "Dens(k) = " << density_extent << std::endl;

		const scalarT& result = scatter_extent + density_extent;

		if (isDebugMode) std::cout << "Scat(k) + Dens(k) = " << result << std::endl;

		if (isDebugMode) std::cout << "Calculation ends." << std::endl;

		return result;

	}




	void set_DebugMode() {isDebugMode = true;}


















private:



	bool isDebugMode = false;



// ------------------------------------------------------------------------------------





	template<typename scalarT>
	scalarT calc_ScatterExtent(const VecX<scalarT>& variance_norms, const scalarT& sum_of_squared_error_norm_of_all, const std::vector<std::vector<int>>& cluster_indices_set) const {

		scalarT result = 0.0;
		const int& cluster_number = variance_norms.size();

		for (int i_cluster = 0; i_cluster < cluster_number; ++i_cluster) {
			const scalarT& cluster_size = static_cast<scalarT>(cluster_indices_set[i_cluster].size());

			const scalarT& variance_norm = variance_norms[i_cluster];
			result += variance_norm * cluster_size / sum_of_squared_error_norm_of_all;
		}

		result /= static_cast<scalarT>(cluster_number);

		return result;
	}





	template<typename scalarT>
	scalarT calc_DensityExtent(const MatX<scalarT>& all_data, const std::vector<std::vector<int>>& cluster_indices_set, const MatX<scalarT>& centroids, const scalarT& average_standard_deviation) const {

		const int& cluster_number = cluster_indices_set.size();

		scalarT result = 0.0;

		for (int i_cluster = 0; i_cluster < cluster_number - 1; ++i_cluster) {
			const std::vector<int>& cluster_indices_of_i = cluster_indices_set[i_cluster];
			const VecX<scalarT>& centroid_of_i = centroids.row(i_cluster);

			const int& num_of_data_close2centroid_i = count_ClosestDatum2Centroid(all_data, cluster_indices_of_i, centroid_of_i, average_standard_deviation);

			for (int j_cluster = i_cluster + 1; j_cluster < cluster_number; ++j_cluster) {
				const std::vector<int>& cluster_indices_of_j = cluster_indices_set[j_cluster];
				const VecX<scalarT>& centroid_of_j = centroids.row(j_cluster);

				const int& num_of_data_close2centroid_j = count_ClosestDatum2Centroid(all_data, cluster_indices_of_j, centroid_of_j, average_standard_deviation);


				const std::vector<int>& cluster_indices_of_i_j = CafeInLess::arithmetic::combine_TwoVectors(cluster_indices_of_i, cluster_indices_of_j);
				const VecX<scalarT>& mid_of_ij_centroids = (centroid_of_i + centroid_of_j) / static_cast<scalarT>(2.);

				const scalarT& num_of_data_close2mid_of_ij = static_cast<scalarT>(count_ClosestDatum2Centroid(all_data, cluster_indices_of_i_j, mid_of_ij_centroids, average_standard_deviation));

				const scalarT& max_num_between_ij = static_cast<scalarT>(std::max(num_of_data_close2centroid_i, num_of_data_close2centroid_j));

				result += num_of_data_close2mid_of_ij / max_num_between_ij;
			}
		}

		const scalarT& f_cluster_number = static_cast<scalarT>(cluster_number);
		result *= (2 / (f_cluster_number * (f_cluster_number - 1)));

		return result;

	}







// -----------------------------------------------------------------------------------------
// calculating centroids


	template<typename scalarT>
	MatX<scalarT> calc_Centroids(const MatX<scalarT>& all_data, const std::vector<std::vector<int>>& cluster_indices_set) const {

		MatX<scalarT> result(static_cast<int>(cluster_indices_set.size()), all_data.cols());

		for (int i_cluster = 0; i_cluster < static_cast<int>(cluster_indices_set.size()); ++i_cluster) {
			const VecX<scalarT>& cluster_centroid = calc_ClusterCentroid(all_data, cluster_indices_set[i_cluster]);
			for (int i_datum = 0; i_datum < all_data.cols(); ++i_datum)
				result(i_cluster, i_datum) = cluster_centroid[i_datum];
		}

		return result;
	}




	template<typename scalarT>
	VecX<scalarT> calc_CentroidOfAllData(const MatX<scalarT>& all_data) const {

		VecX<scalarT> result(all_data.cols(), 1);

		for (int index = 0; index < all_data.rows(); ++index)
			result += all_data.row(index);

		result /= static_cast<scalarT>(all_data.rows());

		return result;
	}



	template<typename scalarT>
	VecX<scalarT> calc_ClusterCentroid(const MatX<scalarT>& all_data, const std::vector<int>& cluster_indices) const {

		const scalarT& f_cluster_size = static_cast<scalarT>(cluster_indices.size());
		VecX<scalarT> retval(all_data.cols(), 1);

		for (const int& cluster_index : cluster_indices)
			retval += all_data.row(cluster_index);
		retval /= f_cluster_size;

		return retval;
	}



// ---------------------------------------------------------------------------------
// calculating variances, sse, and standard deviation




	template<typename scalarT>
	scalarT calc_SumSquaredErrorOfAllData(const MatX<scalarT>& all_data) const {

		const VecX<scalarT>& centroid_of_all = calc_CentroidOfAllData(all_data);

		VecX<scalarT> sum_of_squared_errors(centroid_of_all.size(), 1);

		for (int index = 0; index < all_data.rows(); ++index) {
			const VecX<scalarT>& data_in_frame = all_data.row(index);

			for (int i_datum = 0; i_datum < data_in_frame.size(); ++i_datum) {
				const scalarT& datum_difference_from_center = data_in_frame[i_datum] - centroid_of_all[i_datum];
				sum_of_squared_errors[i_datum] += datum_difference_from_center * datum_difference_from_center;
			}
		}

		scalarT result = 0.0;

		for (int i_datum = 0; i_datum < sum_of_squared_errors.size(); ++i_datum) {
			result += sum_of_squared_errors[i_datum] * sum_of_squared_errors[i_datum];
		}

		result = std::sqrt(result);

		return result;
	}




	template<typename scalarT>
	VecX<scalarT> calc_VarianceNorms(const MatX<scalarT>& all_data, const std::vector<std::vector<int>>& cluster_indices_set, const MatX<scalarT>& centroids) const {

		const int& cluster_number = cluster_indices_set.size();

		VecX<scalarT> result(cluster_number, 1);

		for (int i_cluster = 0; i_cluster < cluster_number; ++i_cluster) {
			result[i_cluster] = calc_ClusterVarianceNorm(all_data, cluster_indices_set[i_cluster], centroids.row(i_cluster));
		}

		return result;
	}




	template<typename scalarT>
	scalarT calc_ClusterVarianceNorm(const MatX<scalarT>& all_data, const std::vector<int>& cluster_indices, const VecX<scalarT>& centroid) const {

		const scalarT& cluster_size = static_cast<scalarT>(cluster_indices.size());

		VecX<scalarT> variances(centroid.size(), 1);

		for (const int& cluster_index : cluster_indices) {
			const VecX<scalarT>& data_in_frame = all_data.row(cluster_index);
			for (int i_datum = 0; i_datum < data_in_frame.size(); ++i_datum) {
				const scalarT& datum_difference_from_center = data_in_frame[i_datum] - centroid[i_datum];
				variances[i_datum] += datum_difference_from_center * datum_difference_from_center;
			}
		}

		scalarT result = 0.0;
		for (int i_datum = 0; i_datum < centroid.size(); ++i_datum) {
			variances[i_datum] /= cluster_size;
			result += variances[i_datum] * variances[i_datum];
		}

		result = std::sqrt(result);


		return result;
	}






	template<typename scalarT>
	scalarT calc_AverageStandardDeviation(const VecX<scalarT>& variances_of_clusters) const {

		scalarT result = 0.0;

		const scalarT& cluster_number = static_cast<scalarT>(variances_of_clusters.size());

		for (int i_cluster = 0; i_cluster < variances_of_clusters.size(); ++i_cluster) {
			result += variances_of_clusters[i_cluster];
		}

		result = std::sqrt(result) / cluster_number;

		return result;
	}



// ---------------------------------------------------------------------------------------
// count data number close to its centroid




	template<typename scalarT>
	int count_ClosestDatum2Centroid(const MatX<scalarT>& all_data, const std::vector<int>& cluster_indices, const VecX<scalarT>& centroid, const scalarT& average_standard_deviation) const {

		int result = 0;

		for (const int& cluster_index : cluster_indices) {
			const VecX<scalarT>& data_in_frame = all_data.row(cluster_index);
			if (is_near_Centroid(data_in_frame, centroid, average_standard_deviation)) {
				++result;
			}
		}

		return result;
	}




	template<typename scalarT>
	bool is_near_Centroid(const VecX<scalarT>& data_in_frame, const VecX<scalarT>& centroid, const scalarT& average_standard_deviation) const {
		const scalarT& distance_from_centroid = calc_Distance(data_in_frame, centroid);
		const bool& result = (distance_from_centroid <= average_standard_deviation);

		return result;
	}




// -----------------------------------------------------------------------------------------


	template<typename scalarT>
	scalarT calc_Distance(const VecX<scalarT>& lhs, const VecX<scalarT>& rhs) const;

};


template<>
template<typename scalarT>
inline scalarT S_Dens_bw<CafeInLess::analysis::S_DBW_DIST_L2>::calc_Distance(const VecX<scalarT>& lhs, const VecX<scalarT>& rhs) const {
	const VecX<scalarT>& relative_diff = lhs - rhs;
	const scalarT retval = makers::distance(relative_diff);
	return retval;
}

}


#endif /* CAFEINLESS_S_DENS_BW_HPP */
