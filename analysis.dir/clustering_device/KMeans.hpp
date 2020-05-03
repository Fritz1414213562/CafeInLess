#ifndef CAFEINLESS_K_MEANS_HPP
#define CAFEINLESS_K_MEANS_HPP
#include<CafeInLess/analysis.dir/clustering_device/KMeans_Flag.hpp>
#include<CafeInLess/analysis.dir/clustering_device/ClusteringDevice_Base.hpp>
#include<coffee-makers/Containers/Containers.hpp>
#include<vector>
#include<random>
#include<string>
#include<limits>



namespace CafeInLess::analysis {


template<KMeans_Flag DIST_FLAG, typename scalarT>
class KMeans : ClusteringDeviceBase<KMeans<DIST_FLAG, scalarT>> {

public:

//	using scalar_type = scalarT;

private:

using MatX = makers::Matrix<scalarT, makers::Variable, makers::Variable>;
using VecX = makers::Vector<scalarT, makers::Variable>;

public:

	KMeans(const int& cluster_number_k, const int& iteration_num, const int& n_seed):
	cluster_number(cluster_number_k),
	MAX_ITERATION_NUM(iteration_num),
	random_seed(n_seed) {
		if ((cluster_number < 1) || (cluster_number > MAX_CLUSTER_NUMBER))
			throw std::invalid_argument("too much or less cluster number.");
	}

	KMeans(const int& cluster_number_k, const int& iteration_num, const int& n_seed, const scalarT& cutoff):
	cluster_number(cluster_number_k),
	MAX_ITERATION_NUM(iteration_num),
	random_seed(n_seed),
	CUTOFF_SIMILARITY(cutoff) {
		if ((cluster_number < 1) || (cluster_number > MAX_CLUSTER_NUMBER))
			throw std::invalid_argument("too much or less cluster number.");
	}

	~KMeans() = default;




	std::vector<std::vector<int>> run(const MatX& all_data) const {

		if (isDebugMode) {
			std::cout << "k-means clustering starts." << std::endl;
			std::cout << "Initialization Methods: ++" << std::endl;
			if (DIST_FLAG == CafeInLess::analysis::KMEANS_DIST_L2) {
				std::cout << "Distance Calculation Methods: Euclidean Distance" << std::endl;
			}
		}



		if (isDebugMode) std::cout << "Initial Clusters" << std::endl;
		const std::vector<std::vector<int>>& initial_cluster_indices_set = init_Cluster(all_data);
		if (isDebugMode) {
			for (const std::vector<int>& initial_cluster_indices : initial_cluster_indices_set)
				std::cout << initial_cluster_indices.size() << std::endl;
		}

		if (isDebugMode) std::cout << "Initial Centroids" << std::endl;
		MatX centroids = calc_Centroids(all_data, initial_cluster_indices_set);

		if (isDebugMode) std::cout << "Iteration starts." << std::endl;
		bool is_to_succeed = false;
		int iterated_num = 1;
		std::vector<std::vector<int>> result(cluster_number);
	//	std::cout << MAX_ITERATION_NUM << std::endl;

		for (int i_updated = 0; i_updated < MAX_ITERATION_NUM; ++i_updated) {
			const std::vector<std::vector<int>>& clustered_indices_set = classify(all_data, centroids);
			const MatX& updated_centroids = calc_Centroids(all_data, clustered_indices_set);

			const VecX& centroids_diff = calc_CentroidsDiff(centroids, updated_centroids);
			if (isDebugMode) {
				for (int idx = 0; idx < centroids_diff.size(); ++idx)
					std::cout << centroids_diff[idx] << std::endl;
			}

			if (is_same_centroid(centroids_diff)) {
				result = clustered_indices_set;
				is_to_succeed = true;
				break;
			}
			else {
				centroids = updated_centroids;
			}
			++iterated_num;
		}

		if (isDebugMode) std::cout << "Iterated number is " << iterated_num << std::endl;
		if (!is_to_succeed)
			throw std::runtime_error("The iteration is insufficient.");

		return result;
	}





	void set_DebugMode() {isDebugMode = true;}






private:

	const int cluster_number;
	const int MAX_ITERATION_NUM = 100;
	const int random_seed = 100000000;
	const scalarT& CUTOFF_SIMILARITY = std::numeric_limits<scalarT>::min();

	constexpr static int MAX_CLUSTER_NUMBER = 100;

	bool isDebugMode = false;


	template<KMeans_Flag Dist_Method, KMeans_Flag Dummy_Flag = DUMMY_FLAG>
	struct DistanceCalculator;


// --------------------------------------------------------------------------------------


	std::vector<std::vector<int>> init_Cluster(const MatX& all_data) const {

		const std::vector<int>& ini_cluster_nuclear_indices = choose_DistanceIndices(all_data);
		if (isDebugMode) std::cout << "Initial centroid seeds are chosen." << std::endl;
		if (isDebugMode) {
			for (const int& initial_cluster_index : ini_cluster_nuclear_indices)
				std::cout << initial_cluster_index << std::endl;
		}

		MatX cluster_nuclears(cluster_number, all_data.cols());
		for (int i_cluster = 0; i_cluster < cluster_number; ++i_cluster) {
			const int& cluster_nuclear_idx = ini_cluster_nuclear_indices[i_cluster];
			for (int i_datum = 0; i_datum < all_data.cols(); ++i_datum)
				cluster_nuclears(i_cluster, i_datum) = all_data(cluster_nuclear_idx, i_datum);
		}

		if (isDebugMode) {
			VecX cluster_0 = cluster_nuclears.row(0);

			for (int i_cluster = 1; i_cluster < cluster_number; ++i_cluster) {
				VecX cluster_i = cluster_nuclears.row(i_cluster);
				std::cout << calc_SquareDistance(cluster_0, cluster_i) << std::endl;
			}
		}

		return classify(all_data, cluster_nuclears);
	}




	std::vector<int> choose_DistanceIndices(const MatX& all_data) const {

		std::mt19937_64 mt_engine_for_uniform_dist(random_seed);
		std::mt19937_64 mt_engine_for_piecewise_dist(random_seed + 1);
		std::uniform_int_distribution<> uni_int_distribution(0, all_data.rows() - 1);
		std::vector<int> result(cluster_number);

		if (isDebugMode) std::cout << 0 << std::endl;
		result[0] = uni_int_distribution(mt_engine_for_uniform_dist);

		for (int i_cluster = 1; i_cluster < cluster_number; ++i_cluster) {
			if (isDebugMode) std::cout << i_cluster << std::endl;
			scalarT total_dist = 0.0;
			std::vector<scalarT> prob_intervals(all_data.rows() + 1, 0);
			std::vector<scalarT> prob_densities(all_data.rows());
			for (int datum_index = 0; datum_index < all_data.rows(); ++datum_index) {
			//	if (isDebugMode) std::cout << "TEST" << std::endl;
				const VecX& datum = all_data.row(datum_index);
				scalarT dist = 0.0;
				for (int chosen_cluster_id = 0; chosen_cluster_id < i_cluster; ++chosen_cluster_id) {
					const VecX& chosen_datum = all_data.row(result[chosen_cluster_id]);
					const scalarT& data_dist = calc_SquareDistance(datum, chosen_datum);
					dist += static_cast<scalarT>(data_dist);
					total_dist += static_cast<scalarT>(data_dist);
				}
				prob_intervals[datum_index + 1] = static_cast<scalarT>(datum_index + 1);
				prob_densities[datum_index] = dist;
				++datum_index;
			}

			for (int prob_dens_idx = 0; prob_dens_idx < static_cast<int>(prob_densities.size()); ++prob_dens_idx) {
				prob_densities[prob_dens_idx] /= total_dist;
			}

			std::piecewise_constant_distribution<> pw_const_dist(
				prob_intervals.begin(),
				prob_intervals.end(),
				prob_densities.begin()
			);

			const int& chosen_index = std::floor(pw_const_dist(mt_engine_for_piecewise_dist));
			result[i_cluster] = chosen_index;

		}

		return result;
	}




	std::vector<std::vector<int>> classify(const MatX& all_data, const MatX& cluster_nuclears) const {

		std::vector<std::vector<int>> result(cluster_number);

		for (int data_index = 0; data_index < all_data.rows(); ++data_index) {
			const VecX& datum = all_data.row(data_index);
			const int& closest_cluster_id = search_ClosestID(datum, cluster_nuclears);
			result[closest_cluster_id].push_back(data_index);
		}
		return result;
	}



	int search_ClosestID(const VecX& datum, const MatX& cluster_nuclears) const {

		int result = 0;
		scalarT minimal_distance = std::numeric_limits<scalarT>::max();

		for (int i_cluster = 0; i_cluster < cluster_number; ++i_cluster) {
			const VecX& cluster_nuclear = cluster_nuclears.row(i_cluster);
			scalarT dist = calc_SquareDistance(datum, cluster_nuclear);
			if (dist < minimal_distance) {
				minimal_distance = dist;
				result = i_cluster;
			}
		}

		return result;
	}



	MatX calc_Centroids(const MatX& all_data, const std::vector<std::vector<int>>& clustered_indices_set) const {

		MatX result(cluster_number, all_data.cols());

		for (int i_cluster = 0; i_cluster < cluster_number; ++i_cluster) {
			const VecX& cluster_centroid = calc_ClusterCentroid(all_data, clustered_indices_set[i_cluster]);
			for (int i_datum = 0; i_datum < all_data.cols(); ++i_datum)
				result(i_cluster, i_datum) = cluster_centroid[i_datum];
		}

		return result;
	}



	VecX calc_ClusterCentroid(const MatX& all_data, const std::vector<int>& clustered_indices) const {

		if (clustered_indices.size() <= 0)
			throw std::runtime_error(
			"A cluster in which has no component exists. Due to the data scattering, such a cluster arise. Please try clustering again with another random seed.");

		VecX result(all_data.cols(), 1);

		for (const int& clustered_index : clustered_indices)
			result += all_data.row(clustered_index);
		result /= static_cast<scalarT>(clustered_indices.size());

		return result;
	}



	VecX calc_CentroidsDiff(const MatX& old_centroids, const MatX& new_centroids) const {

		VecX result(cluster_number, 1);

		for (int i_cluster = 0; i_cluster < cluster_number; ++i_cluster) {
			result[i_cluster] = calc_SquareDistance(old_centroids.row(i_cluster), new_centroids.row(i_cluster));
		}

		return result;
	}


	bool is_same_centroid(const VecX& centroids_difference) const {
		bool retval = true;

		for (int i_cluster = 0; i_cluster < cluster_number; ++i_cluster)
			retval = retval && (centroids_difference[i_cluster] <= CUTOFF_SIMILARITY);
		return retval;
	}


	
	scalarT calc_SquareDistance(const VecX& lhs, const VecX& rhs) const;




// ------------------------------------------------------------------------------------

};


// -- namespace CafeInLess::analysis -----------------------------------------------------



template<KMeans_Flag DIST_FLAG, typename scalarT>
template<KMeans_Flag Dist_Method>
struct KMeans<DIST_FLAG, scalarT>::DistanceCalculator<KMeans_Flag::KMEANS_DIST_L2, Dist_Method> {

	using vector_type = makers::Vector<scalarT, makers::Variable>;

	static scalarT compute(const vector_type& lhs, const vector_type& rhs) {
		const vector_type& relative_diff = lhs - rhs;
		const scalarT& retval = makers::dot(relative_diff, relative_diff);
		return retval;
	}
};


template<KMeans_Flag DIST_FLAG, typename scalarT>
scalarT KMeans<DIST_FLAG, scalarT>::calc_SquareDistance(const makers::Vector<scalarT, makers::Variable>& lhs, const makers::Vector<scalarT, makers::Variable>& rhs) const {
	return DistanceCalculator<DIST_FLAG>::compute(lhs, rhs);
}


//template<typename scalarT>
//inline scalarT KMeans_C_DIST_L2<scalarT>::calc_SquareDistance(const makers::Vector<scalarT, makers::Variable>& lhs, const makers::Vector<scalarT, makers::Variable>& rhs) const {
//	const VecX& relative_diff = lhs - rhs;
//	const scalarT& retval = makers::dot(relative_diff, relative_diff);
//	return retval;
//}



}

#endif /* CAFEINLESS_K_MEANS_HPP */
