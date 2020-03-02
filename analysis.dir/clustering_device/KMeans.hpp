#ifndef K_MEANS_HPP
#define K_MEANS_HPP
#include<CafeInLess/analysis.dir/clustering_device/KMeans_Flag.hpp>
#include<CafeInLess/analysis.dir/clustering_device/ClusteringDevice_interface.hpp>
//#include<KMeans_Type.hpp>
#include<CafeInLess/IO.dir/ErrorMessage.hpp>
#include<CafeInLess/IO.dir/StandardOutput.hpp>
#include<CafeInLess/util.dir/arithmetic>
#include<vector>
#include<random>
#include<string>
#include<limits>



namespace CafeInLess::analysis {


template<KMeans_Flag DIST_FLAG>
class KMeans : ClusteringDeviceInterface<KMeans<DIST_FLAG>> {

public:

	KMeans(const std::size_t& cluster_number_k, const std::size_t& iteration_num, const std::size_t& n_seed):
	cluster_number(cluster_number_k),
	MAX_ITERATION_NUM(iteration_num),
	random_seed(n_seed) {
		if ((cluster_number < 1) || (cluster_number > MAX_CLUSTER_NUMBER)) {
			eout("too much or less cluster number, " + std::to_string(cluster_number));
		}
	}

	KMeans(const std::size_t& cluster_number_k, const std::size_t& iteration_num, const std::size_t& n_seed, const float& cutoff):
	cluster_number(cluster_number_k),
	MAX_ITERATION_NUM(iteration_num),
	random_seed(n_seed),
	CUTOFF_SIMILARITY(cutoff) {
		if ((cluster_number < 1) || (cluster_number > MAX_CLUSTER_NUMBER)) {
			eout("too much or less cluster number, " + std::to_string(cluster_number));
		}
	}

	~KMeans() = default;




	std::vector<std::vector<std::size_t>> run(const std::vector<std::vector<float>>& all_data) {

		if (isDebugMode) {
			sout("k-means clustering starts.");
			sout("Initialization Methods: ++");
			if (DIST_FLAG == CafeInLess::analysis::KMEANS_DIST_L2) {
				sout("Distance Calculation Methods: Euclidean Distance");
			}
		}



		if (isDebugMode) sout("Initial Clusters");
		const std::vector<std::vector<std::size_t>>& initial_cluster_indices_set = init_Cluster(all_data);
	//	for (const std::vector<std::size_t>& indices : initial_cluster_indices_set) {
	//		sout(indices.size());
	//	}

		if (isDebugMode) sout("Initial Centroids");
		std::vector<std::vector<float>> centroids = calc_Centroids(all_data, initial_cluster_indices_set);

		if (isDebugMode) sout("Iteration starts.");
		bool is_to_succeed = false;
		std::size_t iterated_num = 1;
		std::vector<std::vector<std::size_t>> result(cluster_number);

		for (std::size_t i_updated = 0; i_updated < MAX_ITERATION_NUM; ++i_updated) {
			const std::vector<std::vector<std::size_t>>& clustered_indices_set = classify_AllDataIntoClusters(all_data, centroids);
//			for (const std::vector<std::size_t>& indices : clustered_indices_set) {
//				sout(indices.size());
//			}
			const std::vector<std::vector<float>>& updated_centroids = calc_Centroids(all_data, clustered_indices_set);

			const std::vector<float>& centroids_dist = calc_DifferenceBWOldNewCentroids(centroids, updated_centroids);
			if (CafeInLess::arithmetic::is_similar20vec(centroids_dist, CUTOFF_SIMILARITY)) {
				result = clustered_indices_set;
				is_to_succeed = true;
				break;
			}
			else {
				centroids = updated_centroids;
			}
			++iterated_num;
		}

		if (!is_to_succeed) {
			eout("Iteration is finished, but not succeeded. You should cycle more or change cutoff similarity");
		}
		else if (isDebugMode) sout("Iterated number is " + std::to_string(iterated_num));

		return result;
	}





	void set_DebugMode() {isDebugMode = true;}






private:

	const std::size_t cluster_number;
	const std::size_t& MAX_ITERATION_NUM = 100;
	const std::size_t& random_seed = 100000000;
	const float& CUTOFF_SIMILARITY = std::numeric_limits<float>::min();

	const std::size_t& MAX_CLUSTER_NUMBER = 100;

	bool isDebugMode = false;



	CafeInLess::IO::Error_Output eout = CafeInLess::IO::Error_Output();

	CafeInLess::IO::Standard_Output sout = CafeInLess::IO::Standard_Output();


// --------------------------------------------------------------------------------------


	const std::vector<std::vector<std::size_t>> init_Cluster(const std::vector<std::vector<float>>& all_data) {

		const std::vector<std::size_t>& ini_cluster_nuclear_indices = choose_DistanceIndices(all_data);
//		for (const std::size_t& chosen_id : ini_cluster_nuclear_indices) {
//			sout(chosen_id);
//		}

	//	sout("Chosen");
		std::vector<std::vector<float>> cluster_nuclears;
		for (const std::size_t& ini_cluster_nuclear_index : ini_cluster_nuclear_indices) {
			cluster_nuclears.push_back(all_data[ini_cluster_nuclear_index]);
		}
		return classify_AllDataIntoClusters(all_data, cluster_nuclears);
	}




	const std::vector<std::size_t> choose_DistanceIndices(const std::vector<std::vector<float>>& all_data) {

		std::mt19937_64 mt_engine_for_uniform_dist(random_seed);
		std::mt19937_64 mt_engine_for_piecewise_dist(random_seed + 1);
		std::uniform_int_distribution<> uni_int_distribution(0, all_data.size() - 1);
		std::vector<std::size_t> result(cluster_number);

		result[0] = uni_int_distribution(mt_engine_for_uniform_dist);

		for (std::size_t i_cluster = 1; i_cluster < cluster_number; ++i_cluster) {
			double total_dist = 0.0;
			std::vector<double> prob_intervals(all_data.size() + 1, 0);
			std::vector<double> prob_densities(all_data.size());
			std::size_t datum_index = 0;
			for (const std::vector<float>& datum : all_data) {
				double dist = 0.0;
				for (std::size_t chosen_cluster_id = 0; chosen_cluster_id < i_cluster; ++chosen_cluster_id) {
					const std::vector<float> chosen_datum = all_data[result[chosen_cluster_id]];
			//		sout(chosen_datum.size());
			//		sout(datum.size());
			//		sout("Refer");
			//		const float& data_dist_bw_chosen_and_datum = calc_SquareDistance(datum, all_data[result[chosen_cluster_id]]);
					const float& data_dist_bw_chosen_and_datum = calc_SquareDistance(datum, chosen_datum);
					dist += static_cast<double>(data_dist_bw_chosen_and_datum);
					total_dist += static_cast<double>(data_dist_bw_chosen_and_datum);
				}
				prob_intervals[datum_index + 1] = static_cast<double>(datum_index + 1);
				prob_densities[datum_index] = dist;
				++datum_index;
			}

			for (std::size_t prob_dens_idx = 0; prob_dens_idx < prob_densities.size(); ++prob_dens_idx) {
				prob_densities[prob_dens_idx] /= total_dist;
			}

			std::piecewise_constant_distribution<> pw_const_dist(
				prob_intervals.begin(),
				prob_intervals.end(),
				prob_densities.begin()
			);

			const std::size_t& chosen_index = std::floor(pw_const_dist(mt_engine_for_piecewise_dist));
			result[i_cluster] = chosen_index;

		}

		return result;
	}




	const std::vector<std::vector<std::size_t>> classify_AllDataIntoClusters(const std::vector<std::vector<float>>& all_data, const std::vector<std::vector<float>>& cluster_nuclears) {

		std::vector<std::vector<std::size_t>> result(cluster_number);

		for (std::size_t data_index = 0; data_index < all_data.size(); ++data_index) {
			const std::vector<float>& datum = all_data[data_index];
			const std::size_t& closest_cluster_id = calc_ClosestClusterID(datum, cluster_nuclears);
			result[closest_cluster_id].push_back(data_index);
		}
		return result;
	}



	const std::size_t calc_ClosestClusterID(const std::vector<float>& datum, const std::vector<std::vector<float>>& cluster_nuclears) {

		std::size_t result = 0;
		float minimal_distance = std::numeric_limits<float>::max();

		for (std::size_t i_cluster = 0; i_cluster < cluster_number; ++i_cluster) {
		//	sout(datum.size(), cluster_nuclears[i_cluster].size());
			float dist = calc_SquareDistance(datum, cluster_nuclears[i_cluster]);
			if (dist < minimal_distance) {
				minimal_distance = dist;
				result = i_cluster;
			}
		}

		return result;
	}



	const std::vector<std::vector<float>> calc_Centroids(const std::vector<std::vector<float>>& all_data, const std::vector<std::vector<std::size_t>>& clustered_indices_set) {

		std::vector<std::vector<float>> result(cluster_number);

		for (std::size_t i_cluster = 0; i_cluster < cluster_number; ++i_cluster) {
			result[i_cluster] = calc_ClusterCentroid(all_data, clustered_indices_set[i_cluster]);
		}
		return result;
	}



	const std::vector<float> calc_ClusterCentroid(const std::vector<std::vector<float>>& all_data, const std::vector<std::size_t>& clustered_indices) {

		std::vector<std::vector<float>> data_in_cluster;
		for (std::size_t index = 0; index < clustered_indices.size(); ++index) {
			const std::size_t& clustered_index = clustered_indices[index];
			data_in_cluster.push_back(all_data[clustered_index]);
		}
//		sout("data_in_cluster size ", data_in_cluster.size());

		return calc_GeometricCentroid(data_in_cluster);

	}



	const std::vector<float> calc_DifferenceBWOldNewCentroids(const std::vector<std::vector<float>>& old_centroids, const std::vector<std::vector<float>>& new_centroids) {

		std::vector<float> result(cluster_number);

		for (std::size_t i_cluster = 0; i_cluster < cluster_number; ++i_cluster) {
			result[i_cluster] = calc_SquareDistance(old_centroids[i_cluster], new_centroids[i_cluster]);
		}

		return result;
	}


	
	const float calc_SquareDistance(const std::vector<float>& lhs, const std::vector<float>& rhs);



//	const float calc_EuclideanSquareDistance(const std::vector<float>& lhs, const std::vector<float>& rhs) {
//		bool is_same_size = (lhs.size() == rhs.size());
//		if (!is_same_size) {
//			eout("Not same size, ", "lhs " + std::to_string(lhs.size()), "" + std::to_string(rhs.size()));
//		}
//		const std::size_t& iterate_num = lhs.size();
//		float result = 0.0;
//
//		for (std::size_t index = 0; index < iterate_num; ++index) {
//			result += (lhs[index] - rhs[index]) * (lhs[index] - rhs[index]);
//		}
//		return result;
//	}


	const std::vector<float> calc_GeometricCentroid(const std::vector<std::vector<float>>& data_in_cluster) {
		const std::size_t& data_number = data_in_cluster.size();
		const std::size_t& data_size = data_in_cluster.front().size();

		std::vector<float> result(data_size, 0);
		for (std::size_t data_index = 0; data_index < data_number; ++data_index) {
			for (std::size_t index = 0; index < data_size; ++index) {
				result[index] += data_in_cluster[data_index][index];
			}
		}
		for (std::size_t index = 0; index < data_size; ++index) {
			result[index] /= static_cast<float>(data_number);
		}

		return result;
	}

// ------------------------------------------------------------------------------------

};


// -- namespace cafemol::analysis -----------------------------------------------------

template<>
inline const float KMeans<KMeans_Flag::KMEANS_DIST_L2>::calc_SquareDistance(const std::vector<float>& lhs, const std::vector<float>& rhs) {
	return CafeInLess::arithmetic::calc_SquareEuclideanDistance(lhs, rhs);
}



}

#endif /* K_MEANS_HPP */
