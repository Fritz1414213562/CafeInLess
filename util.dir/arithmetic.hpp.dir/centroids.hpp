#ifndef CAFEINLESS_CENTROIDS_HPP
#define CAFEINLESS_CENTROIDS_HPP
#include<vector>
#include<string>
#include<cmath>


namespace CafeInLess::arithmetic {

template<typename realT>
inline const std::vector<realT> calc_GeometricCentroid(const std::vector<std::vector<realT>>& mat_shaped_data) {
	
	const std::size_t mat_col_size = mat_shaped_data[0].size();
	std::vector<realT> result(mat_col_size, 0);

	realT data_number = 0;
	for (const std::vector<realT>& mat_row_vector : mat_shaped_data) {
		for (std::size_t i_col = 0; i_col < mat_col_size; ++i_col) {
			result[i_col] += mat_row_vector[i_col];
		}
		++data_number;
	}
	for (std::size_t i_col = 0; i_col < mat_col_size; ++i_col) {
		result[i_col] /= data_number;
	}

	return result;
}


template<typename realT>
inline const std::vector<realT> calc_PartialGeometricCentroid(const std::vector<std::vector<realT>>& mat_shaped_data, const std::vector<std::size_t>& target_indices) {

	const std::size_t& mat_column_size = mat_shaped_data[target_indices[0]].size();
	std::vector<realT> result(mat_column_size, 0);

	const realT& data_number = static_cast<realT>(target_indices.size());

	for (const std::size_t& target_index : target_indices) {
		for (std::size_t i_col = 0; i_col < mat_column_size; ++i_col) {
			result[i_col] += mat_shaped_data[target_index][i_col];
		}
	}
	for (std::size_t i_col = 0; i_col < mat_column_size; ++i_col) {
		result[i_col] /= data_number;
	}

	return result;
}


template<typename realT>
inline const std::vector<realT> calc_MidPoint(const std::vector<realT>& lhs, const std::vector<realT>& rhs) {

	std::vector<realT> result(lhs.size(), 0);

	for (std::size_t i_result = 0; i_result < lhs.size(); ++i_result) {
		result[i_result] = (lhs[i_result] + rhs[i_result]) / 2.0;
	}

	return result;
}


}


#endif /* CAFEINLESS_CENTROIDS_HPP */
