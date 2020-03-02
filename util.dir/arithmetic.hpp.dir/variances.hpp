#ifndef CAFEINLESS_VARIANCES_HPP
#define CAFEINLESS_VARIANCES_HPP
#include<vector>
#include<CafeInLess/util.dir/arithmetic.hpp.dir/centroids.hpp>


namespace CafeInLess::arithmetic {


template<realT>
inline const std::vector<realT> calc_VarianceVector(const std::vector<std::vector<realT>>& mat_shape_data) {

	const std::vector<realT>& centroid = calc_GeometricCentroid(mat_shape_data);

	std::vector<realT> variances(centroid.size(), 0);
	const realT& data_size = static_cast<realT>(mat_shape_data.size());

	for (const std::vector<realT>& data_in_frame : mat_shape_data) {
		for (std::size_t i_datum = 0; i_datum < data_in_frame.size(); ++i_datum) {
			const realT& datum_difference_from_center = data_in_frame[i_datum] - centroid[i_datum];
			variances[i_datum] += datum_difference_from_center * datum_difference_from_center;
		}
	}

	realT result = 0.0;
	for (std::size_t i_datum = 0; i_datum < variances.size(); ++i_datum) {
		variances[i_datum] /= data_size;
	}

	return variances;
}


}

#endif /* CAFEINLESS_VARIANCES_HPP */
