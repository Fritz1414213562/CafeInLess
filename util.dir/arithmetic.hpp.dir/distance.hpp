#ifndef CAFEINLESS_DISTANCE_HPP
#define CAFEINLESS_DISTANCE_HPP
#include<vector>
#include<cmath>



namespace CafeInLess::arithmetic {


template<typename realT>
inline const realT calc_SquareEuclideanDistance(const std::vector<realT>& lhs, const std::vector<realT>& rhs) {

	realT result = 0;
	for (std::size_t i_lhs = 0; i_lhs < lhs.size(); ++i_lhs) {
		result += (lhs[i_lhs] - rhs[i_lhs]) * (lhs[i_lhs] - rhs[i_lhs]);
	}
	return result;
}



template<typename realT>
inline const realT calc_EuclideanDistance(const std::vector<realT>& lhs, const std::vector<realT>& rhs) {

	realT result = 0;
	for (std::size_t i_lhs = 0; i_lhs < lhs.size(); ++i_lhs) {
		result += (lhs[i_lhs] - rhs[i_lhs]) * (lhs[i_lhs] - rhs[i_lhs]);
	}
	result = std::sqrt(result);

	return result;
}


template<typename realT>
inline const realT calc_ManhattanDistance(const std::vector<realT>& lhs, const std::vector<realT>& rhs) {
	realT result = 0;
	for (std::size_t i_lhs = 0; i_lhs < lhs.size(); ++i_lhs) {
		const realT& difference = lhs[i_lhs] - rhs[i_lhs];
		const realT& diff_sign = static_cast<realT>((difference > 0) - (difference < 0));
		result += diff_sign * difference;
	}
	return result;
}

}



#endif /* CAFEINLESS_DISTANCE_HPP */
