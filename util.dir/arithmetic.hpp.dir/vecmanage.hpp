#ifndef CAFEINLESS_VECMANAGE_HPP
#define CAFEINLESS_VECMANAGE_HPP
#include<vector>



namespace CafeInLess::arithmetic {

template<typename realT>
bool is_similar20vec(const std::vector<realT>& vec, const realT& cutoff) {
	bool result = true;

	for (const realT& component : vec) {
		result = (result && (component <= cutoff));
	}
	return result;
}



template<typename Type>
inline const std::vector<Type> combine_TwoVectors(const std::vector<Type>& lhs, const std::vector<Type>& rhs) {
	std::vector<Type> result(lhs.size() + rhs.size(), 0);

	for (std::size_t i_result = 0; i_result < lhs.size(); ++i_result) {
		result[i_result] = lhs[i_result];
	}
	for (std::size_t i_result = lhs.size(); i_result < lhs.size() + rhs.size(); ++i_result) {
		result[i_result] = rhs[i_result - lhs.size()];
	}

	return result;
}


}


#endif /* CAFEINLESS_VECMANAGE_HPP */
