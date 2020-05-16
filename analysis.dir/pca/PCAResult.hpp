#ifndef CAFEINLESS_PCA_RESULT_HPP
#define CAFEINLESS_PCA_RESULT_HPP
#include<coffee-makers/Containers/Containers.hpp>


namespace CafeInLess {

template<typename scalarT>
struct PCAResult {

	using scalar_type = scalarT;

private:

	using MatX = makers::Matrix<scalarT, makers::Variable, makers::Variable>;
	using VecX = makers::Matrix<scalarT, makers::Variable>;




public:

	PCAResult(const MatX& eigen_vectors, const VecX& eigen_values, const int& pc_num) :
		_axes(pc_num, pc_num),
		_contributions(pc_num, 1) {

		if (eigen_vectors.rows() != eigen_vectors.cols())
			throw std::invalid_argument(
			"The matrix of eigen vectors should be square-shape (N x N).");
		else if (eigen_vectors.rows() != _contributions.size())
			throw std::invalid_argument(
			"The size of eigen vector should be the same as the number of eigen values.")
	}


private:

	MatX _axes;
	VecX _contributions;

};
}


#endif /* CAFEINLESS_PCA_RESULT_HPP */
