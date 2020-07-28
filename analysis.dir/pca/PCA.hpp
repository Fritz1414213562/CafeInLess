#ifndef CAFEINLESS_PCA_HPP
#define CAFEINLESS_PCA_HPP
#include<coffee-makers/Containers/Containers.hpp>
#include<coffee-makers/Solver/Solver.hpp>
#include<coffee-makers/CompileTimeCalculation/matrix_traits.hpp>
#include<coffee-makers/utility/utility.hpp>
#include<vector>
#include<array>
#include<string>
#include<iostream>
#include<fstream>


namespace CafeInLess::analysis {


template<typename MatT>
class PCA {

	static_assert(makers::is_square_matrix<MatT>::value,
	"Matrix shape should be square.");

	using scalarT = typename MatT::scalar_type;
	constexpr static int _dim = MatT::Row_CompileTime;
	using VecT = makers::Vector<scalarT, _dim>;


public:

	PCA() : _dim_runtime(_dim) {
		static_assert(makers::is_fixed<MatT>::value,
		"In the case that template matrix type is variable, the use of a default constructor of PCA class is forbidden.");
	}
	PCA(const int& dim_runtime) : _solver(dim_runtime), _dim_runtime(dim_runtime) {
		if (dim_runtime <= 0)
			throw std::invalid_argument("InvalidArgumentError: Dimension size is not positive.");
	}

	~PCA() = default;

	template<typename matrixT>
	void run(const matrixT& target);

	void turn_onDebugMode() {is_DebugMode = true;}

	VecT contributions() const;
	MatT axes() const;

private: // implement

	template<typename matrixT>
	MatT covariance(const matrixT& mat) const;
	template<typename matrixT>
	VecT mean(const matrixT& mat) const;


private: // member

	makers::JacobiEigenSolver<MatT> _solver;
	bool is_DebugMode = false;
	int _dim_runtime;

};




template<typename MatT>
template<typename matrixT>
void PCA<MatT>::run(const matrixT& target) {

	if (_dim_runtime != target.cols())
		throw std::invalid_argument(
		"The column size of input data is not consistent with the original dimension for PCA");

	if (is_DebugMode)
		std::cout << "Calculating Covariance Matrix" << std::endl;
	
	const MatT& covariance_matrix = covariance(target);
	if (is_DebugMode) {
		std::cout << ".... Done" << std::endl;
		std::cout << "Solving eigenvalue problem" << std::endl;
	}

	_solver.solve(covariance_matrix);

	if (is_DebugMode) {
		std::cout << ".... Done" << std::endl;
	}
}


template<typename MatT>
typename PCA<MatT>::VecT PCA<MatT>::contributions() const {
	const VecT& eigenvalues = _solver.eigen_values();
	const scalarT& sum_of_eigenvalues = makers::sum(eigenvalues);
	return eigenvalues / sum_of_eigenvalues;
}


template<typename MatT>
MatT PCA<MatT>::axes() const {
	return _solver.eigen_vectors();
}



template<typename MatT>
template<typename matrixT>
MatT PCA<MatT>::covariance(const matrixT& mat) const {

	const VecT& center = mean(mat);
	MatT retval(_dim_runtime, _dim_runtime);

	for (int i_row = 0; i_row < mat.rows(); ++i_row) {
		const VecT& deviation = mat.row(i_row) - center;
		retval += deviation * deviation.transpose(); 
		if (is_DebugMode) {
			if (i_row % 100000 == 0) std::cout << "=";
		}
	}

	return retval;
}


template<typename MatT>
template<typename matrixT>
typename PCA<MatT>::VecT PCA<MatT>::mean(const matrixT& mat) const {
	return makers::row_sum(mat) / static_cast<scalarT>(mat.rows());
}

}


#endif /* CAFEINLESS_PCA_HPP */
