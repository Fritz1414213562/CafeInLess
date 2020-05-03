#ifndef CAFEINLESS_EIGEN_DECOMPOSITION_RESULT_HPP
#define CAFEINLESS_EIGEN_DECOMPOSITION_RESULT_HPP
#include<coffee-makers/Containers/Containers.hpp>



template<typename scalarT>
struct EigenDecompositionResult {

	using scalar_type = scalarT;



private:

	makers::Matrix<scalarT, makers::Variable, makers::Variable> _eigen_vectors;
	makers::Vector<scalarT, makers::Variable> _eigen_values;

};


#endif /* CAFEINLESS_EIGEN_DECOMPOSITION_RESULT_HPP */
