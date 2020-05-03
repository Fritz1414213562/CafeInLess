#ifndef PCA_HPP
#define PCA_HPP
#include<coffee-makers/Containers/Containers.hpp>
#include<coffee-makers/Solver/Solver.hpp>
#include<vector>
#include<array>
#include<string>


namespace CafeInLess::analysis {


class PCA {

public:

	PCA() = default;
	~PCA() = default;


	template<typename scalarT>
	using MatX = makers::Matrix<scalarT, makers::Variable, makers::Variable>;
	template<typename scalarT>
	using VecX = makers::Matrix<scalarT, makers::Variable, makers::Variable>;

	template<typename scalarT>
	void run(const MatX<scalarT>& all_data) {
		const VecX<scalarT>& data_average = calc_DataAverage(all_data);

		const MatX<scalarT>& average_cross_covariance_matrix = calc_AverageCrossCovarianceMatrix(all_data, data_average);

		makers::JacobiEigenSolver<MatX<scalarT>> eigen_solver(data_average.size());

		eigen_solver.solve(average_cross_covariance_matrix);
		contribution_rates_sum = eigen_solver.eigenvalues().sum();
		contribution_rates = eigen_solver.eigenvalues() / contribution_rates_sum;
		pca_axis_set = eigen_solver.eigenvectors();

		is_run = true;
	}



	const double& contrib_sum() {return contribution_rates_sum;}
	double contrib_sum() {return contribution_rates_sum;}
	double contrib_sum(const std::size_t component_number) {
		double result = 0.0;
		for (std::size_t i_cmp = 0; i_cmp < component_number; ++i_cmp) {
			result += contribution_rates[i_cmp];
		}
		return result;
	}




	const Eigen::VectorXd& contrib_rates() {return contribution_rates;}
	Eigen::VectorXd contrib_rates() {return contribution_rates;}
	Eigen::VectorXd contrib_rates(const std::size_t& pc_end_id) {
		Eigen::VectorXd result(pc_end_id);
		for (std::size_t pc_id = 0; pc_id < pc_end_id; ++pc_id) {
			result[pc_id] = contribution_rates[pc_id];
		}
		return result;
	}
	Eigen::VectorXd contrib_rates(const std::size_t& pc_begin_id, const std::size_t& pc_end_id) {
		Eigen::VectorXd result(pc_end_id - pc_begin_id);
		for (std::size_t pc_id = pc_begin_id; pc_id < pc_end_id; ++pc_id) {
			result[pc_id - pc_begin_id] = contribution_rates[pc_id];
		}
		return result;
	}




	const double& contrib_rate(const std::size_t& pc_id) {return contribution_rates[pc_id];}
	double contrib_rate(const std::size_t& pc_id) {return contribution_rates[pc_id];}



	const Eigen::MatrixXd& pca_axes() {return pca_axis_set;}
	Eigen::MatrixXd pca_axes() {return pca_axis_set;}
	Eigen::MatrixXd pca_axes(const std::size_t& pc_end_id) {
		Eigen::MatrixXd result(pca_axis_set.rows(), pc_end_id);

		for (std::size_t pc_id = 0; pc_id < pc_end_id; ++pc_id) {
			result.col(pc_id) = pca_axis_set.col(pc_id);
		}
		return result;
	}
	Eigen::MatrixXd pca_axes(const std::size_t& pc_begin_id, const std::size_t& pc_end_id) {
		Eigen::MatrixXd result(pca_axis_set.rows(), pc_end_id - pc_begin_id);

		for (std::size_t pc_id = pc_begin_id; pc_id < pc_end_id; ++pc_id) {
			result.col(pc_id - pc_begin_id) = pca_axis_set.col(pc_id);
		}
		return result;
	}





protected:


	const Eigen::VectorXd calc_DataAverage(const Eigen::MatrixXd& all_data) {
		return all_data.colwise().sum().transpose() / all_data.rows();
	}

	void sum_DataTo(Eigen::VectorXd& sum, const Eigen::MatrixXd& one_data_from_one_file) {
		sum += one_data_from_one_file.colwise().sum().transpose();
	}





	const Eigen::MatrixXd calc_AverageCrossCovarianceMatrix(const Eigen::MatrixXd& all_data, const Eigen::VectorXd& data_average) {
		const std::size_t data_size = all_data.cols();
		const std::size_t data_number = all_data.rows();
		Eigen::MatrixXd result(data_size, data_size);

		for (std::size_t i_frame = 0; i_frame < data_number; ++i_frame) {
			result += calc_CrossCovarianceMatrix(all_data.row(i_frame).transpose(), data_average);
		}
		result /= static_cast<double>(data_number);

		return result;
	}

	
	void sum_CrossCovarianceMatrixTo(Eigen::MatrixXd& sum, const Eigen::MatrixXd& one_data_from_one_file, const Eigen::VectorXd& data_average) {
		const std::size_t data_number = one_data_from_one_file.rows();

		for (std::size_t i_frame = 0; i_frame < data_number; ++i_frame) {
			sum += calc_CrossCovarianceMatrix(one_data_from_one_file.row(i_frame).transpose(), data_average);
		}
	}



	const Eigen::MatrixXd calc_CrossCovarianceMatrix(const Eigen::VectorXd& one_frame, const Eigen::VectorXd& data_average) {
		return (one_frame - data_average) * ((one_frame - data_average).transpose());
	}




	CafeInLess::IO::Standard_Output sout = CafeInLess::IO::Standard_Output();
	CafeInLess::IO::Error_Output eout = CafeInLess::IO::Error_Output();

	bool is_run = false;


private:
	// analysis result
	double contribution_rates_sum;
	Eigen::VectorXd contribution_rates;
	Eigen::MatrixXd pca_axis_set;

	// block_size
	const std::size_t BLOCK_SIZE = 90;


};

}


#endif /* PCA_HPP */
