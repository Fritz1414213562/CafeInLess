#ifndef CAFEINLESS_MATRIX_FILE_WRITER_HPP
#define CAFEINLESS_MATRIX_FILE_WRITER_HPP
#include<coffee-makers/Containers/Containers.hpp>
#include<iostream>
#include<fstream>
#include<vector>
#include<array>
#include<string>
#include<cstring>


namespace CafeInLess::IO {

class MatrixFileWriter {

public:

	MatrixFileWriter() = default;
	MatrixFileWriter(const std::string& outputfile_name) :
		output_name(outputfile_name),
		output_file(output_name, std::ios::out | std::ios::binary) {}

	void open(const std::string& outputfile_name) {
		output_name = outputfile_name;
		open();
	}

	void close() {if (output_file.is_open()) output_file.close();}


	template<typename matrixT>
	void write_MatrixData(const matrixT& mat) {

		using scalarT = typename matrixT::scalar_type;

		std::cout << "Output matrix data to " << output_name << std::endl;

		// dump matrix shape
		write_MatrixShape(mat.rows(), mat.cols());
		std::cout << "Row:Col = " << mat.rows() << ":" << mat.cols() << std::endl;

		// dump matrix data
		std::int64_t mat_data_size = sizeof(scalarT) * mat.size();
//		std::cout << mat_data_size << std::endl;
		write_BinaryOf(mat_data_size);
		for (int idx = 0; idx < mat.size(); ++idx) {
			scalarT copied_cmp = mat[idx];
			write_BinaryOf(copied_cmp);
		}
		write_BinaryOf(mat_data_size);
//		std::cout << mat_data_size << std::endl;
	}


private:

	void open() {
		close();
		output_file.open(output_name, std::ios::out | std::ios::binary);
		if (!output_file.is_open())
			throw std::runtime_error("FileOpenError: The file, " + output_name + " could not be found.");
	}
	

	template<typename T>
	void write_BinaryOf(T& value) {
		output_file.write(reinterpret_cast<char*>(&value), sizeof(T));
	}


	void write_MatrixShape(int row_size, int col_size) {
		std::int64_t block_size = sizeof(int) * 2;
		write_BinaryOf(block_size);
		write_BinaryOf(row_size);
		write_BinaryOf(col_size);
		write_BinaryOf(block_size);
	}


private:

	std::string output_name;
	std::ofstream output_file;

};

}


#endif // CAFEINLESS_MATRIX_FILE_WRITER_HPP
