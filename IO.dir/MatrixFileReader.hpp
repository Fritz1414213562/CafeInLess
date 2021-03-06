#ifndef CAFEINLESS_MATRIX_FILE_READER_HPP
#define CAFEINLESS_MATRIX_FILE_READER_HPP
#include<CafeInLess/IO.dir/MatrixFileParser.hpp>
#include<coffee-makers/Containers/Containers.hpp>
#include<string>
#include<array>

namespace CafeInLess::IO {

class MatrixFileReader : public MatrixFileParser {


	template<typename scalarT>
	using MatX = makers::Matrix<scalarT, makers::Variable, makers::Variable>;

public:

	MatrixFileReader() : MatrixFileParser() {}
	MatrixFileReader(const std::string& inputfile_name) : MatrixFileParser(inputfile_name) {}
	~MatrixFileReader() = default;


	template<typename scalarT>
	auto read_Matrix() -> MatX<scalarT> {
//	auto read_Matrix() -> decltype(std::declval<MatX<scalarT>&>()) {
		open_File();

		std::array<int, 2> mat_shape = read_MatShape();
		MatX<scalarT> result(mat_shape[0], mat_shape[1]);
		std::string s_mat_data = read_Block();
		close_File();

		for (int i_datum = 0; i_datum < result.size(); ++i_datum) {
			std::size_t pos_of_datum = i_datum * sizeof(scalarT);
			result[i_datum] = read_Binary_as<scalarT>(&s_mat_data.at(pos_of_datum));
		}

		return result;
	}



	template<typename scalarT>
	auto read_Matrix(const std::string& input_name) -> MatX<scalarT> {
		open_File(input_name);

		std::array<int, 2> mat_shape = read_MatShape();
		const int& mat_size = mat_shape[0] * mat_shape[1];
//		std::cout << "Row Size -> " << mat_shape[0] << std::endl;
//		std::cout << "Column Size -> " << mat_shape[1] << std::endl << std::endl;

		MatX<scalarT> result(mat_shape[0], mat_shape[1]);
		std::string s_mat_data = read_Block();
		close_File();
//		std::cout << input_name << " closed" << std::endl;

		for (int i_datum = 0; i_datum < mat_size; ++i_datum) {
			std::size_t pos_of_datum = i_datum * sizeof(scalarT);
			result[i_datum] = read_Binary_as<scalarT>(&s_mat_data.at(pos_of_datum));
		}


		return result;
	}

};

}

#endif /* CAFEINLESS_MATRIX_FILE_READER_HPP */
