#ifndef CAFEINLESS_MATRIX_FILE_PARSER_HPP
#define CAFEINLESS_MATRIX_FILE_PARSER_HPP
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<cstring>
#include<array>

namespace CafeInLess::IO {

class MatrixFileParser {

protected:

	MatrixFileParser() = default;
	MatrixFileParser(const std::string& inputfile_name) : input_name(inputfile_name) {
		open_File();
	}

	~MatrixFileParser() = default;

	void open_File(const std::string& inputfile_name) {
		input_name = inputfile_name;
		open_File();
	}

	void close_File() {if (input_file.is_open()) input_file.close();}


	void open_File() {
		close_File();
		input_file.open(input_name, std::ios::in | std::ios::binary);
		if (!input_file.is_open()) 
			throw std::runtime_error("FileOpenError: The file, " + input_name + " could not be found.");
	}

	template<typename T>
	T read_Binary_as(const char *char_ptr) {
		T result;
		std::memcpy(std::addressof(result), char_ptr, sizeof(T));
		return result;
	}


	std::string read_Block() {
		std::int64_t block_size;
		std::vector<char> buffer;

		input_file.read(reinterpret_cast<char*>(&block_size), sizeof(block_size));
		buffer.resize(block_size);
		input_file.read(buffer.data(), block_size);

		std::int64_t check_block_size;
		input_file.read(reinterpret_cast<char*>(&check_block_size), sizeof(check_block_size));
	//	std::cout << check_block_size << std::endl;

		if (block_size != check_block_size)
			throw std::runtime_error(
			"Left Block size is not consistent with to Right one");
		std::string result(buffer.begin(), buffer.end());

		return result;
	}


	std::array<int, 2> read_MatShape() {
		std::string mat_shape_block = read_Block();
		int row_size = read_Binary_as<int>(&mat_shape_block.at(0));
		int col_size = read_Binary_as<int>(&mat_shape_block.at(sizeof(int)));
	//	std::cout << row_size << " " << col_size << std::endl;
		std::array<int, 2> result = {row_size, col_size};
		return result;
	}


	std::string input_name;
	std::ifstream input_file;



};
}

#endif /* CAFEINLESS_MATRIX_FILE_PARSER_HPP */
