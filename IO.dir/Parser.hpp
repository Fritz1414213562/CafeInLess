#ifndef PARSER_BASE_HPP
#define PARSER_BASE_HPP
#include<CafeInLess/IO.dir/ErrorMessage.hpp>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>

namespace CafeInLess::IO {

class FileParser {

public:

	FileParser() = default;
	FileParser(const std::string& inputfile_name) : input_name(inputfile_name) {
		open_File();
	}

	~FileParser() = default;

	void open_File(const std::string& inputfile_name) {
		input_name = inputfile_name;
		open_File();
	}

	void close_File() {if (input_file.is_open()) input_file.close();}

protected:

	CafeInLess::IO::Error_Output eout = CafeInLess::IO::Error_Output();
	std::ifstream input_file;

	void open_File() {
		close_File();
		input_file.open(input_name);
		if (!input_file.is_open()) eout("The file,", input_name, "could not be found.");
	}


private:

//	const char comment_out_char = '*';
//	const char delimiter = ' ';
	std::string input_name;

};
}

#endif /* PARSER_BASE_HPP */
