#ifndef CAFEINLESS_FILE_MANAGE_HPP
#define CAFEINLESS_FILE_MANAGE_HPP
#include<string>
#include<array>
#include<list>

namespace CafeInLess::filemanage {


template<std::size_t file_num>
inline std::array<std::string, file_num> align_Suffix(const std::array<std::string, file_num>& filenames, const std::array<std::string, file_num>& suffix_names) {

	std::array<std::string, file_num> result;

	std::array<std::string, file_num> answer_suffixes;
	for (std::size_t file_index = 0; file_index < file_num; ++file_index) {
		answer_suffixes[file_index] = "." + suffix_names[file_index];
	}

	for (std::size_t file_index = 0; file_index < file_num; ++file_index) {
		const std::string& filename = filenames[file_index];

		std::list<char> buffer_list;
	//	std::string::iterator read_point = filename.rbegin();
		auto read_point = filename.rbegin();

		// search the suffix of filename
		for (auto itr = filename.rbegin(); itr != filename.rend(); ++itr) {

			const char& suffix_str = *itr;

			if (suffix_str == '.') {
				buffer_list.push_front(suffix_str);
				break;
			}
			else {
				buffer_list.push_front(suffix_str);
				++read_point;
			}
		}
		if (read_point == filename.rend()) {
			std::cerr << "Error: Define the file format." << std::endl;
			std::cerr << filename << " <-" << std::endl;
			std::exit(1);
		}

		const std::string suffix(buffer_list.begin(), buffer_list.end());

		// search the index of answer suffix consistent with filename suffix

		std::size_t result_index = 0;
		for (const std::string& answer_suffix : answer_suffixes) {
			if (answer_suffix != suffix) {
				++result_index;
				continue;
			}
			else break;
		}

		if (result_index >= file_num) continue;
		else if (!result[result_index].empty()) {
			std::cerr << "Something wrong! The index is overlapped." << std::endl;
			std::cerr << "Probably, Two file exists specified as the same format." << std::endl;
			std::exit(1);
		}
		else {
			result[result_index] = filename;
		}
	}

	// check whether the empty exists.
	std::size_t check_idx = 0;
	for (const std::string& result_filename : result) {
		if (result_filename.empty()) {
			std::cerr << "Error: The strange file exists." << std::endl;
			std::cerr << "Please use this file format." << std::endl;
			std::cerr << " -> " << suffix_names[check_idx] << std::endl;
			std::exit(1);
		}
		++check_idx;
	}

	return result;

}


}

#endif /* CAFEINLESS_FILE_MANAGE_HPP */
