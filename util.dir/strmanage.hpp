#ifndef CAFEINLESS_STR_MANAGE_HPP
#define CAFEINLESS_STR_MANAGE_HPP
#include<string>
#include<vector>

namespace CafeInLess::strmanage {


std::vector<std::string> split_String(const std::string& line, const char& delimiter, const char& ignore_char) {

	std::vector<std::string> result;
	std::string buffer;

	for (const char& Character : line) {
		if (((Character == delimiter) || (Character == ignore_char)) && (buffer.empty())) continue;
		else if ((Character == delimiter) && (!buffer.empty())) {
			result.push_back(buffer);
			buffer.clear();
		}
		else {
			buffer.push_back(Character);
		}
	}
	if (!buffer.empty()) {
		result.push_back(buffer);
	}
	return result;
}

}


#endif /* CAFEINLESS_STR_MANAGE_HPP */
