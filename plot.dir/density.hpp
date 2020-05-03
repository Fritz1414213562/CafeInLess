#ifndef CAFEINLESS_DATA_DENSITY_HPP
#define CAFEINLESS_DATA_DENSITY_HPP
#include<vector>
#include<utility>
#include<limits>
#include<cmath>
#include<memory>


namespace CafeInLess::plot {


template<typename probT>
class HistogramConverter {

	static_assert(std::is_floating_point<probT>::value, "Template paramaters of Histogram Converter limit a floating point");

public:

	HistogramConverter() = default; 

	~HistogramConverter() = default;

	template<typename T>
	std::vector<std::pair<T, probT>> compute(const std::vector<T>& data, const T& data_width) const;


private:


//	template<typename T>
//	struct Histogram;


};

//template<typename probT>
//template<typename T>
//std::vector<std::pair<T, probT>> HistogramConverter<probT>::compute(const std::vector<T>& data, const T& data_width) const {
//	std::unique_ptr<Histogram<T>> device_ptr = std::make_unique<Histogram<T>>();
//	return device_ptr->compute(data, data_width);
//}


template<typename probT>
template<typename T>
std::vector<std::pair<T, probT>> HistogramConverter<probT>::compute(const std::vector<T>& data, const T& data_width) const {

	T min_datum = std::numeric_limits<T>::max();
	T max_datum = std::numeric_limits<T>::min();
	probT data_number = static_cast<probT>(data.size());

	for (std::size_t i_datum = 0; i_datum < data.size(); ++i_datum) {
		const T& datum = data[i_datum];
		if (datum < min_datum) min_datum = datum;
		if (datum > max_datum) max_datum = datum;
	}

	if (min_datum > max_datum) throw std::invalid_argument("Data searching is Failed");

	const std::size_t& bin_number = std::floor((max_datum - min_datum) / data_width) + 1;

	std::vector<std::pair<T, probT>> retval(bin_number, std::make_pair(0., 0.));

	for (std::size_t i_bin = 0; i_bin < bin_number; ++i_bin) {
		retval[i_bin].first = min_datum + data_width * i_bin;
	}

	for (const T& datum : data) {
		std::size_t bin_index = std::floor((datum - min_datum) / data_width);
		++retval[bin_index].second;
	}

	for (std::size_t i_bin = 0; i_bin < bin_number; ++i_bin) {
		retval[i_bin].second /= data_number;
	}

	return retval;
}
//
//
//template<typename probT>
//template<typename T>
//struct HistogramConverter<probT>::Histogram {
//
//public:
//
//	std::vector<std::pair<T, probT>> operator()(const std::vector<T>& data, const T& data_width) const {
//		return compute(data, data_width);
//	}
//
//	std::vector<std::pair<T, probT>> compute(const std::vector<T>& data, const T& data_width) const {
//
//		T min_datum = std::numeric_limits<T>::max();
//		T max_datum = std::numeric_limits<T>::min();
//		probT data_number = static_cast<probT>(data.size());
//
//		for (std::size_t i_datum = 0; i_datum < data.size(); ++i_datum) {
//			const T& datum = data[i_datum];
//			if (datum < min_datum) min_datum = datum;
//			if (datum > max_datum) max_datum = datum;
//		}
//
//		if (min_datum > max_datum) throw std::invalid_argument("Data searching is Failed");
//
//		const std::size_t& bin_number = std::ceil((max_datum - min_datum) / data_width);
//
//		std::vector<std::pair<T, probT>> retval(bin_number, std::make_pair(0., 0.));
//
//		for (std::size_t i_bin = 0; i_bin < bin_number; ++i_bin) {
//			retval[i_bin].first = min_datum + data_width * i_bin;
//		}
//
//		for (const T& datum : data) {
//			std::size_t bin_index = std::floor((datum - min_datum) / data_width);
//			++retval[bin_index].second;
//		}
//
//		for (std::size_t i_bin = 0; i_bin < bin_number; ++i_bin) {
//			retval[i_bin].second /= data_number;
//		}
//
//		return retval;
//	}
//
//};


//template<typename probT>
//template<typename T>
//struct HistogramConverter<probT>::Histogram<int, T> {
//
//public:
//
//	std::vector<std::pair<int, probT>> compute(const std::vector<int>& data, const int& data_width) const {
//
//		int min_datum = std::numeric_limits<int>::max();
//		int max_datum = std::numeric_limits<int>::min();
//		probT data_number = static_cast<probT>(data.size());
//
//		for (std::size_t i_datum = 0; i_datum < data.size(); ++i_datum) {
//			const int& datum = data[i_datum];
//			if (datum < min_datum) min_datum = datum;
//			if (datum > max_datum) max_datum = datum;
//		}
//
//		if (min_datum > max_datum) throw std::invalid_argument("Data searching is Failed");
//
//		const std::size_t& bin_number = std::ceil((max_datum - min_datum) / data_width);
//
//		std::vector<std::pair<int, probT>> retval(bin_number, std::make_pair(0., 0.));
//
//		for (std::size_t i_bin = 0; i_bin < bin_number; ++i_bin)
//			retval[i_bin].first = min_datum + data_width * i_bin;
//
//		for (const int& datum : data) {
//			std::size_t bin_index = std::floor((datum - min_datum) / data_width);
//			++retval[bin_index].second;
//		}
//
//		for (std::size_t i_bin = 0; i_bin < bin_number; ++i_bin)
//			retval[i_bin].second /= data_number;
//
//
//		return retval;
//	}
//
//
//	std::vector<std::pair<int, probT>> compute(const std::vector<int>& data) const {
//		int min_datum = std::numeric_limits<int>::max();
//		int max_datum = std::numeric_limits<int>::min();
//		probT data_number = static_cast<probT>(data.size());
//
//		for (std::size_t i_datum = 0; i_datum < data.size(); ++i_datum) {
//			const int& datum = data[i_datum];
//			if (datum < min_datum) min_datum = datum;
//			if (datum > max_datum) max_datum = datum;
//		}
//		if (min_datum > max_datum) throw std::invalid_argument("Data searching is failed");
//
//		const std::size_t& bin_number = max_datum - min_datum + 1;
//
//		std::vector<std::pair<int, probT>> retval(bin_number, std::make_pair(0, 0.));
//
//		for (std::size_t i_bin = 0; i_bin < bin_number; ++i_bin)
//			retval[i_bin].first = min_datum + i_bin;
//		
//		for (const int& datum : data)
//			++retval[datum - min_datum].second;
//
//		for (std::size_t i_bin = 0; i_bin < bin_number; ++i_bin)
//			retval[i_bin].second /= data_number;
//
//
//		return retval;
//	}
//
//};
//


//template<typename probT>
//template<typename T>
//std::vector<std::pair<T, probT>> HistogramConverter<probT>::compute(const std::vector<T>& data, const T& data_width) const {
//	return Histogram<T>(data, data_width);
//}



}

#endif /* CAFEINLESS_ARITHMETIC_DENSITY_HPP */
