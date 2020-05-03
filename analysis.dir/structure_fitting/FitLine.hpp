#ifndef CAFEINLESS_FITTING_FITLINE_HPP
#define CAFEINLESS_FITTING_FITLINE_HPP
#include<CafeInLess/analysis.dir/structure_fitting/Fitting_Base.hpp>
#include<coffee-makers/Containers.hpp>
#include<coffee-makers/Solver.hpp>
#include<vector>


namespace CafeInLess::analysis {

template<typename realT>
class FitLine : public Fitting_Base<FitLine<realT>> {


public:

	FitLine() = default;
	~FitLine() = default;


	void fit(const std::vector<makers::FixedVector<realT, 3>>& points);

	
	makers::FixedVector<realT, 3> center_of_mass() const {return _origin;}

	makers::FixedVector<realT, 3> direction() const {return _direction;}


private:

	makers::FixedVector<realT, 3> _origin;
	makers::FixedVector<realT, 3> _direction;

};



template<typename realT>
void FitLine<realT>::fit(const std::vector<makers::FixedVector<realT, 3>>& points) {

	FixedVector<realT, 3> sum(0., 0., 0.);
	for (auto iter = points.cbegin(); iter != points.cend(); ++iter)
		sum += *iter;
	
	_origin = sum * (1. / static_cast<realT>(points.size()));

	FixedMatrix
}


}


#endif /* CAFEINLESS_FITTING_FITLINE_HPP */
