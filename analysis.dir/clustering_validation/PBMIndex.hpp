#ifndef CAFEINLESS_PBM_INDEX_HPP
#define CAFEINLESS_PBM_INDEX_HPP
#include<CafeInLess/analysis.dir/clustering_validation/Validation_Base.hpp>
#include<CafeInLess/analysis.dir/clustering_validation/PBMIndex_Flag.hpp>
#include<CafeInLess/IO.dir/StandardOutput.hpp>
#include<CafeInLess/IO.dir/ErrorMessage.hpp>
#include<CafeInLess/util.dir/arithmetic>
#include<vector>
#include<string>


namespace CafeInLess::analysis {

template<PBM_Index_Flag Dist_Flag>
class PBMIndex : ClusteringValidationBase<PBMIndex<Dist_Flag>> {
};


}


#endif /* CAFEINLESS_PBM_INDEX_HPP */
