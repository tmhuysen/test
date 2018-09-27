#ifndef GQCG_COMMON_HPP
#define GQCG_COMMON_HPP


#include <cstdlib>
#include <vector>
#include <HamiltonianParameters/HamiltonianParameters.hpp>
#include <FockSpace/FockSpace.hpp>
#include "Operator/OneElectronOperator.hpp"
#include "Operator/TwoElectronOperator.hpp"
#include "AOBasis.hpp"



namespace GQCG {


typedef std::vector<size_t> Vectoru;
typedef std::vector<Vectoru> Matrixu;
typedef size_t address_pair[2];
typedef std::shared_ptr<GQCG::AOBasis> AOBasis_sptr;
typedef std::shared_ptr<GQCG::HamiltonianParameters> HamiltonianParameters_sptr;
typedef std::shared_ptr<GQCG::FockSpace> FockSpace_sptr;
typedef std::shared_ptr<GQCG::OneElectronOperator> OneElectronOperator_sptr;
typedef std::shared_ptr<GQCG::TwoElectronOperator> TwoElectronOperator_sptr;


}  // namespace GQCG


#endif  // GQCG_COMMON_HPP
