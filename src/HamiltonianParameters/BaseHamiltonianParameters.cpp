#include "HamiltonianParameters/BaseHamiltonianParameters.hpp"


namespace GQCG {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param ao_basis_sptr
 */
BaseHamiltonianParameters::BaseHamiltonianParameters(std::shared_ptr<GQCG::AOBasis> ao_basis_sptr) :
    ao_basis_sptr (std::move(ao_basis_sptr))
{}



/*
 *  DESTRUCTOR
 */

/**
 *  Provide a pure virtual destructor to make the class abstract
 */
BaseHamiltonianParameters::~BaseHamiltonianParameters() {}



}  // namespace GQCG
