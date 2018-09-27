#ifndef GQCG_BASEHAMILTONIANPARAMETERS_HPP
#define GQCG_BASEHAMILTONIANPARAMETERS_HPP


#include <memory>
#include "common.hpp"
#include "AOBasis.hpp"


namespace GQCG {


class BaseHamiltonianParameters {
protected:
    AOBasis_sptr ao_basis_sptr;  // the initial atomic orbitals

public:
    // CONSTRUCTOR
    /**
     *  Constructor based on a given @param ao_basis_sptr
     */
    explicit BaseHamiltonianParameters(AOBasis_sptr ao_basis_sptr);


    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseHamiltonianParameters() = 0;
};


}  // namespace GQCG


#endif  // GQCG_BASEHAMILTONIANPARAMETERS_HPP
