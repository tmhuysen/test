#include <HamiltonianBuilder/RestrictedHamiltonianBuilder.hpp>


namespace GQCG {

/*
 *  PROTECTED CONSTRUCTORS
 */

/**
 *  Protected constructor given a @param HamiltonianParameters
 */
RestrictedHamiltonianBuilder::RestrictedHamiltonianBuilder(HamiltonianParameters hamiltonian_parameters) :
hamiltonian_parameters(hamiltonian_parameters),
h(hamiltonian_parameters.h),
g(hamiltonian_parameters.g) {}



/*
 *  DESTRUCTOR
 */

/**
 *  Provide a pure virtual destructor to make the class abstract
 */
RestrictedHamiltonianBuilder::~RestrictedHamiltonianBuilder() {}


}  // namespace GQCG
