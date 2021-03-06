#include "HamiltonianParameters/HamiltonianParameters.hpp"

#include "miscellaneous.hpp"


namespace GQCG {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param ao_basis_ptr, overlap @param S, one-electron operator @param h, two-electron
 *  operator @param g and a transformation matrix between the current molecular orbitals and the atomic orbitals
 *  @param C
 */
HamiltonianParameters::HamiltonianParameters(std::shared_ptr<GQCG::AOBasis> ao_basis_sptr, const GQCG::OneElectronOperator& S, const GQCG::OneElectronOperator& h, const GQCG::TwoElectronOperator& g, const Eigen::MatrixXd& C) :
    BaseHamiltonianParameters(std::move(ao_basis_sptr)),
    S (S),
    h (h),
    g (g),
    C (C)
{
    // Check if the dimensions of all matrix representations are compatible
    auto K = this->ao_basis_sptr->number_of_basis_functions;

    if ((S.dim != K) || (h.dim != K) || (g.dim != K) || (C.cols() != K) || (C.rows() != K)) {
        throw std::invalid_argument("The dimensions of the operators and coefficient matrix are incompatible.");
    }
}



/*
 *  PUBLIC METHODS
 */

/**
 *  Given a transformation matrix @param T that links the new molecular orbital basis to the old molecular orbital basis,
 *  in the sense that
 *       b' = b T ,
 *  in which the molecular orbitals are collected as elements of a row vector b, transform
 *      - the one-electron interaction operator (i.e. the core Hamiltonian)
 *      - the two-electron interaction operator
 *
 *  Furthermore
 *      - @member S now gives the overlap matrix in the new molecular orbital basis
 *      - @member C is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
 */
void HamiltonianParameters::transform(const Eigen::MatrixXd& T) {

    this->S.transform(T);

    this->h.transform(T);
    this->g.transform(T);

    this->C = this->C * T;  // use the correct transformation formula for subsequent transformations
}


/**
 *  Given a unitary rotation matrix @param U that links the new molecular orbital basis to the old molecular orbital basis,
 *  in the sense that
 *       b' = b U ,
 *  in which the molecular orbitals are collected as elements of a row vector b, transform
 *      - the one-electron interaction operator (i.e. the core Hamiltonian)
 *      - the two-electron interaction operator
 *
 *  Furthermore, @member C is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
 */
void HamiltonianParameters::rotate(const Eigen::MatrixXd& U) {

    // A rotation leaves the overlap matrix invariant, so we don't have to transform it

    this->h.rotate(U);
    this->g.rotate(U);

    this->C = this->C * U;
}


/**
 *  Given @param jacobi_rotation_parameters that represent a unitary rotation matrix @param U (using a (cos, sin, -sin, cos) definition for the Jacobi rotation matrix) that links the new molecular orbital basis to the old molecular orbital basis,
 *  in the sense that
 *       b' = b U ,
 *  in which the molecular orbitals are collected as elements of a row vector b, transform
 *      - the one-electron interaction operator (i.e. the core Hamiltonian)
 *      - the two-electron interaction operator
 *
 *  Furthermore @member C is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
 */
void HamiltonianParameters::rotate(const GQCG::JacobiRotationParameters& jacobi_rotation_parameters) {

    // A rotation leaves the overlap matrix invariant, so we don't have to transform it

    this->h.rotate(jacobi_rotation_parameters);
    this->g.rotate(jacobi_rotation_parameters);


    // Create a Jacobi rotation matrix to transform the coefficient matrix with
    size_t K = this->h.dim;  // number of spatial orbitals
    auto J = GQCG::jacobiRotationMatrix(jacobi_rotation_parameters, K);
    this->C = this->C * J;
}



}  // namespace GQCG
