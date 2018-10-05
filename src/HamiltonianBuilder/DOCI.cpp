#include "HamiltonianBuilder/DOCI.hpp"



namespace GQCG {


/*
 *  CONSTRUCTORS
 */
DOCI::DOCI(HamiltonianParameters hamiltonian_parameters, FockSpace fock_space):
    RestrictedHamiltonianBuilder(hamiltonian_parameters),
    fock_space(fock_space),
    dim(fock_space.get_dimension())

{
    auto K = this->hamiltonian_parameters->ao_basis_sptr->number_of_basis_functions;
    if(K != this->fock_space->K){
        throw std::invalid_argument("Basis functions of the Fock space and AObasis are incompatible.");
    }
    //this->diagonal = calculateDiagonal();
}



/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @return Hamiltonian matrix as an Eigen::MatrixXd
 */
Eigen::MatrixXd DOCI::constructHamiltonian() {
    Eigen::MatrixXd result_matrix = Eigen::MatrixXd::Zero(dim,dim);
    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one.
    ONV onv = fock_space->get_ONV(0);  // spin string with address 0

    for (size_t I = 0; I < this->dim; I++) {  // I loops over all the addresses of the spin strings
        if(I>0){
            this->fock_space->setNext(onv);
        }
        // Diagonal contribution
        result_matrix(I, I) += this->diagonal(I);

        // Off-diagonal contribution
        for (size_t e1 = 0; e1 < this->fock_space->N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupied_orbital(e1);
            for (size_t q = 0; q < p; q++) {  // q loops over SOs
                if (!onv.isOccupied(q)) {  // if q not in I
                    onv.annihilate(p);
                    onv.create(q);

                    size_t J = fock_space->get_address(onv);  // J is the address of a string that couples to I

                    // The loops are p->K and q<p. So, we should normally multiply by a factor 2 (since the summand is symmetric)
                    // However, we are setting both of the symmetric indices of Hamiltonian, so no factor 2 is required
                    result_matrix(I, J) += g->get(p, q, p, q);
                    result_matrix(J, I) += g->get(p, q, p, q);

                    onv.annihilate(q);  // reset the spin string after previous creation
                    onv.create(p);  // reset the spin string after previous annihilation

                }  // q < p loop
            }
        }  // p loop
    }  // address (I) loop

    return result_matrix;
}


/**
 *  @return the action of the DOCI Hamiltonian of the coefficient vector @param x.
 */
Eigen::VectorXd DOCI::matrixVectorProduct(const Eigen::VectorXd& x) {

    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one.
    ONV onv = fock_space->get_ONV(0);  // spin string with address 0


    // Diagonal contributions
    Eigen::VectorXd matvec = this->diagonal.cwiseProduct(x);


    // Off-diagonal contributions
    for (size_t I = 0; I < this->dim; I++) {  // I loops over all the addresses of the spin strings
        for (size_t e1 = 0; e1 < this->fock_space->N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupied_orbital(e1);
            for (size_t q = 0; q < p; q++) {  // q loops over SOs
                if (!onv.isOccupied(q)) {  // if q not in I

                    onv.annihilate(p);
                    onv.create(q);

                    size_t J = fock_space->get_address(onv);  // J is the address of a string that couples to I

                    // The loops are p->K and q<p. So, we should normally multiply by a factor 2 (since the summand is symmetric)
                    // However, we are setting both of the symmetric indices of Hamiltonian, so no factor 2 is required
                    matvec(I) += g->get(p, q, p, q)* x(J);
                    matvec(J) += g->get(p, q, p, q)* x(I);

                    onv.annihilate(q);  // reset the spin string after previous creation
                    onv.create(p);  // reset the spin string after previous annihilation

                }  // q < p loop
            }
        }  // p loop
    }  // address (I) loop

    return matvec;
}


/**
 *  @return the diagonal of the matrix representation of the DOCI Hamiltonian.
 */
Eigen::VectorXd DOCI::calculateDiagonal() {

    Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(this->dim);
    // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one.
    ONV onv = this->fock_space->get_ONV(0);  // onv with address 0

    for (size_t I = 0; I < this->dim; I++) {  // I loops over addresses of spin strings
        if(I>0){
            this->fock_space->setNext(onv);
        }
        for (size_t e1 = 0; e1 < this->fock_space->N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupied_orbital(e1);
            std::cout<<this->fock_space->K;
            std::cout<<p;
            this->diagonal(I) += 2 * this->h->get(p,p) + this->g->get(p,p,p,p);
            for (size_t e2 = 0; e2 < e1; e2++) {  // e2 (electron 2) loops over the (number of) electrons
                // Since we are doing a restricted summation q<p (and thus e2<e1), we should multiply by 2 since the summand argument is symmetric.
                size_t q = onv.get_occupied_orbital(e2);
                this->diagonal(I) += 2 * (2*this->g->get(p,p,q,q) - this->g->get(p,q,q,p));
            }  // e2 loop
        } // e1 loop
    }  // address (I) loop
    return diagonal;
}



}  // namespace ci
