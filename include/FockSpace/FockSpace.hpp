#ifndef GQCG_FOCKSPACE_HPP
#define GQCG_FOCKSPACE_HPP


#include "BaseFockSpace.hpp"



namespace GQCG {


class FockSpace: public BaseFockSpace {
private:
    const size_t N;  // number of electrons
    const size_t dim;  // dimension of the Fock space
    Matrixu vertex_weights;



public:
    // CONSTRUCTORS
    /**
     *  Constructor given a @param K, N
     */
    explicit FockSpace(size_t K, size_t N);


    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    ~FockSpace() override = default;


    // PUBLIC METHODS
    /**
     *  @return ONV with the corresponding address in the considered space
     */
    ONV get_ONV(size_t address) override;

    /**
     *  @return weights as size_t from the vertex_weight matrix associated with the ONVs in the Fock space
     */
    size_t get_vertex_weights(size_t p, size_t m) const { return this->vertex_weights[p][m];}
};


}  // namespace GQCG


#endif //GQCG_FOCKSPACE_HPP
