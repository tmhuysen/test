#ifndef GQCG_SELECTEDFOCKSPACE_HPP
#define GQCG_SELECTEDFOCKSPACE_HPP


#include "BaseFockSpace.hpp"



namespace GQCG {


class SelectedFockSpace: public GQCG::BaseFockSpace {
private:
    const size_t N_A;  // number of alpha_electrons
    const size_t N_B;  // number of beta_electrons
    const size_t full_dim;  // dimension of the Fock space
    const size_t considered_dim;  // dimension of the Fock space
    /*
     *  UNFINISHED
     */



public:
    // CONSTRUCTORS
    /**
     *  Constructor given a @param K, N_B, N_A, address_list
     */
    explicit SelectedFockSpace(size_t K, size_t N_B, size_t N_A, std::vector<address_pair> address_list);


    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    ~SelectedFockSpace() override = default;


    // PUBLIC METHODS
    /**
     *  @return ONV with the corresponding address in the considered space
     */
    ONV get_ONV(size_t address) override;
};


}  // namespace GQCG


#endif //GQCG_SELECTEDFOCKSPACE_HPP
