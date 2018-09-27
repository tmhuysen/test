#ifndef GQCG_BASEFOCKSPACE_HPP
#define GQCG_BASEFOCKSPACE_HPP


#include "common.hpp"
#include "ONV.hpp"



namespace GQCG {


class BaseFockSpace {
protected:
    const size_t K;  // number spatial orbitals


    // PROTECTED CONSTRUCTORS
    /**
     *  Protected constructor given a @param K, N and dim
     */
    explicit BaseFockSpace(size_t K);



public:
    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseFockSpace() = 0;


    // PURE VIRTUAL PUBLIC METHODS
    /**
     *  @return ONV with the corresponding address in the considered space
     */
    virtual ONV get_ONV(size_t address) = 0;
};


}  // namespace GQCG


#endif //GQCG_BASEFOCKSPACE_HPP
