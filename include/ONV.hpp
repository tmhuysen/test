#ifndef GQCG_ONV_HPP
#define GQCG_ONV_HPP


#include "common.hpp"



namespace GQCG {


class ONV {
private:
    const size_t K;  // number of spatial orbitals
    const size_t N;  // number of electrons
    size_t unsigned_representation;  // unsigned unsigned_representation
    size_t_sptr occupation_indexes;  // the occupied orbital electron indexes


    // PRIVATE METHODS
    /**
     *  Extracts the positions of the set bits from the representation and places them in an array
     */
    void representationToArray(size_t l);

    /**
     *  Tests whether the assigned amount of electrons N and the amount of set bits in the representation match
     */
    void is_compatible(size_t l){;



public:
    // CONSTRUCTORS
    /**
     *  Constructor from a @param K orbitals, N electrons and a representation for the ONV
     */
    ONV(size_t K, size_t N, size_t unsigned_representation);


    // OPERATORS
    /**
     *  Overloading of operator<< for a GQCG::ONV to be used with streams
     */
    std::ostream& operator<<(std::ostream& os, const GQCG::ONV& onv);


    // GETTERS & SETTERS
    void set_representation(size_t unsigned_representation);

    /**
     *  @return occupied orbital based on the electron index
     */
    size_t get_occupied_orbital(size_t electron_index);


    // FRIEND CLASSES
    friend class FockSpace;
};


}  // namespace GQCG

#endif //GQCG_ONV_HPP
