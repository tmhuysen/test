#ifndef GQCG_ONV_HPP
#define GQCG_ONV_HPP


#include "common.hpp"



namespace GQCG {


class ONV {
private:
    const size_t K;  // number of spatial orbitals
    const size_t N;  // number of electrons
    size_t unsigned_representation;  // unsigned representation
    size_t_sptr occupation_indexes;  // the occupied orbital electron indexes


    // PRIVATE METHODS
    /**
     *  Extracts the positions of the set bits from the representation and places them in an array
     */
    void update();

    /**
     *  Tests whether the assigned amount of electrons N and the amount of set bits in the representation match
     */
    void is_compatible(size_t l);



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
    friend std::ostream& operator<<(std::ostream& os, const GQCG::ONV& onv);


    // GETTERS & SETTERS
    void set_representation(size_t unsigned_representation);
    size_t get_urepresentation(){ return unsigned_representation;}

    /**
     *  @return occupied orbital based on the electron index
     */
    size_t get_occupied_orbital(size_t electron_index);


    // PUBLIC METHODS
    /**
     *  @return if the index is occupied (i.e. 1) for the @param p-th spatial orbital, starting from 0
     *  @param p is the lexical index (i.e. read from right to left)
     */
    bool isOccupied(size_t p){
        if (p > this->K-1) {
            throw std::invalid_argument("The index is out of the bitset bounds");
        }

        size_t operator_string = 1U << p;
        return this->unsigned_representation & operator_string;
    }


    /**
     *  performs an annihilation @param p
     *  !!! IMPORTANT: does not update the occupation array if required call "update()" !!!
     */
    void annihilate(size_t p){
        size_t operator_string = 1U << p;
        this->unsigned_representation &= ~operator_string;
    }


    /**
     *  performs a creation @param p
     *  !!! IMPORTANT: does not update the occupation array if required call "update()" !!!
     */
    void create(size_t p){
        size_t operator_string = 1U << p;
        this->unsigned_representation |= operator_string;
    }


    // FRIEND CLASSES
    friend class FockSpace;
};


}  // namespace GQCG

#endif //GQCG_ONV_HPP
