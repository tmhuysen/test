#include "ONV.hpp"


namespace GQCG {


/*
 *  PRIVATE METHODS
 */

/**
 *  Extracts the positions of the set bits from the representation and places them in an array
 */
void ONV::representationToArray(size_t l){
    int i = 0;
    while(l != 0){
        this->occupation_indexes.get()[i] = __builtin_ctzl(l);
        i++;
        l ^= (l & -l);
    }
}

/**
 *  Tests whether the assigned amount of electrons N and the amount of set bits in the representation match
 */
void ONV::is_compatible(size_t l){
    int i = 0;
    while(l != 0){
        if(i == this->N){
            throw std::invalid_argument("representation and electron count are not compatible");
        }
        i++;
        l ^= (l & -l);
    }
}



/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor from a @param K orbitals, N electrons and a representation for the ONV
 */
ONV::ONV(size_t K, size_t N, size_t representation): K(K), N(N), representation(representation){
    occupation_indexes = size_t_sptr(new size_t[N]);
    is_compatible(this->representation);  // throws error if the constructor parameters are not compatible;
    representation_to_array(this->representation);
}



/*
 *  OPERATORS
 */

/**
 *  @return occupied orbital based on the electron index
 */
size_t ONV::get_occupied_orbital(size_t electron_index) { return occupation_indexes.get()[electron_index]; }

/**
 *  Overloading of operator<< for a GQCG::ONV to be used with streams
 */
std::ostream& operator<<(std::ostream& os, const GQCG::ONV& onv) {
    const char* begin = reinterpret_cast<const char*>(&representation);
    const char* end = beg + K;
    while(begin != end) {
        os << std::bitset<CHAR_BIT>(*begin++);
    }
    return os;
}



/*
 *  SETTERS & GETTERS
 */

/**
 *  @set to a new representation
 */
void ONV::set_representation(size_t unsigned_representation) {
    this->unsigned_representation = unsigned_representation;
    representation_to_array(unsigned_representation);
}

/**
 *  @return occupied orbital based on the electron index
 */
size_t ONV::get_occupied_orbital(size_t electron_index) { return occupation_indexes.get()[electron_index]; }


}  // namespace GQCG
