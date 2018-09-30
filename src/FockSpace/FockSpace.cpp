#include "FockSpace/FockSpace.hpp"



namespace GQCG {


/*
 * PRIVATE METHODS
 */

/**
 *  In-place permute the unsigned representation of the @param ONV, giving the next bitstring permutation in reverse lexical ordering.
 *
 *      Examples:
 *          011 -> 101
 *          101 -> 110
 */
size_t ulongNextPermutation(size_t representation) {

    // t gets this->representation's least significant 0 bits set to 1
    unsigned long t = representation | (representation - 1UL);

    // Next set to 1 the most significant bit to change,
    // set to 0 the least significant ones, and add the necessary 1 bits.
    return (t + 1UL) | (((~t & (t+1UL)) - 1UL) >> (__builtin_ctzl(representation) + 1UL));
}




/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param K and @param N on which the dimension of the fock space is based.
 */
FockSpace::FockSpace(size_t K, size_t N) :
        BaseFockSpace(size_t K)
{
    auto dim_double = boost::math::binomial_coefficient<double>(static_cast<unsigned>(K), static_cast<unsigned>(N_P));
    this->dim = boost::numeric::converter<double, size_t>::convert(dim_double);
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return ONV with the corresponding address in the considered space
 */
FockSpace::ONV get_ONV(size_t address){
    size_t representation;
    if (this->N == 0) {
        representation = 0;
    }

    else {
        representation = 0;
        size_t m = this->N;  // counts the number of electrons in the spin string up to orbital p

        for (size_t p = K; p > 0; p--) {  // p is an orbital index
            size_t weight = get_vertex_weights(p-1, m);

            if (weight <= address) {  // the algorithm can move diagonally, so we found an occupied orbital
                address -= weight;
                representation |= (static_cast<T>(1) << (p - 1));  // set the (p-1)th bit: see (https://stackoverflow.com/a/47990)

                m--;  // since we found an occupied orbital, we have one electron less
                if (m == 0) {
                    break;
                }
            }
        }
    }

    return ONV(K,N,representation)

}


/**
 *  sets @param ONV to the next ONV in the space
 *  performs the ulongNextPermutation() function
 *  and updates the corresponding occupation indexes
 */
void FockSpace::setNext(ONV &onv) {
    onv.set(ulongNextPermutation(onv.unsigned_representation));
}



}  // namespace GQCG
