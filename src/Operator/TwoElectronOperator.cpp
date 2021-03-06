#include "Operator/TwoElectronOperator.hpp"

#include <iostream>

#include "miscellaneous.hpp"


namespace GQCG {


/**
 *  Constructor based on a given @param tensor
 */
TwoElectronOperator::TwoElectronOperator(const Eigen::Tensor<double, 4>& tensor) :
    BaseOperator(tensor.dimensions()[0]),
    tensor (tensor)
{
    // Check if the given tensor is 'square'
    auto dims = tensor.dimensions();
    if ((dims[0] != dims[1]) || (dims[1] != dims[2]) || (dims[2] != dims[3]) ) {
        throw std::invalid_argument("The given tensor should have equal dimensions in every rank.");
    }
}



/*
 *  PUBLIC METHODS
 */

/**
 *  Transform the matrix representation of a two-electron operator using the transformation matrix @param T
 *
 *  Note that the transformation matrix @param T is used as
 *      b' = b T ,
 *  in which the basis functions are collected as elements of a row vector b
 */
void TwoElectronOperator::transform(const Eigen::MatrixXd& T) {

    // Since we're only getting T as a matrix, we should make the appropriate tensor to perform contractions
    // For the const Eigen::MatrixXd& argument, we need the const double in the template
    //      For more info, see: https://stackoverflow.com/questions/45283468/eigen-const-tensormap
    Eigen::TensorMap<Eigen::Tensor<const double, 2>> T_tensor (T.data(), T.rows(), T.cols());


    // We will have to do four single contractions, so we specify the contraction indices
    // Eigen3 does not document its tensor contraction clearly, so see the accepted answer on stackoverflow (https://stackoverflow.com/a/47558349/7930415):
    //      Eigen3 does not accept a way to specify the output axes: instead, it retains the order from left to right of the axes that survive the contraction.
    //      This means that, in order to get the right ordering of the axes, we will have to swap axes

    // g(T U V W)  T^*(V R) -> a(T U R W) but we get a(T U W R)
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair1 = {Eigen::IndexPair<int>(2, 0)};
    Eigen::array<int, 4> shuffle_1 {0, 1, 3, 2};

    // a(T U R W)  T(W S) -> b(T U R S) and we get b(T U R S), so no shuffle is needed
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair2 = {Eigen::IndexPair<int>(3, 0)};

    // T(U Q)  b(T U R S) -> c(T Q R S) but we get c(Q T R S)
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair3 = {Eigen::IndexPair<int>(0, 1)};
    Eigen::array<int, 4> shuffle_3 {1, 0, 2, 3};

    // T^*(T P)  c(T Q R S) -> g'(P Q R S) and we get g_SO(p q r s), so no shuffle is needed
    Eigen::array<Eigen::IndexPair<int>, 1> contraction_pair4 = {Eigen::IndexPair<int>(0, 0)};


    // Calculate the contractions. We write this as one large contraction to
    //  1) avoid storing intermediate contractions
    //  2) let Eigen3 figure out some optimizations
    Eigen::Tensor<double, 4> g_transformed = T_tensor.conjugate().contract(T_tensor.contract(this->tensor.contract(T_tensor.conjugate(), contraction_pair1).shuffle(shuffle_1).contract(T_tensor, contraction_pair2), contraction_pair3).shuffle(shuffle_3), contraction_pair4);

    this->tensor = g_transformed;
}


/**
 *  Rotate the matrix representation of a two-electron operator using a unitary rotation matrix @param U
 *
 *  Note that the rotation matrix @param U is used as
 *      b' = b U ,
 *  in which the basis functions are collected as elements of a row vector b.
 */
void TwoElectronOperator::rotate(const Eigen::MatrixXd& U) {

    // Check if the given matrix is actually unitary
    if (!U.isUnitary(1.0e-12)) {
        throw std::invalid_argument("The given matrix is not unitary.");
    }

    this->transform(U);
}


/**
 *  Rotate the matrix representation of a two-electron operator using the unitary Jacobi rotation matrix U constructed from the @param jacobi_rotation_parameters
 *
 *  Note that
 *      - the rotation matrix @param U is used as
 *          b' = b U ,
 *        in which the basis functions are collected as elements of a row vector b.
 *      - we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
 */
void TwoElectronOperator::rotate(const GQCG::JacobiRotationParameters& jacobi_rotation_parameters) {

    /**
     *  While waiting for an analogous Eigen::Tensor Jacobi module, we implement this rotation by constructing a
     *  Jacobi rotation matrix and then doing a rotation with it
     */

    auto dim = static_cast<size_t>(this->tensor.dimension(0));  // .dimension() returns a long
    Eigen::MatrixXd J = jacobiRotationMatrix(jacobi_rotation_parameters, dim);

    this->rotate(J);
}


}  // namespace GQCG
