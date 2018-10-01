#ifndef GQCG_TWOELECTRONOPERATOR_HPP
#define GQCG_TWOELECTRONOPERATOR_HPP


#include <unsupported/Eigen/CXX11/Tensor>

#include "BaseOperator.hpp"


namespace GQCG {


/**
 *  A class that holds the matrix representation of a two-electron operator in an orbital basis
 */
class TwoElectronOperator : public GQCG::BaseOperator {
private:
    Eigen::Tensor<double, 4> tensor;  // the matrix unsigned_representation of the two-electron operator


public:
    // GETTERS
    Eigen::Tensor<double, 4> get_matrix_representation() const { return this->tensor; }
    double get(size_t index1,size_t index2,size_t index3,size_t index4){ return this->tensor(index1,index2,index3,index4);}


    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param tensor
     */
    explicit TwoElectronOperator(const Eigen::Tensor<double, 4>& tensor);


    // PUBLIC METHODS
    /**
     *  Transform the matrix representation of a two-electron operator using the transformation matrix @param T
     *
     *  Note that the transformation matrix @param T is used as
     *      b' = b T ,
     *  in which the basis functions are collected as elements of a row vector b
     */
    void transform(const Eigen::MatrixXd& T) override;

    /**
     *  Rotate the matrix representation of a two-electron operator using a unitary rotation matrix @param U
     *
     *  Note that the rotation matrix @param U is used as
     *      b' = b U ,
     *  in which the basis functions are collected as elements of a row vector b.
     */
    void rotate(const Eigen::MatrixXd& U) override;

    /**
     *  Rotate the matrix representation of a two-electron operator using the unitary Jacobi rotation matrix U constructed from the @param jacobi_rotation_parameters
     *
     *  Note that
     *      - the rotation matrix @param U is used as
     *          b' = b U ,
     *        in which the basis functions are collected as elements of a row vector b.
     *      - we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    void rotate(const GQCG::JacobiRotationParameters& jacobi_rotation_parameters) override;


    // FRIEND CLASSES
    friend class HamiltonianParameters;
};

typedef std::shared_ptr<GQCG::TwoElectronOperator> TwoElectronOperator_sptr;

}  // namespace GQCG


#endif  // GQCG_TWOELECTRONOPERATOR_HPP
