#define BOOST_TEST_MODULE "TwoElectronOperator"


#include "Operator/TwoElectronOperator.hpp"

#include <cpputil.hpp>

#include "miscellaneous.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( TwoElectronOperator_constructor ) {

    // Check a correct constructor
    Eigen::Tensor<double, 4> tensor (3, 3, 3, 3);
    GQCG::TwoElectronOperator O (tensor);


    // Check a faulty constructor
    Eigen::Tensor<double, 4> tensor2 (3, 3, 3, 2);
    BOOST_CHECK_THROW(GQCG::TwoElectronOperator O2 (tensor2), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( TwoElectronOperator_getters ) {

    Eigen::Tensor<double, 4> tensor (3, 3, 3, 3);
    GQCG::TwoElectronOperator O (tensor);

    O.get_matrix_representation();
}


BOOST_AUTO_TEST_CASE ( TwoElectronOperator_transform_trivial ) {

    // Let's test a trivial transformation: i.e. with T being a unit matrix
    Eigen::Tensor<double, 4> g (3, 3, 3, 3);
    g.setRandom();
    GQCG::TwoElectronOperator G (g);

    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(3, 3);
    G.transform(T);

    BOOST_CHECK(cpputil::linalg::areEqual(g, G.get_matrix_representation(), 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( TwoElectronOperator_transform_olsens ) {

    // We can find a reference algorithm in the olsens module from Ayer's lab
    Eigen::Tensor<double, 4> g_transformed_ref (2, 2, 2, 2);
    cpputil::io::readArrayFromFile("../tests/data/rotated_two_electron_integrals_olsens.data", g_transformed_ref);

    // Set an example transformation matrix and two-electron integrals tensor
    Eigen::MatrixXd T (2, 2);
    T << 1, 2, 3, 4;

    Eigen::Tensor<double, 4> g (2, 2, 2, 2);
    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 2; k++) {
                for (size_t l = 0; l < 2; l++) {
                    g(i, j, k, l) = l + 2*k + 4*j + 8*i;
                }
            }
        }
    }
    GQCG::TwoElectronOperator G (g);
    G.transform(T);

    BOOST_CHECK(cpputil::linalg::areEqual(G.get_matrix_representation(), g_transformed_ref, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( TwoElectronOperator_rotate_throws ) {

    // Create a random TwoElectronOperator
    size_t dim = 3;
    Eigen::Tensor<double, 4> g (dim, dim, dim, dim);
    g.setRandom();
    GQCG::TwoElectronOperator G (g);


    // Check if a non-unitary matrix as transformation matrix causes a throw
    Eigen::MatrixXd U (Eigen::MatrixXd::Random(dim, dim));
    BOOST_CHECK_THROW(G.rotate(U), std::invalid_argument);


    // Check if a unitary matrix as transformation matrix is accepted
    G.rotate(Eigen::MatrixXd::Identity(dim, dim));
}


BOOST_AUTO_TEST_CASE ( TwoElectronOperator_rotate_JacobiRotationParameters ) {

    // Create a random TwoElectronOperator
    size_t dim = 5;
    Eigen::Tensor<double, 4> g (dim, dim, dim, dim);
    g.setRandom();
    GQCG::TwoElectronOperator G1 (g);
    GQCG::TwoElectronOperator G2 (g);


    // Check that using a Jacobi transformation (rotation) matrix as U is equal to the custom transformation (rotation)
    // with custom JacobiRotationParameters
    GQCG::JacobiRotationParameters jacobi_rotation_parameters (4, 2, 56.81);

    auto U = GQCG::jacobiRotationMatrix(jacobi_rotation_parameters, dim);


    G1.rotate(jacobi_rotation_parameters);
    G2.rotate(U);


    BOOST_CHECK(cpputil::linalg::areEqual(G1.get_matrix_representation(), G2.get_matrix_representation(), 1.0e-12));
}
