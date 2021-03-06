#define BOOST_TEST_MODULE "LibintCommunicator"


#include "LibintCommunicator.hpp"

#include <cpputil.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( atoms_interface ) {

    std::vector<GQCG::Atom> gqcg_atoms = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };

    std::vector<libint2::Atom> ref_libint_atoms = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };


    // Use the Libint interface to obtain a std::vector<libint2::Atom> from the GQCG ones
    auto test_libint_atoms = GQCG::LibintCommunicator::get().interface(gqcg_atoms);


    /**
     *  Implement a function object to compare (libint_atom) == (libint_atom)
     */
    struct LibintAtomEqual {
        double tolerance;
        explicit LibintAtomEqual(double tolerance) : tolerance (tolerance) {};

        bool operator()(const libint2::Atom& lhs, const libint2::Atom& rhs) {
            return (lhs.atomic_number == rhs.atomic_number) &&
                   (std::abs(lhs.x - rhs.x) < tolerance) &&
                   (std::abs(lhs.y - rhs.y) < tolerance) &&
                   (std::abs(lhs.z - rhs.z) < tolerance);
        }
    };


    // Check if the interfacing between the Atom types works
    BOOST_CHECK((ref_libint_atoms.size() == test_libint_atoms.size()) &&
                std::equal(ref_libint_atoms.begin(), ref_libint_atoms.end(), test_libint_atoms.begin(), LibintAtomEqual(1.0e-08)));
}


BOOST_AUTO_TEST_CASE( Szabo_integrals_h2_sto3g ) {

    // We will follow the example in Szabo, section 3.5.2, where it is stated that R = 1.4 a.u. = 0.740848 Angstrom
    GQCG::Molecule h2 ("../tests/data/h2_szabo.xyz");
    GQCG::AOBasis basis (h2, "STO-3G");
    BOOST_CHECK_EQUAL(basis.get_number_of_basis_functions(), 2);


    // Calculate some integrals
    auto S = GQCG::LibintCommunicator::get().calculateOneElectronIntegrals(libint2::Operator::overlap, basis);

    auto T = GQCG::LibintCommunicator::get().calculateOneElectronIntegrals(libint2::Operator::kinetic, basis);
    auto V = GQCG::LibintCommunicator::get().calculateOneElectronIntegrals(libint2::Operator::nuclear, basis);
    auto H_core = GQCG::OneElectronOperator(T.get_matrix_representation() + V.get_matrix_representation());

    auto g = GQCG::LibintCommunicator::get().calculateTwoElectronIntegrals(libint2::Operator::coulomb, basis);


    // Fill in the reference values from Szabo
    Eigen::MatrixXd ref_S (2, 2);
    ref_S << 1.0,    0.6593,
            0.6593, 1.0;

    Eigen::MatrixXd ref_T (2, 2);
    ref_T << 0.7600, 0.2365,
            0.2365, 0.7600;

    Eigen::MatrixXd ref_H_core (2, 2);
    ref_H_core << -1.1204, -0.9584,
            -0.9584, -1.1204;

    BOOST_CHECK(S.get_matrix_representation().isApprox(ref_S, 1.0e-04));
    BOOST_CHECK(T.get_matrix_representation().isApprox(ref_T, 1.0e-04));
    BOOST_CHECK(H_core.get_matrix_representation().isApprox(ref_H_core, 1.0e-04));


    // The two-electron integrals in Szabo are given in chemist's notation, so this confirms that the LibintCommunicator gives them in chemist's notation as well
    BOOST_CHECK(std::abs(g.get_matrix_representation()(0,0,0,0) - 0.7746) < 1.0e-04);
    BOOST_CHECK(std::abs(g.get_matrix_representation()(0,0,0,0) - g.get_matrix_representation()(1,1,1,1)) < 1.0e-12);

    BOOST_CHECK(std::abs(g.get_matrix_representation()(0,0,1,1) - 0.5697) < 1.0e-04);

    BOOST_CHECK(std::abs(g.get_matrix_representation()(1,0,0,0) - 0.4441) < 1.0e-04);
    BOOST_CHECK(std::abs(g.get_matrix_representation()(1,0,0,0) - g.get_matrix_representation()(1,1,1,0)) < 1.0e-12);

    BOOST_CHECK(std::abs(g.get_matrix_representation()(1,0,1,0) - 0.2970) < 1.0e-04);
}


BOOST_AUTO_TEST_CASE( HORTON_integrals_h2o_sto3g ) {

    // Set up a basis
    GQCG::Molecule water ("../tests/data/h2o.xyz");
    GQCG::AOBasis basis (water, "STO-3G");
    auto nbf = basis.get_number_of_basis_functions();


    // Calculate some integrals
    auto S = GQCG::LibintCommunicator::get().calculateOneElectronIntegrals(libint2::Operator::overlap, basis);
    auto T = GQCG::LibintCommunicator::get().calculateOneElectronIntegrals(libint2::Operator::kinetic, basis);
    auto V = GQCG::LibintCommunicator::get().calculateOneElectronIntegrals(libint2::Operator::nuclear, basis);

    auto g = GQCG::LibintCommunicator::get().calculateTwoElectronIntegrals(libint2::Operator::coulomb, basis);


    // Read in reference data from HORTON
    Eigen::MatrixXd ref_S (nbf, nbf);
    Eigen::MatrixXd ref_T (nbf, nbf);
    Eigen::MatrixXd ref_V (nbf, nbf);
    Eigen::Tensor<double, 4> ref_g (nbf, nbf, nbf, nbf);

    cpputil::io::readArrayFromFile("../tests/data/h2o_sto-3g_overlap_horton.data", ref_S);
    cpputil::io::readArrayFromFile("../tests/data/h2o_sto-3g_kinetic_horton.data", ref_T);
    cpputil::io::readArrayFromFile("../tests/data/h2o_sto-3g_nuclear_horton.data", ref_V);
    cpputil::io::readArrayFromFile("../tests/data/h2o_sto-3g_coulomb_horton.data", ref_g);


    // Check if the calculated integrals are equal to those of HORTON
    BOOST_CHECK(S.get_matrix_representation().isApprox(ref_S, 1.0e-08));
    BOOST_CHECK(T.get_matrix_representation().isApprox(ref_T, 1.0e-08));
    BOOST_CHECK(V.get_matrix_representation().isApprox(ref_V, 1.0e-08));
    BOOST_CHECK(cpputil::linalg::areEqual(g.get_matrix_representation(), ref_g, 1.0e-06));
}


