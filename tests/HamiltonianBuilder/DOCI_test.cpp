               #define BOOST_TEST_MODULE "DOCI"


#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( DOCI_constructor ) {
    // Create an AOBasis
    GQCG::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis_ptr = std::make_shared<GQCG::AOBasis>(water, "STO-3G");


    // Create One- and TwoElectronOperators (and a transformation matrix) with compatible dimensions
    size_t K = ao_basis_ptr->get_number_of_basis_functions();
    GQCG::OneElectronOperator S (Eigen::MatrixXd::Random(K, K));
    GQCG::OneElectronOperator H_core (Eigen::MatrixXd::Random(K, K));

    Eigen::Tensor<double, 4> g_tensor (K, K, K, K);
    g_tensor.setRandom();
    GQCG::TwoElectronOperator g (g_tensor);



    Eigen::MatrixXd C = Eigen::MatrixXd::Random(K, K);

    GQCG::HamiltonianParameters random_hamiltonian_parameters (ao_basis_ptr, S, H_core, g, C);
    GQCG::FockSpace fock_space(K,3);
    std::cout<<" K : "<<K<<" ";
    // Check if a correct constructor works and public methods don't fail
    BOOST_CHECK_NO_THROW(GQCG::DOCI random_doci(random_hamiltonian_parameters,fock_space));
    // Check if faulty constructor parameters throw an error
    GQCG::FockSpace fock_space_2(K+1,3);
    BOOST_CHECK_THROW(GQCG::DOCI random_doci_2(random_hamiltonian_parameters,fock_space_2),std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( DOCI_calculate_diagonal ) {
    // Create an AOBasis
    GQCG::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis_ptr = std::make_shared<GQCG::AOBasis>(water, "STO-3G");


    // Create One- and TwoElectronOperators (and a transformation matrix) with compatible dimensions
    size_t K = ao_basis_ptr->get_number_of_basis_functions();
    GQCG::OneElectronOperator S (Eigen::MatrixXd::Random(K, K));
    GQCG::OneElectronOperator H_core (Eigen::MatrixXd::Random(K, K));

    Eigen::Tensor<double, 4> g_tensor (K, K, K, K);
    g_tensor.setRandom();
    GQCG::TwoElectronOperator g (g_tensor);

    Eigen::MatrixXd C = Eigen::MatrixXd::Random(K, K);

    // Check if a correct constructor works and public methods don't fail
    GQCG::HamiltonianParameters random_hamiltonian_parameters (ao_basis_ptr, S, H_core, g, C);
    GQCG::FockSpace fock_space(K,3);
    GQCG::DOCI random_doci(random_hamiltonian_parameters,fock_space);
    Eigen::VectorXd x = random_doci.calculateDiagonal();
    BOOST_CHECK_NO_THROW(random_doci.constructHamiltonian());
    BOOST_CHECK_NO_THROW(random_doci.matrixVectorProduct(x));
}

/*
BOOST_AUTO_TEST_CASE ( h2_sto3g_szabo_plain ) {

   // In this test case, we will follow section 3.5.2 in Szabo.
   double ref_total_energy = -1.1167;


   // Create a Molecule and an AOBasis
   GQCG::Molecule h2 ("../tests/data/h2_szabo.xyz");
   auto ao_basis = std::make_shared<GQCG::AOBasis>(h2, "STO-3G");

   // Create the molecular Hamiltonian parameters for this molecule and basis
   auto mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);


   // Create a plain RHF SCF solver and solve the SCF equations
   GQCG::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2);
   plain_scf_solver.solve();
   auto rhf = plain_scf_solver.get_solution();


   // Check the total energy
   double total_energy = rhf.get_electronic_energy() + h2.calculateInternuclearRepulsionEnergy();
   BOOST_CHECK(std::abs(total_energy - ref_total_energy) < 1.0e-04);
}
*/

// dim = 120
BOOST_AUTO_TEST_CASE ( DOCI_beh_cation_klaas_dense ) {

   // Klaas' reference DOCI energy for BeH+ (obtained through Caitlin)
   double reference_doci_energy = -14.8782216937;


   // Do a DOCI calculation based on a given FCIDUMP file
   auto ham_par = GQCG::readFCIDUMPFile("../tests/data/beh_cation_631g_caitlin.FCIDUMP");
   GQCG::FockSpace fock_space (16,2);
   GQCG::DOCI doci(ham_par,fock_space);
   Eigen::MatrixXd ham = doci.constructHamiltonian();
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(ham);
   auto doci_eigenvalue = eigenSolver.eigenvalues()(0);
   // Calculate the total energy
   double internuclear_repulsion_energy = 1.5900757460937498e+00;  // this comes straight out of the FCIDUMP file
   double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

   std::cout<<eigenSolver.eigenvalues()<<std::endl;
   BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}

/*
// dim = 120
BOOST_AUTO_TEST_CASE ( DOCI_lih_klaas_dense ) {

   // Klaas' reference DOCI energy for LiH (obtained through Caitlin)
   double reference_doci_energy = -8.0029560313;


   // Do a DOCI calculation based on a given FCIDUMP file
   libwint::SOBasis so_basis ("../tests/reference_data/lih_631g_caitlin.FCIDUMP", 16);  // 16 SOs
   ci::DOCI doci (so_basis, 4);  // 4 electrons

   // Specify solver options and solve the eigenvalue problem
   numopt::eigenproblem::DenseSolverOptions dense_options;
   doci.solve(&dense_options);


   // Calculate the total energy
   double internuclear_repulsion_energy = 9.6074293445896852e-01;  // this comes straight out of the FCIDUMP file
   double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;


   BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


// dim = 816
BOOST_AUTO_TEST_CASE ( DOCI_li2_klaas_dense ) {

   // Klaas' reference DOCI energy for Li2
   double reference_doci_energy = -15.1153976060;


   // Do a DOCI calculation based on a given FCIDUMP file
   libwint::SOBasis so_basis ("../tests/reference_data/li2_321g_klaas.FCIDUMP", 18);  // 18 SOs
   ci::DOCI doci (so_basis, 6);  // 6 electrons

   // Specify solver options and solve the eigenvalue problem
   numopt::eigenproblem::DenseSolverOptions dense_options;
   doci.solve(&dense_options);


   // Calculate the total energy
   double internuclear_repulsion_energy = 3.0036546888874875e+00;  // this comes straight out of the FCIDUMP file
   double test_doci_energy = doci.get_eigenvalue() + internuclear_repulsion_energy;

   BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}
 */