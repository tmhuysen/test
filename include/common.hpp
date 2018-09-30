#ifndef GQCG_COMMON_HPP
#define GQCG_COMMON_HPP


#include <cstdlib>
#include <vector>
#include <Eigen/Dense>


namespace GQCG {


typedef std::vector<size_t> Vectoru;
typedef std::vector<Vectoru> Matrixu;
typedef size_t address_pair[2];
typedef std::shared_ptr<size_t> size_t_sptr;
typedef std::shared_ptr<Eigen::MatrixXd> eigen_matd_sptr;


}  // namespace GQCG


#endif  // GQCG_COMMON_HPP
