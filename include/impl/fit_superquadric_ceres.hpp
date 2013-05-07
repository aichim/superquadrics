#pragma once

#include <pcl/common/centroid.h>
#include <pcl/common/pca.h>


#include "fit_superquadric_ceres.h"
#include "superquadric_formulas.h"



template<typename PointT, typename MatScalar>
SuperquadricFittingCeres<PointT, MatScalar>::SuperquadricFittingCeres ()
  : pre_align_ (true)
  , pre_align_axis_ (2)
{

}


template<typename PointT, typename MatScalar> void
SuperquadricFittingCeres<PointT, MatScalar>::preAlign (Eigen::Matrix<MatScalar, 4, 4> &transformation_prealign,
                                                       Eigen::Matrix<MatScalar, 3, 1> &variances)
{
  /// Compute the centroid
  Eigen::Vector4d centroid;
  pcl::compute3DCentroid (*input_, centroid);
  Eigen::Matrix<MatScalar, 4, 4> transformation_centroid (Eigen::Matrix<MatScalar, 4, 4>::Identity ());
  transformation_centroid (0, 3) = - centroid (0);
  transformation_centroid (1, 3) = - centroid (1);
  transformation_centroid (2, 3) = - centroid (2);

  std::cout << "centroid: " << centroid << std::endl;
  std::cout << "from " << input_->size () << " points" << std::endl;

  /// Compute the PCA
  pcl::PCA<PointT> pca;
  pca.setInputCloud (input_);
  Eigen::Vector3f eigenvalues = pca.getEigenValues ();
  Eigen::Matrix3f eigenvectors = pca.getEigenVectors ();

  std::cout << "eigenvectors:\n" << eigenvectors << std::endl;

  /// Align the first PCA axis with the prealign axis
  Eigen::Vector3f vec_aux = eigenvectors.col (0);
  eigenvectors.col (0) = eigenvectors.col (pre_align_axis_);
  eigenvectors.col (pre_align_axis_) = vec_aux;

  float aux_ev = eigenvalues (0);
  eigenvalues (0) = eigenvalues (pre_align_axis_);
  eigenvalues (pre_align_axis_) = aux_ev;


  Eigen::Matrix<MatScalar, 4, 4> transformation_pca (Eigen::Matrix<MatScalar, 4, 4>::Identity ());
  transformation_pca (0, 0) = eigenvectors (0, 0);
  transformation_pca (1, 0) = eigenvectors (1, 0);
  transformation_pca (2, 0) = eigenvectors (2, 0);

  transformation_pca (0, 1) = eigenvectors (0, 1);
  transformation_pca (1, 1) = eigenvectors (1, 1);
  transformation_pca (2, 1) = eigenvectors (2, 1);

  transformation_pca (0, 2) = eigenvectors (0, 2);
  transformation_pca (1, 2) = eigenvectors (1, 2);
  transformation_pca (2, 2) = eigenvectors (2, 2);

  transformation_prealign = transformation_pca * transformation_centroid;

  std::cout << "pre-align transformation:\n" << transformation_prealign << std::endl;

  /// Set the variances
  eigenvalues /= static_cast<float> (input_->size ());
  variances (0) = sqrt (eigenvalues (0));
  variances (1) = sqrt (eigenvalues (1));
  variances (2) = sqrt (eigenvalues (2));

  std::cout << "variances:\n" << variances << std::endl << eigenvalues << std::endl;
}



template<typename PointT, typename MatScalar> double
SuperquadricFittingCeres<PointT, MatScalar>::fit (SuperquadricParams<MatScalar> &parameters)
{
  Eigen::Matrix<MatScalar, 4, 4> transformation_prealign (Eigen::Matrix<MatScalar, 4, 4>::Identity ());
  Eigen::Matrix<MatScalar, 3, 1> variances;
  variances (0) = variances (1) = variances (2) = static_cast <MatScalar> (1.);

  if (pre_align_)
  {
    preAlign (transformation_prealign, variances);
    input_prealigned_.reset (new Cloud ());
    pcl::transformPointCloud (*input_, *input_prealigned_, transformation_prealign);
  }


  ceres::Problem problem;
  ceres::CostFunction *cost_function = new ceres::AutoDiffCostFunction<SuperquadricCostFunctor, 1, 1> (new SuperquadricCostFunctor);

  double x = 5.;
  for (size_t p_i = 0; p_i < input_prealigned_->size (); ++p_i)
  {
    problem.AddResidualBlock (cost_function, NULL, &x);
  }


  ceres::Solver::Options options;
  options.minimizer_type = ceres::TRUST_REGION;
  options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  options.minimizer_progress_to_stdout = true;

  ceres::Solver::Summary summary;
  ceres::Solve (options, &problem, &summary);

  std::cout << summary.BriefReport () << std::endl;
  std::cout << "x : " << x << std::endl;
}
