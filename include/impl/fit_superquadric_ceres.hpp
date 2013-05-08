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

  double xvec[11];
  xvec[0] = xvec[1] = 1.;
  xvec[2] = variances (0) * 3.;
  xvec[3] = variances (1) * 3.;
  xvec[4] = variances (2) * 3.;
  xvec[5] = xvec[6] = xvec[7] = xvec[8] = xvec[9] = xvec[10] = 0.;

  for (size_t p_i = 0; p_i < input_prealigned_->size (); ++p_i)
  {
    PointT &point = (*input_prealigned_)[p_i];
    ceres::CostFunction *cost_function = new ceres::AutoDiffCostFunction<SuperquadricCostFunctor, 1, 11> (new SuperquadricCostFunctor (point));
//    ceres::CostFunction *cost_function = new ceres::NumericDiffCostFunction<SuperquadricCostFunctor, ceres::CENTRAL, 1, 11> (new SuperquadricCostFunctor (point));
    problem.AddResidualBlock (cost_function, NULL, xvec);
  }


  ceres::Solver::Options options;
  options.minimizer_type = ceres::TRUST_REGION;
  options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  options.minimizer_progress_to_stdout = true;

  ceres::Solver::Summary summary;
  ceres::Solve (options, &problem, &summary);

  std::cout << summary.BriefReport () << std::endl;

  printf ("x = ");
  for (size_t i = 0; i < 11; ++i)
    printf ("%f ", xvec[i]);
  printf ("\n");
  //  std::cout << "x : " << x << std::endl;



  Eigen::Matrix<MatScalar, 4, 4> &transformation = parameters.transform;
  transformation.setZero ();
  transformation (0, 3) = xvec[5];
  transformation (1, 3) = xvec[6];
  transformation (2, 3) = xvec[7];
  transformation (3, 3) = 1.;

  double angle_x = xvec[8],
      angle_y = xvec[9],
      angle_z = xvec[10];
  double aux_a = cos (angle_x),
      aux_b = sin (angle_x),
      aux_c = cos (angle_y),
      aux_d = sin (angle_y),
      aux_e = cos (angle_z),
      aux_f = sin (angle_z),
      aux_ad = aux_a * aux_d,
      aux_bd = aux_b * aux_d;

  transformation (0, 0) = aux_c * aux_e;
  transformation (0, 1) = -aux_c * aux_f;
  transformation (0, 2) = -aux_d;
  transformation (1, 0) = -aux_bd * aux_e + aux_a * aux_f;
  transformation (1, 1) = aux_bd * aux_f + aux_a * aux_e;
  transformation (1, 2) = -aux_b * aux_c;
  transformation (2, 0) = aux_ad * aux_e + aux_b * aux_f;
  transformation (2, 1) = -aux_ad * aux_f + aux_b * aux_e;
  transformation (2, 2) = aux_a * aux_c;



  parameters.e1 = xvec[0];
  parameters.e2 = xvec[1];
  parameters.a = xvec[2];
  parameters.b = xvec[3];
  parameters.c = xvec[4];
  parameters.transform = Eigen::Matrix<MatScalar, 4, 4> (transformation) * transformation_prealign;

  MatScalar final_error = computeSuperQuadricError<PointT, MatScalar> (input_,
                                                                       xvec[0], xvec[1], xvec[2], xvec[3], xvec[4],
                                                                       transformation);


  return (final_error);
}



template <typename PointT, typename MatScalar>
template <typename T> bool
SuperquadricFittingCeres<PointT, MatScalar>::SuperquadricCostFunctor::operator () (const T* const xvec, T* residual) const
{
  T e1 = xvec[0],
    e2 = xvec[1],
    a = xvec[2],
    b = xvec[3],
    c = xvec[4];
  Eigen::Matrix<T, 4, 4> transformation;
  transformation.setZero ();
  transformation (0, 3) = xvec[5];
  transformation (1, 3) = xvec[6];
  transformation (2, 3) = xvec[7];
  transformation (3, 3) = T (1.);

  T angle_x = xvec[8],
      angle_y = xvec[9],
      angle_z = xvec[10];
  T aux_a = cos (angle_x),
      aux_b = sin (angle_x),
      aux_c = cos (angle_y),
      aux_d = sin (angle_y),
      aux_e = cos (angle_z),
      aux_f = sin (angle_z),
      aux_ad = aux_a * aux_d,
      aux_bd = aux_b * aux_d;

  transformation (0, 0) = aux_c * aux_e;
  transformation (0, 1) = -aux_c * aux_f;
  transformation (0, 2) = -aux_d;
  transformation (1, 0) = -aux_bd * aux_e + aux_a * aux_f;
  transformation (1, 1) = aux_bd * aux_f + aux_a * aux_e;
  transformation (1, 2) = -aux_b * aux_c;
  transformation (2, 0) = aux_ad * aux_e + aux_b * aux_f;
  transformation (2, 1) = -aux_ad * aux_f + aux_b * aux_e;
  transformation (2, 2) = aux_a * aux_c;


  Eigen::Matrix<T, 4, 1> xyz (T (point_.x), T (point_.y), T (point_.z), T (1.));
  Eigen::Matrix<T, 4, 1> xyz_tr = transformation * xyz;

//  std::cout << "xyz_tr: " << xyz_tr[0] << " " << xyz_tr[1] << " " << xyz_tr[2] << std::endl;

  //    std::cout << xyz << " " << xyz_tr << std::endl;
  //    std::cout << transformation << std::endl;


  //    double term_1 = pow (fabs(xyz_tr[0] / a), 2./e2);
  //    double term_2 = pow (fabs(xyz_tr[1] / b), 2./e2);
  //    double term_3 = pow (fabs(xyz_tr[2] / c), 2./e1);
  //    double superellipsoid_f = pow (fabs(term_1 + term_2), e2/e1) + term_3;

  T op = (Eigen::Matrix<T, 3, 1> (xvec[5], xvec[6], xvec[7]) -
      Eigen::Matrix<T, 3, 1> (xyz_tr[0], xyz_tr[1], xyz_tr[2])).norm ();

//  std::cout << "op: " << op << std::endl;

//  std::cout << "params xvec: " << xvec[0] << " " << xvec[1] << " " << xvec[2] << " " << xvec[3] << " " << xvec[4] << std::endl;
//  std::cout << "params before: " << e1 << " " << e2 << " " << a << " " << b << " " << c << std::endl;
  residual[0] = op *superquadric_function<T> (xyz_tr[0], xyz_tr[1], xyz_tr[2], e1, e2, a, b, c);


//  std::cout << "residual (" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << ") = " << residual[0] << std::endl;

  //    double op = Eigen::Matrix<MatScalar, 3, 1> (xyz_tr[0], xyz_tr[1], xyz_tr[2]).norm ();


  //    fvec[i] = /*op */ (pow (superellipsoid_f, e1 / 2.) - 1.) * pow (a*b*c, 0.25);
  //    PCL_INFO ("fvec[%ld] = %f\n", i, fvec[i]);


  return (true);
}
