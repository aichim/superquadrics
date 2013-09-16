#pragma once

#include <pcl/common/centroid.h>
#include <pcl/common/pca.h>

#include <ceres/ceres.h>

#include "superquadric_rigid_registration.h"
#include "superquadric_formulas.h"


////////////////////////////////////////////////////////////////////////////////
template<typename PointT, typename MatScalar>
sq::SuperquadricRigidRegistration<PointT, MatScalar>::SuperquadricRigidRegistration ()
{
  init_transform_ = Eigen::Matrix<MatScalar, 4, 4>::Identity ();
}


////////////////////////////////////////////////////////////////////////////////
template<typename PointT, typename MatScalar> double
sq::SuperquadricRigidRegistration<PointT, MatScalar>::fit (Eigen::Matrix<MatScalar, 4, 4> &transform)
{
  ceres::Problem problem;

  Eigen::Matrix<MatScalar, 4, 4> init_transform_inverse = init_transform_.inverse ();

  double xvec[6];
  /// TODO put in the init transform
  xvec[0] = init_transform_inverse (0, 3);
  xvec[1] = init_transform_inverse (1, 3);
  xvec[2] = init_transform_inverse (2, 3);
  Eigen::Matrix<MatScalar, 3, 3> rot_mat = init_transform_inverse.block (0, 0, 3, 3);
  Eigen::Matrix<MatScalar, 3, 1> euler_angles = rot_mat.eulerAngles (2, 0, 2);
  xvec[3] = euler_angles (0, 0);
  xvec[4] = euler_angles (1, 0);
  xvec[5] = euler_angles (2, 0);

  for (size_t p_i = 0; p_i < input_->size (); ++p_i)
  {
    const PointT &point = (*input_)[p_i];
    ceres::CostFunction *cost_function = new ceres::AutoDiffCostFunction<SuperquadricCostFunctor, 1, 6> (new SuperquadricCostFunctor (point, e1_, e2_, a_, b_, c_));
    problem.AddResidualBlock (cost_function, NULL, xvec);
  }

  ceres::Solver::Options options;
  options.minimizer_type = ceres::TRUST_REGION;
  options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
//  options.minimizer_progress_to_stdout = true;
  options.num_threads = 8;

  ceres::Solver::Summary summary;
  ceres::Solve (options, &problem, &summary);

  /// If we did not converge, return infinite error
  if (summary.termination_type != ceres::FUNCTION_TOLERANCE)
  {
    PCL_ERROR ("Did not converge.\n");
    return (std::numeric_limits<double>::infinity ());
  }

  std::cout << summary.BriefReport () << std::endl;

  printf ("x = ");
  for (size_t i = 0; i < 6; ++i)
    printf ("%f ", xvec[i]);
  printf ("\n");


  transform.setZero ();
  transform (0, 3) = xvec[0];
  transform (1, 3) = xvec[1];
  transform (2, 3) = xvec[2];
  transform (3, 3) = 1.;
  transform.block (0, 0, 3, 3) = Eigen::AngleAxis<MatScalar> (xvec[3], Eigen::Matrix<MatScalar, 3, 1>::UnitZ ()) *
                                 Eigen::AngleAxis<MatScalar> (xvec[4], Eigen::Matrix<MatScalar, 3, 1>::UnitX ()) *
                                 Eigen::AngleAxis<MatScalar> (xvec[5], Eigen::Matrix<MatScalar, 3, 1>::UnitZ ()).matrix ();


  MatScalar final_error = computeSuperQuadricError<PointT, MatScalar> (input_,
                                                                       e1_, e2_, a_, b_, c_,
                                                                       transform);

  return (final_error);
}


////////////////////////////////////////////////////////////////////////////////
template <typename PointT, typename MatScalar>
template <typename T> bool
sq::SuperquadricRigidRegistration<PointT, MatScalar>::SuperquadricCostFunctor::operator () (const T* const xvec, T* residual) const
{
  Eigen::Matrix<T, 4, 4> transformation;
  transformation.setZero ();
  transformation (0, 3) = xvec[0];
  transformation (1, 3) = xvec[1];
  transformation (2, 3) = xvec[2];
  transformation (3, 3) = T (1.);
  transformation.block (0, 0, 3, 3) = Eigen::AngleAxis<T> (xvec[3], Eigen::Matrix<T, 3, 1>::UnitZ ()) *
                                      Eigen::AngleAxis<T> (xvec[4], Eigen::Matrix<T, 3, 1>::UnitX ()) *
                                      Eigen::AngleAxis<T> (xvec[5], Eigen::Matrix<T, 3, 1>::UnitZ ()).matrix ();


  Eigen::Matrix<T, 4, 1> xyz (T (point_.x), T (point_.y), T (point_.z), T (1.));
  Eigen::Matrix<T, 4, 1> xyz_tr = transformation * xyz;

  residual[0] = superquadric_function<T> (xyz_tr[0], xyz_tr[1], xyz_tr[2], T (e1_), T (e2_), T (a_), T (b_), T (c_));


  return (true);
}
