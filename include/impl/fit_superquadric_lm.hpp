#pragma once

#include "fit_superquadric_lm.h"
#include "superquadric_formulas.h"
#include <pcl/common/centroid.h>
#include <unsupported/Eigen/NonLinearOptimization>
#include <pcl/common/pca.h>


////////////////////////////////////////////////////////////////////////////////
template<typename PointT, typename MatScalar>
sq::SuperquadricFittingLM<PointT, MatScalar>::SuperquadricFittingLM ()
  : pre_align_ (true)
  , pre_align_axis_ (2)
{
}


////////////////////////////////////////////////////////////////////////////////
template<typename PointT, typename MatScalar> void
sq::SuperquadricFittingLM<PointT, MatScalar>::preAlign (Eigen::Matrix<MatScalar, 4, 4> &transformation_prealign,
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
  transformation_pca (1, 0) = eigenvectors (0, 1);
  transformation_pca (2, 0) = eigenvectors (0, 2);

  transformation_pca (0, 1) = eigenvectors (1, 0);
  transformation_pca (1, 1) = eigenvectors (1, 1);
  transformation_pca (2, 1) = eigenvectors (1, 2);

  transformation_pca (0, 2) = eigenvectors (2, 0);
  transformation_pca (1, 2) = eigenvectors (2, 1);
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


////////////////////////////////////////////////////////////////////////////////
template<typename PointT, typename MatScalar> double
sq::SuperquadricFittingLM<PointT, MatScalar>::fit (SuperquadricParameters<MatScalar> &parameters)
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

  // Initialize the parameters based on the input cloud
  int n_unknowns = 11;
  /// e1, e2, a, b, c
  VectorX xvec (n_unknowns);
  xvec[0] = xvec[1] = 1.;
  xvec[2] = variances (0) * 3.;
  xvec[3] = variances (1) * 3.;
  xvec[4] = variances (2) * 3.;
  xvec[5] = xvec[6] = xvec[7] = xvec[8] = xvec[9] = xvec[10] = 0.;

  std::cout << "initial parameters:\n" << xvec << std::endl;

  /// TODO take care of the case with indices

  OptimizationFunctor functor (input_prealigned_->size (), this);

  //  Eigen::NumericalDiff<OptimizationFunctor> num_diff (functor);
  //  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<OptimizationFunctor>, MatScalar> lm (num_diff);
  Eigen::LevenbergMarquardt<OptimizationFunctor, MatScalar> lm (functor);
  int info = lm.minimize (xvec);

  std::cout << xvec << std::endl;


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


////////////////////////////////////////////////////////////////////////////////
template<typename PointT, typename MatScalar> int
sq::SuperquadricFittingLM<PointT, MatScalar>::OptimizationFunctor::operator () (const VectorX &xvec,
                                                                                VectorX &fvec) const
{
  const Cloud &points = *estimator_->input_prealigned_;

  MatScalar e1 = xvec[0],
      e2 = xvec[1],
      a = xvec[2],
      b = xvec[3],
      c = xvec[4];
  Eigen::Matrix<MatScalar, 4, 4> transformation;
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

  for (int i = 0; i < values (); ++i)
  {
    Eigen::Matrix<MatScalar, 4, 1> xyz (points[i].x, points[i].y, points[i].z, 1.);
    Eigen::Matrix<MatScalar, 4, 1> xyz_tr = transformation * xyz;

    double op = Eigen::Matrix<MatScalar, 3, 1> (xyz_tr[0], xyz_tr[1], xyz_tr[2]).norm ();

    fvec[i] = op *superquadric_function (xyz_tr[0], xyz_tr[1], xyz_tr[2], e1, e2, a, b, c);
  }
  return (0);
}


////////////////////////////////////////////////////////////////////////////////
template<typename PointT, typename MatScalar> int
sq::SuperquadricFittingLM<PointT, MatScalar>::OptimizationFunctor::df (const VectorX &xvec,
                                                                       Eigen::Matrix<MatScalar, Eigen::Dynamic, Eigen::Dynamic> &fjac) const
{
  const Cloud &points = *estimator_->input_prealigned_;
  double e1 = xvec[0],
      e2 = xvec[1],
      a = xvec[2],
      b = xvec[3],
      c = xvec[4],
      tx = xvec[5],
      ty = xvec[6],
      tz = xvec[7],
      ax = xvec[8],
      ay = xvec[9],
      az = xvec[10];


  for (int i = 0; i < values (); ++i)
  {
    MatScalar x = points[i].x,
        y = points[i].y,
        z = points[i].z;

    superquadric_derivative<MatScalar> (x, y, z,
                                        xvec[0], xvec[1], xvec[2], xvec[3], xvec[4], xvec[5], xvec[6], xvec[7], xvec[8], xvec[9], xvec[10],
                                        fjac (i, 0), fjac (i, 1), fjac (i, 2), fjac (i, 3), fjac (i, 4), fjac (i, 5), fjac (i, 6), fjac (i, 7), fjac (i, 8), fjac (i, 9), fjac (i, 10));
  }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
template<typename PointT, typename MatScalar> int
sq::SuperquadricFittingLM<PointT, MatScalar>::OptimizationFunctorWithIndices::operator() (const VectorX &xvec,
                                                                                          VectorX &fvec) const
{
  /// TODO fill up the indices version
  const Cloud &points = *estimator_->input_prealigned_;
  const std::vector<int> &indices = *estimator_->indices_;

  for (int i = 0; i < values (); ++i)
  {
    MatScalar x = points[indices[i]].x,
        y = points[indices[i]].y,
        z = points[indices[i]].z;

    MatScalar e1 = xvec[0],
        e2 = xvec[1],
        a = xvec[2],
        b = xvec[3],
        c = xvec[4];

    fvec[i] = pow (pow (x / a, 2 / e2) + pow (y / b, 2 / e2), e2 / e1) + pow (z / c, 2 / e1);
  }
  return (0);
}


