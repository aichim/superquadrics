#pragma once

#include "fit_superquadric_lm.h"
#include "superquadric_formulas.h"
#include <pcl/common/centroid.h>
#include <unsupported/Eigen/NonLinearOptimization>


template<typename PointT, typename MatScalar>
SuperquadricFittingLM<PointT, MatScalar>::SuperquadricFittingLM ()
{
}

template<typename PointT, typename MatScalar> double
SuperquadricFittingLM<PointT, MatScalar>::fit (VectorX &parameters,
                                               Eigen::Matrix<MatScalar, 4, 4> &transformation) const
{
  int n_unknowns = 11;
  /// e1, e2, a, b, c
  VectorX xvec (n_unknowns);
  xvec[0] = xvec[1] = xvec[2] = xvec[3] = xvec[4] = 1.;
  xvec[5] = xvec[6] = xvec[7] = xvec[8] = xvec[9] = xvec[10] = 0.;

  // Initialize the parameters based on the input cloud
  // The center of the superquadric should be the centroid of the cloud
  Eigen::Vector4d centroid;
  pcl::compute3DCentroid (*input_, centroid);
  xvec[5] = -centroid[0];
  xvec[6] = -centroid[1];
  xvec[7] = -centroid[2];


  // The scales are the dimensions of the cloud on each axis.
  Eigen::Vector3f min (Eigen::Vector3f (1.f, 1.f, 1.f) * std::numeric_limits<float>::max ()),
      max (Eigen::Vector3f (1.f, 1.f, 1.f) * std::numeric_limits<float>::min ());
  for (size_t i = 0; i < input_->size (); ++i)
  {
    const Eigen::Vector3f &p = (*input_)[i].getVector3fMap ();

    for (size_t j = 0; j < 3; ++j)
    {
      if (p[j] < min[j])
        min[j] = p[j];

      if (p[j] > max[j])
        max[j] = p[j];
    }
  }

  xvec[2] = (max[0] - min[0]) / 2.;
  xvec[3] = (max[1] - min[1]) / 2.;
  xvec[4] = (max[2] - min[2]) / 2.;

  std::cout << "initial parameters:\n" << xvec << std::endl;

  /// TODO take care of the case with indices

  OptimizationFunctor functor (input_->size (), this);

  /// TODO replace numerical diff with actual df - use Maple
//  Eigen::NumericalDiff<OptimizationFunctor> num_diff (functor);
//  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<OptimizationFunctor>, MatScalar> lm (num_diff);
  Eigen::LevenbergMarquardt<OptimizationFunctor, MatScalar> lm (functor);
  /*int info = */lm.minimize (xvec);

  std::cout << xvec << std::endl;

  parameters = xvec;







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

  MatScalar final_error = computeSuperQuadricError<PointT, MatScalar> (input_,
                                                                       xvec[0], xvec[1], xvec[2], xvec[3], xvec[4],
                                                                       transformation);
  return (final_error);
}


template<typename PointT, typename MatScalar> int
SuperquadricFittingLM<PointT, MatScalar>::OptimizationFunctor::operator () (const VectorX &xvec,
                                                                            VectorX &fvec) const
{
  const Cloud &points = *estimator_->input_;

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

//    std::cout << xyz << " " << xyz_tr << std::endl;
//    std::cout << transformation << std::endl;


//    double term_1 = pow (fabs(xyz_tr[0] / a), 2./e2);
//    double term_2 = pow (fabs(xyz_tr[1] / b), 2./e2);
//    double term_3 = pow (fabs(xyz_tr[2] / c), 2./e1);
//    double superellipsoid_f = pow (fabs(term_1 + term_2), e2/e1) + term_3;

    double op = (Eigen::Matrix<MatScalar, 3, 1> (xvec[5], xvec[6], xvec[7]) -
        Eigen::Matrix<MatScalar, 3, 1> (xyz_tr[0], xyz_tr[1], xyz_tr[2])).norm ();

    fvec[i] = op *superquadric_function (xyz_tr[0], xyz_tr[1], xyz_tr[2],
                                     e1, e2, a, b, c);

//    double op = Eigen::Matrix<MatScalar, 3, 1> (xyz_tr[0], xyz_tr[1], xyz_tr[2]).norm ();


//    fvec[i] = /*op */ (pow (superellipsoid_f, e1 / 2.) - 1.) * pow (a*b*c, 0.25);
//    PCL_INFO ("fvec[%ld] = %f\n", i, fvec[i]);
  }
  return (0);
}


template<typename PointT, typename MatScalar> int
SuperquadricFittingLM<PointT, MatScalar>::OptimizationFunctor::df (const VectorX &xvec,
                                                                   Eigen::Matrix<MatScalar, Eigen::Dynamic, Eigen::Dynamic> &fjac) const
{
//  PCL_INFO ("%d %d %d %d\n", xvec.rows (), xvec.cols (), fjac.rows (), fjac.cols ());

  const Cloud &points = *estimator_->input_;
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




template<typename PointT, typename MatScalar> int
SuperquadricFittingLM<PointT, MatScalar>::OptimizationFunctorWithIndices::operator() (const VectorX &xvec,
                                                                                      VectorX &fvec) const
{
  const Cloud &points = *estimator_->input_;
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


