#pragma once

#include "fit_superquadric_lm.h"

#include <unsupported/Eigen/NonLinearOptimization>


template<typename PointT, typename MatScalar>
SuperquadricFittingLM<PointT, MatScalar>::SuperquadricFittingLM ()
{

}

template<typename PointT, typename MatScalar> void
SuperquadricFittingLM<PointT, MatScalar>::fit (VectorX &parameters,
                                               Eigen::Matrix<MatScalar, 4, 4> &transformation) const
{
  /// TODO 11 to include the transformation
  int n_unknowns = 5;
  /// e1, e2, a, b, c
  VectorX xvec (n_unknowns);
  xvec[0] = 1.;
  xvec[1] = 1.;
  xvec[2] = 1.;
  xvec[3] = 1.;
  xvec[4] = 1.;
//  xvec[0] = 1.0;//0.8;
//  xvec[1] = 0.8;//1.0;
//  xvec[2] = 0.3;
//  xvec[3] = 0.3;
//  xvec[4] = 0.4;

  /// TODO take care of the case with indices

  OptimizationFunctor functor (input_->size (), this);

  /// TODO replace numerical diff with actual df - use Maple
  Eigen::NumericalDiff<OptimizationFunctor> num_diff (functor);
  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<OptimizationFunctor>, MatScalar> lm (num_diff);
  int info = lm.minimize (xvec);

  std::cout << xvec << std::endl;

  parameters = xvec;
}


template<typename PointT, typename MatScalar> int
SuperquadricFittingLM<PointT, MatScalar>::OptimizationFunctor::operator () (const VectorX &xvec,
                                                                            VectorX &fvec) const
{
  const Cloud &points = *estimator_->input_;

  for (int i = 0; i < values (); ++i)
  {
    MatScalar x = points[i].x,
        y = points[i].y,
        z = points[i].z;

    MatScalar e1 = xvec[0],
        e2 = xvec[1],
        a = xvec[2],
        b = xvec[3],
        c = xvec[4];

    PCL_INFO ("x %f  y %f  z %f\n", x, y, z);
    PCL_INFO ("e1 %f  e2 %f  a %f  b %f  c %f\n",
              e1, e2, a, b, c);

    double term_1 = pow (fabs(x/a), 2./e2);
    double term_2 = pow (fabs(y/b), 2./e2);
    double term_3 = pow (fabs(z/c), 2./e1);
    double superellipsoid_f = pow (fabs(term_1 + term_2), e2/e1) + term_3;

    fvec[i] = (pow (superellipsoid_f, e1 / 2.) - 1.) * pow (a*b*c, 0.25);
//    fvec[i] = superellipsoid_f - 1.0;
    PCL_INFO ("fvec[%d] = %f\n", i, fvec[i]);
  }
  return (0);
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
